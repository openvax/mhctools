"""Named sequence-context cleavage assessments with lossless epitope overlays.

Input records and output manifests are JSON-compatible mappings. All intervals
are zero-based and half-open; bond b divides sequence[:b] and sequence[b:].
Scenario labels describe the user's question, not a learned cell-type model.
"""

from copy import deepcopy
from dataclasses import asdict, replace
import html
import json
from pathlib import Path

from .cleavage import CleavageInput, CleavageModel, CleavageResult
from .peptidases import cleavage_models, get_cleavage_model
from .substrate_reference import PeptidaseSubstrateReference


SCENARIO_COMPARTMENTS = {
    "tumor": ("cytosol", "er"),
    "apc": ("cytosol", "er", "endosome"),
    "extracellular": ("extracellular", "serum", "plasma"),
}

_COVERAGE_GAPS = {
    "tumor": [
        "No cell-specific enzyme abundance, IFN response, MHC protection or processing kinetics.",
        "Cytosolic peptidase rules are partial; THOP1/NLN entries are exact-source observations.",
    ],
    "apc": [
        "ITCell B/S internal and H initial-trimming profiles are assay-scoped; no CatL or legumain/AEP predictor (#470).",
        "IRAP is an exact-substrate reference, not prediction on new vaccine sequences.",
        "No uptake, endosomal escape, pH/activation or primary-DC presentation calibration.",
    ],
    "extracellular": [
        "Enzyme-location annotations do not establish activity or calibration in serum/plasma.",
        "No competing-reaction kinetics, intact-peptide half-life or systemic clearance prediction.",
    ],
}


def _json_copy(value):
    # This also refuses NaN, custom objects and non-JSON evidence before inference.
    return json.loads(json.dumps(value, allow_nan=False))


def _identifier(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("%s must be a nonempty string" % name)
    return value


def _interval(value, length):
    start, end = value["start"], value["end"]
    if type(start) is not int or type(end) is not int or not 0 <= start < end <= length:
        raise ValueError("Intervals must be zero-based, half-open and inside the sequence")
    return start, end


def normalize_cleavage_input(value):
    """Validate a named full sequence or a named peptide with explicit flanks.

    Extra JSON source/evidence fields are retained. ``source_start`` refers to
    the start of the entire supplied sequence, including flanks. A native
    window has unknown molecular termini. A construct's terminal chemistry
    must be supplied to enable chemistry-dependent models.
    """
    value = _json_copy(value)
    _identifier(value.get("id"), "input id")
    scope = value.setdefault("scope", "native_window")
    if scope not in ("native_window", "protein", "construct"):
        raise ValueError("scope must be native_window, protein or construct")
    if "sequence" not in value:
        peptide = CleavageInput(value["peptide"]).sequence
        n_flank, c_flank = value["n_flank"], value["c_flank"]
        if not isinstance(n_flank, str) or not isinstance(c_flank, str):
            raise ValueError("Flanks must be strings, including explicit empty strings")
        value["sequence"] = n_flank + peptide + c_flank
        value.setdefault("epitopes", [{
            "id": value["id"], "start": len(n_flank),
            "end": len(n_flank) + len(peptide), "sequence": peptide}])
    elif any(key in value for key in ("peptide", "n_flank", "c_flank")):
        # Normalized peptide-mode records retain their original fields.
        if value.get("n_flank", "") + value.get("peptide", "") + value.get("c_flank", "") != value["sequence"]:
            raise ValueError("Sequence contradicts supplied peptide and flanks")
    value.setdefault("n_term", "unknown")
    value.setdefault("c_term", "unknown")
    value.setdefault("source_start", 0)
    value.setdefault("source_id", value["id"])
    parent = _peptide(value)
    if scope == "native_window" and (parent.n_term != "unknown" or parent.c_term != "unknown"):
        raise ValueError("Native windows do not establish molecular termini")
    epitope_ids = []
    for epitope in value.setdefault("epitopes", []):
        epitope_ids.append(_identifier(epitope.get("id"), "epitope id"))
        start, end = _interval(epitope, len(parent.sequence))
        expected = parent.sequence[start:end]
        if epitope.setdefault("sequence", expected) != expected:
            raise ValueError("Epitope sequence contradicts its interval")
    if len(epitope_ids) != len(set(epitope_ids)):
        raise ValueError("Epitope IDs must be unique within an input")
    fragment_ids = []
    for fragment in value.setdefault("fragments", []):
        fragment_ids.append(_identifier(fragment.get("id"), "fragment id"))
        _identifier(fragment.get("assumption"), "fragment production assumption")
        start, end = _interval(fragment, len(parent.sequence))
        parent.fragment(start, end, n_term=fragment["n_term"], c_term=fragment["c_term"])
    if len(fragment_ids) != len(set(fragment_ids)):
        raise ValueError("Fragment IDs must be unique within an input")
    return value


def _peptide(value):
    return CleavageInput(**{key: value[key] for key in (
        "sequence", "n_term", "c_term", "source_id", "source_start")})


def _normalize_scenarios(scenarios):
    result = _json_copy(list(scenarios))
    names = []
    compartments = {c for m in cleavage_models() for c in m.compartments}
    for scenario in result:
        names.append(_identifier(scenario.get("id"), "scenario id"))
        context = scenario["context"]
        if context not in SCENARIO_COMPARTMENTS:
            raise ValueError("Scenario context must be tumor, apc or extracellular")
        selected = scenario.setdefault("compartments", list(SCENARIO_COMPARTMENTS[context]))
        if not isinstance(selected, list) or not selected or set(selected) - compartments:
            raise ValueError("Scenario requires known compartments")
        models = scenario.get("models")
        if not isinstance(models, list) or not models or any(not isinstance(m, str) for m in models):
            raise ValueError("Each scenario requires an explicit nonempty model list")
        if len(set(models)) != len(models):
            raise ValueError("Duplicate model in scenario")
        states = scenario.setdefault("enzyme_states", {})
        if not isinstance(states, dict) or any(v not in ("active", "inactive", "zymogen", "unknown") for v in states.values()):
            raise ValueError("Invalid enzyme_states")
    if len(names) != len(set(names)):
        raise ValueError("Scenario IDs must be unique")
    return result


def cleavage_overlays(value, result, *, fragment=None):
    """Map native evidence to epitope interiors and separate N/C boundaries.

    No threshold or aggregate probability is applied. A fragment boundary is
    not an experimentally produced terminus; its production assumption stays
    on every overlay. Unreturned bonds remain explicit unassessed entries.
    """
    offset = result.peptide.source_start - value["source_start"]
    by_bond = {site["bond"] + offset: site for site in result.to_dict()["sites"]}

    def observation(bond):
        row = {"bond": bond, "source_bond": value["source_start"] + bond}
        if bond in (0, len(value["sequence"])):
            return dict(row, status="sequence_endpoint")
        if bond in by_bond:
            return dict(by_bond[bond], **row)
        reason = result.unsupported_reason or "Model did not assess this bond"
        if fragment and not fragment["start"] < bond < fragment["end"]:
            reason = "Bond is outside this conditional fragment or at its endpoint"
        return dict(row, status="unassessed", reason=reason)

    return [{
        "epitope": deepcopy(epitope),
        "conditional_on": deepcopy(fragment),
        "substrate_observation": result.substrate_observation,
        "n_boundary": observation(epitope["start"]),
        "c_boundary": observation(epitope["end"]),
        "internal": [observation(b) for b in range(epitope["start"] + 1, epitope["end"])],
    } for epitope in value["epitopes"]]


def predict_cleavage_batch(inputs, scenarios, *, predictors=None, reference_panels=(),
                           raise_on_error=False):
    """Assess named occurrences under explicitly selected biological scenarios.

    Parameters
    ----------
    inputs : iterable of dict
        Named sequences/peptide occurrences accepted by normalize_cleavage_input.
    scenarios : iterable of dict
        Named tumor/APC/extracellular questions with explicit model lists.
        Compartments filter enzyme locations, not experimental validation.
    predictors : mapping, optional
        Already constructed canonical predictors keyed by exact model name.
        Useful for user-managed assets and source-observation adapters.
    reference_panels : iterable of dict
        Experimental source catalogs with ``model`` and ``cases`` fields,
        using the PeptidaseSubstrateReference contract. One named panel per
        assay/condition; exact sequence/chemistry lookup never extrapolates.
        Panels are retained in the output so imported observations round-trip.
    raise_on_error : bool
        Raise backend errors instead of retaining failed assessment records.

    Returns
    -------
    dict
        Versioned JSON manifest retaining inputs, scenarios, model results,
        conditional fragments, epitope overlays and declared coverage gaps.
    """
    inputs = [normalize_cleavage_input(value) for value in inputs]
    if len({v["id"] for v in inputs}) != len(inputs):
        raise ValueError("Input IDs must be unique")
    scenarios = _normalize_scenarios(scenarios)
    predictors = {} if predictors is None else dict(predictors)
    reference_panels = _json_copy(list(reference_panels))
    catalog = {m.name: m for m in cleavage_models(include_optional=True)}
    for panel in reference_panels:
        predictor = PeptidaseSubstrateReference(panel["model"], panel["cases"])
        name = predictor.model.name
        if name in catalog or name in predictors:
            raise ValueError("Imported reference panel duplicates model name %r" % name)
        predictors[name] = predictor
    catalog.update({name: p.model for name, p in predictors.items()})
    # Validate model selection before running any optional backend.
    for scenario in scenarios:
        for name in scenario["models"]:
            if name not in catalog:
                raise ValueError("Unknown cleavage model %r" % name)
            if not set(catalog[name].compartments) & set(scenario["compartments"]):
                raise ValueError("Model %s is outside the scenario compartments" % name)
        selected_enzymes = {catalog[name].enzyme for name in scenario["models"]}
        if set(scenario["enzyme_states"]) - selected_enzymes:
            raise ValueError("Enzyme state supplied without selecting that enzyme")
    groups, tasks, assessments = {}, [], []
    for scenario in scenarios:
        for value in inputs:
            parent = _peptide(value)
            for fragment in [None] + value["fragments"]:
                peptide = parent if fragment is None else parent.fragment(
                    fragment["start"], fragment["end"],
                    n_term=fragment["n_term"], c_term=fragment["c_term"])
                for name in scenario["models"]:
                    state = scenario["enzyme_states"].get(catalog[name].enzyme)
                    model_input = replace(peptide, source_id=None, source_start=0)
                    model_key = (name, state)
                    groups.setdefault(model_key, {}).setdefault(model_input, None)
                    row = dict(input_id=value["id"], scenario_id=scenario["id"],
                               model=name, fragment=deepcopy(fragment), result=None,
                               requested_model=asdict(catalog[name]),
                               overlays=[], error=None,
                               context_application="Compartments filter locations; scenario labels and "
                               "extra annotations do not change model parameters or establish assay calibration")
                    tasks.append((row, value, peptide, model_key, model_input))
    for (name, state), pending in groups.items():
        try:
            if name in predictors:
                if state is not None:
                    raise ValueError("Configure activation on a supplied predictor directly")
                predictor = predictors[name]
            else:
                predictor = get_cleavage_model(name, enzyme_state=state)
            predict_many = getattr(predictor, "predict_many", None)
            if predict_many is not None:
                returned = tuple(predict_many(tuple(pending)))
                if len(returned) != len(pending):
                    raise ValueError("Batch predictor returned the wrong number of results")
                pending.update(zip(pending, returned))
            else:
                for peptide in pending:
                    try:
                        pending[peptide] = predictor.predict(peptide)
                    except Exception as error:
                        if raise_on_error:
                            raise
                        pending[peptide] = error
        except Exception as error:
            if raise_on_error:
                raise
            pending.update({p: error for p in pending})
    for row, value, peptide, model_key, model_input in tasks:
        try:
            prediction = groups[model_key][model_input]
            if isinstance(prediction, Exception):
                raise prediction
            if prediction.peptide != model_input or prediction.model.name != row["model"]:
                raise ValueError("Predictor changed input or model identity")
            result = replace(prediction, peptide=peptide)
            row.update(status="unsupported" if result.unsupported_reason else "assessed",
                       result=result.to_dict(),
                       overlays=cleavage_overlays(value, result, fragment=row["fragment"]))
        except Exception as error:
            if raise_on_error:
                raise
            row.update(status="failed", error="%s: %s" % (type(error).__name__, error))
        assessments.append(row)
    return _json_copy(dict(schema_version=1, inputs=inputs, scenarios=scenarios,
                          reference_panels=reference_panels,
                          assessments=assessments, coverage_gaps={
                              s["id"]: _COVERAGE_GAPS[s["context"]] for s in scenarios}))


def load_cleavage_batch(path):
    """Read a saved assessment without inference, validating result coordinates."""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if value.get("schema_version") != 1:
        raise ValueError("Unsupported cleavage batch schema version")
    inputs = {v["id"]: normalize_cleavage_input(v) for v in value["inputs"]}
    scenarios = {s["id"]: s for s in _normalize_scenarios(value["scenarios"])}
    for panel in value.get("reference_panels", ()):
        PeptidaseSubstrateReference(panel["model"], panel["cases"])
    if len(inputs) != len(value["inputs"]):
        raise ValueError("Duplicate input IDs")
    keys = set()
    for row in value["assessments"]:
        source = inputs[row["input_id"]]
        scenario = scenarios[row["scenario_id"]]
        if row["model"] not in scenario["models"]:
            raise ValueError("Assessment model was not requested")
        if CleavageModel(**row["requested_model"]).name != row["model"]:
            raise ValueError("Requested model metadata contradicts its name")
        fragment = row["fragment"]
        if fragment is not None and fragment not in source["fragments"]:
            raise ValueError("Unknown conditional fragment")
        key = (source["id"], scenario["id"], row["model"], fragment["id"] if fragment else None)
        if key in keys:
            raise ValueError("Duplicate assessment")
        keys.add(key)
        if row["status"] == "failed":
            if not row["error"] or row["result"] is not None or row["overlays"]:
                raise ValueError("Failed assessment must retain an error without results")
            continue
        result = CleavageResult.from_dict(row["result"])
        expected = _peptide(source)
        if fragment:
            expected = expected.fragment(fragment["start"], fragment["end"],
                                         n_term=fragment["n_term"], c_term=fragment["c_term"])
        if result.peptide != expected or result.model.name != row["model"]:
            raise ValueError("Saved result contradicts input/model identity")
        status = "unsupported" if result.unsupported_reason else "assessed"
        if row["status"] != status or row["error"] is not None:
            raise ValueError("Saved assessment status contradicts evidence")
        if row["overlays"] != _json_copy(cleavage_overlays(source, result, fragment=fragment)):
            raise ValueError("Saved overlays contradict cleavage evidence")
    expected_count = sum(len(s["models"]) * sum(1 + len(v["fragments"]) for v in inputs.values())
                         for s in scenarios.values())
    if len(keys) != expected_count:
        raise ValueError("Saved assessment is missing requested results")
    return value


def write_cleavage_batch(report, path, *, html_path=None):
    """Save full precision JSON and optionally a self-contained evidence report."""
    if html_path is not None and Path(path).resolve() == Path(html_path).resolve():
        raise ValueError("JSON and HTML output paths must be distinct")
    Path(path).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if html_path is None:
        return
    def escape(value):
        return html.escape(str(value))
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'>",
             "<title>Cleavage assessment</title><style>body{font:16px system-ui;max-width:1100px;"
             "margin:40px auto;padding:0 20px;color:#192a36}table{border-collapse:collapse;width:100%}"
             "th,td{border:1px solid #cbd5db;padding:8px;text-align:left}pre{white-space:pre-wrap;"
             "overflow-wrap:anywhere}summary{cursor:pointer}h2{margin-top:2em}</style>",
             "<h1>Cleavage evidence by sequence and scenario</h1>",
             "<p>Native scores, recognition rules and source observations are separate evidence. "
             "Unassessed bonds are not evidence of resistance. No combined survival or presentation "
             "probability is calculated. Intervals are zero-based, half-open.</p>"]
    for scenario in report["scenarios"]:
        parts.append("<h2>%s (%s)</h2><ul>%s</ul>" % (
            escape(scenario["id"]), escape(scenario["context"]),
            "".join("<li>%s</li>" % escape(gap) for gap in report["coverage_gaps"][scenario["id"]])))
        for row in report["assessments"]:
            if row["scenario_id"] != scenario["id"]:
                continue
            parts.append("<h3>%s · %s · %s</h3>" % (
                escape(row["input_id"]), escape(row["model"]), escape(row["status"])))
            if row["fragment"]:
                parts.append("<p><strong>Conditional fragment:</strong> %s</p>" % escape(row["fragment"]))
            if row["error"]:
                parts.append("<p>%s</p>" % escape(row["error"]))
            for overlay in row["overlays"]:
                parts.append("<h4>Epitope %s [%d, %d)</h4><table><tr><th>Region</th>"
                             "<th>Bond</th><th>Evidence</th><th>Native score</th></tr>" % (
                                 escape(overlay["epitope"]["id"]), overlay["epitope"]["start"],
                                 overlay["epitope"]["end"]))
                for region, sites in (("N boundary", [overlay["n_boundary"]]),
                                      ("Internal", overlay["internal"]),
                                      ("C boundary", [overlay["c_boundary"]])):
                    for site in sites:
                        parts.append("<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td></tr>" % (
                            region, site["bond"], escape(site["status"]),
                            escape(site.get("score") if site.get("score") is not None else "—")))
                parts.append("</table>")
            parts.append("<details><summary>Full evidence, conditions and provenance</summary><pre>%s</pre>"
                         "</details>" % escape(json.dumps(row, indent=2)))
    parts.append("<h2>Original inputs and imported evidence</h2><pre>%s</pre>" %
                 escape(json.dumps(report["inputs"], indent=2)))
    additional = {key: value for key, value in report.items()
                  if key not in ("inputs", "scenarios", "assessments", "coverage_gaps")}
    parts.append("<details><summary>Additional assay and validation evidence</summary><pre>%s</pre>"
                 "</details></html>" % escape(json.dumps(additional, indent=2)))
    Path(html_path).write_text("\n".join(parts), encoding="utf-8")
