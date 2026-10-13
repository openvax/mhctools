"""Dose-normalized, per-core assessments of supplied vaccine windows.

This is a batch adapter over the trajectory model, not a sequence-to-human-PK
predictor. Every core occurrence is retained; marginal runs are not a joint
multi-epitope model and must not be added as fractions of administered molecules.
"""

from copy import deepcopy
from dataclasses import asdict
import math

from Bio import __version__ as biopython_version
from Bio.SeqUtils import molecular_weight
from mhcgnomes import parse

from .peptide_input import PeptideInput, sequence_only_chemistry_error
from .serum_contributions import EnzymeCutRate
from .serum_degradation import DegradationTarget
from .vaccine_exposure import summarize_core_exposure
from .vaccine_trajectory import (
    VaccineCleavageRates, VaccineDose, VaccineTrajectoryInput,
    simulate_vaccine_trajectory,
)
from .vaccine_trajectory_report import summarize_vaccine_trajectory, _number, _text


def fragment_cleavage_table(sequence, assessments):
    """Validate an explicit fragment table and return a trajectory callback.

    Coordinates refer to the original input. Cut bonds in each channel refer
    to the current fragment. Unlisted fragments remain unassessed.
    """
    table = {}
    for item in assessments:
        start, end = item["start"], item["end"]
        if (type(start) is not int or type(end) is not int or
                not 0 <= start < end <= len(sequence) or
                item["sequence"] != sequence[start:end]):
            raise ValueError("Cleavage assessment must match its original-parent interval")
        key = item["compartment"], start, end
        if key in table:
            raise ValueError("Duplicate fragment cleavage assessment")
        table[key] = VaccineCleavageRates(
            tuple(EnzymeCutRate(**channel) for channel in item["channels"]),
            item["basis"], item["source"])
    return lambda compartment, fragment: table.get((compartment, fragment.start, fragment.end))


def _dose(candidate, peptide):
    if candidate.get("dose") is not None:
        return VaccineDose(**candidate["dose"])
    if candidate["delivery"] != "slp":
        return None
    if sequence_only_chemistry_error(peptide):
        return None
    return VaccineDose(
        amount=1, unit="ug", basis="assumed",
        source="Requested 1 ug comparison; average molecular weight of a free linear unmodified peptide, Biopython %s" % biopython_version,
        molecular_weight_g_per_mol=molecular_weight(peptide.sequence, seq_type="protein"))


def _aggregate(cores, delivery):
    result = {
        "status": "unassessed", "score_weighted_mean_node_entry_per_input": None,
        "definition": "Sum(score contribution * per-core nodal APC arrivals per input molecule) / sum(score contribution)",
        "units": "mean core arrivals per input peptide" if delivery == "slp" else "mean core arrivals per input RNA",
        "any_core_molecule_fraction": None,
        "any_core_status": "not_identified_from_separate_core_marginals",
        "limitations": "A score-weighted descriptor, not a fraction of unique molecules carrying any core, not a joint MHC competition model, and not vaccine efficacy",
    }
    weights = [core.get("score_contribution") for core in cores]
    if not cores or any(w is None or w < 0 for w in weights) or math.fsum(weights) <= 0:
        result["reason"] = "Requires a nonempty set of nonnegative score contributions with positive total"
        return result
    if any(core["status"] != "conditional_scenario" for core in cores):
        result["reason"] = "Every contributing core must be assessed; missing cores are not assigned zero"
        return result
    result.update(
        status="conditional_scenario",
        score_weighted_mean_node_entry_per_input=math.fsum(
            w * core["exposure"]["node_apc_entry"]["copies_per_input_molecule"]
            for w, core in zip(weights, cores)) / math.fsum(weights))
    return result


def assess_vaccine_windows(request):
    """Assess every supplied candidate and scored core under explicit kinetics.

    Parameters
    ----------
    request : dict
        Schema 1 with ``scenario``, ``times_hours`` and ``candidates``. Each
        candidate has ``id``, ``peptide``, ``delivery`` and ``cores`` (original
        window start/end, label, allele, MHC class, optional score contribution).
        ``kinetics`` maps MHC class to rates, a TAP/class-II length bound and a
        fragment cleavage table. Candidate kinetics replace the corresponding
        request-level class profile. Each profile must apply to that delivery.
        SLP dose defaults to 1 ug of the actual free linear peptide. RNA requires
        an explicit dose, using transcript molecular weight for a mass dose.

    Returns
    -------
    dict
        Per-window, per-core summaries and full trajectory audits. Missing
        parameters remain unassessed. No ranking or physiological fit occurs.
    """
    if request.get("schema_version") != 1:
        raise ValueError("Window schema_version must be 1")
    if not isinstance(request.get("scenario"), str) or not request["scenario"].strip():
        raise ValueError("A nonempty scenario label is required")
    times = list(request["times_hours"])
    if not times or any(isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(t) or t < 0 for t in times):
        raise ValueError("times_hours must be a nonempty list of finite nonnegative hours")
    from . import __version__
    output = {
        "schema_version": 1, "mhctools_version": __version__, "scenario": request["scenario"],
        "horizon_hours": max(times), "times_hours": times, "candidates": [],
        "coordinate_system": "zero-based half-open on each full candidate input",
        "input_request": deepcopy(request),
        "interpretation": "Conditional scenarios, not calibrated predictions of human IM bioavailability or circulation half-life",
    }
    seen = set()
    for candidate in request["candidates"]:
        name = candidate["id"]
        if not isinstance(name, str) or not name.strip() or name in seen:
            raise ValueError("Candidate ids must be unique nonempty strings")
        seen.add(name)
        delivery = candidate["delivery"]
        if delivery not in ("slp", "secreted_mrna", "cytosolic_mrna"):
            raise ValueError("Unknown delivery: %s" % delivery)
        raw = candidate["peptide"]
        peptide = PeptideInput.from_dict(raw) if isinstance(raw, dict) else PeptideInput(raw)
        dose = _dose(candidate, peptide)
        chemistry_error = sequence_only_chemistry_error(peptide)
        window = {"id": name, "peptide_input": peptide.to_dict(), "delivery": delivery,
                  "dose": asdict(dose) if dose else None, "initial_copies": dose.copies if dose else None,
                  "metadata": candidate.get("metadata", {}), "cores": []}
        identities = set()
        for core in candidate["cores"]:
            target = DegradationTarget(core["label"], core["start"], core["end"])
            target.validate(peptide.sequence)
            if core.get("sequence", peptide.sequence[target.start:target.end]) != peptide.sequence[target.start:target.end]:
                raise ValueError("Core sequence does not match its candidate interval")
            allele = parse(core["allele"], raise_on_error=False)
            mhc_class = core["mhc_class"]
            if mhc_class not in ("I", "II") or allele is None or not (allele.is_class1 if mhc_class == "I" else allele.is_class2):
                raise ValueError("Each core allele must parse and match its MHC class")
            key = (target.start, target.end, allele.to_string())
            if key in identities:
                raise ValueError("Duplicate core occurrence/allele in candidate " + name)
            identities.add(key)
            score = core.get("score_contribution")
            if score is not None and (isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score)):
                raise ValueError("score_contribution must be finite or absent")
            record = {**core, "allele": allele.to_string(), "sequence": peptide.sequence[target.start:target.end],
                      "status": "unassessed", "missing_inputs": [], "trajectory": None,
                      "exposure": None, "presentation": None}
            profiles = {**request.get("kinetics", {}), **candidate.get("kinetics", {})}
            profile = profiles.get(mhc_class, {})
            gate = "tap_max_length" if mhc_class == "I" else "class_ii_max_length"
            if chemistry_error:
                record["missing_inputs"].append(chemistry_error)
            if dose is None:
                record["missing_inputs"].append("Explicit dose required; RNA mass requires actual transcript molecular weight")
            if profile.get(gate) is None:
                record["missing_inputs"].append("kinetics.%s.%s" % (mhc_class, gate))
            if delivery == "secreted_mrna" and candidate.get("signal_end") is None:
                record["missing_inputs"].append("signal_end for the full secreted construct")
            if not record["missing_inputs"]:
                model = VaccineTrajectoryInput.from_dict({
                    "schema_version": 1, "scenario": request["scenario"], "peptide": peptide.to_dict(),
                    "delivery": delivery, "target": asdict(target), "allele": allele.to_string(),
                    "mhc_class": mhc_class, "dose": asdict(dose), "signal_end": candidate.get("signal_end"),
                    "administration": candidate.get("administration", request.get("administration", "IM")),
                    "organism": candidate.get("organism", request.get("organism")),
                    gate: profile[gate], "rates": profile.get("rates", {}),
                })
                trajectory = simulate_vaccine_trajectory(
                    model, times, fragment_cleavage_table(peptide.sequence, profile.get("cleavage_assessments", [])))
                record.update(status=trajectory["status"], trajectory=trajectory,
                              exposure=summarize_core_exposure(trajectory),
                              presentation=summarize_vaccine_trajectory(trajectory)["presentation"],
                              missing_inputs=trajectory["missing_kinetics"])
            window["cores"].append(record)
        window["aggregate"] = _aggregate(window["cores"], delivery)
        window["status"] = ("no_scored_cores" if not window["cores"] else
                            "conditional_scenario" if all(c["status"] == "conditional_scenario" for c in window["cores"]) else "unassessed")
        output["candidates"].append(window)
    return output


def vaccine_window_rows(result):
    """Flatten to one row per core occurrence/allele; preserve absent estimates."""
    rows = []
    horizon = result["horizon_hours"]
    for window in result["candidates"]:
        for core in window["cores"]:
            exposure = core["exposure"] or {}
            node, blood = exposure.get("node_apc_entry") or {}, exposure.get("blood") or {}
            decline = blood.get("post_peak_half_decline", {})
            interval = decline.get("hours_after_sampled_peak_interval") or [None, None]
            presentation = (core["presentation"] or {}).get("node", {})
            dose = window["dose"] or {}
            initial = window["initial_copies"]
            denominator = initial * horizon if initial is not None and horizon > 0 else None
            slp_fraction = node.get("fraction_of_dosed_peptide_copies")
            rows.append({
                "candidate_id": window["id"], "window_sequence": window["peptide_input"]["sequence"],
                "delivery": window["delivery"], "core_label": core["label"], "core_sequence": core["sequence"],
                "core_start": core["start"], "core_end": core["end"], "allele": core["allele"],
                "score_contribution": core.get("score_contribution"), "status": core["status"],
                "horizon_hours": horizon, "dose_amount": dose.get("amount"), "dose_unit": dose.get("unit"),
                "initial_copies": initial,
                "node_apc_entry_copies": node.get("cumulative_core_bearing_copies"),
                "node_apc_fraction_of_dosed_peptide": slp_fraction,
                "node_apc_percent_of_dosed_peptide": 100 * slp_fraction if slp_fraction is not None else None,
                "node_apc_entries_per_input": node.get("copies_per_input_molecule"),
                "node_pmhc_auc_copy_hours": presentation.get("auc_copy_hours"),
                "node_pmhc_mean_copies_per_input": (
                    presentation["auc_copy_hours"] / denominator
                    if denominator is not None and presentation.get("auc_copy_hours") is not None else None),
                "blood_core_auc_copy_hours": blood.get("auc_copy_hours"),
                "blood_core_mean_copies_per_input": (
                    blood["auc_copy_hours"] / denominator
                    if denominator is not None and blood.get("auc_copy_hours") is not None else None),
                "blood_free_core_copies_at_horizon": blood.get("exact_free_core_copies_at_horizon"),
                "blood_precursor_copies_at_horizon": blood.get("longer_precursor_copies_at_horizon"),
                "blood_half_decline_status": decline.get("status"),
                "blood_half_decline_lower_hours": interval[0], "blood_half_decline_upper_hours": interval[1],
                "intrinsic_circulation_half_life_hours": None,
                "missing_inputs": "; ".join(core["missing_inputs"]),
            })
    return rows


def format_vaccine_window_report(result):
    """Readable endpoint table with dose bases, gaps and score-weighted summaries."""
    lines = ["# Vaccine window exposure", "", _text(result["scenario"]), "", result["interpretation"], "",
             "Observation horizon: %s hours after injection." % _number(result["horizon_hours"]), "",
             "Each core includes the exact free peptide and all longer fragments retaining that exact occurrence. "
             "Arrival counts first presence in a draining-node APC, including migration, endogenous production and cross-dressing.", "",
             "| Window | Core [start, end) / allele | Status | Nodal APC arrivals (% SLP dose) | Arrivals per input molecule | Mean nodal pMHC per input | Mean blood core per input |",
             "| --- | --- | --- | ---: | ---: | ---: | ---: |"]
    for row in vaccine_window_rows(result):
        lines.append("| %s | %s [%s, %s) / %s | %s | %s | %s | %s | %s |" % (
            _text(row["candidate_id"]), row["core_sequence"], row["core_start"], row["core_end"], _text(row["allele"]),
            row["status"], _number(row["node_apc_percent_of_dosed_peptide"]), _number(row["node_apc_entries_per_input"]),
            _number(row["node_pmhc_mean_copies_per_input"]), _number(row["blood_core_mean_copies_per_input"])))
    lines += ["", "Means are copy-hours divided by input molecule count and observation hours; they are dimensionless time-averaged amounts. "
              "At a zero-hour horizon they are undefined. RNA arrivals per input transcript can exceed one. "
              "SLP percentages refer to administered molecule count, not residual fragment mass. "
              "Separate cores/alleles are marginal scenarios; do not sum them as a unique-molecule delivery fraction.", "",
              "The blood curve follows core-bearing material after IM/SC input. Its sampled post-peak half-decline bracket is in the CSV/JSON. "
              "An intrinsic circulation elimination half-life is not identified by this model.", ""]
    for window in result["candidates"]:
        dose = window["dose"]
        lines += ["## " + _text(window["id"]), "", "Sequence: `%s`; delivery: %s; status: %s." % (
            window["peptide_input"]["sequence"], window["delivery"], window["status"]), ""]
        if dose:
            lines += ["Dose: %s %s; %s. %s" % (_number(dose["amount"]), dose["unit"], dose["basis"], _text(dose["source"])), ""]
        aggregate = window["aggregate"]
        lines += ["Score-weighted mean core arrivals per input: %s (%s). %s." % (
            _number(aggregate["score_weighted_mean_node_entry_per_input"]), aggregate["status"], aggregate["limitations"]), ""]
        for core in window["cores"]:
            if core["missing_inputs"]:
                lines += ["- %s / %s: missing %s." % (_text(core["label"]), _text(core["allele"]), _text("; ".join(core["missing_inputs"])))]
        lines += [""]
    return "\n".join(lines)
