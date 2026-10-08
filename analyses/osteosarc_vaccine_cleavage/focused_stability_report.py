"""Export an epitope-centered report from frozen, local cascade evidence.

No external prediction service is used. Rendering dependencies are optional;
install mhctools[vaccine-report]. Inputs and outputs may contain private data.
"""

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import html
import io
import json
import math
from pathlib import Path

from mhctools import (
    CleavageResult, DegradationTarget, PeptideInput, __version__, annotate_target_cleavage,
    predict_cleavage, simulate_target_degradation, summarize_target_degradation,
)


MODEL_LABELS = {
    "dpp4-qpisa": "DPP4", "ace-dipeptidyl": "ACE", "cpn-basic": "CPN",
    "app2-xp": "Aminopeptidase P", "fap-endo-gp": "FAP internal",
    "fap-dipeptidyl": "FAP N-terminal", "mme-hydrophobic": "MME",
    "anpep-ala": "ANPEP", "enpep-acidic": "ENPEP",
    "prep-pro": "PREP", "cpb2-basic": "CPB2 (activation unknown)",
}
SCENARIOS = ("Uniform cuts", "Recognition 3x", "Recognition 10x")
PROCESSING = {
    "pepsickle-in-vivo-human-only": "Pepsickle human",
    "netchop-3.1-20s-3.0": "NetChop 20S",
    "netchop-3.1-cterm-3.0": "NetChop Cterm",
    "netcleave-ii-hla": "NetCleave II",
}
SOURCES = (
    ("Stable products after parent cleavage", "https://doi.org/10.1021/acs.jmedchem.1c00795"),
    ("Blood, plasma and serum give different stability profiles", "https://doi.org/10.1371/journal.pone.0178943"),
    ("PeptiVerse: 130 sequence half-life examples", "https://doi.org/10.1038/s41467-026-74167-w"),
    ("Cavaco composition-equation baseline", "https://doi.org/10.1111/cts.12985"),
    ("DPP4 purified-enzyme depletion assay", "https://doi.org/10.1038/s44320-024-00071-4"),
    ("Constrained cyclotides can persist beyond 24 hours", "https://doi.org/10.1021/ja405108p"),
    ("SLP processing linkers can affect presentation after uptake", "https://doi.org/10.1080/2162402X.2018.1560919"),
    ("Semaglutide: albumin binding and a separate DPP4-resistant substitution", "https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/209637s038lbl.pdf"),
)


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(root):
    """Validate both report and upstream source manifests before reading data."""
    root = root.resolve()
    manifest = json.loads((root / "SHA256SUMS.json").read_text())
    manifest = manifest.get("files", manifest)
    for name, expected in manifest.items():
        if isinstance(expected, dict):
            expected = expected["sha256"]
        path = root / name
        if not path.resolve().is_relative_to(root):
            raise ValueError("Manifest path leaves frozen input directory")
        if sha256(path) != expected:
            raise ValueError("Frozen input checksum mismatch: " + name)


def positive_hours(value):
    """Unavailable estimates remain absent; malformed durations are rejected."""
    if value in (None, ""):
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Half-life must be finite and positive")
    return value


def select_primary_targets(targets):
    """Select one best mapped-mutant target per class, never compare classes."""
    groups = defaultdict(list)
    for target in targets:
        if (target.get("mhc_class") in ("I", "II") and
                target.get("mutant_overlap") is True and
                target.get("target_label", "").startswith("Predicted mutant")):
            rank = float(target["percentile_rank"])
            if math.isfinite(rank) and 0 <= rank <= (2 if target["mhc_class"] == "I" else 5):
                groups[target["sequence_record_id"], target["mhc_class"]].append(target)
    selected = [min(rows, key=lambda t: (
        float(t["percentile_rank"]), t["start"], t["end"], t.get("allele", "")))
        for rows in groups.values()]
    return sorted(selected, key=lambda t: (t["gene"], t["sequence_record_id"], t["mhc_class"]))


def is_candidate(row):
    """Only display evidence flags; native scored models keep their semantics."""
    return (row["status"] in ("matched", "reported") or
            row["model"] == "dpp4-qpisa" and row["status"] == "scored" and row["score"] > 0)


def load_bundle(root):
    """Read without modifying frozen outputs and verify their stored manifest."""
    root = root.resolve()
    verify_manifest(root)
    provenance = json.loads((root / "provenance.json").read_text())
    source = Path(provenance["source"])
    if sha256(source / "SHA256SUMS.json") != provenance["source_manifest_sha256"]:
        raise ValueError("Upstream source manifest differs from frozen cascade provenance")
    verify_manifest(source)
    tables = source / "processing-final" / "tables"
    records = {r["sequence_record_id"]: r for r in read_csv(tables / "vaccine_sequence_inventory.csv")
               if r["sequence_type"] == "synthetic_long_peptide"}
    targets = read_csv(root / "targets.csv")
    for t in targets:
        t["start"], t["end"] = int(t["start"]), int(t["end"])
        raw = t["mutant_overlap"].lower()
        t["mutant_overlap"] = True if raw == "true" else False if raw == "false" else None
        seq = records[t["sequence_record_id"]]["sequence"]
        DegradationTarget(t["target_label"], t["start"], t["end"]).validate(seq)
        if seq[t["start"]:t["end"]] != t["target_sequence"]:
            raise ValueError("Target sequence does not match its exact occurrence")
    estimates = {name: {} for name in ("PeptiVerse", "Cavaco")}
    for row in read_csv(root / "fragment_half_lives.csv"):
        estimates[row["estimator"]][row["sequence"]] = positive_hours(row["half_life_hours"])
    native = json.loads((root / "fragment_cleavage_native.json").read_text())
    annotations = {}
    for t in targets:
        seq = records[t["sequence_record_id"]]["sequence"]
        results = [CleavageResult.from_dict(r) for r in native[seq]]
        # Cheap built-in assessments extend annotation coverage only. They do
        # not retrospectively change the frozen simulation's cut allocations.
        results.extend(predict_cleavage(seq, models=("prep-pro", "cpb2-basic")))
        annotations[t["sequence_record_id"], t["target_label"]] = annotate_target_cleavage(
            DegradationTarget(t["target_label"], t["start"], t["end"]), results)
    curves = defaultdict(list)
    for row in read_csv(root / "curves.csv"):
        key = row["sequence_record_id"], row["target_label"], row["estimator"], row["scenario"]
        curves[key].append({k: float(row[k]) for k in (
            "time_hours", "parent_remaining", "target_in_circulation", "unknown")})
    summaries = {(r["sequence_record_id"], r["target_label"], r["estimator"], r["scenario"]): r
                 for r in read_csv(root / "summary.csv")}
    coverage = {r["sequence_record_id"]: r for r in read_csv(root / "coverage.csv")}
    quantitative = read_csv(tables / "slp_quantitative_bond_scores.csv")
    catalog = {r["model"]: r for r in read_csv(tables / "model_catalog.csv")}
    paths = read_csv(root / "sampled_paths.csv")
    primary = select_primary_targets(targets)
    uniform_summaries = {}
    for t in primary:
        sequence = records[t["sequence_record_id"]]["sequence"]
        target = DegradationTarget(t["target_label"], t["start"], t["end"])
        for estimator in estimates:
            parts = (t["sequence_record_id"], t["start"], t["end"], estimator, "Uniform cuts")
            seed = int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")
            sampled = simulate_target_degradation(
                PeptideInput(sequence, occurrence_id=t["sequence_record_id"]), target, estimates[estimator],
                lambda s: [1] * (len(s) - 1), estimator=estimator, scenario="Uniform cuts",
                n_paths=provenance["sample_size_per_target_estimator_scenario"],
                horizon_hours=provenance["horizon_hours"], seed=seed)
            uniform_summaries[t["sequence_record_id"], t["target_label"], estimator] = (
                summarize_target_degradation(sampled))
    return dict(root=root, source=source, records=records, targets=targets,
                primary=primary, estimates=estimates, uniform_summaries=uniform_summaries,
                annotations=annotations, curves=curves, summaries=summaries,
                coverage=coverage, quantitative=quantitative, catalog=catalog,
                provenance=provenance, paths=paths)


def fmt_hours(value):
    if value in (None, ""):
        return "Unavailable"
    value = float(value)
    return "%d min" % round(value * 60) if value < 1 else "%.1f h" % value


def fmt_retention(summary):
    if summary["median_status"] == "beyond_horizon":
        return "> " + fmt_hours(summary["retention_median_lower_bound_hours"])
    return fmt_hours(summary["retention_median_hours"])


def target_data(bundle, target):
    key = target["sequence_record_id"], target["target_label"]
    sequence = bundle["records"][key[0]]["sequence"]
    coverage = bundle["coverage"][key[0]]
    a, b = coverage["mutant_region_start"], coverage["mutant_region_end"]
    mutant = [int(float(a)), int(float(b))] if a and b else None
    estimates, curves, medians = {}, {}, {}
    for model in ("PeptiVerse", "Cavaco"):
        estimates[model] = dict(parent=bundle["estimates"][model].get(sequence),
                                released_target=bundle["estimates"][model].get(target["target_sequence"]))
        curves[model], medians[model] = {}, {}
        for scenario in SCENARIOS:
            curves[model][scenario] = bundle["curves"][key + (model, scenario)]
            summary = bundle["summaries"][key + (model, scenario)]
            unknown = float(summary["unknown_fraction"])
            medians[model][scenario] = (positive_hours(summary["target_retention_median_hours"])
                                       if unknown == 0 else None)
        medians[model]["Uniform cuts"] = bundle["uniform_summaries"][key + (model,)]["retention_median_hours"]
    scores = []
    for row in bundle["quantitative"]:
        if row["sequence_record_id"] != key[0] or row["assessable"].lower() != "true":
            continue
        bond = int(float(row["bond"]))
        if (target["start"] < bond < target["end"] and
                (row["model"] in PROCESSING or "-pwm-" in row["model"])):
            scores.append(dict(model=row["model"], bond=bond, score=float(row["score"]),
                               units=row["score_units"], context=row["biological_context"]))
    paths = [r for r in bundle["paths"] if r["sequence_record_id"] == key[0] and
             r["target_label"] == key[1] and r["sample_index"] == "0" and
             r["scenario"] in SCENARIOS]
    # Recreate rather than mutate the frozen table dictionaries.
    paths = [dict(r, cut_mechanism_evidence=r["cut_mechanism_evidence"].replace(
        "background / unknown mechanism", "Assumed bond choice; no supporting flag in sampling models")) for r in paths]
    return dict(target=target, sequence=sequence, mutant=mutant, estimates=estimates,
                annotations=bundle["annotations"][key], curves=curves,
                conditional_medians=medians, internal_native_scores=scores,
                uniform_summaries={m: bundle["uniform_summaries"][key + (m,)]
                                   for m in ("PeptiVerse", "Cavaco")},
                calibrated_epitope_half_life_hours=None, sampled_paths=paths)


def sequence_figure(item):
    """One residue-centered map per model; every red line lies BETWEEN letters."""
    import matplotlib.pyplot as plt
    rows = item["annotations"]
    candidates = [r for r in rows if is_candidate(r)]
    models = [m for m in MODEL_LABELS if any(r["model"] == m for r in candidates)]
    if not models:
        models = ["No flag in these models"]
    sequence, target = item["sequence"], item["target"]
    fig, ax = plt.subplots(figsize=(7.3, max(1.35, .27 * len(models) + .75)))
    ax.axvspan(target["start"], target["end"], color="#ffe7a0")
    for i, aa in enumerate(sequence):
        mutant = item["mutant"] and item["mutant"][0] <= i < item["mutant"][1]
        ax.text(i + .5, len(models) + .45, aa, ha="center", va="center", fontsize=13,
                fontweight="bold", color="#733492" if mutant else "#192d3e")
        ax.text(i + .5, len(models) + .05, str(i + 1), ha="center", fontsize=7, color="#657786")
    for y, model in enumerate(models):
        for row in candidates:
            if row["model"] == model:
                ax.plot([row["bond"]] * 2, [y - .28, y + .28], color="#ce2333", ls="--", lw=2)
    ax.set_yticks(range(len(models)), [MODEL_LABELS.get(m, m) for m in models])
    ax.tick_params(axis="y", labelsize=9, length=0)
    ax.set_xlim(0, len(sequence)); ax.set_ylim(-.5, len(models) + .85)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=.15)
    stream = io.BytesIO(); fig.savefig(stream, format="png", dpi=180); plt.close(fig)
    stream.seek(0)
    return stream


def curve_figure(item):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.3, 1.85))
    for model, color in (("PeptiVerse", "#007e87"), ("Cavaco", "#b7600a")):
        rows = item["curves"][model]["Uniform cuts"]
        ax.plot([r["time_hours"] for r in rows], [100 * r["target_in_circulation"] for r in rows],
                color=color, lw=2, label=model + " target retained")
        ax.plot([r["time_hours"] for r in rows], [100 * r["parent_remaining"] for r in rows],
                color=color, lw=1, ls=":", alpha=.6, label=model + " parent intact")
    ax.set_xlim(0, 8); ax.set_ylim(0, 102)
    ax.set_xlabel("Hours since peptide was present in the modeled serum compartment")
    ax.set_ylabel("Remaining (%)"); ax.grid(axis="y", alpha=.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, ncol=2, frameon=False, loc="upper right")
    fig.tight_layout(pad=.4)
    stream = io.BytesIO(); fig.savefig(stream, format="png", dpi=180); plt.close(fig)
    stream.seek(0)
    return stream


def model_table_rows(item):
    rows = []
    for model, label in MODEL_LABELS.items():
        assessed = [r for r in item["annotations"] if r["model"] == model]
        matched = [r for r in assessed if is_candidate(r)]
        unavailable = next((r["reason"] for r in assessed if r["status"] == "unavailable"), None)
        if unavailable:
            rows.append([label, "Unavailable", unavailable]); continue
        if not matched:
            continue
        inside = [r for r in matched if r["target_effect"] == "target_split"]
        outside = [r for r in matched if r["target_effect"] != "target_split"]
        detail = ", ".join(r["bond_label"] for r in inside) or "No flag"
        if model != "dpp4-qpisa" and matched:
            label += " (" + assessed[0]["motif_strictness"] + " motif)"
        other = ", ".join(r["bond_label"] + (" (boundary)" if r["target_effect"] == "boundary_release" else "")
                          for r in outside) or "No flag"
        if model == "dpp4-qpisa" and matched:
            signal = 100 * (1 - 2 ** (-matched[0]["score"]))
            label += " (%.0f%% assay loss)" % signal
        rows.append([label, detail, other])
    return rows


def build_pdf(bundle, output, title):
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("Small", parent=styles["BodyText"], fontSize=8, leading=10.5))
    styles.add(ParagraphStyle("Caption", parent=styles["BodyText"], fontSize=7.5, leading=9))
    styles["BodyText"].fontSize = 9; styles["BodyText"].leading = 12
    styles["Heading1"].fontSize = 18; styles["Heading1"].leading = 22
    story = []
    def para(text, style="BodyText"):
        return Paragraph(text, styles[style])
    def table(data, widths):
        data = [[para(html.escape(str(v)), "Small") for v in row] for row in data]
        obj = Table(data, colWidths=widths, hAlign="LEFT", repeatRows=1)
        obj.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eff4")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"), ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, -1), .3, colors.HexColor("#dce4ea"))]))
        return obj
    primary, all_targets = bundle["primary"], bundle["targets"]
    story += [para(html.escape(title), "Heading1"), Spacer(1, 12),
              para("<b>Primary question: how long does the selected mutant epitope remain intact, including in shorter fragments?</b>"),
              Spacer(1, 12), para("<b>Calibrated epitope-loss half-life: not established.</b> Parent and released-target sequence estimates are available. Successive-cut trajectories remain illustrative because enzyme-specific serum rates and cut probabilities are unavailable."),
              Spacer(1, 12), para("%d top predicted mutant targets across %d of %d vaccine constructs. Class I and class II are selected separately. %d source/predicted annotations remain in the audit; disclosed targets appear as references." % (
                  len(primary), len({t["sequence_record_id"] for t in primary}), len(bundle["records"]), len(all_targets))),
              Spacer(1, 16), para("How to read a target page", "Heading2"),
              para("<b>Gold:</b> exact tracked target. <b>Purple letters:</b> source-mapped mutation. <b>Dashed red lines:</b> candidate extracellular bonds, separated by named model. An internal cut splits the target; flank or boundary cuts can retain it. A flag is not a cleavage probability."),
              Spacer(1, 10), para("<b>Half-life table:</b> parent disappearance and a hypothetical already-released free target are different sequence endpoints. Each estimator is kept separate. The conditional time-to-50%-loss uses equal weights at all bonds; it is an assumption, not a biological calibration."),
              Spacer(1, 10), para("<b>Coverage:</b> no flag does not mean protection. Missing coefficients, model length limits and unknown activation remain visible. Enzyme access, concentration, binding, chemistry and formulation are not established for these constructs."),
              Spacer(1, 10), para("<b>Supporting detail:</b> arbitrary 3x/10x weighting, all native scores, processing after uptake, source targets and assumed unflagged cuts are retained in the interactive report and audit files. They do not provide measured circulation half-lives."), PageBreak()]
    overview = [["Selected target", "PV: parent / free target", "Cavaco: parent / free target"]]
    for t in primary:
        item = target_data(bundle, t)
        overview.append([t["gene"] + " / " + t["mhc_class"] + " / " + t["target_sequence"],
                         " / ".join(fmt_hours(item["estimates"]["PeptiVerse"][k]) for k in ("parent", "released_target")),
                         " / ".join(fmt_hours(item["estimates"]["Cavaco"][k]) for k in ("parent", "released_target"))])
    story += [para("Selected targets: sequence estimates", "Heading1"),
              para("Free-target values assume that sequence is already present with free termini. They are not total target survival times after injection."),
              Spacer(1, 12), table(overview, [224, 152, 152]), PageBreak()]
    for t in primary:
        item = target_data(bundle, t)
        story += [para(html.escape(t["gene"] + " / mutant class " + t["mhc_class"]), "Heading1"),
                  para(html.escape(t["sequence_record_id"].split(":")[-1] + " | " + t["allele"] +
                                   " | rank " + format(float(t["percentile_rank"]), ".3g") + "% | " + t["kind"]), "Small"),
                  para("<b>" + html.escape(t["target_sequence"]) + "</b> | residues %d-%d | mutation included" % (t["start"] + 1, t["end"])),
                  Spacer(1, 6)]
        data = [["Separate estimator", "Parent half-life", "Free-target half-life", "Conditional 50% target loss"]]
        for model in ("PeptiVerse", "Cavaco"):
            data.append([model, fmt_hours(item["estimates"][model]["parent"]),
                         fmt_hours(item["estimates"][model]["released_target"]),
                         fmt_retention(item["uniform_summaries"][model])])
        story += [table(data, [100, 120, 140, 168]),
                  para("Last column: illustrative equal-weight cuts, zero clearance/uptake. Neither estimator is validated on these vaccine fragments. Class-II mutant flanks, when required, are included in the tracked span.", "Small"),
                  Spacer(1, 6)]
        stream = sequence_figure(item)
        image = Image(stream, width=528, height=528 * Image(stream).imageHeight / Image(stream).imageWidth)
        story += [image, para("Candidate bonds if the enzyme has access; labels identify evidence, not activity or serum rates.", "Small"),
                  table([["Model", "Would split target", "Would retain target / coverage"]] + model_table_rows(item), [100, 195, 233]),
                  Spacer(1, 4), Image(curve_figure(item), width=528, height=119),
                  para("Illustrative equal-weight cuts, including unflagged bonds. Solid: target in any retained fragment; dotted: parent intact. No clearance, formulation release or uptake is inferred. No flag does not imply protection; full model coverage is in the interactive report.", "Caption"), PageBreak()]
    refs = [t for t in all_targets if not t["target_label"].startswith("Predicted mutant")]
    story += [para("Disclosed target references", "Heading1"),
              para("Source annotations are retained as references. They are not automatically the best mapped-mutant predictions or experimentally confirmed presentation."), Spacer(1, 10),
              table([["Construct", "Reference", "Target / positions"]] + [
                  [t["gene"] + " / " + t["sequence_record_id"].split(":")[-1], t["target_label"],
                   t["target_sequence"] + " / %d-%d" % (t["start"] + 1, t["end"])] for t in refs], [165, 138, 225]), PageBreak()]
    processing_rows = [["Selected target", "After-uptake model", "Maximum internal native score"]]
    for t in primary:
        item = target_data(bundle, t)
        relevant_models = ("netcleave-ii-hla",) if t["mhc_class"] == "II" else tuple(
            m for m in PROCESSING if m != "netcleave-ii-hla")
        for model in relevant_models:
            scores = [r for r in item["internal_native_scores"] if r["model"] == model]
            if scores:
                peak = max(scores, key=lambda r: r["score"])
                b, seq = peak["bond"], item["sequence"]
                signal = "%s%d|%s%d: %.3f (%s)" % (seq[b - 1], b, seq[b], b + 1,
                                                       peak["score"], peak["units"])
            else:
                signal = "No assessable internal score; no protection inferred"
            processing_rows.append([t["gene"] + " / class " + t["mhc_class"] + " / " +
                                    t["sequence_record_id"].split(":")[-1], PROCESSING[model], signal])
    story += [para("Target-internal processing scores after uptake", "Heading1"),
              para("These scores do not set the extracellular clock. A maximum inside a target is supporting evidence, not a calibrated cleavage probability or proof of lost presentation. Pepsickle human is the primary human processing model; NetChop 20S and Cterm have different training endpoints. NetCleave-II is a class-II proxy without an enzyme assignment."),
              Spacer(1, 10), table(processing_rows, [160, 130, 238]), PageBreak()]
    selected = {t["sequence_record_id"] for t in primary}
    story += [para("Construct coverage", "Heading1"),
              para("Missing target or mutation mapping is missing information, not protection."), Spacer(1, 10),
              table([["Construct", "Primary mutant target", "Mutation mapping / limitation"]] + [
                  [r["gene"] + " / " + key.split(":")[-1], "Available" if key in selected else "Not resolved",
                   bundle["coverage"][key]["mutation_mapping"]] for key, r in bundle["records"].items()], [155, 110, 263]), PageBreak()]
    story += [para("Evidence, assumptions and next measurements", "Heading1"),
              para("<b>Measure the target-bearing products.</b> An LC-MS time course should quantify intact parent and products retaining the selected mutant target, with peptide chemistry and serum/plasma preparation recorded. Site predictions alone do not establish rates."), Spacer(1, 10),
              para("<b>Keep compartments separate.</b> Extracellular cleavage can affect delivery before uptake. Proteasomal, lysosomal and MHC-ligand processing scores describe later mechanisms; they are not extracellular rate inputs. NetChop Cterm and 20S have different training endpoints. Human Pepsickle is the primary human processing model."), Spacer(1, 10),
              para("<b>Long lifetime is possible.</b> Constrained structures, binding and chemical modifications can produce lifetimes over 24 hours. The 24-hour simulation horizon is not a biological cap. A missing native duration or unknown path must not become a zero or a protected target."), Spacer(1, 10),
              para("<b>Frozen simulation assumptions.</b> Each fragment has an exponential next-cut clock derived from one sequence half-life estimate. Equal-weight cuts allocate this total hazard across all bonds. The optional 3x/10x scenarios favor flagged bonds by arbitrary multipliers; overlapping flags do not stack. An unflagged sampled cut is an assumed cut without supporting enzyme prediction."), Spacer(1, 10),
              para("<b>Endpoint limitations.</b> A cut strictly inside the exact selected target destroys that sequence. That does not establish loss of every alternative epitope, T-cell recognition or productive MHC presentation. Serum estimates are not patient circulation PK; both half-life estimators have unestablished vaccine-fragment accuracy."), Spacer(1, 14)]
    for label, url in SOURCES:
        story += [para('<a href="%s">%s</a>' % (url, html.escape(label)), "Small"), Spacer(1, 6)]
    def footer(canvas, doc):
        canvas.setFont("Helvetica", 7); canvas.setFillColor(colors.HexColor("#657786"))
        canvas.drawString(42, 25, "Local sequence estimates | Conditional simulations | No calibrated circulation half-life")
        canvas.drawRightString(570, 25, str(doc.page))
    SimpleDocTemplate(str(output), pagesize=(612, 792), rightMargin=42, leftMargin=42,
                      topMargin=38, bottomMargin=42, title=title).build(story, onFirstPage=footer, onLaterPages=footer)


def export(root, output, title, overwrite=False):
    root, output = root.resolve(), output.resolve()
    if output == root or output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError("Output must be separate from frozen inputs")
    if output.exists() and not overwrite:
        raise ValueError("Output already exists; choose a new directory")
    bundle = load_bundle(root)
    output.mkdir(parents=True, exist_ok=overwrite)
    data = dict(title=title, generated=datetime.now(timezone.utc).isoformat(), mhctools_version=__version__,
                source=str(root), input_manifest_sha256=sha256(root / "SHA256SUMS.json"),
                primary=[target_data(bundle, t) for t in bundle["primary"]],
                reference_targets=[t for t in bundle["targets"] if not t["target_label"].startswith("Predicted mutant")],
                coverage=list(bundle["coverage"].values()), sources=SOURCES,
                calibrated_epitope_half_life_hours=None)
    (output / "focused_evidence.json").write_text(json.dumps(data, indent=2, allow_nan=False))
    flat = [dict(sequence_record_id=key[0], target_label=key[1], **row)
            for key, rows in bundle["annotations"].items() for row in rows]
    with (output / "target_cut_annotations.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]) if flat else [])
        writer.writeheader(); writer.writerows(flat)
    build_pdf(bundle, output / "epitope-stability-report.pdf", title)
    template = Path(__file__).with_name("focused_stability_template.html").read_text()
    payload = json.dumps(data, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    (output / "epitope-stability.html").write_text(template.replace("__DATA__", payload))
    (output / "provenance.json").write_text(json.dumps(dict(
        **{k: data[k] for k in ("generated", "mhctools_version", "source", "input_manifest_sha256")},
        exporter_sha256=sha256(Path(__file__)),
        template_sha256=sha256(Path(__file__).with_name("focused_stability_template.html")),
        new_learned_predictor_inference=False, extra_annotation_rules=["prep-pro", "cpb2-basic"],
        frozen_simulation_unchanged=True, uniform_summaries_recomputed=True,
        sampler_sha256=sha256(Path(simulate_target_degradation.__code__.co_filename)),
        primary_target_count=len(data["primary"]),
        primary_construct_count=len({t["target"]["sequence_record_id"] for t in data["primary"]})), indent=2))
    print(output / "epitope-stability-report.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="Selected mutant epitope stability")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    export(args.input, args.output, args.title, args.overwrite)


if __name__ == "__main__":
    main()
