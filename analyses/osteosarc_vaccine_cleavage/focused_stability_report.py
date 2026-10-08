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
    "prep-pro": "PREP", "cpb2-basic": "CPB2",
}
SCENARIOS = ("Uniform cuts", "Recognition 3x", "Recognition 10x")
PROCESSING = {
    "pepsickle-in-vivo-human-only": "Pepsickle human",
    "pepsickle-in-vivo-all-mammal": "Pepsickle all-mammal (comparison)",
    "netchop-3.1-20s-3.0": "NetChop 20S",
    "netchop-3.1-cterm-3.0": "NetChop Cterm",
    "netcleave-i-hla": "NetCleave I",
    "netcleave-ii-hla": "NetCleave II",
    **{"itcell-cat%s-%s" % (enzyme, time): "Cathepsin %s / %s-min matrix" % (enzyme.upper(), time)
       for enzyme in ("b", "h", "s") for time in (15, 60, 240)},
    "eramer-step": "ERAP1 / ERAMER",
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
    ("ACE versus CPN contribution changes with substrate concentration", "https://doi.org/10.1152/ajpheart.2000.278.4.H1069"),
    ("DPP4 cleavage of particular hormones in human blood specimens", "https://doi.org/10.1371/journal.pone.0134427"),
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
    model_summary = read_csv(tables / "slp_model_summary.csv")
    processing_motifs = [r for r in read_csv(tables / "slp_motif_assessments.csv")
                         if set(r["compartments"].split(";")) & {"cytosol", "er", "endosome"}]
    mmp_paths = [name for name in json.loads((source / "SHA256SUMS.json").read_text())
                 if name.endswith("mmp_window_scores.csv")]
    if len(mmp_paths) > 1:
        raise ValueError("Ambiguous frozen MMP evidence")
    mmp = read_csv(source / mmp_paths[0]) if mmp_paths else []
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
                model_summary=model_summary, processing_motifs=processing_motifs, mmp=mmp,
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


def tracked_description(target):
    if target["mhc_class"] == "I":
        return "Predicted minimal class-I epitope"
    if target["binding_core"] == target["target_sequence"]:
        return "Predicted minimal class-II binding core"
    return "Class-II binding core plus mutant flank"


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
    processing = []
    processing_sites = []
    for row in bundle["quantitative"]:
        if (row["sequence_record_id"] == key[0] and row["model"] in PROCESSING and
                row["assessable"].lower() == "true"):
            bond = int(float(row["bond"]))
            processing_sites.append(dict(model=row["model"], bond=bond,
                                         score=float(row["score"]), units=row["score_units"],
                                         display_flag=row["above_display_threshold"].lower() == "true"))
    for model, label in PROCESSING.items():
        native_scores = [r for r in scores if r["model"] == model]
        peak = max(native_scores, key=lambda r: r["score"]) if native_scores else None
        full_scores = [r for r in processing_sites if r["model"] == model]
        full_peak = max(full_scores, key=lambda r: r["score"]) if full_scores else None
        summary = next((r for r in bundle["model_summary"]
                        if r["sequence_record_id"] == key[0] and r["model"] == model), {})
        processing.append(dict(model=model, label=label, peak=peak, full_peak=full_peak,
                               reason=summary.get("unsupported_reason") or
                               ("" if peak else "No assessable bond inside tracked sequence")))
    motifs = []
    for row in bundle["processing_motifs"]:
        if row["sequence_record_id"] != key[0]:
            continue
        row = dict(row)
        raw = row["bond"]
        bond = int(float(raw)) if raw and math.isfinite(float(raw)) else None
        row["bond"] = bond
        row["target_effect"] = (None if bond is None else
                                "target_split" if target["start"] < bond < target["end"] else
                                "boundary_release" if bond in (target["start"], target["end"]) else "flank_trim")
        motifs.append(row)
    return dict(target=target, sequence=sequence, mutant=mutant, estimates=estimates,
                protein_change=bundle["records"][key[0]]["protein_change"],
                tracked_description=tracked_description(target),
                annotations=bundle["annotations"][key], curves=curves,
                conditional_medians=medians, internal_native_scores=scores,
                processing_predictions=processing, processing_sites=processing_sites,
                processing_motifs=motifs,
                mmp_windows=[r for r in bundle["mmp"] if r["sequence_record_id"] == key[0]],
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
    for y, model in enumerate(models):
        for row in candidates:
            if row["model"] == model:
                ax.plot([row["bond"]] * 2, [y - .16, y + .16], color="#ce2333", lw=2)
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
                color=color, lw=2, label=model + " epitope sequence present")
        ax.plot([r["time_hours"] for r in rows], [100 * r["parent_remaining"] for r in rows],
                color=color, lw=1, ls=":", alpha=.6, label=model + " parent intact")
    ax.set_xlim(0, 8); ax.set_ylim(0, 102)
    ax.set_xlabel("Hours after starting with full vaccine peptides")
    ax.set_ylabel("Copies remaining (%)"); ax.grid(axis="y", alpha=.2)
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
        if not matched:
            continue
        inside = [r for r in matched if r["target_effect"] == "target_split"]
        outside = [r for r in matched if r["target_effect"] != "target_split"]
        detail = ", ".join(r["bond_label"] for r in inside) or "No flag"
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
        data = [[v if isinstance(v, Paragraph) else para(html.escape(str(v)), "Small")
                 for v in row] for row in data]
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
              para("<b>Gold:</b> tracked epitope sequence. <b>Purple letters:</b> mutation. <b>Short red marks:</b> candidate extracellular cuts, separated by enzyme. A cut inside gold splits the tracked sequence; outside or at the boundary can leave it intact."),
              Spacer(1, 10), para("<b>What is tracked?</b> Class I: the predicted minimal epitope. Class II: the predicted binding core, extended only when a flank is needed to keep the mutation. The exact tracked sequence and class-II core appear on each page."),
              Spacer(1, 10), para("<b>Three different times:</b> parent half-life estimates disappearance of the full vaccine peptide. Released-sequence half-life estimates disappearance of the epitope/span already alone with free termini. Illustrative time to 50% epitope loss starts with the full vaccine peptide and follows cuts until the tracked sequence is split."),
              Spacer(1, 10), para("<b>Curves:</b> imagine 100 starting vaccine-peptide copies. Solid lines count copies still carrying the entire tracked sequence, including in shorter fragments. Dotted lines count full vaccine peptides still intact. Colors are separate half-life estimators. Cut locations are assumed; the curves are not measured circulation survival."),
              Spacer(1, 10), para("<b>Coverage:</b> no flag does not mean protection. Missing coefficients, model length limits and unknown activation remain visible. Enzyme access, concentration, binding, chemistry and formulation are not established for these constructs."),
              Spacer(1, 10), para("<b>Antigen processing:</b> all frozen processing models, including cathepsins, ERAP1 and intracellular rules, are available in a visible interactive section and the PDF processing appendix. They describe events after uptake and do not set the circulation clock."), PageBreak()]
    guide = json.loads(Path(__file__).with_name("focused_enzyme_guide.json").read_text())
    for section in ("Before uptake", "After uptake", "Inflammation and tissue exposure"):
        story += [para("Enzyme guide: " + section.lower(), "Heading1"),
                  para("Importance depends on the substrate and actual enzyme exposure. This guide gives established roles, not a ranked activity score for the vaccine."), Spacer(1, 10),
                  table([["Enzyme", "Where it acts", "What it cuts", "Circulation / epitope relevance"]] +
                        [[para('<a href="%s">%s</a>' % (r["references"][0], html.escape(r["name"])), "Small")] +
                         [r[k] for k in ("location", "action", "relevance")]
                         for r in guide if r["section"] == section], [83, 108, 122, 215]), Spacer(1, 12)]
        if section == "Inflammation and tissue exposure":
            story += [para("Can we rank DPP4, FAP and ACE?", "Heading2"),
                      para("Not from these scores. Human-plasma experiments found ACE dominated dilute bradykinin degradation while CPN dominated at higher concentrations. DPP4-mediated loss is established for particular hormones. Those substrate-specific findings cannot supply relative vaccine cut rates."), Spacer(1, 8),
                      para("Measure parent and epitope-bearing products with selective inhibitor conditions, then fit fragment-specific rates. Purified-enzyme catalytic efficiency and active concentration can support a kinetic model, but membrane access, inhibitors and competition must also be represented. Reassess each new fragment; a global enzyme multiplier is not justified.")]
        story.append(PageBreak())
    overview = [["Selected target", "PV: parent / free target", "Cavaco: parent / free target"]]
    for t in primary:
        item = target_data(bundle, t)
        overview.append([t["gene"] + " " + item["protein_change"] + " / " + t["mhc_class"] + " / " + t["target_sequence"],
                         " / ".join(fmt_hours(item["estimates"]["PeptiVerse"][k]) for k in ("parent", "released_target")),
                         " / ".join(fmt_hours(item["estimates"]["Cavaco"][k]) for k in ("parent", "released_target"))])
    story += [para("Selected targets: sequence estimates", "Heading1"),
              para("Free-target values assume that sequence is already present with free termini. They are not total target survival times after injection."),
              Spacer(1, 12), table(overview, [224, 152, 152]), PageBreak()]
    for t in primary:
        item = target_data(bundle, t)
        story += [para(html.escape(t["gene"] + " " + item["protein_change"] + " / class " + t["mhc_class"]), "Heading1"),
                  para(html.escape(t["sequence_record_id"].split(":")[-1] + " | " + t["allele"] +
                                   " | rank " + format(float(t["percentile_rank"]), ".3g") + "% | " + t["kind"]), "Small"),
                  para("<b>" + html.escape(t["target_sequence"]) + "</b> | " + item["tracked_description"]),
                  Spacer(1, 6)]
        if t["mhc_class"] == "II" and t["binding_core"] != t["target_sequence"]:
            story += [para("Binding core: " + html.escape(t["binding_core"]) + "; gold also includes the mapped mutant flank.", "Small")]
        data = [["Estimator", "Full vaccine peptide half-life", "Released sequence alone half-life", "Illustrative time to 50% epitope loss"]]
        for model in ("PeptiVerse", "Cavaco"):
            data.append([model, fmt_hours(item["estimates"][model]["parent"]),
                         fmt_hours(item["estimates"][model]["released_target"]),
                         fmt_retention(item["uniform_summaries"][model])])
        story += [table(data, [100, 120, 140, 168]),
                  para("Released sequence: starts alone with free termini. Last column: starts in the full vaccine peptide and counts survival through flank cuts; equal-weight cuts are an assumption.", "Small"),
                  Spacer(1, 6)]
        stream = sequence_figure(item)
        image = Image(stream, width=528, height=528 * Image(stream).imageHeight / Image(stream).imageWidth)
        story += [image, para("Red marks: possible cuts if the enzyme has access. Gold: entire tracked sequence. No mark does not establish protection.", "Small"),
                  table([["Enzyme", "Would split epitope sequence", "Would leave sequence intact"]] + model_table_rows(item), [100, 195, 233]),
                  para("Out of 100 starting vaccine copies", "Heading3"),
                  para("Solid: copies still carrying <b>" + html.escape(t["target_sequence"]) + "</b>, including trimmed fragments. Dotted: full vaccine peptide still intact.", "Small"),
                  Spacer(1, 4), Image(curve_figure(item), width=528, height=119),
                  para("Equal weights at all bonds; zero clearance/uptake. Illustrative, not measured circulation survival. Processing: appendix.", "Caption"), PageBreak()]
    refs = [t for t in all_targets if not t["target_label"].startswith("Predicted mutant")]
    story += [para("Disclosed target references", "Heading1"),
              para("Source annotations are retained as references. They are not automatically the best mapped-mutant predictions or experimentally confirmed presentation."), Spacer(1, 10),
              table([["Construct", "Reference", "Target / positions"]] + [
                  [t["gene"] + " / " + t["sequence_record_id"].split(":")[-1], t["target_label"],
                   t["target_sequence"] + " / %d-%d" % (t["start"] + 1, t["end"])] for t in refs], [165, 138, 225]), PageBreak()]
    processing_rows = [["Selected target", "After-uptake model", "Maximum internal native score"]]
    for t in primary:
        item = target_data(bundle, t)
        for prediction in item["processing_predictions"]:
            model, peak = prediction["model"], prediction["peak"]
            if peak:
                b, seq = peak["bond"], item["sequence"]
                units = peak["units"].replace("native ", "").replace("sum of log2 profile/background ratios", "log2 profile/background sum")
                signal = "%s%d|%s%d: %.3f (%s)" % (seq[b - 1], b, seq[b], b + 1,
                                                       peak["score"], units)
            else:
                signal = prediction["reason"]
            processing_rows.append([t["gene"] + " " + item["protein_change"] + " / " + t["mhc_class"] +
                                    " / construct " + t["sequence_record_id"].rsplit("-", 1)[-1], PROCESSING[model], signal])
    story += [para("Target-internal processing scores after uptake", "Heading1"),
              para("All frozen quantitative processing models are restored here: Pepsickle human and all-mammal, NetChop 20S/Cterm, NetCleave I/II, cathepsins B/H/S at each assay matrix time, and ERAP1/ERAMER. Each maximum is taken strictly inside the tracked sequence. N-terminal trimming models may have only a flank/boundary site, so absence of an internal score is not protection. Complete bond scores and qualitative intracellular rules are in the interactive report."),
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
                enzyme_guide=json.loads(Path(__file__).with_name("focused_enzyme_guide.json").read_text()),
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
        enzyme_guide_sha256=sha256(Path(__file__).with_name("focused_enzyme_guide.json")),
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
