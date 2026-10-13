"""Paired, explicitly synthetic SLP/cytosolic-RNA scenarios for disclosed SLPs.

Run from the repository root with ``PYTHONPATH=examples python
analyses/osteosarc_vaccine_cleavage/paired_exposure.py --output-dir NEW_DIR``.
The input ligand ranks select target occurrences; they are never converted
into pharmacokinetic rates or interpreted as observed epitopes.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics

from mhctools import assess_vaccine_windows, vaccine_window_rows
from vaccine_sensitivity_scenarios import synthetic_kinetics


FROZEN = Path("analyses/osteosarc_vaccine_cleavage/results/2026-10-06T142010-273003-0400/tables")
LOSSES = (0.1, 1.0, 10.0)
TIMES = (0, 0.5, 1, 2, 4, 8, 12, 24, 48)


def _records(inventory_path, ligand_path):
    with inventory_path.open(newline="") as stream:
        inventory = [row for row in csv.DictReader(stream)
                     if row["sequence_type"] == "synthetic_long_peptide"]
    with ligand_path.open(newline="") as stream:
        ligands = [row for row in csv.DictReader(stream) if row["mhc_class"] == "I"]
    by_id = {}
    for ligand in ligands:
        by_id.setdefault(ligand["sequence_record_id"], []).append(ligand)
    if len(inventory) != 40 or len({row["sequence_record_id"] for row in inventory}) != 40:
        raise ValueError("Frozen input must contain 40 distinct disclosed SLP records")
    selected = []
    for row in inventory:
        choices = by_id.get(row["sequence_record_id"], [])
        if not choices:
            raise ValueError("No class-I ligand for " + row["sequence_record_id"])
        minimal = row["minimal_epitope"]
        exact_minimal = [x for x in choices if minimal and x["peptide"] == minimal]
        pool = exact_minimal or choices
        ligand = min(pool, key=lambda x: (float(x["percentile_rank"]),
                                             -float(x["score"]), x["allele"],
                                             int(x["start"]), int(x["end"])))
        # Frozen ligand table uses one-based inclusive coordinates; the model
        # uses zero-based half-open coordinates. Check the actual residues.
        start, end = int(ligand["start"]) - 1, int(ligand["end"])
        if row["sequence"][start:end] != ligand["peptide"]:
            raise ValueError("Ligand coordinates do not match " + row["sequence_record_id"])
        selected.append((row, ligand, start, end, bool(exact_minimal)))
    return selected


def _candidate(row, ligand, start, end, disclosed_core, delivery, loss):
    peptide_id = row["sequence_record_id"]
    # Translation initiates with methionine. Keep the published SLP unchanged,
    # and model a retained initiator Met in the counterfactual RNA product.
    initiator_added = delivery == "cytosolic_mrna" and not row["sequence"].startswith("M")
    peptide = ("M" if initiator_added else "") + row["sequence"]
    offset = int(initiator_added)
    if peptide[start + offset:end + offset] != ligand["peptide"]:
        raise ValueError("Target core does not match translated antigen for " + peptide_id)
    candidate = {
        "id": peptide_id + ":" + delivery,
        "peptide": peptide,
        "delivery": delivery,
        "cores": [{
            "label": "disclosed minimal epitope" if disclosed_core else "top frozen MHCflurry class-I ligand",
            "start": start + offset, "end": end + offset, "sequence": ligand["peptide"],
            "allele": ligand["allele"], "mhc_class": "I",
        }],
        "metadata": {
            "sequence_record_id": peptide_id, "gene": row["gene"],
            "vaccines": row["vaccines"], "source_url": row["source_url"],
            "sequence_sha256": row["sequence_sha256"],
            "disclosed_slp_sequence": row["sequence"],
            "initiator_methionine_added": initiator_added,
            "ligand_percentile_rank": float(ligand["percentile_rank"]),
            "ligand_score": float(ligand["score"]),
            "disclosed_minimal_epitope": row["minimal_epitope"] or None,
            "target_basis": "disclosed_minimal_epitope" if disclosed_core else "predicted_ligand",
        },
    }
    if delivery == "cytosolic_mrna":
        candidate["dose"] = {
            "amount": 1, "unit": "copies", "basis": "assumed",
            "source": "One hypothetical input transcript encoding the displayed SLP with a retained N-terminal initiator Met; no transcript mass or formulation is known",
        }
    candidate["kinetics"] = synthetic_kinetics(candidate, loss)
    return candidate


def _scenario(selected, loss):
    candidates = [_candidate(row, ligand, start, end, disclosed_core, delivery, loss)
                  for row, ligand, start, end, disclosed_core in selected
                  for delivery in ("slp", "cytosolic_mrna")]
    request = {
        "schema_version": 1, "administration": "IM", "organism": "human scenario",
        "scenario": "Synthetic sensitivity: all unproductive loss/clearance rates = %g/h; not physiological calibration" % loss,
        "times_hours": TIMES, "candidates": candidates,
    }
    result = assess_vaccine_windows(request)
    if any(window["status"] != "conditional_scenario" for window in result["candidates"]):
        missing = [(window["id"], window["cores"][0]["missing_inputs"])
                   for window in result["candidates"] if window["status"] != "conditional_scenario"]
        raise ValueError("Unassessed candidates: %r" % missing)
    return result


def _flat(result, loss):
    rows = {}
    for window, summary in zip(result["candidates"], vaccine_window_rows(result)):
        core = window["cores"][0]
        meta = window["metadata"]
        entry = core["exposure"]["node_apc_entry"]
        curves = core["trajectory"]["curves"]
        retained = []
        for curve in curves:
            supplied = curve["antigen_copies_supplied"]
            retained.append((100 * curve["target_retained_copies"] / supplied)
                            if supplied > 0 else None)
        row = {
            "sequence_record_id": meta["sequence_record_id"],
            "gene": meta["gene"], "vaccines": meta["vaccines"],
            "source_url": meta["source_url"],
            "sequence_sha256": meta["sequence_sha256"],
            "disclosed_slp_sequence": meta["disclosed_slp_sequence"],
            "modeled_antigen_sequence": window["peptide_input"]["sequence"],
            "initiator_methionine_added": meta["initiator_methionine_added"],
            "delivery": window["delivery"], "loss_per_hour": loss,
            "core": core["sequence"], "core_start_zero_based": core["start"],
            "core_end_exclusive": core["end"], "allele": core["allele"],
            "mhcflurry_percentile_rank": meta["ligand_percentile_rank"],
            "disclosed_minimal_epitope": meta["disclosed_minimal_epitope"],
            "target_basis": meta["target_basis"],
            "node_apc_arrivals_per_input": entry["copies_per_input_molecule"],
            "node_apc_percent_slp_dose": summary["node_apc_percent_of_dosed_peptide"],
            "node_apc_percent_antigen_produced": (
                100 * entry["fraction_of_antigen_supplied_by_horizon"]
                if entry["fraction_of_antigen_supplied_by_horizon"] is not None else None),
            "mean_node_pmhc_per_input": summary["node_pmhc_mean_copies_per_input"],
            "mean_blood_core_per_input": summary["blood_core_mean_copies_per_input"],
            "blood_post_peak_half_decline_hours": core["exposure"]["blood"]["post_peak_half_decline"]["hours_after_sampled_peak_interval"],
            "retained_percent_at_48h": retained[-1],
            "retained_percent_curve": retained,
        }
        rows[(meta["sequence_record_id"], window["delivery"])] = row
    return rows


def _number(value):
    return "—" if value is None else "%.3g" % value


def _report(selected, scenarios, inputs, figure_name):
    central = scenarios[1.0]
    lines = [
        "# Sid's disclosed SLPs: paired route scenarios", "",
        "This is a **48-hour calculation under synthetic rates**, not a measured fraction",
        "of Sid's dose, a clinically calibrated prediction, or a vaccine recommendation.",
        "The 40 publicly disclosed SLP sequences are from [osteosarc.com/vaccines](https://osteosarc.com/vaccines/).",
        "For each, the target is the disclosed minimal epitope if it occurs in",
        "that SLP and has an exact frozen MHCflurry row; otherwise it is the",
        "lowest-percentile class-I ligand in that table. These are candidate",
        "ligands, **not observed presented epitopes**. Each disclosed sequence is",
        "modeled as a free SLP; the hypothetical mRNA product has a retained",
        "N-terminal initiator methionine, with the core coordinate shifted one",
        "residue. [Eukaryotic translation starts with methionine](https://www.ncbi.nlm.nih.gov/books/NBK9849/),",
        "but that residue is [removed from many proteins](https://www.ncbi.nlm.nih.gov/books/NBK9843/).",
        "This model does not account for such removal or its downstream effects.", "",
        "[Log-scale core-retention time courses](%s) summarize all 40 peptides." % figure_name,
        "The SLP panel is the percent of dosed peptide molecules with this exact",
        "core still present in any modeled fragment or pMHC state. The RNA panel",
        "is the percent of antigen copies produced by that time still retaining",
        "the core; at time zero its denominator is zero and the point is omitted.",
        "These are per-core marginals, not an 'any of several cores' union.", "",
        "## What the two routes mean", "",
        "- **SLP:** 1 µg of each unmodified free peptide, using its sequence-derived",
        "  molecular weight. Nodal APC arrival is a percent of input peptide copies.",
        "- **mRNA:** one input transcript encoding the SLP plus an initiator Met. RNA",
        "  arrival is antigen copies per input transcript, and also a percent of",
        "  antigen copies produced by 48 hours. The real transcript, UTRs, signal",
        "  peptide, LNP, and mass are unknown. Do not read it as an equal-µg comparison.",
        "  This cytosolic route includes direct APC expression, antigen transfer",
        "  and preformed pMHC transfer; secretion would require a different construct.",
        "- The existing 19 disclosed mRNA context sequences are a different",
        "  set of constructs and are outside this matched 40-SLP counterfactual.",
        "- All transport, translation, processing and loading rates are declared",
        "  synthetic assumptions from `examples/vaccine_sensitivity_scenarios.py`.",
        "  The loss/clearance scenarios are jointly 0.1, 1, and 10 per hour.",
        "  They are **not** physiological bounds. MHCflurry scores only select",
        "  the displayed target; they do not set any rate or pMHC loading affinity.", "",
        "## Per-peptide results", "",
        "The three-number ranges below run from 0.1 to 10 /h joint-loss scenarios;",
        "their endpoints are sensitivity results, not confidence intervals.",
        "Node APC arrival counts first core-bearing entry. Mean pMHC is surface",
        "copy-hours per input molecule divided by 48 hours; it is not a response",
        "probability. RNA and SLP means have different input-molecule denominators.", "",
        "| Gene / record | Chosen core / allele (rank %, basis) | SLP APC arrivals: central % (low–high loss) | RNA APC arrivals/transcript: central (low–high loss) | RNA APC arrivals / antigen made: central % | Mean node pMHC/input: SLP ; RNA |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for record, _, _, _, _ in selected:
        key = record["sequence_record_id"]
        slp, rna = (central[(key, route)] for route in ("slp", "cytosolic_mrna"))
        slp_low, slp_high = (scenarios[loss][(key, "slp")]["node_apc_percent_slp_dose"] for loss in (10.0, 0.1))
        rna_low, rna_high = (scenarios[loss][(key, "cytosolic_mrna")]["node_apc_arrivals_per_input"] for loss in (10.0, 0.1))
        lines.append("| %s / `%s` | `%s` / %s (%.2f, %s) | %s (%s–%s) | %s (%s–%s) | %s | %s ; %s |" % (
            record["gene"], key, slp["core"], slp["allele"], slp["mhcflurry_percentile_rank"],
            "disclosed" if slp["target_basis"] == "disclosed_minimal_epitope" else "predicted",
            _number(slp["node_apc_percent_slp_dose"]), _number(slp_low), _number(slp_high),
            _number(rna["node_apc_arrivals_per_input"]), _number(rna_low), _number(rna_high),
            _number(rna["node_apc_percent_antigen_produced"]),
            _number(slp["mean_node_pmhc_per_input"]), _number(rna["mean_node_pmhc_per_input"])))
    slp_values = [central[(row["sequence_record_id"], "slp")]["node_apc_percent_slp_dose"]
                  for row, _, _, _, _ in selected]
    rna_values = [central[(row["sequence_record_id"], "cytosolic_mrna")]["node_apc_arrivals_per_input"]
                  for row, _, _, _, _ in selected]
    lines.extend([
        "", "Across the 40 selected targets, the central-scenario median is",
        "**%s%%** SLP dose arriving in nodal APCs, and **%s** antigen arrivals" % (
            _number(statistics.median(slp_values)), _number(statistics.median(rna_values))),
        "per input transcript. The shared assumed rates, rather than measured",
        "sequence-specific kinetics, drive these values. SLP APC arrival is",
        "identical to displayed precision across all rows by construction; small",
        "RNA differences arise from the assumed route and target geometry.",
        "This is **not a validated ranking of Sid's peptides**.", "",
        "The CSV contains each scenario's core-retention percentages, blood",
        "exposure and sampled post-peak half-decline bracket. This blood decline",
        "is not an intrinsic circulation half-life while input continues.",
        "The model cannot derive a T-cell response or vaccine efficacy from pMHC",
        "alone. For empirical gaps and pathway evidence, see the",
        "[calibration audit](../../../../docs/research/vaccine-calibration.md) and",
        "[mRNA presentation research](../../../../docs/evidence/mrna-antigen-presentation.md).", "",
        "## Provenance", "",
        "Generated at %s UTC with mhctools %s." % (
            datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
            __import__("mhctools").__version__),
    ])
    for label, info in inputs.items():
        lines.append("- `%s`: SHA-256 `%s`" % (label, info))
    lines.extend([
        "- Script: `analyses/osteosarc_vaccine_cleavage/paired_exposure.py`.", "",
        "From the repository root, install `.[vaccine-report]` and reproduce",
        "into a **new** directory with:", "",
        "```sh",
        "PYTHONPATH=examples python analyses/osteosarc_vaccine_cleavage/paired_exposure.py \\",
        "  --output-dir analyses/osteosarc_vaccine_cleavage/results/NEW_RUN",
        "```", "",
    ])
    return "\n".join(lines)


def _plot(selected, scenarios, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, route in zip(axes, ("slp", "cytosolic_mrna")):
        for loss, color in zip(LOSSES, ("#317a4d", "#2673b8", "#ba5537")):
            curves = [scenarios[loss][(record["sequence_record_id"], route)]["retained_percent_curve"]
                      for record, _, _, _, _ in selected]
            medians = [statistics.median(x[i] for x in curves if x[i] is not None)
                       if any(x[i] is not None for x in curves) else None
                       for i in range(len(TIMES))]
            points = [(t, value) for t, value in zip(TIMES, medians)
                      if value is not None and value > 0]
            ax.plot([p[0] for p in points], [p[1] for p in points], marker="o",
                    markersize=3, color=color, label="joint loss %g/h" % loss)
        ax.set_yscale("log")
        ax.set_xlim(0, 48)
        ax.set_ylim(1e-22, 110)
        ax.grid(True, which="both", alpha=0.2)
        ax.set_title("1 µg free SLP" if route == "slp" else "cytosolic mRNA (per antigen made)")
        ax.set_ylabel("% input peptide copies retaining core" if route == "slp"
                      else "% produced antigen copies retaining core")
        ax.set_xlabel("Hours after IM input")
        ax.legend(fontsize=8)
    fig.suptitle("Median of 40 selected targets · synthetic sensitivity scenarios")
    fig.text(0.5, 0.01, "Log axis stops at 10⁻²²%; lower/zero modeled values are off-scale."
             " Shared rates make the 40 curves nearly coincide.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path, help="New directory; never overwrites a run")
    parser.add_argument("--inventory", type=Path, default=FROZEN / "vaccine_sequence_inventory.csv")
    parser.add_argument("--ligands", type=Path, default=FROZEN / "slp_mhc_ligand_predictions.csv")
    args = parser.parse_args(argv)
    if args.output_dir.exists():
        parser.error("Output directory already exists")
    selected = _records(args.inventory, args.ligands)
    scenarios = {loss: _flat(_scenario(selected, loss), loss) for loss in LOSSES}
    args.output_dir.mkdir(parents=True)
    figure = "core-retention-log.pdf"
    _plot(selected, scenarios, args.output_dir / figure)
    rows = []
    for loss in LOSSES:
        for row in scenarios[loss].values():
            rows.append({**row, "retained_percent_curve": json.dumps(row["retained_percent_curve"]),
                         "blood_post_peak_half_decline_hours": json.dumps(row["blood_post_peak_half_decline_hours"])})
    with (args.output_dir / "paired-exposure.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    inputs = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in (args.inventory, args.ligands)}
    (args.output_dir / "REPORT.md").write_text(_report(selected, scenarios, inputs, figure))
    print(args.output_dir / "REPORT.md")


if __name__ == "__main__":
    main()
