#!/usr/bin/env python3
"""Reproduce public reference kinetics and an explicitly synthetic cut scenario."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from mhctools import (
    DegradationTarget, empirical_terminal_cut_rates, enzyme_removal_effects,
    serum_contribution_evidence, serum_reference_kinetics,
)


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit_library(path, evidence):
    import pandas as pd

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != evidence["table_s1_sha256"]:
        raise ValueError("Table S1 hash differs from the curated source snapshot")
    table = pd.read_excel(path, sheet_name="Cleavage Sites")
    columns = {"1": "Bond1-2", "2": "Bond2-3", "3": "Bond3-4", "4": "Bond4-5",
               "5": "Bond5-6", "10": "Bond10-11", "12": "Bond12-13", "13": "Bond 13-14"}
    counts = {bond: int(table[column].count()) for bond, column in columns.items()}
    if counts != evidence["boundary_context_counts_by_bond"]:
        raise ValueError("Boundary counts differ from the published source worksheet")
    return dict(table_sha256=digest, source_boundary_counts=counts, matched=True)


def run(output, library_table=None, n_paths=20000):
    output.mkdir(parents=True, exist_ok=True)
    evidence = serum_contribution_evidence()
    reference_rows, calculations = [], []
    for concentration in (5.0, 10.0, 70.0):
        for substrate in evidence["kinetic_reference"]["substrates"]:
            result = serum_reference_kinetics(
                substrate, concentration_um=concentration, serum_fraction=0.1)
            calculations.append(result)
            for enzyme in result["enzymes"]:
                reference_rows.append(dict(
                    substrate_id=substrate, concentration_um=concentration,
                    serum_fraction=0.1, enzyme=enzyme["enzyme"], product=enzyme["product"],
                    local_rate_per_hour=enzyme["local_rate_per_hour"],
                    share_of_reported_channels=enzyme["share_of_reported_channels"],
                    source=result["source"]))

    def rates(sequence):
        return empirical_terminal_cut_rates(
            sequence, total_rate_per_hour=math.log(2) / .5, allow_transfer=True)

    synthetic = enzyme_removal_effects(
        "AGDEFGHIKLNPQR", DegradationTarget("synthetic central target FGHIKL", 4, 10),
        rates, enzymes=("n_mono", "n_di", "c_mono", "c_di"),
        times_hours=[i / 4 for i in range(33)], n_paths=n_paths, seed=42,
        scenario="Transferred pig-plasma terminal-only boundary spectrum; assumed 0.5 h next-cut half-life per fragment")
    payload = dict(
        evidence=evidence, reference_calculations=calculations, synthetic_scenario=synthetic,
        library_audit=audit_library(library_table, evidence["plasma_library"]) if library_table else None,
        limitations="Reference kinetics are assay/substrate-specific. Synthetic trajectories are not vaccine predictions. No universal named-enzyme serum weights or calibrated human epitope half-life.")
    (output / "serum-model.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    write_csv(output / "reference-kinetics.csv", reference_rows)
    write_csv(output / "synthetic-target-survival.csv", synthetic["baseline"])
    write_csv(output / "synthetic-removal-effects.csv", synthetic["enzyme_removal_effects"])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return payload
    rows = synthetic["baseline"]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot([r["time_hours"] for r in rows], [100 * r["parent_remaining"] for r in rows],
            color="#606b78", linestyle="--", linewidth=2.5, label="Full parent remaining")
    ax.plot([r["time_hours"] for r in rows], [100 * r["target_in_circulation"] for r in rows],
            color="#1269b0", linewidth=3, label="Complete target retained, including fragments")
    ax.set(xlabel="Time (hours)", ylabel="Starting copies retained (%)", ylim=(0, 102), xlim=(0, 8))
    ax.set_title("Successive terminal cuts can spare a central target", loc="left", fontsize=15, pad=20)
    ax.grid(alpha=.2)
    ax.legend(frameon=False, fontsize=10)
    fig.text(.1, .015, "Synthetic example: assumed 0.5 h next-cut clock; transferred pig-plasma boundary counts.\n"
             "Terminal-only scenario omits internal cuts. This is not a human vaccine half-life prediction.",
             fontsize=9, color="#48515c")
    fig.tight_layout(rect=(0, .09, 1, 1))
    fig.savefig(output / "synthetic-target-survival.png", dpi=160)
    plt.close(fig)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--library-table", type=Path)
    parser.add_argument("--n-paths", type=int, default=20000)
    args = parser.parse_args()
    run(args.output, args.library_table, args.n_paths)
    print(args.output)


if __name__ == "__main__":
    main()
