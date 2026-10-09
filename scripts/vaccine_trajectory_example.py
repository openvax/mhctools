"""Reproducible synthetic SLP/secreted-mRNA scenarios, not vaccine forecasts."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from mhctools import (
    DegradationTarget, EnzymeCutRate, VaccineCleavageRates, VaccineRate,
    VaccineTrajectoryInput, simulate_vaccine_trajectory, vaccine_route_steps,
)


def example(delivery):
    """Construct a fully labeled toy scenario; no empirical rate transfer."""
    offset = 4 if delivery == "secreted_mrna" else 0
    peptide = "MLLL"[:offset] + "ACDEFGHIKLMNPQRSTVWY"
    rates = {step.name: VaccineRate(0, "disabled", "Excluded branch in synthetic demonstration")
             for step in vaccine_route_steps(delivery, "I")}
    values = dict(lymph_entry=.8, vascular_absorption=.25, lymph_transit=1.2, node_outflow=.12,
                  local_antigen_uptake=.18, node_antigen_uptake=.7, systemic_antigen_uptake=.02,
                  interstitium_clearance=.1, afferent_lymph_clearance=.03,
                  node_fluid_clearance=.03, blood_clearance=1.1)
    for region in ("local", "node", "systemic"):
        values.update({region + "_cross_escape": .35, region + "_tap": 1.5,
                       region + "_er_loading": 2, region + "_er_loss": .1,
                       region + "_endosome_loss": .1, region + "_cytosol_loss": .1,
                       region + "_surface_export": 1, region + "_loaded_loss": .1,
                       region + "_surface_loss": .04})
    for phase in ("endosome", "cytosol", "er", "loaded", "surface"):
        values["migration_" + phase] = .06
    if offset:
        values.update(carrier_drainage=.25, producer_rna_uptake=.35, local_rna_uptake=.08,
                      node_rna_uptake=.4, rna_carrier_loss=.2, node_rna_carrier_loss=.15)
        for cell in ("producer", "local", "node"):
            values.update({cell + "_rna_escape": .5, cell + "_endosomal_rna_loss": .5,
                           cell + "_rna_decay": .12, cell + "_translation": 2,
                           cell + "_signal_entry": 5, cell + "_secretion": .8,
                           cell + "_secretory_loss": .12, cell + "_failed_entry": .2})
            if cell != "producer":
                values[cell + "_erad"] = .04
        for phase in ("rna_endosome", "rna", "protein", "secretory_er"):
            values["migration_" + phase] = .06
    rates.update({name: VaccineRate(value, "assumed", "Synthetic demonstration value; not estimated physiology")
                  for name, value in values.items()})
    model = VaccineTrajectoryInput(
        peptide, DegradationTarget("synthetic chosen target", 4+offset, 13+offset),
        "HLA-A*02:01", "I", delivery,
        "Synthetic rates for sensitivity demonstration; not Sid or patient predictions",
        rates, signal_end=offset or None, tap_max_length=16)

    def cleavage(compartment, fragment):
        length = len(fragment.sequence)
        channels = []

        def add(bond, rate, mechanism):
            if 0 < bond < length:
                channels.append(EnzymeCutRate(bond, mechanism, rate, "assumed",
                                             "Synthetic absolute hazard; no native predictor conversion"))

        local_end = model.target.end - fragment.start
        if compartment.endswith("_er"):
            add(1, 1.2 if fragment.start < model.target.start else .15, "assumed ER amino-terminal trimming")
        else:
            base = 1 if compartment == "blood" else .3 if compartment in (
                "interstitium", "afferent_lymph", "node_fluid") else .4
            add(1, base, "assumed N-terminal trimming")
            add(length-1, base/2, "assumed C-terminal trimming")
            add(model.target.start + 4 - fragment.start, .05, "assumed destructive endoproteolysis")
            if compartment.endswith("_cytosol") or compartment.endswith("_endosome"):
                add(local_end, 4 if compartment.endswith("_cytosol") else 1,
                    "assumed target C-terminal release")
        return VaccineCleavageRates(tuple(channels), "assumed", "Entire cut pattern is synthetic; enzymes are not identified")

    return model, cleavage


def write_examples(output):
    output.mkdir(parents=True, exist_ok=False)
    results = []
    times = [i/4 for i in range(193)]
    for name, delivery in (("slp", "slp"), ("mrna", "secreted_mrna")):
        model, cleavage = example(delivery)
        result = simulate_vaccine_trajectory(model, times, cleavage)
        if result["status"] != "conditional_scenario":
            raise RuntimeError("Synthetic example is missing required kinetics")
        payload = {"schema_version": 1, "peptide": model.peptide, "target": asdict(model.target),
                   "allele": model.allele, "mhc_class": model.mhc_class, "delivery": delivery,
                   "scenario": model.scenario, "signal_end": model.signal_end,
                   "tap_max_length": model.tap_max_length, "times_hours": times,
                   "rates": {key: asdict(rate) for key, rate in model.rates.items()},
                   "cleavage_assessments": result["cuts"]}
        for filename, value in ((name + "-input.json", payload), (name + "-trajectory.json", result)):
            (output / filename).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
        results.append(result)
    (output / "README.txt").write_text(
        "SYNTHETIC SCENARIOS, NOT SID OR HUMAN VACCINE FORECASTS\n\n"
        "Local IM/SC free SLP versus secretion-tagged mRNA-LNP. Each numerical rate\n"
        "is an assumed demonstration value; no physiological rates are estimated.\n"
        "SLP output is per input peptide; mRNA output is per input transcript.\n"
        "Their raw yields do not compare efficacy at an equal administered dose.\n"
        "The full translated mRNA example and its signal boundary are synthetic.\n"
        "Every JSON includes rate/cut sources, missing inputs and model exclusions.\n")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New directory for synthetic inputs/results")
    args = parser.parse_args()
    results = write_examples(args.output)
    for result in results:
        print(result["delivery"], "states:", len(result["states"]),
              "max balance error:", max(abs(row["antigen_balance_error"]) for row in result["curves"]))
    print(args.output)


if __name__ == "__main__":
    main()
