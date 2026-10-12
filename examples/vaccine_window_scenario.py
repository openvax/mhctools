"""Write a synthetic two-window/two-core input; no human kinetic calibration."""

import json
import sys

from mhctools import vaccine_route_steps


def example_request():
    antigen = "AACDEFGHIKLMNPQRSTVWYACDEFGHI"
    candidates = []
    for index, (start, clearance) in enumerate(((0, .5), (2, 1.5))):
        sequence = antigen[start:start + 25]
        rates = {step.name: {"value_per_hour": 0, "basis": "disabled", "source": "Synthetic example: exclude unused pathway"}
                 for step in vaccine_route_steps("slp", "I")}
        for name, value in {"lymph_entry": 1, "lymph_transit": 1, "node_antigen_uptake": 1,
                            "vascular_absorption": .5, "blood_clearance": .5,
                            "interstitium_clearance": clearance}.items():
            rates[name] = {"value_per_hour": value, "basis": "assumed", "source": "Synthetic sensitivity example, not a measured or sequence-predicted rate"}
        cuts = [{"compartment": compartment, "start": 0, "end": len(sequence), "sequence": sequence,
                 "channels": [], "basis": "assumed", "source": "Synthetic example: exclude cleavage"}
                for compartment in ("interstitium", "afferent_lymph", "node_fluid", "blood", "node_endosome")]
        candidates.append({"id": "window-%s" % (index + 1), "peptide": sequence, "delivery": "slp",
                           "metadata": {"source_start": start, "note": "Different clearance rates are illustrative assumptions, not inferred from these flanks"},
                           "cores": [{"label": "core-%s" % (i + 1), "start": a-start, "end": b-start,
                                      "allele": "HLA-A*02:01", "mhc_class": "I", "score_contribution": score,
                                      "score_source": "Arbitrary demonstration weight; not a binding prediction"}
                                     for i, (a, b, score) in enumerate(((3, 12, 1), (11, 20, 2)))],
                           "kinetics": {"I": {"tap_max_length": 25, "rates": rates, "cleavage_assessments": cuts}}})
    return {"schema_version": 1, "scenario": "Synthetic sensitivity example; NOT human IM vaccine estimates",
            "administration": "IM", "times_hours": [0, .5, 1, 2, 4, 8, 24, 48], "candidates": candidates}


if __name__ == "__main__":
    json.dump(example_request(), sys.stdout, indent=2, allow_nan=False)
    print()
