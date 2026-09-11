"""Evaluate supplied measurements and predictions without merging assay scales."""

import argparse
import json
from pathlib import Path

from mhctools._resources import load_json_resource
from mhctools.benchmark import (
    AssayMeasurement, BenchmarkPrediction, ModelLineage, evaluate_benchmark,
    model_lineage_inventory, predict_cleavage_measurements,
)

#: Single source of truth mapping a --reference-cleavage choice to its packaged
#: data file, so the two can never drift apart.
REFERENCE_PANELS = {
    "starter": "cleavage_reference.json",
    "serum": "serum_cleavage_reference.json",
    "intracellular": "intracellular_cleavage_reference.json",
}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="mhctools benchmark")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", help="JSON containing measurements and predictions or --model inputs")
    mode.add_argument("--lineage-inventory", action="store_true")
    mode.add_argument("--reference-cleavage", nargs="?", const="starter", choices=tuple(REFERENCE_PANELS),
                      help="Run a source-linked reproduction panel (default: starter)")
    parser.add_argument("--model", action="append", help="Run an exact cleavage model against site records")
    parser.add_argument("--evaluation", choices=("external_validation", "reproduction"), default="external_validation")
    parser.add_argument("--out")
    args = parser.parse_args(argv)
    try:
        if args.lineage_inventory:
            if args.model:
                parser.error("--model cannot be combined with --lineage-inventory")
            result = model_lineage_inventory()
        else:
            if args.reference_cleavage:
                filename = REFERENCE_PANELS[args.reference_cleavage]
                data = load_json_resource(filename)
                evaluation = "reproduction"
            else:
                data = json.loads(Path(args.input).read_text())
                evaluation = args.evaluation
            measurements = [AssayMeasurement(**m) for m in data["measurements"]]
            if args.model and data.get("predictions"):
                parser.error("Choose supplied predictions or --model, not both")
            if args.model or args.reference_cleavage:
                predictions = predict_cleavage_measurements(measurements, args.model or data["models"])
            else:
                predictions = [BenchmarkPrediction(**p) for p in data.get("predictions", [])]
            lineages = [ModelLineage(**x) for x in data.get("lineages", [])]
            result = evaluate_benchmark(measurements, predictions, lineages, evaluation, data.get("requested_domains", ()))
            if args.reference_cleavage:
                result["dataset_notice"] = data["notice"]
        output = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.out:
            Path(args.out).write_text(output)
        else:
            print(output, end="")
    except (ValueError, TypeError, KeyError, OSError) as error:
        parser.error(str(error))
