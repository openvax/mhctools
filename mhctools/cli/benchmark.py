"""Evaluate supplied measurements and predictions without merging assay scales."""

import argparse
import json
from pathlib import Path
from importlib.resources import files

from mhctools.benchmark import (
    AssayMeasurement, BenchmarkPrediction, ModelLineage, evaluate_benchmark,
    model_lineage_inventory, predict_cleavage_measurements,
)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="mhctools benchmark")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", help="JSON containing measurements and predictions or --model inputs")
    mode.add_argument("--lineage-inventory", action="store_true")
    mode.add_argument("--reference-cleavage", nargs="?", const="starter", choices=("starter", "serum"),
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
                filename = "cleavage_reference.json" if args.reference_cleavage == "starter" else "serum_cleavage_reference.json"
                data = json.loads(files("mhctools").joinpath("data/" + filename).read_text())
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
