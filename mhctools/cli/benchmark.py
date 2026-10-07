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
    parser = argparse.ArgumentParser(
        prog="mhctools benchmark",
        description="Evaluate supplied measurements and predictions without "
                    "merging assay scales. Writes a JSON report.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input", help="JSON containing measurements and predictions or --model inputs")
    mode.add_argument("--lineage-inventory", action="store_true",
                      help="Report the curated model/dataset relationship "
                           "inventory and exit")
    mode.add_argument("--reference-cleavage", nargs="?", const="starter", choices=tuple(REFERENCE_PANELS),
                      help="Run a source-linked reproduction panel (default: starter)")
    mode.add_argument("--peptiverse-cpp-metadata", help="Audit the pinned local CPP split CSV")
    parser.add_argument("--predict-cpp", action="store_true",
                        help="Opt into offline native CPP source-label reproduction")
    parser.add_argument("--source-id", action="append",
                        help="Select a validation source ID for a partial CPP prediction cohort (repeatable)")
    parser.add_argument("--model", action="append", help="Run an exact cleavage model against site records")
    parser.add_argument("--evaluation", choices=("external_validation", "reproduction"), default="external_validation",
                       help="How to label the run for --input data "
                            "(default: %(default)s). --reference-cleavage "
                            "and --peptiverse-cpp-metadata are always reproduction.")
    parser.add_argument("--out", help="Write JSON to this path instead of stdout")
    args = parser.parse_args(argv)
    try:
        if (args.predict_cpp or args.source_id) and not args.peptiverse_cpp_metadata:
            parser.error("--predict-cpp and --source-id require --peptiverse-cpp-metadata")
        if args.peptiverse_cpp_metadata:
            if args.model:
                parser.error("--model cannot be combined with --peptiverse-cpp-metadata")
            from mhctools.peptiverse_cpp_benchmark import evaluate_cpp_metadata
            result = evaluate_cpp_metadata(args.peptiverse_cpp_metadata,
                                          predict=args.predict_cpp, source_ids=args.source_id)
        elif args.lineage_inventory:
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
