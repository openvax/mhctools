"""JSON reporting for per-bond cleavage evidence."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from mhctools.cleavage import CleavageInput
from mhctools.peptidases import cleavage_models, predict_cleavage


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mhctools cleavage",
        description="Report peptidase motif evidence and native model scores per bond.")
    parser.add_argument("--list-models", action="store_true")
    parser.add_argument("--sequence", action="append", help="Canonical linear L-peptide; repeat for multiple inputs")
    parser.add_argument("--model", action="append", help="Exact model name; repeat to select multiple models")
    parser.add_argument("--compartment", help="Filter enzyme locations: serum, plasma, extracellular, cytosol, er")
    parser.add_argument("--n-term", choices=("free", "acetylated", "unknown"), default="free")
    parser.add_argument("--c-term", choices=("free", "amidated", "unknown"), default="free")
    parser.add_argument("--source-id")
    parser.add_argument("--source-start", type=int, default=0)
    parser.add_argument("--out", help="Write JSON to this path instead of stdout")
    args = parser.parse_args(argv)
    if args.list_models:
        if args.sequence or args.model or args.compartment:
            parser.error("--list-models cannot be combined with prediction options")
        result = {"models": [asdict(m) for m in cleavage_models(include_optional=True)]}
    else:
        if not args.sequence:
            parser.error("Provide --sequence or --list-models")
        try:
            result = {"results": [
                prediction.to_dict()
                for seq in args.sequence
                for prediction in predict_cleavage(CleavageInput(
                    seq, args.n_term, args.c_term, args.source_id, args.source_start),
                    models=args.model, compartment=args.compartment)]}
        except (ValueError, TypeError, OSError, ImportError) as error:
            parser.error(str(error))
    result["schema_version"] = 1
    output = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output, end="")
