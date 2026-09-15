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
        description="Report peptidase motif evidence and native model scores "
                    "per bond, as JSON.")
    parser.add_argument("--list-models", action="store_true",
                        help="List the curated cleavage models with their "
                             "enzymes, compartments and evidence, then exit")
    parser.add_argument("--sequence", action="append", help="Canonical linear L-peptide; repeat for multiple inputs")
    parser.add_argument("--model", action="append",
                        help="Exact model name (see --list-models); repeat to "
                             "select multiple models. Defaults to every "
                             "always-available model for --compartment; "
                             "models needing a fetched asset (eramer-step) "
                             "must be named explicitly.")
    parser.add_argument("--compartment", help="Filter enzyme locations: serum, plasma, extracellular, cytosol, er, endosome")
    parser.add_argument("--n-term", choices=("free", "acetylated", "unknown"), default="free",
                        help="N-terminal chemistry of the input peptide "
                             "(default: %(default)s)")
    parser.add_argument("--c-term", choices=("free", "amidated", "unknown"), default="free",
                        help="C-terminal chemistry of the input peptide "
                             "(default: %(default)s)")
    parser.add_argument("--source-id",
                        help="Identifier of the protein the peptide came from, "
                             "echoed back in the report")
    parser.add_argument("--source-start", type=int, default=0,
                        help="Offset of the peptide within --source-id, used "
                             "to report source coordinates "
                             "(default: %(default)s)")
    parser.add_argument("--enzyme-state", action="append", default=[],
                        metavar="ENZYME=STATE",
                        help="Declare one enzyme's state, e.g. "
                             "'CPB2=active'; STATE is active, zymogen, "
                             "inactive or unknown. Repeat per enzyme.")
    parser.add_argument("--out", help="Write JSON to this path instead of stdout")
    args = parser.parse_args(argv)
    if args.list_models:
        if args.sequence or args.model or args.compartment or args.enzyme_state:
            parser.error("--list-models cannot be combined with prediction options")
        result = {"models": [asdict(m) for m in cleavage_models(include_optional=True)]}
    else:
        if not args.sequence:
            parser.error("Provide --sequence or --list-models")
        try:
            states = {}
            for value in args.enzyme_state:
                enzyme, state = value.split("=", 1)
                if enzyme in states:
                    raise ValueError("Duplicate enzyme state for %s" % enzyme)
                states[enzyme] = state
            result = {"results": [
                prediction.to_dict()
                for seq in args.sequence
                for prediction in predict_cleavage(CleavageInput(
                    seq, args.n_term, args.c_term, args.source_id, args.source_start),
                    models=args.model, compartment=args.compartment, enzyme_states=states)]}
        except Exception as error:
            # Deliberately broad: this is the CLI's user-facing error
            # boundary. Resolving or running a named model can load an
            # external asset (e.g. eramer-step's PWM workbook), and any
            # failure there should exit cleanly with parser.error() rather
            # than a raw traceback, regardless of the exception type raised
            # by that asset's own loader.
            parser.error(str(error))
    result["schema_version"] = 1
    output = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output, end="")
