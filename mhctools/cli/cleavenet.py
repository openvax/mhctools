# Licensed under the Apache License, Version 2.0 (the "License");
# https://www.apache.org/licenses/LICENSE-2.0

"""Native whole-substrate MMP evidence, deliberately separate from bond tracks."""

from argparse import ArgumentParser
import json

from ..cleavenet import CleaveNet
from ..peptide_input import PeptideInput


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--peptides", nargs="+", required=True)
    parser.add_argument("--windows", action="store_true", help="Scan ten-residue windows")
    parser.add_argument("--source-name")
    parser.add_argument("--source-start", type=int, default=0)
    parser.add_argument("--cleavenet-home")
    parser.add_argument("--cleavenet-python")
    args = parser.parse_args(argv)
    inputs = [PeptideInput(sequence, source_sequence_name=args.source_name,
                           source_start=args.source_start) for sequence in args.peptides]
    model = CleaveNet(cleavenet_home=args.cleavenet_home,
                      cleavenet_python=args.cleavenet_python)
    results = ([result for item in inputs for result in model.predict_windows(item)]
               if args.windows else model.predict(inputs))
    print(json.dumps({"results": [result.to_dict() for result in results]},
                     indent=2, allow_nan=False))
    return 0
