# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""CSV prediction command for the MixTCRpred wrapper."""

from argparse import ArgumentParser

import pandas as pd


def main(args_list=None):
    parser = ArgumentParser(
        prog="mhctools mixtcrpred",
        description=(
            "Score a paired-alpha/beta TCR table with one pMHC-specific "
            "MixTCRpred model."),
    )
    parser.add_argument("--model", required=True)
    parser.add_argument("--input", required=True, help="Input TCR CSV")
    parser.add_argument("--out", required=True, help="Output annotated CSV")
    parser.add_argument(
        "--mixtcrpred-path",
        help="Licensed upstream checkout (default: MIXTCRPRED_HOME/cache)")
    parser.add_argument(
        "--mixtcrpred-python",
        help="Python with the MixTCRpred optional dependencies")
    parser.add_argument("--batch-size", type=int, default=256)
    args = parser.parse_args(args_list)

    # Keep heavyweight/optional predictor imports out of ordinary CLI startup.
    from ..mixtcrpred import MixTCRpred
    try:
        predictor = MixTCRpred(
            model=args.model,
            mixtcrpred_path=args.mixtcrpred_path,
            mixtcrpred_python=args.mixtcrpred_python,
            batch_size=args.batch_size,
        )
        dataframe = pd.read_csv(args.input, keep_default_na=False)
        output = predictor.annotate_dataframe(dataframe)
    except (FileNotFoundError, RuntimeError, TypeError, ValueError) as error:
        parser.error(str(error))
    output.to_csv(args.out, index=False)
    print("Wrote: %s (%d rows, %d columns)" % (
        args.out, len(output), len(output.columns)))
