# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""``mhctools predict-table`` — annotate a CSV of peptides with predictor scores.

Thin I/O wrapper around :func:`mhctools.annotate.annotate_table`: reads a
CSV (optionally compressed), runs the requested predictors, and writes the
input table back with one appended score column per predictor.
"""

from argparse import ArgumentParser

import pandas as pd

from ..annotate import (
    annotate_table,
    parse_annotation_spec,
    output_field_tokens,
)
from ..pred import best_direction
from ..logging import get_logger


logger = get_logger(__name__)


def make_arg_parser():
    parser = ArgumentParser(
        prog="mhctools predict-table",
        description=(
            "Annotate a CSV table of peptides (and optional alleles) with "
            "predictor score columns, preserving all input columns."))
    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV (compression inferred from extension: .csv, .bz2, .gz).")
    parser.add_argument(
        "--out",
        required=True,
        help="Output CSV path (compression inferred from extension).")
    parser.add_argument(
        "--peptide-column",
        default="peptide",
        help="Column holding the peptide sequence (default: peptide).")
    parser.add_argument(
        "--alleles-column",
        default=None,
        help=(
            "Column holding the row's allele(s). Multiple alleles per cell "
            "may be separated by whitespace, comma, or semicolon. Omit for "
            "allele-free predictors."))
    parser.add_argument(
        "--predictor",
        dest="predictors",
        action="append",
        required=True,
        metavar="NAME:OUTPUT_COLUMN:FIELD",
        help=(
            "Predictor spec, repeatable. OUTPUT_COLUMN and FIELD are optional "
            "(default column NAME_FIELD, default field 'affinity'). FIELD is "
            "one of: %s. E.g. 'netmhcpan42-ba:netmhcpan4.2.ba:affinity'."
            % ", ".join(output_field_tokens())))
    parser.add_argument(
        "--predictor-info",
        default=None,
        help=(
            "Optional path to write a sidecar CSV describing each output "
            "column: predictor spec, output_column, score_field, "
            "higher_is_better."))
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite output columns if they already exist in the input.")
    return parser


def _write_predictor_info(path, raw_specs, specs):
    rows = []
    for raw, spec in zip(raw_specs, specs):
        rows.append({
            "predictor": raw,
            "output_column": spec.output_column,
            "score_field": spec.field,
            "higher_is_better":
                best_direction(spec.kind, spec.prediction_field) == "max",
        })
    pd.DataFrame(rows).to_csv(path, index=False)
    print("Wrote predictor info: %s" % path)


def main(args_list=None):
    args = make_arg_parser().parse_args(args_list)

    df = pd.read_csv(args.input)
    logger.info("Read %d rows, %d columns from %s",
                len(df), len(df.columns), args.input)

    specs = [parse_annotation_spec(token) for token in args.predictors]

    annotated = annotate_table(
        df,
        specs,
        peptide_column=args.peptide_column,
        allele_column=args.alleles_column,
        overwrite=args.overwrite)

    annotated.to_csv(args.out, index=False)
    print("Wrote: %s (%d rows, %d columns)"
          % (args.out, len(annotated), len(annotated.columns)))

    if args.predictor_info:
        _write_predictor_info(args.predictor_info, args.predictors, specs)


if __name__ == "__main__":
    main()
