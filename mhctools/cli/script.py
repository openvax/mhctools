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

import sys

from pyensembl.fasta import parse_fasta_dictionary

from .args import (
    _flatten_predictor_names,
    make_mhc_arg_parser,
    predictors_from_args,
)
from .errors import CLI_ERROR_TYPES, cli_error_message


arg_parser = make_mhc_arg_parser(
    prog="mhctools",
    description=("Predict MHC ligands from protein sequences."))

def add_input_args(arg_parser):
    input_group = arg_parser.add_argument_group("Inputs")
    input_source = input_group.add_mutually_exclusive_group()
    input_source.add_argument(
        "--sequence",
        action="extend",
        nargs="+",
        help=(
            "Peptide sequences; may be repeated"))
    input_group.add_argument(
        "--extract-subsequences",
        default=False,
        action="store_true",
        help=(
            "Extract subsequences from peptides supplied by --sequence or "
            "--input-peptides-file, lengths specified by "
            "--mhc-peptide-lengths argument."))
    input_source.add_argument(
        "--input-peptides-file",
        help="Path to file with one peptide per line; blank lines are ignored")
    input_source.add_argument(
        "--input-fasta-file",
        help="Path to FASTA file which contains protein sequences")
    return input_group

def add_filter_args(parser):
    filter_group = parser.add_argument_group(
        "Filtering",
        description="Filter predictions before output. "
        "All active filters must pass for a row to be kept.")
    filter_group.add_argument(
        "--max-affinity",
        type=float,
        default=None,
        help="Keep predictions with affinity (IC50 nM) <= this value. "
             "E.g. --max-affinity 500")
    filter_group.add_argument(
        "--max-percentile-rank",
        type=float,
        default=None,
        help="Keep predictions with percentile rank <= this value. "
             "E.g. --max-percentile-rank 2")
    filter_group.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="Keep predictions with normalized score >= this value. "
             "E.g. --min-score 0.5")
    return filter_group

def add_output_args(parser):
    output_group = parser.add_argument_group("Outputs")
    output_group.add_argument(
        "--output-csv",
        default=None,
        help=(
            "Write the prediction table to this CSV path. affinity is IC50 "
            "in nM; percentile_rank is a percentile; score is "
            "predictor-specific."))
    return output_group

add_input_args(arg_parser)
add_filter_args(arg_parser)
add_output_args(arg_parser)

def parse_args(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    return arg_parser.parse_args(args_list)

def _run_single_predictor(predictor, args):
    # The legacy prediction CLI is built on the BindingPrediction model
    # (predict_peptides / predict_subsequences). New-model-only predictors
    # (e.g. bigmhc, calis, deeptap, eramer, netchop, pepsickle) implement only
    # predict(); route the user to the predict-table subcommand rather than
    # failing later with an opaque AttributeError.
    if not hasattr(predictor, "predict_peptides"):
        raise ValueError(
            "%s does not support this command (it implements the new "
            "prediction model only). Use `mhctools predict-table` instead."
            % type(predictor).__name__)
    if args.input_fasta_file:
        input_dictionary = parse_fasta_dictionary(args.input_fasta_file)
        if not input_dictionary:
            raise ValueError(
                "No sequences could be parsed from fasta file: %s" % (
                    args.input_fasta_file))
        # Capitalize sequences
        input_dictionary = dict(
            (key, value.upper()) for (key, value) in input_dictionary.items())
        return predictor.predict_subsequences(input_dictionary)
    elif args.sequence:
        if args.extract_subsequences:
            return predictor.predict_subsequences(args.sequence)
        else:
            return predictor.predict_peptides(args.sequence)
    elif args.input_peptides_file:
        with open(args.input_peptides_file) as f:
            peptides = [line.strip() for line in f if line.strip()]
        if not peptides:
            raise ValueError(
                "No peptide sequences found in file: %s" % (
                    args.input_peptides_file,))
        if args.extract_subsequences:
            return predictor.predict_subsequences(peptides)
        else:
            return predictor.predict_peptides(peptides)
    else:
        raise ValueError(
            ("No input sequences provided, "
             "use --sequence, --input-fasta-file, or input-peptides-file"))


def run_predictor(args):
    from mhctools.binding_prediction_collection import BindingPredictionCollection
    predictors = predictors_from_args(args)
    predictor_names = _flatten_predictor_names(args)
    has_source_coordinates = bool(
        args.input_fasta_file or args.extract_subsequences)
    all_predictions = []
    for predictor_name, predictor in zip(predictor_names, predictors):
        results = _run_single_predictor(predictor, args)
        for prediction in results:
            updates = {"prediction_method_name": predictor_name}
            if not has_source_coordinates:
                # Names such as PEPLIST / Sequence and row-number offsets are
                # scratch-file details emitted by some wrapped tools. Plain
                # peptide input has no source coordinates, so expose the same
                # representation regardless of predictor.
                updates.update(source_sequence_name="", offset=0)
            all_predictions.append(prediction.clone_with_updates(**updates))
    if has_source_coordinates:
        all_predictions.sort(key=lambda prediction: (
            prediction.source_sequence_name or "",
            prediction.offset,
        ))
    return BindingPredictionCollection(all_predictions)

def apply_filters(df, args):
    """Apply --max-affinity, --max-percentile-rank, --min-score filters."""
    if args.max_affinity is not None:
        df = df[df["affinity"].isna() | (df["affinity"] <= args.max_affinity)]
    if args.max_percentile_rank is not None:
        df = df[df["percentile_rank"].isna() | (df["percentile_rank"] <= args.max_percentile_rank)]
    if args.min_score is not None:
        df = df[df["score"].isna() | (df["score"] >= args.min_score)]
    return df


def format_predictions(df):
    """Render a prediction DataFrame as an aligned, index-free text table.

    Parameters
    ----------
    df : pandas.DataFrame
        Prediction table from
        :meth:`BindingPredictionCollection.to_dataframe`.

    Returns
    -------
    str
        Every row and column of ``df``, with floats trimmed to six
        significant digits, or a short notice when ``df`` is empty.
    """
    if len(df) == 0:
        return "No predictions."
    return df.to_string(index=False, float_format=lambda value: "%.6g" % value)


def main(args_list=None):
    """
    Script to make pMHC binding predictions from amino acid sequences.

    Usage example:
        mhctools
            --sequence SFFPIQQQQQAAALLLI \
            --sequence SILQQQAQAQQAQAASSSC \
            --extract-subsequences \
            --mhc-predictor netmhc \
            --mhc-alleles HLA-A0201 H2-Db \
            --mhc-predictor netmhc \
            --output-csv epitope.csv

    The ``predict-table`` subcommand annotates an existing CSV of peptides
    (and optional alleles) with predictor score columns; see
    ``mhctools predict-table --help``.
    """
    if args_list is None:
        args_list = sys.argv[1:]
    if args_list and args_list[0] == "ls":
        from .artifacts import ls_main
        return ls_main(args_list[1:])
    if args_list and args_list[0] == "fetch":
        from .artifacts import fetch_main
        return fetch_main(args_list[1:])
    if args_list and args_list[0] in ("predictors", "integrations"):
        from .integrations import predictors_main
        return predictors_main(args_list[1:], command_name=args_list[0])
    if args_list and args_list[0] == "benchmark":
        from .benchmark import main as benchmark_main
        return benchmark_main(args_list[1:])
    if args_list and args_list[0] == "cleavage":
        from .cleavage import main as cleavage_main
        return cleavage_main(args_list[1:])
    if args_list and args_list[0] == "predict-table":
        from .annotate_table import main as annotate_table_main
        return annotate_table_main(args_list[1:])
    if args_list and args_list[0] == "mixtcrpred":
        from .mixtcrpred import main as mixtcrpred_main
        return mixtcrpred_main(args_list[1:])

    args = parse_args(args_list)
    try:
        binding_predictions = run_predictor(args)
        df = binding_predictions.to_dataframe()
        n_before = len(df)
        df = apply_filters(df, args)
        n_after = len(df)
        if n_before != n_after:
            # Keep the note off stdout so the table stays machine-readable.
            print("Filtered %d -> %d predictions" % (n_before, n_after),
                  file=sys.stderr)
        if args.output_csv:
            # Don't also dump the table to a terminal the user redirected to a
            # file; a proteome-wide scan is millions of rows.
            df.to_csv(args.output_csv, index=False, float_format="%.6g")
            print("Wrote: %s (%d rows, %d columns)" % (
                args.output_csv, len(df), len(df.columns)))
        else:
            print(format_predictions(df))
    except CLI_ERROR_TYPES as error:
        arg_parser.error(cli_error_message(error))
