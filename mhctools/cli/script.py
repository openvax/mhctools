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

from argparse import RawDescriptionHelpFormatter
from collections import namedtuple
import csv
from importlib import import_module
from io import StringIO
import math
import os
import sys

from pyensembl.fasta import parse_fasta_dictionary

from .args import (
    _flatten_predictor_names,
    make_mhc_arg_parser,
    predictors_from_args,
)
from .errors import CLI_ERROR_TYPES, cli_error_message


class Subcommand(namedtuple(
        "Subcommand", ["name", "help", "module", "entry_point", "aliases"])):
    """One `mhctools <command>` entry.

    The module is imported on dispatch, not at startup: several subcommands
    pull heavy optional runtimes that a plain prediction run must not pay for.
    """

    def load(self):
        module = import_module(self.module, package=__package__)
        return getattr(module, self.entry_point)


def _subcommand(name, help, module, entry_point, aliases=()):
    return Subcommand(name, help, module, entry_point, aliases)


# The single source of truth for `mhctools <command>`: dispatch reads it and
# --help renders it, so a new subcommand cannot appear in one and not the
# other (#352).
#
# These are dispatched before argparse rather than through add_subparsers for
# two reasons. add_mhc_args declares --mhc-predictor required=True, so a
# subparser on this parser makes every subcommand fail on the missing flag;
# relaxing that is #444's decision and changes behavior for topiary, which
# embeds add_mhc_args. And argparse's REMAINDER does not capture a
# subcommand's own leading options, so `ls --all` would become "unrecognized
# arguments: --all" unless every subcommand exposed its parser for
# registration here.
SUBCOMMANDS = (
    _subcommand(
        "ls", "List model weights and other predictor artifacts",
        ".artifacts", "ls_main"),
    _subcommand(
        "fetch", "Fetch artifacts required by a predictor wrapper",
        ".artifacts", "fetch_main"),
    _subcommand(
        "predictors",
        "Report artifact location, launchability and reference reproduction",
        ".integrations", "predictors_main", aliases=("integrations",)),
    _subcommand(
        "predict-table",
        "Annotate a CSV of peptides (and optional alleles) with score columns",
        ".annotate_table", "main"),
    _subcommand(
        "benchmark",
        "Evaluate supplied measurements and predictions; writes a JSON report",
        ".benchmark", "main"),
    _subcommand(
        "cleavage",
        "Report peptidase motif evidence and native model scores per bond",
        ".cleavage", "main"),
    _subcommand(
        "cleavenet", "Score whole substrates or ten-residue windows for 18 MMPs",
        ".cleavenet", "main"),
    _subcommand(
        "mixtcrpred",
        "Score a paired alpha/beta TCR table with one MixTCRpred model",
        ".mixtcrpred", "main"),
    _subcommand(
        "vaccine-report",
        "Generate a route-aware peptide cleavage and MHC-window report",
        ".vaccine_report", "main"),
)


def find_subcommand(name):
    """Return the Subcommand invoked by `name`, or None for the legacy form."""
    for subcommand in SUBCOMMANDS:
        if name == subcommand.name or name in subcommand.aliases:
            return subcommand
    return None


def _subcommand_epilog():
    width = max(len(subcommand.name) for subcommand in SUBCOMMANDS)
    lines = ["Subcommands:"]
    for subcommand in SUBCOMMANDS:
        lines.append("  %-*s  %s" % (width, subcommand.name, subcommand.help))
    lines.append("")
    lines.append(
        "Run `mhctools <command> --help` for a subcommand's own options. "
        "The flags above\nbelong to the bare prediction command.")
    return "\n".join(lines)


arg_parser = make_mhc_arg_parser(
    prog="mhctools",
    description=("Predict MHC ligands from protein sequences."),
    epilog=_subcommand_epilog(),
    formatter_class=RawDescriptionHelpFormatter)

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


def _prediction_cell(value):
    """Return one lossless-enough, machine-readable stdout field."""
    if value is None:
        return ""
    try:
        if math.isnan(value):
            return ""
    except TypeError:
        pass
    if isinstance(value, float):
        return "%.6g" % value
    return str(value)


def write_predictions(df, stream):
    """Stream a tab-separated prediction table without a full-table buffer.

    TSV preserves empty strings as actual fields and remains readable at a
    terminal.  Rows are written one at a time, so a proteome-wide result does
    not create a second, formatted copy of the whole table in memory.

    An empty frame writes nothing. "No predictions." is a message for the
    reader, not a row, and this stream is the machine-readable one: emitted
    here it is consumed as data by ``mhctools ... | awk 'NR > 1 {print $3}'``.
    The CLI reports the empty case on stderr instead.
    """
    if len(df) == 0:
        return
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(df.columns)
    for row in df.itertuples(index=False, name=None):
        writer.writerow(_prediction_cell(value) for value in row)


def format_predictions(df):
    """Render predictions as TSV for callers that explicitly need a string.

    Parameters
    ----------
    df : pandas.DataFrame
        Prediction table from
        :meth:`BindingPredictionCollection.to_dataframe`.

    Returns
    -------
    str
        Every row and column of ``df`` as tab-separated fields, with floats
        trimmed to six significant digits. An empty frame renders as the empty
        string. The CLI uses :func:`write_predictions` directly so it streams.
    """
    stream = StringIO()
    write_predictions(df, stream)
    return stream.getvalue().rstrip("\n")


def _silence_closed_stdout():
    """Prevent Python's shutdown flush from reporting a second broken pipe."""
    try:
        stdout_fd = sys.stdout.fileno()
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        try:
            os.dup2(devnull_fd, stdout_fd)
        finally:
            os.close(devnull_fd)
    except (AttributeError, OSError, ValueError):
        pass


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
    if not args_list:
        # A bare `mhctools` is someone looking for the interface, not a
        # malformed prediction run. Show it instead of an argparse error.
        arg_parser.print_help()
        return 0
    subcommand = find_subcommand(args_list[0])
    if subcommand is not None:
        entry_point = subcommand.load()
        if subcommand.aliases:
            # Subcommands with an alias report the name they were invoked by.
            return entry_point(args_list[1:], command_name=args_list[0])
        return entry_point(args_list[1:])

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
        elif len(df) == 0:
            # Same reason the filter summary goes to stderr: stdout carries
            # the table and nothing else, so an empty result leaves it empty.
            print("No predictions.", file=sys.stderr)
        else:
            write_predictions(df, sys.stdout)
    except BrokenPipeError:
        # A downstream consumer such as ``head`` intentionally closed the
        # pipe. This is successful early termination, not malformed usage.
        _silence_closed_stdout()
        return 0
    except CLI_ERROR_TYPES as error:
        arg_parser.error(cli_error_message(error))
