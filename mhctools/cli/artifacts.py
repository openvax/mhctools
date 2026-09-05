# Copyright (c) 2026 Mount Sinai School of Medicine
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

"""CLI commands for discovering and fetching predictor artifacts."""

from argparse import ArgumentParser
import json

from ..artifacts import fetch, list_artifacts


def _print_table(statuses):
    headers = ("NAME", "STATUS", "MANAGER", "VERSION", "FETCHABLE", "PATH")
    rows = [
        (
            status.name,
            status.status,
            status.manager,
            status.version or "-",
            "yes" if status.fetchable else "no",
            status.path or "-",
        )
        for status in statuses
    ]
    widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows))
        for i in range(len(headers))
    ]
    print("  ".join(value.ljust(widths[i]) for i, value in enumerate(headers)))
    for row in rows:
        print("  ".join(value.ljust(widths[i]) for i, value in enumerate(row)))


def ls_main(args_list=None):
    """Run ``mhctools ls``."""
    parser = ArgumentParser(
        prog="mhctools ls",
        description="List model weights and other predictor artifacts.",
    )
    parser.add_argument("name", nargs="*", help="Optional artifact names")
    parser.add_argument(
        "--data-dir",
        help="Override MHCTOOLS_DATA_DIR for mhctools-managed artifacts")
    parser.add_argument(
        "--json", action="store_true", help="Emit machine-readable JSON")
    parser.add_argument(
        "--models", action="store_true",
        help="List the per-pMHC model catalog (mixtcrpred only)")
    parser.add_argument(
        "--downloaded", action="store_true",
        help="With --models, show only locally available weights")
    parser.add_argument(
        "--high-confidence", action="store_true",
        help="With --models, show only upstream high-confidence models")
    args = parser.parse_args(args_list)
    if args.models:
        if args.name != ["mixtcrpred"]:
            parser.error("--models requires exactly: mhctools ls mixtcrpred")
        from ..mixtcrpred import print_model_catalog
        try:
            return print_model_catalog(
                data_dir=args.data_dir,
                downloaded=args.downloaded,
                high_confidence=args.high_confidence,
                json_output=args.json,
            )
        except (OSError, RuntimeError, ValueError) as error:
            parser.error(str(error))
    if args.downloaded or args.high_confidence:
        parser.error("--downloaded/--high-confidence require --models")
    try:
        statuses = list_artifacts(args.name or None, data_dir=args.data_dir)
    except ValueError as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps([status.to_dict() for status in statuses], indent=2))
    else:
        _print_table(statuses)


def fetch_main(args_list=None):
    """Run ``mhctools fetch``."""
    parser = ArgumentParser(
        prog="mhctools fetch",
        description="Fetch artifacts required by a predictor wrapper.",
    )
    parser.add_argument("name", help="Artifact or predictor name")
    parser.add_argument(
        "--version", help="Specific upstream artifact release to fetch")
    parser.add_argument(
        "--data-dir",
        help="Override MHCTOOLS_DATA_DIR for mhctools-managed artifacts")
    parser.add_argument(
        "--accept-license",
        action="store_true",
        help="Confirm acceptance when an upstream license requires it")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--model", dest="models", action="append", metavar="NAME",
        help="Fetch one model weight; repeat for multiple models")
    selection.add_argument(
        "--all-models", action="store_true",
        help="Fetch every available model weight (mixtcrpred only)")
    selection.add_argument(
        "--high-confidence", action="store_true",
        help="Fetch upstream's high-confidence weights (mixtcrpred only)")
    parser.add_argument(
        "--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args(args_list)
    try:
        status = fetch(
            args.name,
            version=args.version,
            data_dir=args.data_dir,
            accept_license=args.accept_license,
            models=args.models,
            all_models=args.all_models,
            high_confidence=args.high_confidence,
        )
    except (RuntimeError, ValueError) as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps(status.to_dict(), indent=2))
    else:
        _print_table([status])
