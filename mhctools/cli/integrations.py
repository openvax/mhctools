# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""CLI for honest optional-integration capability reports."""

from .help import HelpArgumentParser as ArgumentParser
import json
import sys

from ..capabilities import CAPABILITY_LEVELS
from .formatting import format_table


def _display(value):
    if value is None:
        return "not checked"
    return "yes" if value else "no"


def _print_table(statuses):
    headers = ("NAME", "LOCATED", "RUNNABLE", "REPRODUCED", "CAPABILITY", "PATH")
    rows = [
        (
            status.name,
            _display(status.located),
            _display(status.runnable),
            _display(status.reproduced),
            status.capability,
            status.path or "-",
        )
        for status in statuses
    ]
    print(format_table(headers, rows, stream=sys.stdout))
    notes = [status for status in statuses if status.detail]
    if notes:
        print()
        for status in notes:
            print("%s: %s" % (status.name, status.detail))


def make_parser(command_name="predictors"):
    """Build the command parser without running predictions or probing assets."""
    parser = ArgumentParser(
        prog="mhctools %s" % command_name,
        description=(
            "Report artifact location, bounded launchability, and reference "
            "inference reproduction as separate capabilities."),
    )
    parser.add_argument("name", nargs="*", help="Optional integration names")
    parser.add_argument(
        "--check", choices=CAPABILITY_LEVELS, default="runnable",
        help="Highest capability to attempt (default: runnable)")
    parser.add_argument(
        "--strict", action="store_true",
        help="Exit nonzero unless every selected integration meets --check")
    parser.add_argument(
        "--timeout", type=float, default=10,
        help="Seconds allowed for each bounded launch probe (default: 10)")
    parser.add_argument(
        "--data-dir",
        help="Override MHCTOOLS_DATA_DIR for mhctools-managed artifacts")
    parser.add_argument(
        "--json", action="store_true", help="Emit machine-readable JSON")
    return parser


def predictors_main(args_list=None, command_name="predictors"):
    """Run the predictor capability-report command."""
    parser = make_parser(command_name=command_name)
    args = parser.parse_args(args_list)

    from ..integrations import list_integrations
    try:
        statuses = list_integrations(
            args.name or None,
            check=args.check,
            data_dir=args.data_dir,
            timeout=args.timeout,
        )
    except ValueError as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps([status.to_dict() for status in statuses], indent=2))
    else:
        _print_table(statuses)
    if args.strict and not all(status.meets(args.check) for status in statuses):
        return 1
    return 0
