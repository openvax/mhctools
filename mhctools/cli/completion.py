"""Complete the legacy prediction form and lazily dispatched subcommands."""

import os

from .help import HelpArgumentParser


def autocomplete(prediction_parser, subcommands):
    """Complete using the same parsers as execution, loading only one command."""
    if "_ARGCOMPLETE" not in os.environ:
        return

    import argcomplete
    from argcomplete.lexers import split_line

    # Use argcomplete's lexer so quoted paths and a cursor inside the line work.
    line = os.environ.get("COMP_LINE", "")
    point = int(os.environ.get("COMP_POINT", len(line)))
    words = split_line(line, point)[3]
    offset = int(os.environ["_ARGCOMPLETE"])
    command_name = words[offset] if len(words) > offset else None
    selected_command = next((command for command in subcommands
                             if command_name in (command.name, *command.aliases)), None)
    # A subcommand has its own arguments; prediction requirements do not apply.
    parser = HelpArgumentParser(
        parents=[prediction_parser] if selected_command is None else [], add_help=False)
    commands = parser.add_subparsers()
    for command in subcommands:
        selected = command is selected_command
        commands.add_parser(
            command.name, aliases=list(command.aliases), help=command.help,
            add_help=False, parents=[command.make_parser()] if selected else [])
    argcomplete.autocomplete(parser)
