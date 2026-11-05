"""Terminal styling must leave redirected help and data unmodified."""

import io
import re

import pytest

from mhctools.cli.formatting import format_table
from mhctools.cli.help import HelpArgumentParser, color_enabled


class Terminal(io.StringIO):
    def isatty(self):
        return True


def parser():
    result = HelpArgumentParser(prog="mhctools", description="Plain description.")
    result.add_argument("--example", help="Plain explanation.")
    return result


def test_only_printed_terminal_help_has_color(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    p = parser()
    assert "\033" not in p.format_help()
    terminal = Terminal()
    p.print_help(terminal)
    styled = terminal.getvalue()
    assert "\033[36m--example\033[0m" in styled
    assert "\033[1;36moptions:" in styled
    assert "\033[1musage:\033[0m" in styled
    assert re.sub(r"\033\[[0-9;]*m", "", styled) == p.format_help()
    output = io.StringIO()
    p.print_help(output)
    assert output.getvalue() == p.format_help()


@pytest.mark.parametrize("name,value", [("NO_COLOR", ""), ("NO_COLOR", "1"), ("TERM", "dumb")])
def test_terminal_color_opt_out(monkeypatch, name, value):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm")
    monkeypatch.setenv(name, value)
    stream = Terminal()
    assert not color_enabled(stream)
    parser().print_help(stream)
    assert "\033" not in stream.getvalue()
    assert "\033" not in format_table(("NAME",), [("example",)], stream=stream)


def test_table_styling_preserves_values_and_alignment(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm")
    headers = ("NAME", "STATUS", "PATH")
    rows = [("a", "ready", "/models"), ("longer", "missing", "-")]
    plain = format_table(headers, rows)
    colored = format_table(headers, rows, stream=Terminal())
    assert colored.startswith("\033[1;36mNAME")
    assert "\033[33mmissing" in colored
    assert "\033[2m-" in colored
    assert "\033[32m" not in colored
    assert "\n".join(line.rstrip() for line in re.sub(r"\033\[[0-9;]*m", "", colored).splitlines()) == plain
    assert format_table(headers, rows, stream=io.StringIO()) == plain
