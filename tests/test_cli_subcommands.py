# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The subcommand table is the only place `mhctools <command>` is declared.

Dispatch and --help both read SUBCOMMANDS, so these tests pin the property
that made #352 possible: a subcommand that exists but is undiscoverable.
"""

import subprocess
import sys

import pytest

from mhctools.cli.script import (
    SUBCOMMANDS,
    arg_parser,
    main,
)


def test_help_lists_every_subcommand():
    help_text = arg_parser.format_help()
    missing = [
        subcommand.name for subcommand in SUBCOMMANDS
        if subcommand.name not in help_text]
    assert not missing, (
        "subcommand(s) dispatched but absent from --help: %s"
        % ", ".join(missing))


def test_help_describes_every_subcommand():
    help_text = arg_parser.format_help()
    for subcommand in SUBCOMMANDS:
        assert subcommand.help, "%s has no help text" % subcommand.name
        # The listing wraps nothing, so the first words survive verbatim.
        assert " ".join(subcommand.help.split()[:4]) in help_text


def test_bare_invocation_prints_help_and_succeeds(capsys):
    """`mhctools` alone used to exit 2 on a missing --mhc-predictor."""
    exit_code = main([])

    assert exit_code == 0
    captured = capsys.readouterr()
    assert "Subcommands:" in captured.out
    assert "predict-table" in captured.out


@pytest.mark.parametrize("command", [
    name for subcommand in SUBCOMMANDS
    for name in (subcommand.name, *subcommand.aliases)])
def test_subcommand_help_dispatches_to_its_parser(command, capsys):
    """Exercise imports, routing, argument forwarding, and alias help names."""
    with pytest.raises(SystemExit) as exit_info:
        main([command, "--help"])
    assert exit_info.value.code == 0
    captured = capsys.readouterr()
    assert "usage: mhctools %s " % command in captured.out
    assert not captured.err


def test_legacy_prediction_flags_still_produce_output(capsys):
    main([
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--sequence", "SIINFEKL"])
    captured = capsys.readouterr()
    assert "SIINFEKL" in captured.out
    assert "peptide" in captured.out
    assert not captured.err


def test_bare_help_keeps_subcommand_modules_and_heavy_runtimes_lazy():
    # Other tests load subcommands, so check startup in a fresh interpreter.
    result = subprocess.run([
        sys.executable, "-c", """
import sys
from mhctools.cli.script import SUBCOMMANDS, main
assert main([]) == 0
assert not {'torch', 'tensorflow'} & sys.modules.keys()
for subcommand in SUBCOMMANDS:
    assert 'mhctools.cli' + subcommand.module not in sys.modules
"""], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Subcommands:" in result.stdout


def test_subcommand_names_are_unique():
    names = []
    for subcommand in SUBCOMMANDS:
        names.append(subcommand.name)
        names.extend(subcommand.aliases)
    assert len(names) == len(set(names))
