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

import pytest

from mhctools.cli.script import (
    SUBCOMMANDS,
    arg_parser,
    find_subcommand,
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


@pytest.mark.parametrize(
    "subcommand", SUBCOMMANDS, ids=[s.name for s in SUBCOMMANDS])
def test_subcommand_entry_point_resolves(subcommand):
    """Dispatch imports by name, so a typo would only surface at runtime."""
    entry_point = subcommand.load()

    assert callable(entry_point)


@pytest.mark.parametrize(
    "subcommand", SUBCOMMANDS, ids=[s.name for s in SUBCOMMANDS])
def test_find_subcommand_resolves_own_name(subcommand):
    assert find_subcommand(subcommand.name) is subcommand


def test_find_subcommand_resolves_aliases():
    assert find_subcommand("integrations") is find_subcommand("predictors")


def test_find_subcommand_returns_none_for_legacy_flags():
    """The bare prediction command must still fall through to argparse."""
    assert find_subcommand("--mhc-predictor") is None
    assert find_subcommand("--sequence") is None


def test_subcommand_names_are_unique():
    names = []
    for subcommand in SUBCOMMANDS:
        names.append(subcommand.name)
        names.extend(subcommand.aliases)
    assert len(names) == len(set(names))
