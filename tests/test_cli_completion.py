"""Shell-protocol regressions: completion must never initialize predictors."""

import os
from pathlib import Path
import subprocess
import sys

import pytest

from mhctools.cli.script import SUBCOMMANDS


COMPLETE_SCRIPT = """
import sys
import argcomplete
from mhctools.cli.script import SUBCOMMANDS, main
original = argcomplete.autocomplete

def complete(parser):
    def exit_method(code):
        assert not {'pandas', 'numpy', 'pyarrow', 'torch', 'tensorflow',
                    'pyensembl', 'mhcgnomes', 'biotite', 'matplotlib'} & sys.modules.keys()
        assert 'mhctools.artifacts' not in sys.modules
        assert 'mhctools.integrations' not in sys.modules
        selected = sys.argv[1]
        for command in SUBCOMMANDS:
            if command.module != selected:
                assert 'mhctools.cli' + command.module not in sys.modules
        raise SystemExit(code)
    return original(parser, exit_method=exit_method)

argcomplete.autocomplete = complete
main([])
"""


def complete(tmp_path, line, selected="", point=None):
    output = tmp_path / "completions"
    env = dict(os.environ, _ARGCOMPLETE="1", COMP_LINE=line,
               COMP_POINT=str(len(line) if point is None else point),
               _ARGCOMPLETE_STDOUT_FILENAME=str(output), _ARGCOMPLETE_IFS="\013")
    result = subprocess.run([sys.executable, "-c", COMPLETE_SCRIPT, selected],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
    return {word.rstrip() for word in output.read_text().split("\013")}


def test_completes_commands_aliases_and_prediction_flags(tmp_path):
    choices = complete(tmp_path, "mhctools ")
    assert {"--version", "-V", "--sequence", "--mhc-predictor"} <= choices
    for command in SUBCOMMANDS:
        assert {command.name, *command.aliases} <= choices


@pytest.mark.parametrize("command,option,module", [
    ("ls", "--json", ".artifacts"),
    ("fetch", "--accept-license", ".artifacts"),
    ("predictors", "--check", ".integrations"),
    ("integrations", "--check", ".integrations"),
    ("predict-table", "--predictor-info", ".annotate_table"),
    ("benchmark", "--reference-cleavage", ".benchmark"),
    ("cleavage", "--enzyme-state", ".cleavage"),
    ("cleavenet", "--windows", ".cleavenet"),
    ("mixtcrpred", "--batch-size", ".mixtcrpred"),
    ("vaccine-report", "--max-mhc-windows", ".vaccine_report"),
    ("vaccine-trajectory", "--report", ".vaccine_trajectory"),
    ("vaccine-windows", "--table", ".vaccine_windows"),
    ("vaccine-compare", "--input", ".vaccine_compare"),
])
def test_completes_each_commands_options_without_loading_others(tmp_path, command, option, module):
    assert option in complete(tmp_path, f"mhctools {command} --", module)


def test_completes_predictor_names_and_comma_separated_names(tmp_path):
    for prefix in ("", "random,"):
        choices = complete(tmp_path, f"mhctools --mhc-predictor {prefix}netmhcpan4")
        assert {prefix + "netmhcpan4", prefix + "netmhcpan42-ba"} <= choices


def test_completes_selected_command_choices(tmp_path):
    assert complete(tmp_path, "mhctools integrations --check r", ".integrations") == {
        "runnable", "reproduced"}


def test_completion_respects_cursor_before_remaining_arguments(tmp_path):
    prefix = "mhctools fetch --ver"
    assert complete(tmp_path, prefix + " --data-dir /tmp", ".artifacts", len(prefix)) == {"--version"}


@pytest.mark.parametrize("command,module", [
    ("predict-table --input", ".annotate_table"), ("--input-fasta-file", "")])
def test_completion_handles_quoted_paths(tmp_path, command, module):
    source = tmp_path / "peptide inputs.csv"
    source.touch()
    choices = complete(tmp_path, f'mhctools {command} "{tmp_path}/peptide', module)
    assert any("peptide inputs.csv" in choice for choice in choices)


def test_global_completion_marker_is_in_entry_point_prefix():
    from mhctools.cli import script
    assert b"PYTHON_ARGCOMPLETE_OK" in Path(script.__file__).read_bytes()[:1024]


@pytest.mark.parametrize("command", [None, *[c.name for c in SUBCOMMANDS], "integrations"])
def test_help_imports_no_scientific_packages(command):
    args = [] if command is None else [command, "--help"]
    source = f"""
import sys
from mhctools.cli.script import main
try:
    main({args!r})
except SystemExit as e:
    assert e.code == 0
assert not {{'pandas', 'numpy', 'pyarrow', 'torch', 'tensorflow',
            'pyensembl', 'mhcgnomes', 'biotite', 'matplotlib'}} & sys.modules.keys()
"""
    result = subprocess.run([sys.executable, "-c", source], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


def test_package_exports_remain_discoverable_and_resolve_to_original_objects():
    source = """
import sys
import mhctools
assert set(mhctools.__all__) <= set(dir(mhctools))
assert not {'pandas', 'numpy', 'torch', 'pyensembl'} & sys.modules.keys()
from mhctools import BindingPrediction, NetMHCpan, PeptideInput
from mhctools.binding_prediction import BindingPrediction as original_prediction
from mhctools.netmhc_pan import NetMHCpan as original_factory
from mhctools.peptide_input import PeptideInput as original_input
assert BindingPrediction is original_prediction
assert NetMHCpan is original_factory
assert PeptideInput is original_input
assert mhctools.BindingPrediction is BindingPrediction
assert not hasattr(mhctools, 'nonexistent_export')
"""
    result = subprocess.run([sys.executable, "-c", source], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
