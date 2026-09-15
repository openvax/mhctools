"""Regression tests for concise, actionable command-line failures."""

from pathlib import Path

import pandas as pd
import pytest

from mhctools.artifacts import ArtifactStatus
from mhctools.cli.annotate_table import main as predict_table_main
from mhctools.cli.errors import cli_error_message, predictor_error_message


def _run_predict_table(capsys, args):
    with pytest.raises(SystemExit) as raised:
        predict_table_main(args)
    assert raised.value.code == 2
    stderr = capsys.readouterr().err
    assert "Traceback" not in stderr
    return stderr


def test_predict_table_missing_file_has_no_traceback(capsys, tmp_path):
    missing = tmp_path / "missing.csv"
    stderr = _run_predict_table(capsys, [
        "--input", str(missing),
        "--out", str(tmp_path / "out.csv"),
        "--predictor", "random",
    ])
    assert str(missing) in stderr


def test_predict_table_key_error_is_not_repr_quoted(capsys, tmp_path):
    source = tmp_path / "input.csv"
    pd.DataFrame({"peptide": ["SIINFEKL"]}).to_csv(source, index=False)
    stderr = _run_predict_table(capsys, [
        "--input", str(source),
        "--out", str(tmp_path / "out.csv"),
        "--peptide-column", "pep",
        "--predictor", "random",
    ])
    assert "peptide column 'pep' not found" in stderr
    assert '"peptide column' not in stderr


def test_predict_table_unknown_predictor_has_no_traceback(capsys, tmp_path):
    source = tmp_path / "input.csv"
    pd.DataFrame({"peptide": ["SIINFEKL"]}).to_csv(source, index=False)
    stderr = _run_predict_table(capsys, [
        "--input", str(source),
        "--out", str(tmp_path / "out.csv"),
        "--predictor", "nosuchmodel",
    ])
    assert "Unknown predictor 'nosuchmodel'" in stderr


def test_predict_table_collision_names_cli_overwrite_flag(capsys, tmp_path):
    source = tmp_path / "input.csv"
    pd.DataFrame({
        "peptide": ["SIINFEKL"],
        "random_affinity": [1.0],
    }).to_csv(source, index=False)
    stderr = _run_predict_table(capsys, [
        "--input", str(source),
        "--out", str(tmp_path / "out.csv"),
        "--predictor", "random",
    ])
    assert "pass --overwrite" in stderr
    assert "overwrite=True" not in stderr


def test_missing_backend_reuses_artifact_install_guidance(monkeypatch):
    detail = "Install mhctools[pepsickle] to obtain the packaged weights"
    status = ArtifactStatus(
        name="pepsickle",
        status="unavailable",
        manager="pepsickle package",
        version=None,
        path=None,
        fetchable=False,
        detail=detail,
    )
    monkeypatch.setattr(
        "mhctools.artifacts.artifact_status", lambda name: status)
    error = ModuleNotFoundError("No module named 'pepsickle'")
    assert predictor_error_message(error, ["pepsickle::processing"]) == detail


def test_cli_license_message_omits_library_keyword_argument():
    message = cli_error_message(RuntimeError(
        "Review the license, then rerun with --accept-license "
        "(or accept_license=True)."))
    assert message == "Review the license, then rerun with --accept-license."


def test_no_partial_output_is_created_after_input_failure(capsys, tmp_path):
    output = tmp_path / "out.csv"
    _run_predict_table(capsys, [
        "--input", str(tmp_path / "missing.csv"),
        "--out", str(output),
        "--predictor", "random",
    ])
    assert not Path(output).exists()

