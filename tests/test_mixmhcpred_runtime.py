"""Interpreter integrity, per-process isolation and runtime failure handling."""

import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from mhctools.mixmhcpred_runtime import MixMHCpredRuntime, resolve_runtime
from mhctools.cli.args import make_mhc_arg_parser, predictors_from_args


def diagnostics(pandas_version="2.3.3", executable=None):
    return dict(reported_executable=executable or sys.executable,
                python_version="3.11.12", pandas_version=pandas_version,
                numpy_version="2.2.6", scipy_version="1.15.3",
                logomaker_version="0.8.7", matplotlib_version="3.10.3")


def test_unconfigured_runtime_retains_existing_launchers(monkeypatch):
    monkeypatch.delenv("MIXMHCPRED_PYTHON", raising=False)
    assert resolve_runtime() is None


def test_missing_configured_runtime_is_not_replaced_by_host(monkeypatch):
    monkeypatch.setenv("MIXMHCPRED_PYTHON", "/absent/venv/bin/python")
    with pytest.raises(FileNotFoundError, match="separate venv.*MIXMHCPRED_PYTHON"):
        resolve_runtime()


def test_explicit_selection_preserves_venv_symlink_and_beats_environment(tmp_path, monkeypatch):
    interpreter = tmp_path / "venv/bin/python"
    interpreter.parent.mkdir(parents=True)
    interpreter.symlink_to(sys.executable)
    monkeypatch.setenv("MIXMHCPRED_PYTHON", "/invalid")
    calls = []

    def probe(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout=json.dumps(diagnostics(executable=str(interpreter))))

    monkeypatch.setattr(subprocess, "run", probe)
    runtime = resolve_runtime(interpreter)
    assert calls[0][0][0] == str(interpreter)
    assert runtime.python_executable == str(interpreter) != str(interpreter.resolve())
    assert runtime.reported_executable == str(interpreter)
    assert calls[0][1]["timeout"] == 30


@pytest.mark.parametrize("output", ["invalid JSON", "{}", json.dumps(diagnostics("not-a-version"))])
def test_invalid_runtime_diagnostics_raise_actionable_error(output, monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs:
                        SimpleNamespace(returncode=0, stdout=output))
    with pytest.raises(RuntimeError, match="Invalid.*separate venv"):
        resolve_runtime(sys.executable)


def test_pandas3_backend_is_rejected_before_inference(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs:
                        SimpleNamespace(returncode=0, stdout=json.dumps(diagnostics("3.0.3"))))
    with pytest.raises(RuntimeError, match="backend pandas<3.*pandas 3.0.3.*separate venv"):
        resolve_runtime(sys.executable)


def test_missing_backend_dependency_has_install_guidance(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs:
                        SimpleNamespace(returncode=1, stdout="", stderr="No module named 'logomaker'"))
    with pytest.raises(RuntimeError, match="dependencies failed.*logomaker.*separate venv"):
        resolve_runtime(sys.executable)


def test_probe_timeout_has_bounded_actionable_failure(monkeypatch):
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs["timeout"])
    monkeypatch.setattr(subprocess, "run", timeout)
    with pytest.raises(RuntimeError, match="Could not probe.*timed out.*MIXMHCPRED_PYTHON"):
        resolve_runtime(sys.executable, timeout=0.1)


def test_child_only_shim_quotes_path_sanitizes_imports_and_cleans_up(tmp_path, monkeypatch):
    interpreter = tmp_path / "runtime with spaces/py'thon$literal"
    interpreter.parent.mkdir()
    interpreter.symlink_to(sys.executable)
    runtime = MixMHCpredRuntime(python_executable=str(interpreter), **diagnostics())
    monkeypatch.setenv("PYTHONPATH", "/invalid/host/imports")
    monkeypatch.setenv("PYTHONHOME", "/invalid/host/home")
    before = dict(os.environ)
    with pytest.raises(RuntimeError, match="forced failure"):
        with runtime.environment() as environment:
            launcher = Path(environment["PATH"].split(os.pathsep)[0]) / "python3"
            assert launcher.exists()
            assert "PYTHONPATH" not in environment and "PYTHONHOME" not in environment
            output = subprocess.check_output(
                [str(launcher), "-c", "import json,sys;print(json.dumps(dict(executable=sys.executable)))"],
                env=environment, text=True)
            # macOS Python 3.11 can report the base executable through a
            # symlink outside its venv, even though our quoted path executed.
            assert os.path.samefile(json.loads(output)["executable"], sys.executable)
            assert os.environ == before
            raise RuntimeError("forced failure")
    assert not launcher.parent.exists()
    assert os.environ == before


@pytest.mark.parametrize("predictor", ["mixmhcpred", "prime"])
def test_cli_routes_configured_interpreter_to_a_real_runtime_probe(predictor, tmp_path, monkeypatch):
    parser = make_mhc_arg_parser()
    executable = tmp_path / "MixMHCpred"
    executable.write_text("#!/bin/sh\n")
    executable.chmod(0o755)
    args = parser.parse_args(["--mhc-predictor", predictor, "--mhc-alleles", "HLA-A*02:01",
                              "--mhc-predictor-path", str(executable),
                              "--mixmhcpred-python", "/missing/backend/python"])
    model, = predictors_from_args(args)
    if predictor == "prime":
        monkeypatch.setattr(model, "_validate_mixmhcpred", lambda: ("/unused/MixMHCpred", "3.0"))
    with pytest.raises(FileNotFoundError, match="backend Python was not found"):
        model.predict(["SIINFEKL"])


def test_cli_does_not_silently_ignore_runtime_for_other_predictors():
    args = make_mhc_arg_parser().parse_args([
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--mixmhcpred-python", "/unused"])
    with pytest.raises(ValueError, match="requires MixMHCpred or PRIME"):
        predictors_from_args(args)
