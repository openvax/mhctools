# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

import json
import subprocess
import sys
from types import SimpleNamespace

import pytest

from mhctools.artifacts import ArtifactStatus
from mhctools.cli.script import main
from mhctools.integrations import IntegrationStatus, integration_status
from mhctools import integrations


def _artifact(name="example", ready=True):
    return ArtifactStatus(
        name=name,
        status="ready" if ready else "missing",
        manager="test",
        version="1",
        path="/tools/%s" % name if ready else "",
        fetchable=False,
        detail="artifact detail",
    )


def test_located_check_does_not_claim_runtime_capability(monkeypatch):
    monkeypatch.setattr(
        "mhctools.integrations.artifact_status", lambda *args, **kwargs: _artifact())

    status = integration_status("example", check="located")

    assert status.located is True
    assert status.runnable is None
    assert status.reproduced is None
    assert status.capability == "located"
    assert status.meets("located") is True
    assert status.meets("runnable") is False


def test_missing_artifact_blocks_requested_capabilities(monkeypatch):
    monkeypatch.setattr(
        "mhctools.integrations.artifact_status",
        lambda *args, **kwargs: _artifact(ready=False))

    status = integration_status("example", check="reproduced")

    assert status.located is False
    assert status.runnable is False
    assert status.reproduced is False
    assert status.capability == "missing"


def test_unregistered_probe_remains_unknown(monkeypatch):
    monkeypatch.setattr(
        "mhctools.integrations.artifact_status", lambda *args, **kwargs: _artifact())

    status = integration_status("example", check="runnable")

    assert status.located is True
    assert status.runnable is None
    assert status.capability == "located"
    assert "not been established" in status.detail


def test_executable_probe_rejects_false_zero_nested_launcher(monkeypatch):
    monkeypatch.setattr(
        integrations,
        "probe_executable",
        lambda *args, **kwargs: SimpleNamespace(
            runnable=False,
            reason="reported failure marker 'no binaries found'",
        ),
    )
    monkeypatch.setattr(
        integrations, "artifact_status",
        lambda *args, **kwargs: _artifact(name="netmhcpan"),
    )

    status = integration_status("netmhcpan", check="runnable")

    assert status.located is True
    assert status.runnable is False
    assert status.capability == "located"
    assert "no binaries found" in status.detail


def test_calis_reference_probe_is_reproduced():
    status = integration_status("calis", check="reproduced")

    assert status.located is True
    assert status.runnable is True
    assert status.reproduced is True
    assert status.capability == "reproduced"
    assert "0.30484" in status.detail


def test_integrations_cli_json_and_strict_failure(monkeypatch, capsys):
    statuses = [IntegrationStatus(
        name="example",
        located=True,
        runnable=True,
        reproduced=None,
        capability="runnable",
        path="/tools/example",
        detail="launch succeeded",
    )]
    monkeypatch.setattr(
        integrations, "list_integrations", lambda *args, **kwargs: statuses)

    result = main([
        "integrations", "example", "--check", "reproduced", "--strict",
        "--json",
    ])
    payload = json.loads(capsys.readouterr().out)

    assert result == 1
    assert payload[0]["located"] is True
    assert payload[0]["runnable"] is True
    assert payload[0]["reproduced"] is None


def test_integrations_cli_strict_success(monkeypatch, capsys):
    statuses = [IntegrationStatus(
        name="example",
        located=True,
        runnable=True,
        reproduced=True,
        capability="reproduced",
        path="/tools/example",
        detail="reference reproduced",
    )]
    monkeypatch.setattr(
        integrations, "list_integrations", lambda *args, **kwargs: statuses)

    result = main([
        "integrations", "example", "--check", "reproduced", "--strict",
    ])
    output = capsys.readouterr().out

    assert result == 0
    assert "REPRODUCED" in output
    assert "reference reproduced" in output


def test_predictors_is_canonical_command_and_integrations_remains_alias(
        monkeypatch, capsys):
    statuses = [IntegrationStatus(
        name="example",
        located=True,
        runnable=True,
        reproduced=None,
        capability="runnable",
        path="/tools/example",
        detail="launch succeeded",
    )]
    monkeypatch.setattr(
        integrations, "list_integrations", lambda *args, **kwargs: statuses)

    assert main(["predictors", "--check", "runnable", "--strict"]) == 0
    canonical_output = capsys.readouterr().out
    assert main(["integrations", "--check", "runnable", "--strict"]) == 0
    alias_output = capsys.readouterr().out

    assert canonical_output == alias_output


@pytest.mark.parametrize("name", ["deepimmuno", "deeptap", "tlimmuno2", "mixtcrpred"])
@pytest.mark.parametrize("returncode", [0, 1])
def test_sidecar_probe_uses_configured_interpreter_and_environment(
        name, returncode, tmp_path, monkeypatch):
    python = tmp_path / "sidecar-python"
    python.symlink_to(sys.executable)
    variable = name.upper() + "_PYTHON"
    monkeypatch.setenv(variable, str(python))
    monkeypatch.delenv("TF_USE_LEGACY_KERAS", raising=False)
    monkeypatch.setattr(integrations, "artifact_status", lambda *a, **kw: _artifact(name))
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, returncode)

    monkeypatch.setattr(subprocess, "run", run)
    status = integration_status(name, check="reproduced", timeout=7)
    assert len(calls) == 1
    command, options = calls[0]
    assert command[0] == str(python)
    assert options["timeout"] == 7
    if name in ("deepimmuno", "tlimmuno2"):
        assert options["env"]["TF_USE_LEGACY_KERAS"] == "1"
        assert "tf_keras" in command[2]
    elif name == "mixtcrpred":
        assert options["env"]["PYTHONNOUSERSITE"] == "1"
        assert "torchvision" in command[2]
    else:
        assert "pytorch_lightning" in command[2]
    assert status.runnable is (returncode == 0)
    assert status.reproduced is None
    assert str(python) in status.detail


@pytest.mark.parametrize("name", ["deepimmuno", "deeptap", "tlimmuno2", "mixtcrpred"])
def test_invalid_sidecar_is_reported_without_falling_back(name, monkeypatch):
    monkeypatch.setenv(name.upper() + "_PYTHON", "/missing/mhctools-sidecar-python")
    monkeypatch.setattr(integrations, "artifact_status", lambda *a, **kw: _artifact(name))
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("wrong runtime launched"))
    status = integration_status(name)
    assert status.runnable is False
    assert name.upper() + "_PYTHON" in status.detail
    assert "/missing/mhctools-sidecar-python" in status.detail


def test_runtime_timeout_is_reported(monkeypatch):
    monkeypatch.setenv("DEEPTAP_PYTHON", sys.executable)
    monkeypatch.setattr(integrations, "artifact_status", lambda *a, **kw: _artifact("deeptap"))

    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", timeout)
    status = integration_status("deeptap", timeout=0.5)
    assert status.runnable is False
    assert "timed out" in status.detail
    assert status.reproduced is None


@pytest.mark.parametrize("broken", [False, True])
def test_pepsickle_probes_the_interpreters_prediction_selects(tmp_path, monkeypatch, broken):
    from mhctools import Pepsickle
    normal, gb = [tmp_path / name for name in ("normal-python", "gb-python")]
    for path in (normal, gb):
        path.symlink_to(sys.executable)
    monkeypatch.setenv("PEPSICKLE_PYTHON", str(normal))
    monkeypatch.setenv("PEPSICKLE_GB_PYTHON", str(gb))
    normal_predictor = Pepsickle()
    gb_predictor = Pepsickle(model_type="in-vitro", proteasome_type="C", human_only=False)
    monkeypatch.setattr(integrations, "artifact_status", lambda *a, **kw: _artifact("pepsickle"))
    calls = []

    def run(command, **kwargs):
        calls.append((command[0], kwargs["timeout"]))
        return subprocess.CompletedProcess(
            command, int(broken and command[0] == str(gb)), "", "missing package")

    monkeypatch.setattr(subprocess, "run", run)
    status = integration_status("pepsickle", timeout=3)
    assert [path for path, timeout in calls] == [
        normal_predictor.python_executable, gb_predictor.python_executable]
    assert all(0 < timeout <= 3 for path, timeout in calls)
    assert status.runnable is (not broken)
    assert status.reproduced is None
    assert str(gb) in status.detail


def test_pepsickle_shared_interpreter_is_probed_once(monkeypatch):
    monkeypatch.setenv("PEPSICKLE_PYTHON", sys.executable)
    monkeypatch.setenv("PEPSICKLE_GB_PYTHON", sys.executable)
    monkeypatch.setattr(integrations, "artifact_status", lambda *a, **kw: _artifact("pepsickle"))
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(subprocess, "run", run)
    assert integration_status("pepsickle").runnable is True
    assert len(calls) == 1


def test_every_artifact_has_a_probe_or_explicit_gap():
    names = {status.name for status in integrations.list_artifacts()}
    registered = set(integrations._EXECUTABLE_ARGS) | set(integrations._SIDECAR_PROBES) | {
        "weirdo", "calis", "netchop", "mhcflurry", "mhcflurry-affinity", "pepsickle"}
    assert not registered.intersection(integrations._NO_RUNTIME_PROBE)
    assert names == registered | set(integrations._NO_RUNTIME_PROBE)
