# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

import json
from types import SimpleNamespace

from mhctools.artifacts import ArtifactStatus
from mhctools.cli import integrations as integration_cli
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
        integration_cli, "list_integrations", lambda *args, **kwargs: statuses)

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
        integration_cli, "list_integrations", lambda *args, **kwargs: statuses)

    result = main([
        "integrations", "example", "--check", "reproduced", "--strict",
    ])
    output = capsys.readouterr().out

    assert result == 0
    assert "REPRODUCED" in output
    assert "reference reproduced" in output
