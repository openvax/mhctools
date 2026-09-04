# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

from pathlib import Path
from types import SimpleNamespace

import pytest

from mhctools import artifacts
from mhctools.artifacts import ArtifactStatus, artifact_status, fetch, list_artifacts
from mhctools.cli import artifacts as artifact_cli
from mhctools.cli.script import main


def test_list_includes_native_and_packaged_artifacts():
    statuses = {status.name: status for status in list_artifacts()}
    assert set(statuses) == {
        "calis", "mhcflurry", "mhcflurry-affinity", "pepsickle"}
    assert statuses["calis"].manager == "mhctools package"
    assert statuses["calis"].status == "ready"
    assert statuses["calis"].fetchable is False
    assert statuses["pepsickle"].manager == "pepsickle package"
    assert statuses["mhcflurry"].manager == "mhcflurry"
    assert statuses["mhcflurry"].fetchable is True


def test_alias_resolves_to_presentation_artifact():
    assert artifact_status("mhcflurry-presentation").name == "mhcflurry"


def test_unknown_artifact_error_lists_choices():
    with pytest.raises(ValueError, match="Available: calis, mhcflurry"):
        artifact_status("unknown")


def test_fetch_packaged_artifact_is_noop():
    status = fetch("calis")
    assert status.status == "ready"
    assert status.manager == "mhctools package"
    assert Path(status.path).is_file()


def test_fetch_packaged_artifact_rejects_other_version():
    with pytest.raises(ValueError, match="cannot fetch version never"):
        fetch("calis", version="never")


def test_mhcflurry_status_uses_native_manager(monkeypatch, tmp_path):
    models = tmp_path / "models_class1_presentation" / "models"
    models.mkdir(parents=True)
    fake_downloads = SimpleNamespace(
        get_path=lambda name, test_exists=False: str(tmp_path / name),
        get_current_release=lambda: "test-release",
    )
    monkeypatch.setattr(
        artifacts, "_mhcflurry_downloads", lambda: fake_downloads)

    status = artifacts._mhcflurry_presentation_status()
    assert status.status == "ready"
    assert status.version == "test-release"
    assert status.path == str(models)


def test_fetch_mhcflurry_delegates_and_returns_model_path(monkeypatch, tmp_path):
    download_root = tmp_path / "models_class1_presentation"
    model_path = download_root / "models"
    model_path.mkdir(parents=True)
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if "path" in command:
            return SimpleNamespace(stdout=str(download_root) + "\n")
        return SimpleNamespace(stdout=None)

    monkeypatch.setattr(artifacts.shutil, "which", lambda name: "/bin/mhcflurry-downloads")
    monkeypatch.setattr(artifacts.subprocess, "run", fake_run)

    status = fetch("mhcflurry", version="2.2.0")
    assert status.path == str(model_path)
    assert status.version == "2.2.0"
    assert calls[0][0] == [
        "/bin/mhcflurry-downloads", "fetch", "--release", "2.2.0",
        "models_class1_presentation"]
    assert calls[0][1]["env"][
        "MHCFLURRY_DOWNLOADS_CURRENT_RELEASE"] == "2.2.0"
    assert calls[1][0] == [
        "/bin/mhcflurry-downloads", "path", "models_class1_presentation"]


def test_ls_cli_table(monkeypatch, capsys):
    status = ArtifactStatus(
        name="example",
        status="ready",
        manager="package",
        version="1",
        path="/models/example",
        fetchable=False,
        detail="example",
    )
    monkeypatch.setattr(artifact_cli, "list_artifacts", lambda names: [status])
    result = main(["ls"])
    output = capsys.readouterr().out
    assert result is None
    assert "NAME" in output
    assert "MANAGER" in output
    assert "/models/example" in output


def test_ls_cli_json(monkeypatch, capsys):
    status = ArtifactStatus(
        name="example",
        status="missing",
        manager="mhctools",
        version="v1",
        path="",
        fetchable=True,
        detail="example",
    )
    monkeypatch.setattr(artifact_cli, "list_artifacts", lambda names: [status])
    main(["ls", "--json"])
    output = capsys.readouterr().out
    assert '"name": "example"' in output
    assert '"fetchable": true' in output


def test_fetch_cli_dispatch(monkeypatch, capsys):
    status = ArtifactStatus(
        name="example",
        status="ready",
        manager="upstream",
        version="v1",
        path="/models/example",
        fetchable=True,
        detail="example",
    )
    monkeypatch.setattr(
        artifact_cli, "fetch", lambda name, version=None: status)
    result = main(["fetch", "example", "--version", "v1"])
    assert result is None
    assert "upstream" in capsys.readouterr().out
