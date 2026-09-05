# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from mhctools import artifacts
from mhctools.artifacts import ArtifactStatus, artifact_status, fetch, list_artifacts
from mhctools.cli import artifacts as artifact_cli
from mhctools.cli.script import main


def test_list_includes_native_and_packaged_artifacts():
    statuses = {status.name: status for status in list_artifacts()}
    assert set(statuses) == {
        "bigmhc", "calis", "caphla", "deepimmuno", "deeptap", "eramer",
        "mhcflurry", "mhcflurry-affinity", "mixmhc2pred", "mixmhcpred",
        "mixtcrpred",
        "netchop", "netcleave", "netmhc", "netmhccons", "netmhciipan",
        "netmhcpan", "netmhcstabpan", "nettcr", "pepsickle", "prime",
        "tlimmuno2", "tulip"}
    assert statuses["calis"].manager == "mhctools package"
    assert statuses["calis"].status == "ready"
    assert statuses["calis"].fetchable is False
    assert statuses["pepsickle"].manager == "pepsickle package"
    assert statuses["mhcflurry"].manager == "mhcflurry"
    assert statuses["mhcflurry"].fetchable is True
    assert statuses["deeptap"].manager == "mhctools"
    assert statuses["mixmhcpred"].manager == "manual"
    assert statuses["mixmhcpred"].fetchable is False
    assert "identity-bound" in statuses["netmhcpan"].detail
    assert "--accept-license" in statuses["netmhcpan"].detail
    assert statuses["mixtcrpred"].manager == "mhctools"
    assert statuses["mixtcrpred"].fetchable is True


def test_alias_resolves_to_presentation_artifact():
    assert artifact_status("mhcflurry-presentation").name == "mhcflurry"


def test_unknown_artifact_error_lists_choices():
    with pytest.raises(ValueError, match="Available: bigmhc, calis"):
        artifact_status("unknown")


def test_fetch_packaged_artifact_is_noop():
    status = fetch("calis")
    assert status.status == "ready"
    assert status.manager == "mhctools package"
    assert Path(status.path).is_file()


def test_fetch_packaged_artifact_rejects_other_version():
    with pytest.raises(ValueError, match="cannot fetch version never"):
        fetch("calis", version="never")


def test_data_path_precedence(monkeypatch, tmp_path):
    configured = tmp_path / "configured"
    explicit = tmp_path / "explicit"
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(configured))
    assert artifacts.data_path() == configured
    assert artifacts.data_path(explicit) == explicit


def test_managed_status_reports_destination_and_pinned_version(tmp_path):
    status = artifact_status("eramer", data_dir=tmp_path)
    snapshot = artifacts._SNAPSHOTS["eramer"]
    assert status.status == "missing"
    assert status.manager == "mhctools"
    assert status.version == snapshot.revision
    assert status.path == str(
        tmp_path / "artifacts" / "eramer" / snapshot.revision)


def test_managed_status_rejects_missing_provenance(tmp_path):
    target = artifacts.managed_path("eramer", data_dir=tmp_path)
    target.mkdir(parents=True)
    (target / "PWM.xlsx").touch()
    status = artifact_status("eramer", data_dir=tmp_path)
    assert status.status == "missing"
    assert "invalid provenance" in status.detail


def test_managed_status_finds_user_install(monkeypatch, tmp_path):
    (tmp_path / "PWM.xlsx").touch()
    monkeypatch.setenv("ERAMER_HOME", str(tmp_path))
    status = artifact_status("eramer")
    assert status.status == "ready"
    assert status.manager == "user"
    assert status.path == str(tmp_path)


def test_eramer_status_finds_direct_pwm(monkeypatch, tmp_path):
    pwm_path = tmp_path / "custom.xlsx"
    pwm_path.touch()
    monkeypatch.setenv("ERAMER_PWM", str(pwm_path))
    status = artifact_status("eramer")
    assert status.status == "ready"
    assert status.manager == "user"
    assert status.path == str(pwm_path)


def test_license_gated_snapshot_requires_acceptance(tmp_path):
    with pytest.raises(RuntimeError, match="--accept-license"):
        fetch("nettcr", data_dir=tmp_path)
    with pytest.raises(RuntimeError, match="--accept-license"):
        fetch("mixtcrpred", data_dir=tmp_path)


def test_snapshot_rejects_untested_revision(tmp_path):
    with pytest.raises(ValueError, match="tested only at revision"):
        fetch("eramer", version="main", data_dir=tmp_path)


def test_fetches_pinned_sparse_snapshot(monkeypatch, tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    subprocess.run(["git", "init", str(source)], check=True, capture_output=True)
    subprocess.run(
        ["git", "-C", str(source), "config", "user.email", "test@example.com"],
        check=True)
    subprocess.run(
        ["git", "-C", str(source), "config", "user.name", "Test"],
        check=True)
    (source / "models").mkdir()
    (source / "models" / "model.bin").write_bytes(b"weights")
    (source / "ignored").mkdir()
    (source / "ignored" / "large.bin").write_bytes(b"ignored")
    subprocess.run(
        ["git", "-C", str(source), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(source), "commit", "-m", "fixture"],
        check=True,
        capture_output=True)
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    snapshot = artifacts._Snapshot(
        repository=str(source),
        revision=revision,
        sparse_paths=("models",),
        required_paths=("models/model.bin",),
        license_name="MIT",
    )
    monkeypatch.setitem(artifacts._SNAPSHOTS, "fixture", snapshot)

    status = fetch("fixture", data_dir=tmp_path / "data")
    path = Path(status.path)
    assert (path / "models" / "model.bin").read_bytes() == b"weights"
    assert not (path / "ignored").exists()
    assert not (path / ".git").exists()
    manifest = json.loads((path / ".mhctools-artifact.json").read_text())
    assert manifest["revision"] == revision
    assert artifact_status("fixture", data_dir=tmp_path / "data").status == "ready"


def test_git_progress_is_sent_to_stderr(monkeypatch):
    calls = []
    monkeypatch.setattr(artifacts.shutil, "which", lambda name: "/bin/git")
    monkeypatch.setattr(
        artifacts.subprocess,
        "run",
        lambda command, **kwargs: calls.append((command, kwargs)))
    artifacts._run_git(("status",))
    assert calls == [(["/bin/git", "status"], {
        "check": True,
        "stdout": sys.stderr,
    })]


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
    monkeypatch.setattr(
        artifact_cli,
        "list_artifacts",
        lambda names, data_dir=None: [status])
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
    monkeypatch.setattr(
        artifact_cli,
        "list_artifacts",
        lambda names, data_dir=None: [status])
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
        artifact_cli,
        "fetch",
        lambda name, version=None, data_dir=None, accept_license=False,
        models=None, all_models=False, high_confidence=False: status)
    result = main([
        "fetch", "example", "--version", "v1", "--data-dir", "/models",
        "--accept-license"])
    assert result is None
    assert "upstream" in capsys.readouterr().out


def test_fetch_cli_passes_mixtcrpred_model_selection(monkeypatch):
    calls = []
    status = ArtifactStatus(
        name="mixtcrpred", status="ready", manager="mhctools",
        version="v1", path="/models", fetchable=True, detail="example")

    def fake_fetch(name, **kwargs):
        calls.append((name, kwargs))
        return status

    monkeypatch.setattr(artifact_cli, "fetch", fake_fetch)
    main([
        "fetch", "mixtcrpred", "--model", "A0201_GILGFVFTL",
        "--accept-license",
    ])
    assert calls[0][0] == "mixtcrpred"
    assert calls[0][1]["models"] == ["A0201_GILGFVFTL"]
    assert calls[0][1]["accept_license"] is True
