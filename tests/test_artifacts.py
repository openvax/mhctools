# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

import io
import json
from importlib.resources import files
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from mhctools import artifacts
from mhctools.artifacts import ArtifactStatus, artifact_status, fetch, list_artifacts
from mhctools.cli import artifacts as artifact_cli
from mhctools.cli.script import main
from mhctools.netcleave import _find_netcleave_dir
from mhctools.optional_backend import common_checkout_paths
from mhctools.tlimmuno2 import _find_tlimmuno2_home


@pytest.fixture
def no_user_installs(monkeypatch, tmp_path):
    """Hide the developer's own checkouts from the inventory.

    The inventory deliberately reports a user-managed install ahead of a
    managed snapshot, so any assertion about an artifact being "missing" or
    "mhctools"-managed otherwise depends on what the machine running the
    tests happens to have in ``~`` or ``~/code``.
    """
    home = tmp_path / "isolated-home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    # data_path() prefers MHCTOOLS_DATA_DIR over anything derived from HOME,
    # so leaving it set would let a developer's real fetched snapshots answer
    # "missing"/"mhctools" assertions despite the patched home.
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(tmp_path / "isolated-data"))
    for snapshot in artifacts._SNAPSHOTS.values():
        if snapshot.environment_variable:
            monkeypatch.delenv(snapshot.environment_variable, raising=False)
    for definition in artifacts._MANUAL_DIRECTORIES.values():
        monkeypatch.delenv(definition["environment_variable"], raising=False)
    for definition in artifacts._MANUAL_EXECUTABLES.values():
        for variable in definition.get("environment_variables", ()):
            monkeypatch.delenv(variable, raising=False)
    monkeypatch.delenv("ERAMER_PWM", raising=False)
    return home


def test_curated_json_files_are_explicit_package_resources():
    resources = files("mhctools.data")
    assert {
        resource.name for resource in resources.iterdir()
        if resource.name.endswith(".json")
    } == {
        "cleavage_reference.json",
        "dpp4_qpisa.json",
        "intracellular_cleavage_reference.json",
        "intracellular_substrate_evidence.json",
        "model_lineage.json",
        "serum_cleavage_reference.json",
    }


def test_list_includes_native_and_packaged_artifacts(no_user_installs):
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


def test_fetch_is_a_no_op_success_for_any_ready_artifact(monkeypatch):
    """Ready is success in every tier, not only for packaged predictors.

    A manually installed tool used to raise even though the request's goal was
    already met, which made ``fetch`` unusable in a provisioning loop. Who
    manages it stays visible through ``manager`` and ``fetchable``.
    """
    status = ArtifactStatus(
        name="netmhcpan",
        status="ready",
        manager="manual",
        version="unknown",
        path="/tools/netMHCpan",
        fetchable=False,
        detail="Install manually",
    )
    monkeypatch.setattr(artifacts, "artifact_status", lambda *args, **kwargs: status)
    assert fetch("netmhcpan") == status


def test_unfetchable_missing_artifact_names_its_resolution_mechanism():
    """One error shape that says what to install and what the wrapper reads."""
    for name, expected in (
            ("prime", "PRIME_EXECUTABLE"),
            ("netchop", "NETCHOP_HOME"),
            ("mixmhc2pred", "MIXMHC2PRED_EXECUTABLE")):
        status = ArtifactStatus(
            name=name, status="missing", manager="manual", version="",
            path="", fetchable=False, detail="Install %s somehow" % name)
        message = artifacts._unfetchable_message(name, status)
        assert message.startswith(
            "%s is not installed and mhctools cannot fetch it:" % name)
        # Case-sensitive identifiers must survive sentence formatting.
        assert expected in message
        assert "once installed." in message


def test_inventory_covers_wrapper_checkout_paths():
    """``ls``/``fetch`` must not call an install missing that a wrapper runs.

    Each wrapper below resolves conventional checkouts through
    ``common_checkout_paths``, which searches ``~/NAME`` and ``~/code/NAME``.
    When the inventory listed only ``~/NAME``, ``mhctools fetch tlimmuno2``
    exited 2 with "not installed" on a machine whose wrapper was happily
    using ``~/code/TLimmuno2``.
    """
    expected = {
        "caphla": ("CapHLA",),
        "deepimmuno": ("DeepImmuno",),
        "deeptap": ("DeepTAP",),
        "eramer": ("ERAMER",),
        "netcleave": ("NetCleave",),
        "nettcr": ("NetTCR-2.2", "nettcr"),
        "tlimmuno2": ("TLimmuno2",),
    }
    for name, directory_names in expected.items():
        wrapper_paths = set(common_checkout_paths(*directory_names))
        if name in artifacts._SNAPSHOTS:
            listed = artifacts._SNAPSHOTS[name].legacy_paths
        else:
            listed = artifacts._MANUAL_DIRECTORIES[name]["legacy_paths"]
        inventory_paths = {Path(path).expanduser() for path in listed}
        assert wrapper_paths <= inventory_paths, name


def test_data_path_precedence(monkeypatch, tmp_path):
    configured = tmp_path / "configured"
    explicit = tmp_path / "explicit"
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(configured))
    assert artifacts.data_path() == configured
    assert artifacts.data_path(explicit) == explicit


def test_managed_status_reports_destination_and_pinned_version(
        no_user_installs, tmp_path):
    status = artifact_status("eramer", data_dir=tmp_path)
    snapshot = artifacts._SNAPSHOTS["eramer"]
    assert status.status == "missing"
    assert status.manager == "mhctools"
    assert status.version == snapshot.revision
    assert status.path == str(
        tmp_path / "artifacts" / "eramer" / snapshot.revision)


def test_managed_status_rejects_missing_provenance(no_user_installs, tmp_path):
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


@pytest.fixture(params=[
    ("netcleave", "NetCleave", "NetCleave.py", _find_netcleave_dir),
    ("tlimmuno2", "TLimmuno2", "Python/TLimmuno2.py", _find_tlimmuno2_home),
])
def snapshot_discovery(request, no_user_installs):
    name, directory, entrypoint, resolver = request.param
    snapshot = artifacts._SNAPSHOTS[name]
    target = artifacts.managed_path(name)
    _populate_snapshot_assets(target, snapshot)
    (target / ".mhctools-artifact.json").write_text(json.dumps({
        "name": name,
        "repository": snapshot.repository,
        "revision": snapshot.revision,
    }))
    return SimpleNamespace(
        name=name, directory=directory, entrypoint=entrypoint,
        resolver=resolver, snapshot=snapshot, target=target,
        home=no_user_installs,
    )


def _populate_snapshot_assets(root, snapshot):
    for relative in snapshot.required_paths:
        path = root / relative
        if path.suffix:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        else:
            path.mkdir(parents=True, exist_ok=True)


@pytest.mark.parametrize("location", ["", "code"])
@pytest.mark.parametrize("contents", ["empty", "entrypoint-only", "complete"])
def test_automatic_discovery_agrees_with_inventory(
        snapshot_discovery, location, contents):
    install = snapshot_discovery
    checkout = install.home / location / install.directory
    checkout.mkdir(parents=True)
    if contents == "entrypoint-only":
        entrypoint = checkout / install.entrypoint
        entrypoint.parent.mkdir(parents=True, exist_ok=True)
        entrypoint.touch()
    elif contents == "complete":
        _populate_snapshot_assets(checkout, install.snapshot)

    expected = checkout if contents == "complete" else install.target
    status = artifact_status(install.name)
    assert status.status == "ready"
    assert status.manager == ("user" if contents == "complete" else "mhctools")
    assert Path(status.path) == expected.resolve()
    assert Path(install.resolver()).resolve() == expected.resolve()


@pytest.mark.parametrize("override", ["argument", "environment"])
def test_explicit_snapshot_override_is_preserved(
        snapshot_discovery, monkeypatch, override):
    install = snapshot_discovery
    checkout = install.home / "custom"
    entrypoint = checkout / install.entrypoint
    entrypoint.parent.mkdir(parents=True)
    entrypoint.touch()
    if override == "environment":
        monkeypatch.setenv(install.snapshot.environment_variable, str(checkout))
        actual = install.resolver()
    else:
        monkeypatch.setenv(
            install.snapshot.environment_variable, str(install.home / "absent"))
        actual = install.resolver(str(checkout))
    assert Path(actual) == checkout


@pytest.mark.parametrize("override", ["argument", "environment"])
def test_missing_explicit_snapshot_override_does_not_fall_back(
        snapshot_discovery, monkeypatch, override):
    install = snapshot_discovery
    missing = str(install.home / "absent")
    with pytest.raises(FileNotFoundError):
        if override == "environment":
            monkeypatch.setenv(install.snapshot.environment_variable, missing)
            install.resolver()
        else:
            install.resolver(missing)


def test_eramer_status_finds_direct_pwm(monkeypatch, tmp_path):
    pwm_path = tmp_path / "custom.xlsx"
    pwm_path.touch()
    monkeypatch.setenv("ERAMER_PWM", str(pwm_path))
    status = artifact_status("eramer")
    assert status.status == "ready"
    assert status.manager == "user"
    assert status.path == str(pwm_path)


def test_netchop_status_uses_wrapper_installation_resolution(
        monkeypatch, tmp_path):
    executable = tmp_path / "bin" / "netChop"
    executable.parent.mkdir()
    executable.touch()
    monkeypatch.setenv("NETCHOP_HOME", str(tmp_path))
    monkeypatch.delenv("NETMHC_BUNDLE_HOME", raising=False)
    monkeypatch.setattr(artifacts.shutil, "which", lambda name: None)

    status = artifact_status("netchop")

    assert status.status == "ready"
    assert status.path == str(executable)


def test_license_gated_snapshot_requires_acceptance(tmp_path):
    with pytest.raises(RuntimeError, match="--accept-license"):
        fetch("nettcr", data_dir=tmp_path)
    with pytest.raises(RuntimeError, match="--accept-license"):
        fetch("mixtcrpred", data_dir=tmp_path)


def test_unlicensed_snapshot_is_gated_without_claiming_a_license(tmp_path):
    """NetCleave has no license, so the gate must not imply accepting terms."""
    with pytest.raises(RuntimeError) as caught:
        fetch("netcleave", data_dir=tmp_path)
    message = str(caught.value)
    assert "--accept-license" in message
    assert "publishes no license" in message
    assert "cannot grant permission" in message
    # "distributed under the ..." is the wording for a real license and would
    # be a false statement of terms here.
    assert "distributed under" not in message


def test_unlicensed_snapshot_records_absent_license_in_inventory():
    from mhctools.artifacts import artifact_status

    status = artifact_status("netcleave")
    assert status.manager == "mhctools"
    assert status.fetchable is True
    assert "publishes no license" in status.detail


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


def test_mhcflurry_download_progress_is_sent_to_stderr(monkeypatch, tmp_path):
    """``fetch --json`` must emit only JSON, whoever does the downloading.

    MHCflurry's downloader prints a progress table; inheriting stdout put that
    ahead of the JSON document and broke machine-readable consumers.
    """
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, stdout=str(tmp_path))

    monkeypatch.setattr(
        artifacts.shutil, "which", lambda name: "/bin/mhcflurry-downloads")
    monkeypatch.setattr(artifacts.subprocess, "run", fake_run)
    monkeypatch.setattr(artifacts.os.path, "exists", lambda path: True)
    artifacts._fetch_mhcflurry(
        name="mhcflurry",
        download_name="models_class1_presentation",
        relative_path="models",
    )
    fetch_call = next(
        kwargs for command, kwargs in calls if "fetch" in command)
    assert fetch_call["stdout"] is sys.stderr


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


def test_ls_cli_missing_artifact_shows_note_not_destination(monkeypatch, capsys):
    revision = "c7e37a249317704bf96a1e3881a7ece3c3c977a6"
    destination = "/not/yet/installed"
    status = ArtifactStatus(
        name="bigmhc",
        status="missing",
        manager="mhctools",
        version=revision,
        path=destination,
        fetchable=True,
        detail="Run mhctools fetch bigmhc --accept-license",
    )
    monkeypatch.setattr(
        artifact_cli,
        "list_artifacts",
        lambda names, data_dir=None: [status])
    main(["ls"])
    output = capsys.readouterr().out
    assert revision[:12] in output
    assert revision not in output
    assert destination not in output
    assert "bigmhc: Run mhctools fetch bigmhc --accept-license" in output


def test_ls_cli_json_keeps_full_revision_and_destination(monkeypatch, capsys):
    revision = "c7e37a249317704bf96a1e3881a7ece3c3c977a6"
    destination = "/not/yet/installed"
    status = ArtifactStatus(
        name="bigmhc",
        status="missing",
        manager="mhctools",
        version=revision,
        path=destination,
        fetchable=True,
        detail="Fetch it",
    )
    monkeypatch.setattr(
        artifact_cli,
        "list_artifacts",
        lambda names, data_dir=None: [status])
    main(["ls", "--json"])
    output = capsys.readouterr().out
    assert revision in output
    assert destination in output


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


def test_fetch_mhcflurry_is_a_no_op_when_the_models_are_present(monkeypatch):
    """The native tier has to honour "ready is success" like the others.

    Re-running otherwise restarted MHCflurry's downloader, and failed with
    exit 2 whenever ``mhcflurry-downloads`` was off PATH -- a venv invoked by
    absolute path, cron, an IDE -- even though ``ls`` called it ready.
    """
    status = ArtifactStatus(
        name="mhcflurry", status="ready", manager="mhcflurry",
        version="2.2.0", path="/models", fetchable=True,
        detail="Managed by MHCflurry")
    monkeypatch.setattr(
        artifacts, "artifact_status", lambda *args, **kwargs: status)

    def fail(*args, **kwargs):
        raise AssertionError("a ready artifact must not run the downloader")

    monkeypatch.setattr(artifacts.subprocess, "run", fail)
    assert fetch("mhcflurry") == status


def test_progress_forwarding_survives_a_replaced_stderr(monkeypatch, capsys):
    """``sys.stderr`` need not own a file descriptor.

    Passing it straight to the child raised ``io.UnsupportedOperation:
    fileno`` under ``contextlib.redirect_stderr``, pytest's capture, or any
    host that replaces the stream, which broke the Python fetch API.
    """
    monkeypatch.setattr(sys, "stderr", io.StringIO())
    completed = artifacts._run_showing_progress(
        [sys.executable, "-c", "print('downloading')"], check=True)
    assert completed.returncode == 0
    assert "downloading" in sys.stderr.getvalue()


def test_unfetchable_hint_is_a_sentence_without_a_leading_or():
    """An artifact with no environment variable still reads as English."""
    status = ArtifactStatus(
        name="netmhcpan", status="missing", manager="manual", version="",
        path="", fetchable=False, detail="Install NetMHCpan from DTU")
    message = artifacts._unfetchable_message("netmhcpan", status)
    assert "Put netMHCpan on PATH once installed." in message
    assert "Or put" not in message


def test_mhcflurry_downloader_failure_becomes_a_clean_cli_error(
        monkeypatch, tmp_path):
    """A failed download exits 2 with an ``error:`` line, not a traceback."""
    def fake_run(command, **kwargs):
        raise subprocess.CalledProcessError(
            1, command, stderr="no space left on device")

    monkeypatch.setattr(
        artifacts.shutil, "which", lambda name: "/bin/mhcflurry-downloads")
    monkeypatch.setattr(artifacts.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="no space left on device"):
        artifacts._fetch_mhcflurry(
            name="mhcflurry",
            download_name="models_class1_presentation",
            relative_path="models",
            version="2.2.0",
        )


def _stage_losing_fetch(monkeypatch, tmp_path, winner_contents):
    """Drive _fetch_snapshot to the rename with the destination populated.

    git is stubbed out and the staged checkout is built by hand, so the test
    exercises mhctools' own race recovery rather than the network.
    """
    snapshot = artifacts._SNAPSHOTS["eramer"]
    target = artifacts.managed_path("eramer", data_dir=tmp_path)

    def fake_run_git(arguments):
        # The final git call is the checkout. Populate the staged directory
        # and, at the same moment, let the "winner" appear at the
        # destination: _fetch_snapshot returns early if the target already
        # exists when it starts, so the race can only be reproduced by the
        # destination being created mid-fetch.
        if "checkout" in arguments:
            checkout = Path(arguments[arguments.index("-C") + 1])
            checkout.mkdir(parents=True, exist_ok=True)
            (checkout / "PWM.xlsx").write_text("loser", encoding="utf-8")
            (checkout / ".git").mkdir(exist_ok=True)
            target.mkdir(parents=True, exist_ok=True)
            for relative, text in winner_contents.items():
                (target / relative).write_text(text, encoding="utf-8")

    monkeypatch.setattr(artifacts, "_run_git", fake_run_git)
    monkeypatch.setattr(
        artifacts.shutil, "which", lambda name: "/usr/bin/git")
    monkeypatch.setattr(
        artifacts.subprocess, "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, stdout=snapshot.revision))
    return snapshot, target


def test_losing_a_fetch_race_returns_the_winners_snapshot(
        monkeypatch, tmp_path):
    """Recovery has to run on a real ENOTEMPTY, not on FileExistsError.

    Renaming onto a populated directory never raises FileExistsError, so the
    original ``except FileExistsError`` guard was unreachable and the loser
    of a race got a traceback instead of the completed snapshot.
    """
    eramer = artifacts._SNAPSHOTS["eramer"]
    snapshot, target = _stage_losing_fetch(
        monkeypatch, tmp_path,
        {"PWM.xlsx": "winner", ".mhctools-artifact.json": json.dumps({
            "name": "eramer",
            "repository": eramer.repository,
            "revision": eramer.revision,
        })})

    status = artifacts._fetch_snapshot("eramer", data_dir=tmp_path)

    assert status.status == "ready"
    assert status.path == str(target)
    # The winner's content survived; the loser did not overwrite it.
    assert (target / "PWM.xlsx").read_text(encoding="utf-8") == "winner"


def test_losing_a_fetch_race_to_an_invalid_snapshot_still_raises(
        monkeypatch, tmp_path):
    """Recovery must not paper over a destination that is not usable."""
    _stage_losing_fetch(monkeypatch, tmp_path, {"PWM.xlsx": "winner"})
    with pytest.raises(OSError):
        artifacts._fetch_snapshot("eramer", data_dir=tmp_path)


def test_tlimmuno2_is_a_fetchable_unlicensed_snapshot():
    """No published license is NetCleave's situation, not netMHC's.

    Refusing to fetch grants no rights either way, so TLimmuno2 is acquired
    on an explicit acknowledgement rather than left as manual-only.
    """
    snapshot = artifacts._SNAPSHOTS["tlimmuno2"]
    assert snapshot.unlicensed is True
    assert snapshot.acceptance_required is True
    assert snapshot.environment_variable == "TLIMMUNO2_HOME"
    # Cone mode materializes every root-level file, and this repository keeps
    # an 88 MB .RData session dump beside its code.
    assert snapshot.sparse_cone is False
    assert snapshot.sparse_paths == ("/Python/",)
    assert "tlimmuno2" not in artifacts._MANUAL_DIRECTORIES


def test_unlicensed_snapshots_require_acknowledgement_not_acceptance(tmp_path):
    for name in ("netcleave", "tlimmuno2"):
        with pytest.raises(RuntimeError, match="publishes no license"):
            fetch(name, data_dir=tmp_path)


def test_fetchable_wrappers_all_expose_a_fetch_classmethod():
    """Every snapshot-backed wrapper offers the Python fetch() shortcut.

    Asserting this for one wrapper let the invariant lapse twice: NetCleave
    had no fetch() until it was added here, and TLimmuno2 arrived as a
    snapshot without one. The mapping is explicit so a new snapshot entry
    fails until its wrapper is named.
    """
    import importlib

    wrappers = {
        "bigmhc": ("bigmhc", "BigMHC"),
        "caphla": ("caphla", "CapHLA"),
        "deepimmuno": ("deepimmuno", "DeepImmuno"),
        "deeptap": ("deeptap", "DeepTAP"),
        "eramer": ("eramer", "ERAMER"),
        "netcleave": ("netcleave", "NetCleave"),
        "nettcr": ("nettcr", "NetTCR"),
        "mixtcrpred": ("mixtcrpred", "MixTCRpred"),
        "tlimmuno2": ("tlimmuno2", "TLimmuno2"),
        "tulip": ("tulip", "Tulip"),
    }
    assert set(wrappers) == set(artifacts._SNAPSHOTS), (
        "a snapshot was added or removed without naming its wrapper here")
    for name, (module_name, class_name) in wrappers.items():
        module = importlib.import_module("mhctools.%s" % module_name)
        predictor = getattr(module, class_name)
        assert callable(getattr(predictor, "fetch", None)), name


def test_progress_fallback_preserves_the_childs_explanation(monkeypatch):
    """A failure under a replaced stderr must not discard the reason.

    The fallback branch captured only stdout, so CalledProcessError.stderr
    was None and the wrapped message lost the actual cause (disk full,
    network error) on exactly the path the fallback exists to serve.
    """
    monkeypatch.setattr(sys, "stderr", io.StringIO())
    with pytest.raises(subprocess.CalledProcessError) as raised:
        artifacts._run_showing_progress(
            [sys.executable, "-c",
             "import sys; sys.stderr.write('no space left'); sys.exit(1)"],
            check=True)
    assert "no space left" in artifacts._process_failure_detail(raised.value)


def test_mhcflurry_path_failure_is_not_reported_as_a_failed_download(
        monkeypatch):
    """Naming the wrong step sends the user debugging the wrong thing."""
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        if "path" in command:
            raise subprocess.CalledProcessError(
                1, command, stderr="unknown download")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(
        artifacts.shutil, "which", lambda name: "/bin/mhcflurry-downloads")
    monkeypatch.setattr(artifacts.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="could not report its path"):
        artifacts._fetch_mhcflurry(
            name="mhcflurry",
            download_name="models_class1_presentation",
            relative_path="models",
            version="2.2.0",
        )
