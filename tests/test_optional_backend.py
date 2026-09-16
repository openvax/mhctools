# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

import os
from pathlib import Path
import sys

import pytest

from mhctools.optional_backend import (
    BackendSpec,
    backend_inventory,
    common_checkout_paths,
    inspect_artifact,
    probe_executable,
    run_python_sidecar,
    sha256_file,
)


_SPEC = BackendSpec(
    name="fixture",
    endpoint="serum_half_life",
    developed_against="fixture-v1",
    license="MIT",
    serialization="JSON data only",
    entry_point="prediction_only",
    supported_platforms=("linux", "macos"),
    supported_interpreters=("Python 3.9+",),
)


def test_common_checkout_paths_cover_home_and_code(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    assert common_checkout_paths("Tool", "tool-lower") == (
        tmp_path / "Tool",
        tmp_path / "code" / "Tool",
        tmp_path / "tool-lower",
        tmp_path / "code" / "tool-lower",
    )


def test_executable_probe_distinguishes_missing_launch_and_false_zero():
    missing = probe_executable("mhctools-certainly-missing-executable")
    assert not missing.located
    assert not missing.runnable
    assert "not found" in missing.reason

    working = probe_executable(
        sys.executable, args=("-c", "print('ready')"))
    assert working.located
    assert working.runnable

    false_zero = probe_executable(
        sys.executable,
        args=("-c", "print('nested backend failed')"),
        failure_patterns=("backend failed",),
    )
    assert false_zero.located
    assert not false_zero.runnable
    assert "failure marker" in false_zero.reason


def test_executable_probe_bounds_launch_time():
    result = probe_executable(
        sys.executable,
        args=("-c", "import time; time.sleep(2)"),
        timeout=0.1,
    )
    assert not result.runnable
    assert "timed out" in result.reason


def test_inspection_hashes_content_without_loading_it(tmp_path):
    artifact = tmp_path / "model.pickle"
    artifact.write_bytes(b"not actually a pickle")
    expected = sha256_file(artifact)
    identity = inspect_artifact(
        "weights", "model_weights", artifact,
        expected_sha256=expected, serialization="pickle")
    assert identity.status == "verified"
    assert identity.sha256 == expected
    assert identity.size == len(b"not actually a pickle")


def test_missing_mismatched_and_unverified_states_are_distinct(tmp_path):
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"model")
    missing = inspect_artifact("missing", "weights", tmp_path / "missing")
    mismatch = inspect_artifact(
        "mismatch", "weights", artifact, expected_sha256="0" * 64)
    unverified = inspect_artifact("unverified", "weights", artifact)
    assert missing.status == "missing"
    assert mismatch.status == "mismatch"
    assert unverified.status == "unverified"


def test_inventory_identity_ignores_machine_specific_paths(tmp_path):
    first = tmp_path / "first" / "model.bin"
    second = tmp_path / "second" / "model.bin"
    first.parent.mkdir()
    second.parent.mkdir()
    first.write_bytes(b"same")
    second.write_bytes(b"same")
    expected = sha256_file(first)
    left = backend_inventory(_SPEC, [
        inspect_artifact("weights", "model_weights", first, expected)])
    right = backend_inventory(_SPEC, [
        inspect_artifact("weights", "model_weights", second, expected)])
    assert left.identity_sha256 == right.identity_sha256
    assert left.predictor_version == right.predictor_version


def test_verification_metadata_does_not_change_content_identity(tmp_path):
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"same")
    unverified = backend_inventory(_SPEC, [
        inspect_artifact("weights", "model_weights", artifact)])
    verified = backend_inventory(_SPEC, [inspect_artifact(
        "weights", "model_weights", artifact,
        expected_sha256=sha256_file(artifact))])
    assert unverified.status == "unverified"
    assert verified.status == "verified"
    assert unverified.identity_sha256 == verified.identity_sha256


def test_inventory_identity_includes_settings_and_content(tmp_path):
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"one")
    first = backend_inventory(
        _SPEC, [inspect_artifact("weights", "model_weights", artifact)],
        settings={"mode": "one"})
    second = backend_inventory(
        _SPEC, [inspect_artifact("weights", "model_weights", artifact)],
        settings={"mode": "two"})
    artifact.write_bytes(b"two")
    changed = backend_inventory(
        _SPEC, [inspect_artifact("weights", "model_weights", artifact)],
        settings={"mode": "one"})
    assert len({first.identity_sha256, second.identity_sha256,
                changed.identity_sha256}) == 3


def test_capability_distinguishes_location_verification_and_inference(tmp_path):
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"model")
    located = backend_inventory(
        _SPEC, [inspect_artifact("weights", "model_weights", artifact)])
    verified = backend_inventory(_SPEC, [inspect_artifact(
        "weights", "model_weights", artifact,
        expected_sha256=sha256_file(artifact))])
    assert located.capability == "artifacts_located"
    assert verified.capability == "artifacts_verified"
    assert verified.with_inference_reproduced().capability == \
        "inference_reproduced"


def test_unverified_executable_assets_require_explicit_opt_in(tmp_path):
    artifact = tmp_path / "model.pickle"
    artifact.write_bytes(b"untrusted")
    inventory = backend_inventory(_SPEC, [inspect_artifact(
        "weights", "model_weights", artifact, serialization="pickle")])
    with pytest.raises(RuntimeError, match="allow_unverified_assets"):
        inventory.require_usable()
    assert inventory.require_usable(allow_unverified=True) is inventory


def _write_sidecar(path, source):
    Path(path).write_text(source)


def test_sidecar_runs_with_offline_flags_and_no_user_site(tmp_path):
    output = tmp_path / "output.txt"
    sidecar = tmp_path / "sidecar.py"
    _write_sidecar(sidecar, """
import os
import sys
from pathlib import Path
assert os.environ["HF_HUB_OFFLINE"] == "1"
assert os.environ["TRANSFORMERS_OFFLINE"] == "1"
assert os.environ["PYTHONNOUSERSITE"] == "1"
Path(sys.argv[1]).write_text("ok")
""")
    run_python_sidecar(
        "fixture", sys.executable, sidecar, [output],
        environment=os.environ.copy(), timeout=5)
    assert output.read_text() == "ok"


def test_sidecar_executes_exactly_once(tmp_path):
    output = tmp_path / "executions.txt"
    sidecar = tmp_path / "sidecar.py"
    _write_sidecar(sidecar, """
import sys
from pathlib import Path
with Path(sys.argv[1]).open("a") as result:
    result.write("executed\\n")
""")
    run_python_sidecar(
        "fixture", sys.executable, sidecar, [output],
        environment=os.environ.copy(), timeout=5)
    assert output.read_text().splitlines() == ["executed"]


def test_sidecar_blocks_socket_connections(tmp_path):
    sidecar = tmp_path / "sidecar.py"
    _write_sidecar(sidecar, """
import socket
socket.create_connection(("127.0.0.1", 9), timeout=0.1)
""")
    with pytest.raises(RuntimeError, match="Network access disabled"):
        run_python_sidecar(
            "fixture", sys.executable, sidecar, [],
            environment=os.environ.copy(), timeout=5)


def test_sidecar_timeout_is_reported(tmp_path):
    sidecar = tmp_path / "sidecar.py"
    _write_sidecar(sidecar, "import time; time.sleep(5)\n")
    with pytest.raises(RuntimeError, match="timed out after 0.1 seconds"):
        run_python_sidecar(
            "fixture", sys.executable, sidecar, [],
            environment=os.environ.copy(), timeout=0.1)
