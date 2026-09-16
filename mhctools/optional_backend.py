# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Artifact identity and isolated execution for optional predictor backends."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Mapping, Optional, Sequence, Tuple


ARTIFACT_STATUS_VALUES = ("missing", "mismatch", "unverified", "verified")
INFERENCE_STATUS_VALUES = ("not_run", "reproduced")


def common_checkout_paths(*directory_names):
    """Return conventional user checkout locations in deterministic order.

    Optional backends historically checked either ``~/NAME`` or
    ``~/code/NAME`` (and often only one of the two). Keeping that policy here
    gives wrappers and their integration tests one shared definition without
    searching the filesystem broadly or guessing outside the user's home.
    """
    home = Path.home()
    paths = []
    for directory_name in directory_names:
        paths.extend((home / directory_name, home / "code" / directory_name))
    return tuple(paths)


@dataclass(frozen=True)
class ExecutableCapability:
    """Bounded launch probe for an optional command-line backend."""

    program: str
    path: str
    runnable: bool
    returncode: Optional[int]
    reason: str

    @property
    def located(self):
        """Whether the executable was resolved on this machine."""
        return bool(self.path)

    def to_dict(self):
        """Return a JSON-serializable launch-capability record."""
        return asdict(self)


def probe_executable(
        program, args=("--help",), timeout=10, failure_patterns=()):
    """Check that an optional executable can launch, with an honest reason.

    Some legacy launchers incorrectly exit zero after a nested command fails,
    so callers may supply diagnostic *failure_patterns* that invalidate such
    output. This is deliberately a launch probe, not proof of inference.
    """
    path = shutil.which(str(program))
    if not path:
        return ExecutableCapability(
            program=str(program), path="", runnable=False, returncode=None,
            reason="%s was not found on PATH" % program)
    command = [path] + [str(arg) for arg in args]
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return ExecutableCapability(
            program=str(program), path=path, runnable=False, returncode=None,
            reason="%s launch probe timed out after %s seconds" % (
                program, timeout))
    output = "\n".join(
        part.strip() for part in (completed.stdout, completed.stderr)
        if part and part.strip())
    matched = next(
        (pattern for pattern in failure_patterns if pattern in output), None)
    if completed.returncode != 0 or matched:
        detail = output[-500:] if output else "no diagnostic output"
        cause = (
            "reported failure marker %r" % matched if matched
            else "exited with code %d" % completed.returncode)
        return ExecutableCapability(
            program=str(program), path=path, runnable=False,
            returncode=completed.returncode,
            reason="%s %s: %s" % (program, cause, detail))
    return ExecutableCapability(
        program=str(program), path=path, runnable=True,
        returncode=completed.returncode, reason="launch probe succeeded")


@dataclass(frozen=True)
class BackendSpec:
    """Reviewed integration contract for one optional prediction endpoint."""

    name: str
    endpoint: str
    developed_against: str
    license: str
    serialization: str
    entry_point: str
    supported_platforms: Tuple[str, ...]
    supported_interpreters: Tuple[str, ...]
    network_policy: str = "offline"

    def __post_init__(self):
        if self.network_policy != "offline":
            raise ValueError("Optional prediction backends must run offline")
        if self.entry_point != "prediction_only":
            raise ValueError(
                "Optional backends require a noninteractive prediction-only entry point")


@dataclass(frozen=True)
class ArtifactIdentity:
    """Content identity and verification state of one runtime file."""

    name: str
    role: str
    path: str
    sha256: str
    size: int
    expected_sha256: str
    status: str
    serialization: str = "data"

    def __post_init__(self):
        if self.status not in ARTIFACT_STATUS_VALUES:
            raise ValueError(
                "Unknown artifact status %r; expected one of %s"
                % (self.status, ARTIFACT_STATUS_VALUES))

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class BackendInventory:
    """Resolved files, settings, and runtime status for an optional backend."""

    spec: BackendSpec
    artifacts: Tuple[ArtifactIdentity, ...]
    settings: Tuple[Tuple[str, str], ...] = ()
    inference_status: str = "not_run"

    def __post_init__(self):
        if self.inference_status not in INFERENCE_STATUS_VALUES:
            raise ValueError(
                "Unknown inference status %r; expected one of %s"
                % (self.inference_status, INFERENCE_STATUS_VALUES))
        names = [artifact.name for artifact in self.artifacts]
        if len(names) != len(set(names)):
            raise ValueError("Artifact names must be unique: %s" % names)

    @property
    def status(self):
        """Worst artifact state, ordered from missing through verified."""
        states = {artifact.status for artifact in self.artifacts}
        for status in ARTIFACT_STATUS_VALUES:
            if status in states:
                return status
        return "missing"

    @property
    def capability(self):
        """Honest distinction between location, verification, and inference."""
        if self.status in ("missing", "mismatch"):
            return "blocked"
        if self.inference_status == "reproduced":
            return "inference_reproduced"
        if self.status == "verified":
            return "artifacts_verified"
        return "artifacts_located"

    @property
    def identity_sha256(self):
        """Stable identity of content and settings, independent of local paths."""
        payload = {
            "backend": self.spec.name,
            "endpoint": self.spec.endpoint,
            "artifacts": [
                {
                    "name": artifact.name,
                    "role": artifact.role,
                    "sha256": artifact.sha256,
                    "size": artifact.size,
                    "serialization": artifact.serialization,
                }
                for artifact in sorted(self.artifacts, key=lambda item: item.name)
            ],
            "settings": dict(self.settings),
        }
        encoded = json.dumps(
            payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    @property
    def predictor_version(self):
        """Version string carrying both compatibility and loaded-asset identity."""
        return (
            "developed-against=%s;assets-sha256=%s;status=%s"
            % (self.spec.developed_against, self.identity_sha256, self.status)
        )

    def with_inference_reproduced(self):
        """Return a copy recording that an offline prediction completed."""
        return replace(self, inference_status="reproduced")

    def require_usable(self, allow_unverified=False):
        """Raise unless required artifacts are present and acceptably verified."""
        if self.status == "missing":
            missing = [
                artifact.name for artifact in self.artifacts
                if artifact.status == "missing"]
            raise FileNotFoundError(
                "%s is missing required artifacts: %s"
                % (self.spec.name, ", ".join(missing)))
        if self.status in ("mismatch", "unverified") and not allow_unverified:
            affected = [
                artifact.name for artifact in self.artifacts
                if artifact.status in ("mismatch", "unverified")]
            raise RuntimeError(
                "%s has %s artifacts: %s. Review their provenance and pass "
                "allow_unverified_assets=True to use them explicitly; hashes do "
                "not make executable model serialization safe."
                % (self.spec.name, self.status, ", ".join(affected)))
        return self

    def to_dict(self):
        """Return a JSON-serializable capability and provenance record."""
        return {
            "spec": asdict(self.spec),
            "artifacts": [artifact.to_dict() for artifact in self.artifacts],
            "settings": dict(self.settings),
            "artifact_status": self.status,
            "inference_status": self.inference_status,
            "capability": self.capability,
            "identity_sha256": self.identity_sha256,
            "predictor_version": self.predictor_version,
        }


def sha256_file(path, chunk_size=1024 * 1024):
    """Return the SHA-256 digest of a file without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as input_file:
        for chunk in iter(lambda: input_file.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_artifact(
        name, role, path, expected_sha256="", serialization="data"):
    """Resolve and hash one file without importing or deserializing it."""
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        return ArtifactIdentity(
            name=name,
            role=role,
            path=str(resolved),
            sha256="",
            size=0,
            expected_sha256=expected_sha256,
            status="missing",
            serialization=serialization,
        )
    actual = sha256_file(resolved)
    if expected_sha256:
        status = "verified" if actual == expected_sha256 else "mismatch"
    else:
        status = "unverified"
    return ArtifactIdentity(
        name=name,
        role=role,
        path=str(resolved),
        sha256=actual,
        size=resolved.stat().st_size,
        expected_sha256=expected_sha256,
        status=status,
        serialization=serialization,
    )


def backend_inventory(spec, artifacts, settings=None):
    """Build a deterministic inventory from resolved artifact records."""
    return BackendInventory(
        spec=spec,
        artifacts=tuple(artifacts),
        settings=tuple(sorted(
            (str(key), str(value)) for key, value in (settings or {}).items())),
    )


def prediction_cache_key(peptide_input, inventory):
    """Deterministic identity of exact input, model assets, and settings.

    Source occurrence metadata is deliberately absent from
    ``peptide_input.inference_identity_sha256``. Duplicate occurrences can
    reuse inference while remaining distinct records in returned results.
    """
    from .peptide_input import PeptideInput
    if not isinstance(peptide_input, PeptideInput):
        raise TypeError("prediction_cache_key requires PeptideInput")
    if not isinstance(inventory, BackendInventory):
        raise TypeError("prediction_cache_key requires BackendInventory")
    payload = {
        "schema": "mhctools.optional_prediction_cache.v1",
        "input_identity_sha256": peptide_input.inference_identity_sha256,
        "backend_identity_sha256": inventory.identity_sha256,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


_OFFLINE_BOOTSTRAP = """
import runpy
import socket
import sys

def _network_disabled(*args, **kwargs):
    raise RuntimeError("Network access disabled for optional backend inference")

socket.socket.connect = _network_disabled
socket.socket.connect_ex = _network_disabled
socket.create_connection = _network_disabled
socket.getaddrinfo = _network_disabled
sys.argv = sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
""".strip()


def run_python_sidecar(
        backend_name: str,
        python: str,
        sidecar,
        args: Sequence[str],
        cwd=None,
        timeout: Optional[float] = None,
        environment: Optional[Mapping[str, str]] = None):
    """Run a prediction-only Python sidecar with network access disabled."""
    env = dict(os.environ if environment is None else environment)
    env.update({
        "HF_DATASETS_OFFLINE": "1",
        "HF_HUB_OFFLINE": "1",
        "PYTHONNOUSERSITE": "1",
        "TRANSFORMERS_OFFLINE": "1",
        "WANDB_MODE": "offline",
    })
    command = [
        python,
        "-c",
        _OFFLINE_BOOTSTRAP,
        str(Path(sidecar).resolve()),
    ] + [str(value) for value in args]
    try:
        process = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            "%s inference timed out after %s seconds"
            % (backend_name, timeout)) from error
    if process.returncode != 0:
        raise RuntimeError(
            "%s inference failed (exit %d):\n%s"
            % (backend_name, process.returncode,
               (process.stderr or process.stdout).strip()))
    return process
