# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Discover and fetch model artifacts used by mhctools predictors.

Predictors that already manage their own downloads retain ownership of their
cache.  mhctools provides a common entry point and reports the resolved path;
it does not copy those artifacts into a second cache.
"""

from dataclasses import asdict, dataclass
from importlib import metadata, util
import os
from pathlib import Path
import shutil
import subprocess


@dataclass(frozen=True)
class ArtifactStatus:
    """Status of the artifacts required by one predictor wrapper."""

    name: str
    status: str
    manager: str
    version: str
    path: str
    fetchable: bool
    detail: str

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)


def _mhcflurry_downloads():
    """Import MHCflurry's lightweight download configuration lazily."""
    from mhcflurry import downloads
    return downloads


def _mhcflurry_status(name, download_name, relative_path):
    try:
        downloads = _mhcflurry_downloads()
        root = downloads.get_path(download_name, test_exists=False)
        path = os.path.join(root, relative_path)
        version = downloads.get_current_release() or "custom"
        ready = os.path.exists(path)
        detail = "Managed by MHCflurry"
    except (ImportError, ModuleNotFoundError):
        path = ""
        version = ""
        ready = False
        detail = "Install the mhcflurry package before fetching its models"
    return ArtifactStatus(
        name=name,
        status="ready" if ready else "missing",
        manager="mhcflurry",
        version=str(version),
        path=os.path.abspath(path) if path else "",
        fetchable=True,
        detail=detail,
    )


def _mhcflurry_presentation_status():
    return _mhcflurry_status(
        name="mhcflurry",
        download_name="models_class1_presentation",
        relative_path="models",
    )


def _mhcflurry_affinity_status():
    return _mhcflurry_status(
        name="mhcflurry-affinity",
        download_name="models_class1_pan",
        relative_path="models.combined",
    )


def _calis_status():
    from . import __version__
    from . import calis
    return ArtifactStatus(
        name="calis",
        status="ready",
        manager="mhctools package",
        version=__version__,
        path=os.path.abspath(calis.__file__),
        fetchable=False,
        detail="Published model parameters are embedded in mhctools",
    )


def _pepsickle_status():
    spec = util.find_spec("pepsickle")
    if spec is None or spec.origin is None:
        return ArtifactStatus(
            name="pepsickle",
            status="missing",
            manager="pepsickle package",
            version="",
            path="",
            fetchable=False,
            detail="Install mhctools[pepsickle] to obtain the packaged weights",
        )

    package_dir = Path(spec.origin).resolve().parent
    required = (
        package_dir / "model.joblib",
        package_dir / "trained_model_dict.pickle",
    )
    ready = all(path.is_file() for path in required)
    try:
        version = metadata.version("pepsickle")
    except metadata.PackageNotFoundError:
        version = ""
    return ArtifactStatus(
        name="pepsickle",
        status="ready" if ready else "missing",
        manager="pepsickle package",
        version=version,
        path=str(package_dir),
        fetchable=False,
        detail=(
            "Model weights are installed with the pepsickle package"
            if ready else
            "The pepsickle package is present but its model weights are missing"
        ),
    )


_STATUS_FUNCTIONS = {
    "calis": _calis_status,
    "mhcflurry": _mhcflurry_presentation_status,
    "mhcflurry-affinity": _mhcflurry_affinity_status,
    "pepsickle": _pepsickle_status,
}

_ALIASES = {
    "mhcflurry-presentation": "mhcflurry",
}


def _canonical_name(name):
    normalized = name.strip().lower()
    return _ALIASES.get(normalized, normalized)


def artifact_status(name):
    """Return the current status for a named predictor's artifacts."""
    canonical = _canonical_name(name)
    try:
        status_function = _STATUS_FUNCTIONS[canonical]
    except KeyError:
        raise ValueError(
            "Unknown artifact %r. Available: %s" % (
                name, ", ".join(sorted(_STATUS_FUNCTIONS))))
    return status_function()


def list_artifacts(names=None):
    """List known packaged, native, and optionally downloadable artifacts."""
    if names is None:
        canonical_names = sorted(_STATUS_FUNCTIONS)
    else:
        canonical_names = []
        for name in names:
            canonical = _canonical_name(name)
            if canonical not in canonical_names:
                canonical_names.append(canonical)
    return [artifact_status(name) for name in canonical_names]


def _fetch_mhcflurry(name, download_name, relative_path, version=None):
    executable = shutil.which("mhcflurry-downloads")
    if executable is None:
        raise RuntimeError(
            "The mhcflurry-downloads command is unavailable. Reinstall the "
            "mhcflurry package in this Python environment.")

    environment = os.environ.copy()
    command = [executable, "fetch"]
    if version is not None:
        # MHCflurry selects its storage directory at import time using this
        # environment variable. Passing --release alone can otherwise fetch a
        # historical release into the current release's directory.
        environment["MHCFLURRY_DOWNLOADS_CURRENT_RELEASE"] = str(version)
        command.extend(["--release", str(version)])
    command.append(download_name)
    subprocess.run(command, check=True, env=environment)

    path_result = subprocess.run(
        [executable, "path", download_name],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    root = path_result.stdout.strip()
    path = os.path.abspath(os.path.join(root, relative_path))
    if not os.path.exists(path):
        raise RuntimeError(
            "MHCflurry reported a completed download, but the expected path "
            "does not exist: %s" % path)

    resolved_version = str(version) if version is not None else ""
    if not resolved_version:
        resolved_version = str(_mhcflurry_downloads().get_current_release() or "custom")
    return ArtifactStatus(
        name=name,
        status="ready",
        manager="mhcflurry",
        version=resolved_version,
        path=path,
        fetchable=True,
        detail="Managed by MHCflurry",
    )


def fetch(name, version=None):
    """Fetch a predictor's external artifacts and return their status.

    Native download managers retain ownership of their files. Predictors whose
    parameters are already packaged are treated as successful no-ops.
    """
    canonical = _canonical_name(name)
    if canonical == "mhcflurry":
        return _fetch_mhcflurry(
            name="mhcflurry",
            download_name="models_class1_presentation",
            relative_path="models",
            version=version,
        )
    if canonical == "mhcflurry-affinity":
        return _fetch_mhcflurry(
            name="mhcflurry-affinity",
            download_name="models_class1_pan",
            relative_path="models.combined",
            version=version,
        )

    status = artifact_status(canonical)
    if status.status == "ready":
        if version is not None and version != status.version:
            raise ValueError(
                "%s is managed by %s at version %s; mhctools cannot fetch "
                "version %s" % (
                    canonical, status.manager, status.version, version))
        return status
    raise RuntimeError(status.detail)
