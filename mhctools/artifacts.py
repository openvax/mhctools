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
cache. For upstream repositories that have no download manager, mhctools can
install a pinned, minimal git snapshot in its data directory. Other tools are
inventory-only because their licenses or distribution mechanisms require a
manual installation.
"""

from dataclasses import asdict, dataclass
from importlib import metadata, util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from platformdirs import user_data_path


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


@dataclass(frozen=True)
class _Snapshot:
    """A tested upstream git snapshot that mhctools can acquire."""

    repository: str
    revision: str
    required_paths: tuple
    sparse_paths: tuple = ()
    sparse_cone: bool = True
    required_globs: tuple = ()
    license_name: str = ""
    license_url: str = ""
    acceptance_required: bool = False
    environment_variable: str = ""
    legacy_paths: tuple = ()


# Revisions are deliberately immutable. Updating one is a reviewed mhctools
# change, which gives wrappers a stable upstream layout and reproducible model
# identity instead of silently following a moving default branch.
_SNAPSHOTS = {
    "bigmhc": _Snapshot(
        repository="https://github.com/KarchinLab/bigmhc.git",
        revision="c7e37a249317704bf96a1e3881a7ece3c3c977a6",
        sparse_paths=("data", "models", "src"),
        required_paths=("src/bigmhc.py", "data/pseudoseqs.csv", "models"),
        license_name="BigMHC Academic License",
        license_url=(
            "https://github.com/KarchinLab/bigmhc/blob/"
            "c7e37a249317704bf96a1e3881a7ece3c3c977a6/LICENSE"),
        acceptance_required=True,
        environment_variable="BIGMHC_DIR",
        legacy_paths=("~/bigmhc",),
    ),
    "caphla": _Snapshot(
        repository="https://github.com/changyunjian/CapHLA.git",
        revision="33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831",
        sparse_paths=("params",),
        required_paths=(
            "BA_model.py", "EL_model.py", "HLA_library.csv", "params"),
        required_globs=("params/ba_fold*.params", "params/el_fold*.params"),
        license_name="MIT",
        license_url=(
            "https://github.com/changyunjian/CapHLA/blob/"
            "33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831/LICENSE"),
        environment_variable="CAPHLA_HOME",
        legacy_paths=("~/CapHLA",),
    ),
    "deepimmuno": _Snapshot(
        repository="https://github.com/frankligy/DeepImmuno.git",
        revision="df42ac5b6bddfe531268335e2dcb496559cd488b",
        sparse_paths=("data", "models"),
        required_paths=("deepimmuno-cnn.py", "data", "models"),
        license_name="MIT",
        license_url=(
            "https://github.com/frankligy/DeepImmuno/blob/"
            "df42ac5b6bddfe531268335e2dcb496559cd488b/LICENSE"),
        environment_variable="DEEPIMMUNO_HOME",
        legacy_paths=("~/DeepImmuno",),
    ),
    "deeptap": _Snapshot(
        repository="https://github.com/zjupgx/DeepTAP.git",
        revision="d2dad5bdea146ecc245304f509e97ae8137f94fd",
        sparse_paths=("data", "model"),
        required_paths=("deeptap.py", "data", "model"),
        license_name="Apache-2.0",
        license_url=(
            "https://github.com/zjupgx/DeepTAP/blob/"
            "d2dad5bdea146ecc245304f509e97ae8137f94fd/LICENSE"),
        environment_variable="DEEPTAP_HOME",
        legacy_paths=("~/DeepTAP",),
    ),
    "eramer": _Snapshot(
        repository="https://github.com/aalokaily/ERAMER.git",
        revision="7745d5cf72d99bcda1c73f26ca746d025b46b7f3",
        required_paths=("PWM.xlsx",),
        license_name="GPL-3.0",
        license_url=(
            "https://github.com/aalokaily/ERAMER/blob/"
            "7745d5cf72d99bcda1c73f26ca746d025b46b7f3/LICENSE"),
        environment_variable="ERAMER_HOME",
        legacy_paths=("~/ERAMER",),
    ),
    "nettcr": _Snapshot(
        repository="https://github.com/mnielLab/NetTCR-2.2.git",
        revision="7cead3fe6dcb539ff8e2d9121586dafca1e059c2",
        sparse_paths=(
            "/README.md",
            "/academic_software_license_agreement.pdf",
            "/models/nettcr_2_2_pan/checkpoint/*.tflite",
        ),
        sparse_cone=False,
        required_paths=("models/nettcr_2_2_pan/checkpoint",),
        required_globs=("models/nettcr_2_2_pan/checkpoint/*.tflite",),
        license_name="NetTCR Academic Software License Agreement",
        license_url=(
            "https://github.com/mnielLab/NetTCR-2.2/blob/"
            "7cead3fe6dcb539ff8e2d9121586dafca1e059c2/"
            "academic_software_license_agreement.pdf"),
        acceptance_required=True,
        environment_variable="NETTCR_DIR",
        legacy_paths=("~/NetTCR-2.2", "~/code/NetTCR-2.2"),
    ),
    "mixtcrpred": _Snapshot(
        repository="https://github.com/GfellerLab/MixTCRpred.git",
        revision="acd6f57444bde675840890207c74ca3b0c7ffac2",
        sparse_paths=("src", "pretrained_models"),
        required_paths=(
            "LICENSE.md",
            "MixTCRpred.py",
            "src/dataloaders.py",
            "src/models.py",
            "src/utils.py",
            "pretrained_models/anchors_perc_rank.pickle",
            "pretrained_models/info_models.csv",
            # These two checkpoints ship in the upstream repository. Keep
            # them visible alongside optional Zenodo downloads rather than
            # copying them into a second cache.
            "pretrained_models/model_A0201_ELAGIGILTV.ckpt",
            "pretrained_models/model_A0201_GILGFVFTL.ckpt",
        ),
        license_name="MixTCRpred Academic Non-Commercial License",
        license_url=(
            "https://github.com/GfellerLab/MixTCRpred/blob/"
            "acd6f57444bde675840890207c74ca3b0c7ffac2/LICENSE.md"),
        acceptance_required=True,
        environment_variable="MIXTCRPRED_HOME",
        legacy_paths=("~/MixTCRpred", "~/code/MixTCRpred"),
    ),
    "tulip": _Snapshot(
        repository="https://github.com/barthelemymp/TULIP-TCR.git",
        revision="798fab97a3b13d08dcbfc381ea643e8dc14297c2",
        sparse_paths=(
            "aatok", "configs", "mhctok", "model_weights", "src"),
        required_paths=(
            "predict.py", "aatok", "configs", "mhctok", "model_weights",
            "src"),
        license_name="GPL-3.0",
        license_url=(
            "https://github.com/barthelemymp/TULIP-TCR/blob/"
            "798fab97a3b13d08dcbfc381ea643e8dc14297c2/LICENSE"),
        environment_variable="TULIP_HOME",
        legacy_paths=("~/TULIP-TCR", "~/code/TULIP-TCR"),
    ),
}


# These distributions cannot be fetched safely and unattended. They still
# belong in the inventory so ``mhctools ls`` reports the executable/checkout
# actually used by their wrappers and identifies who manages it.
_MANUAL_EXECUTABLES = {
    "mixmhc2pred": {
        "environment_variables": ("MIXMHC2PRED_EXECUTABLE",),
        "executables": ("MixMHC2pred", "MixMHC2pred_unix"),
        "detail": "Install MixMHC2pred under its upstream license",
    },
    "mixmhcpred": {
        "environment_variables": ("MIXMHCPRED_PATH",),
        "executables": ("MixMHCpred",),
        "detail": "Install MixMHCpred under its academic/non-commercial license",
    },
    "netchop": {
        "executables": ("netChop",),
        "detail": (
            "Install NetChop from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "netmhc": {
        "executables": ("netMHC",),
        "detail": (
            "Install NetMHC from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "netmhccons": {
        "executables": ("netMHCcons",),
        "detail": (
            "Install NetMHCcons from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "netmhciipan": {
        "executables": ("netMHCIIpan",),
        "detail": (
            "Install NetMHCIIpan from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "netmhcpan": {
        "executables": ("netMHCpan",),
        "detail": (
            "Install NetMHCpan from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "netmhcstabpan": {
        "executables": ("netMHCstabpan",),
        "detail": (
            "Install NetMHCstabpan from DTU Health Tech; its identity-bound "
            "request form and emailed private link cannot be replaced by "
            "--accept-license"),
    },
    "prime": {
        "environment_variables": ("PRIME_EXECUTABLE",),
        "executables": ("PRIME",),
        "detail": "Install PRIME under its academic/non-commercial license",
    },
}

_MANUAL_DIRECTORIES = {
    "netcleave": {
        "environment_variable": "NETCLEAVE_DIR",
        "legacy_paths": ("~/NetCleave", "~/code/NetCleave"),
        "required_path": "NetCleave.py",
        "detail": (
            "Install NetCleave manually; its repository has no license file"),
    },
    "tlimmuno2": {
        "environment_variable": "TLIMMUNO2_HOME",
        "legacy_paths": ("~/TLimmuno2",),
        "required_path": "Python/TLimmuno2.py",
        "detail": (
            "Install TLimmuno2 manually; its repository has no license file"),
    },
}


def data_path(data_dir=None):
    """Return the base directory used for mhctools-managed artifacts.

    The explicit argument takes precedence over ``MHCTOOLS_DATA_DIR`` and the
    platform-native user data directory.
    """
    if data_dir is not None:
        return Path(data_dir).expanduser().resolve()
    configured = os.environ.get("MHCTOOLS_DATA_DIR")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(user_data_path("mhctools", appauthor=False)).resolve()


def managed_path(name, data_dir=None):
    """Return the pinned installation path for an mhctools-managed snapshot."""
    canonical = _canonical_name(name)
    try:
        snapshot = _SNAPSHOTS[canonical]
    except KeyError:
        raise ValueError("%s is not an mhctools-managed artifact" % name)
    return data_path(data_dir) / "artifacts" / canonical / snapshot.revision


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
    "tulip-tcr": "tulip",
}


def _canonical_name(name):
    normalized = name.strip().lower()
    return _ALIASES.get(normalized, normalized)


def _snapshot_is_valid(path, snapshot):
    path = Path(path)
    return (
        all((path / relative).exists() for relative in snapshot.required_paths)
        and all(list(path.glob(pattern)) for pattern in snapshot.required_globs)
    )


def _managed_snapshot_is_valid(name, path, snapshot):
    if not _snapshot_is_valid(path, snapshot):
        return False
    try:
        with (Path(path) / ".mhctools-artifact.json").open() as input_file:
            manifest = json.load(input_file)
    except (OSError, ValueError):
        return False
    return (
        manifest.get("name") == name
        and manifest.get("repository") == snapshot.repository
        and manifest.get("revision") == snapshot.revision
    )


def _user_snapshot_path(name, snapshot):
    if name == "eramer":
        pwm_path = os.environ.get("ERAMER_PWM")
        if pwm_path and Path(pwm_path).expanduser().is_file():
            return Path(pwm_path).expanduser().resolve()
    if snapshot.environment_variable:
        configured = os.environ.get(snapshot.environment_variable)
        if configured and _snapshot_is_valid(configured, snapshot):
            return Path(configured).expanduser().resolve()
    for legacy in snapshot.legacy_paths:
        candidate = Path(legacy).expanduser()
        if _snapshot_is_valid(candidate, snapshot):
            return candidate.resolve()
    return None


def _snapshot_status(name, data_dir=None):
    snapshot = _SNAPSHOTS[name]
    user_path = _user_snapshot_path(name, snapshot)
    if user_path is not None:
        return ArtifactStatus(
            name=name,
            status="ready",
            manager="user",
            version="unknown",
            path=str(user_path),
            fetchable=True,
            detail="Using an existing user-managed upstream installation",
        )

    path = managed_path(name, data_dir=data_dir)
    ready = _managed_snapshot_is_valid(name, path, snapshot)
    detail = "Pinned upstream snapshot managed by mhctools"
    if path.exists() and not ready:
        detail = "Managed destination is incomplete or has invalid provenance"
    if snapshot.acceptance_required:
        detail += "; initial fetch requires explicit acceptance of %s (%s)" % (
            snapshot.license_name, snapshot.license_url)
    if name == "mixtcrpred" and ready:
        try:
            from .mixtcrpred import model_catalog
            models = model_catalog(mixtcrpred_path=path)
            downloaded = sum(model.status == "ready" for model in models)
            detail += "; %d/%d pMHC model weights available" % (
                downloaded, len(models))
        except (OSError, RuntimeError, ValueError):
            pass
    return ArtifactStatus(
        name=name,
        status="ready" if ready else "missing",
        manager="mhctools",
        version=snapshot.revision,
        path=str(path),
        fetchable=True,
        detail=detail,
    )


def _manual_executable_status(name):
    definition = _MANUAL_EXECUTABLES[name]
    path = ""
    for variable in definition.get("environment_variables", ()):
        configured = os.environ.get(variable)
        if configured:
            candidate = Path(configured).expanduser()
            if candidate.is_dir() and name == "mixmhcpred":
                candidate = candidate / "MixMHCpred"
            if candidate.is_file():
                path = str(candidate)
                break
    if not path:
        for executable in definition["executables"]:
            path = shutil.which(executable) or ""
            if path:
                break
    return ArtifactStatus(
        name=name,
        status="ready" if path else "missing",
        manager="manual",
        version="unknown" if path else "",
        path=os.path.abspath(path) if path else "",
        fetchable=False,
        detail=definition["detail"],
    )


def _manual_directory_status(name):
    definition = _MANUAL_DIRECTORIES[name]
    candidates = []
    configured = os.environ.get(definition["environment_variable"])
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.extend(Path(path).expanduser() for path in definition["legacy_paths"])
    root = next((
        path.resolve() for path in candidates
        if (path / definition["required_path"]).is_file()), None)
    return ArtifactStatus(
        name=name,
        status="ready" if root else "missing",
        manager="manual",
        version="unknown" if root else "",
        path=str(root) if root else "",
        fetchable=False,
        detail=definition["detail"],
    )


def _all_names():
    return set(_STATUS_FUNCTIONS) | set(_SNAPSHOTS) | set(
        _MANUAL_EXECUTABLES) | set(_MANUAL_DIRECTORIES)


def artifact_status(name, data_dir=None):
    """Return the current status for a named predictor's artifacts."""
    canonical = _canonical_name(name)
    if canonical in _SNAPSHOTS:
        return _snapshot_status(canonical, data_dir=data_dir)
    if canonical in _MANUAL_EXECUTABLES:
        return _manual_executable_status(canonical)
    if canonical in _MANUAL_DIRECTORIES:
        return _manual_directory_status(canonical)
    try:
        status_function = _STATUS_FUNCTIONS[canonical]
    except KeyError:
        raise ValueError(
            "Unknown artifact %r. Available: %s" % (
                name, ", ".join(sorted(_all_names()))))
    return status_function()


def list_artifacts(names=None, data_dir=None):
    """List known packaged, native, managed, and manual artifacts."""
    if names is None:
        canonical_names = sorted(_all_names())
    else:
        canonical_names = []
        for name in names:
            canonical = _canonical_name(name)
            if canonical not in canonical_names:
                canonical_names.append(canonical)
    return [
        artifact_status(name, data_dir=data_dir) for name in canonical_names]


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
        resolved_version = str(
            _mhcflurry_downloads().get_current_release() or "custom")
    return ArtifactStatus(
        name=name,
        status="ready",
        manager="mhcflurry",
        version=resolved_version,
        path=path,
        fetchable=True,
        detail="Managed by MHCflurry",
    )


def _resolve_revision(name, snapshot, version):
    if version is None or version == "default":
        return snapshot.revision
    requested = str(version).lower()
    if snapshot.revision.lower().startswith(requested):
        return snapshot.revision
    raise ValueError(
        "%s is tested only at revision %s, not %s" % (
            name, snapshot.revision, version))


def _run_git(arguments):
    executable = shutil.which("git")
    if executable is None:
        raise RuntimeError("git is required to fetch upstream model artifacts")
    # Keep command progress visible without corrupting ``fetch --json`` stdout.
    subprocess.run(
        [executable] + list(arguments), check=True, stdout=sys.stderr)


def _fetch_snapshot(name, version=None, data_dir=None, accept_license=False):
    snapshot = _SNAPSHOTS[name]
    revision = _resolve_revision(name, snapshot, version)
    target = managed_path(name, data_dir=data_dir)
    if target.exists():
        if _managed_snapshot_is_valid(name, target, snapshot):
            return ArtifactStatus(
                name=name,
                status="ready",
                manager="mhctools",
                version=revision,
                path=str(target),
                fetchable=True,
                detail="Pinned upstream snapshot managed by mhctools",
            )
        raise RuntimeError(
            "Artifact directory exists but is incomplete: %s. Move it aside "
            "and fetch again." % target)

    if snapshot.acceptance_required and not accept_license:
        raise RuntimeError(
            "%s is distributed under the %s. Review %s, then rerun with "
            "--accept-license (or accept_license=True)." % (
                name, snapshot.license_name, snapshot.license_url))

    target.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(tempfile.mkdtemp(
        prefix=".%s-" % name, dir=str(target.parent)))
    checkout = temporary_root / "checkout"
    try:
        _run_git(("init", str(checkout)))
        _run_git(("-C", str(checkout), "remote", "add", "origin",
                  snapshot.repository))
        if snapshot.sparse_paths:
            mode = "--cone" if snapshot.sparse_cone else "--no-cone"
            _run_git((
                "-C", str(checkout), "sparse-checkout", "init", mode))
            _run_git((
                "-C", str(checkout), "sparse-checkout", "set", mode)
                + snapshot.sparse_paths)
        _run_git(("-C", str(checkout), "fetch", "--depth", "1",
                  "--filter=blob:none", "origin", revision))
        _run_git(("-C", str(checkout), "checkout", "--detach", "FETCH_HEAD"))
        executable = shutil.which("git")
        result = subprocess.run(
            [executable, "-C", str(checkout), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        actual_revision = result.stdout.strip()
        if actual_revision != revision:
            raise RuntimeError(
                "Expected revision %s but git checked out %s" % (
                    revision, actual_revision))
        if not _snapshot_is_valid(checkout, snapshot):
            raise RuntimeError(
                "%s revision %s lacks files required by its mhctools wrapper"
                % (name, revision))

        manifest = {
            "license": snapshot.license_name,
            "license_accepted": bool(
                snapshot.acceptance_required and accept_license),
            "license_url": snapshot.license_url,
            "name": name,
            "repository": snapshot.repository,
            "revision": revision,
            "sparse_paths": list(snapshot.sparse_paths),
        }
        with (checkout / ".mhctools-artifact.json").open("w") as output:
            json.dump(manifest, output, indent=2, sort_keys=True)
            output.write("\n")
        # The checked-out content is immutable and identified by the manifest;
        # discard git objects so large model blobs are not stored twice.
        shutil.rmtree(checkout / ".git")
        try:
            checkout.rename(target)
        except FileExistsError:
            if not _managed_snapshot_is_valid(name, target, snapshot):
                raise
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "Could not fetch %s revision %s with git" % (name, revision)) from error
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)

    return ArtifactStatus(
        name=name,
        status="ready",
        manager="mhctools",
        version=revision,
        path=str(target),
        fetchable=True,
        detail="Pinned upstream snapshot managed by mhctools",
    )


def fetch(
        name,
        version=None,
        data_dir=None,
        accept_license=False,
        models=None,
        all_models=False,
        high_confidence=False):
    """Fetch a predictor's external artifacts and return their status.

    Native download managers retain ownership of their files. Tested upstream
    snapshots are installed under :func:`data_path`. Predictors whose
    parameters are already packaged are successful no-ops.
    """
    canonical = _canonical_name(name)
    model_selection = bool(models or all_models or high_confidence)
    if model_selection and canonical != "mixtcrpred":
        raise ValueError(
            "Model selection is supported only for the mixtcrpred artifact")
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
    if canonical in _SNAPSHOTS:
        current = _snapshot_status(canonical, data_dir=data_dir)
        if data_dir is None and version is None and current.status == "ready":
            status = current
        else:
            status = _fetch_snapshot(
                canonical,
                version=version,
                data_dir=data_dir,
                accept_license=accept_license,
            )
        if canonical == "mixtcrpred" and model_selection:
            from .mixtcrpred import fetch_models
            fetch_models(
                mixtcrpred_path=status.path,
                models=models,
                all_models=all_models,
                high_confidence=high_confidence,
            )
            return _snapshot_status(canonical, data_dir=data_dir)
        return status

    status = artifact_status(canonical, data_dir=data_dir)
    if status.status == "ready":
        if version is not None and version != status.version:
            raise ValueError(
                "%s is managed by %s at version %s; mhctools cannot fetch "
                "version %s" % (
                    canonical, status.manager, status.version, version))
        return status
    raise RuntimeError(status.detail)
