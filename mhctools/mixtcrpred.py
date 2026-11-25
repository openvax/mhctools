# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0

"""Wrapper for the pMHC-specific MixTCRpred TCR-binding models.

MixTCRpred has one transformer checkpoint per fixed peptide/MHC target. The
upstream program is academic/non-commercial software and is never vendored by
mhctools. ``mhctools fetch mixtcrpred --accept-license`` obtains a pinned copy
from its original repository; additional CC-BY-4.0 checkpoints are downloaded
individually from the immutable Zenodo record 7930623.

Inference runs upstream's model and input-encoding classes in a subprocess.
This keeps its runtime dependencies and license boundary separate from the
Apache-licensed mhctools process while preserving upstream V-gene correction,
V-derived CDR1/2 encoding, raw binding score, and calibrated percentile rank.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import hashlib
from importlib.metadata import version as package_version
import json
import os
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.error import HTTPError
from urllib.request import urlopen as _default_urlopen

from datacache import OfflineError, fetch_bytes, fetch_file, require_online
from mhcgnomes import Allele, Pair, ParseError, parse
import pandas as pd

from .artifacts import resolve_artifact_path
from .process_helpers import python_imports_available, resolve_sidecar_python
from .allele_normalization import (
    AlleleIdentityError,
    PARSER_OPTIONS,
    native_allele_identity,
    normalize_allele_name_or_raw,
)
from .pred import (
    COLUMNS,
    Kind,
    MeasurementContext,
    PeptideResult,
    Prediction,
    intern_measurement_context,
)
from .unsupported_residue import UnsupportedResidueError
from .tcr import TCR

urlopen = _default_urlopen


logger = logging.getLogger(__name__)

UPSTREAM_REVISION = "acd6f57444bde675840890207c74ca3b0c7ffac2"
ZENODO_RECORD = "7930623"
ZENODO_API_URL = "https://zenodo.org/api/records/%s" % ZENODO_RECORD
ZENODO_RETRY_DELAYS = (2, 5)
_RETRYABLE_HTTP_STATUS = frozenset((408, 429, 500, 502, 503, 504))
_MODEL_MANIFEST = ".mhctools-mixtcrpred-models.json"
_STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYX-")


@dataclass(frozen=True)
class MixTCRpredMHCTarget:
    """Classified MHC target identity from one MixTCRpred catalog row.

    Attributes
    ----------
    upstream_mhc : str
        Exact ``MHC`` string as written in ``pretrained_models/info_models.csv``.
    canonical_target : str
        Lossless mhcgnomes canonical spelling for the parsed target (for
        example ``HLA-A*02:01``, ``H2-q``, ``HLA-B*08``, ``HLA-DQA1*01``, or
        ``HLA-DQ2``). Non-specific targets are never relabeled with a
        representative allele.
    target_kind : str
        Structured target category: ``"allele"``, ``"pair"``, or ``"raw"`` for
        specific targets; ``"haplotype"``, ``"allele_group"``, ``"serotype"``,
        ``"typing_group"``, ``"invalid_pair"``, ``"gene"``, etc. for
        non-specific targets.
    is_specific_allele : bool
        ``True`` only when the target resolves to a canonical single allele or
        explicit class-II alpha/beta pair under strict normalization.
    reason : str
        Diagnostic reason for exclusion when ``is_specific_allele`` is
        ``False``; empty string for specific targets.
    mhcgnomes_version : str
        Installed ``mhcgnomes`` version that produced this classification.
    """

    upstream_mhc: str
    canonical_target: str
    target_kind: str
    is_specific_allele: bool
    reason: str = ""
    mhcgnomes_version: str = ""

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class MixTCRpredModel:
    """One single-allele pMHC-specific model in the upstream MixTCRpred catalog."""

    name: str
    peptide: str
    allele: str
    upstream_allele: str
    mhc_class: str
    host_species: str
    origin: str
    training_tcrs: int
    auc_5fold: float
    high_confidence: bool
    status: str
    manager: str
    path: str
    version: str = "zenodo.%s" % ZENODO_RECORD

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class MixTCRpredNonSpecificModel:
    """One MixTCRpred catalog model whose MHC target is not a single allele.

    Covers catalog entries specified only as a haplotype (such as ``H2-q``),
    an allele group (``HLA-B*08``, ``HLA-DQA1*01``), or a serotype
    (``HLA-DQ2``). These entries are kept strictly separate from single-allele
    :class:`MixTCRpredModel` records and deliberately expose ``mhc_target`` and
    ``mhc_target_kind`` instead of a single-allele ``allele`` attribute.
    """

    name: str
    peptide: str
    mhc_target: str
    mhc_target_kind: str
    upstream_mhc: str
    mhc_class: str
    host_species: str
    origin: str
    training_tcrs: int
    auc_5fold: float
    high_confidence: bool
    status: str
    manager: str
    path: str
    reason: str
    upstream_revision: str = UPSTREAM_REVISION
    mhcgnomes_version: str = ""
    version: str = "zenodo.%s" % ZENODO_RECORD

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class MixTCRpredCatalog:
    """Partitioned MixTCRpred catalog separating specific and non-specific targets."""

    models: tuple[MixTCRpredModel, ...]
    non_specific_models: tuple[MixTCRpredNonSpecificModel, ...]

    def to_dict(self):
        """Return a JSON-serializable representation."""
        return {
            "models": [model.to_dict() for model in self.models],
            "non_specific_models": [
                model.to_dict() for model in self.non_specific_models
            ],
        }


def classify_catalog_mhc(raw_mhc):
    """Classify one MixTCRpred catalog MHC string without guessing alleles.

    Specific single alleles and explicit class-II alpha/beta pairs resolve
    through :func:`mhctools.allele_normalization.native_allele_identity`.
    Non-specific targets (haplotypes, allele groups, serotypes, typing groups,
    genes, or invalid pairs) are preserved with their canonical mhcgnomes
    target name, structured category, and rejection reason rather than
    choosing a representative allele.
    """
    upstream = str(raw_mhc).strip()
    normalized = upstream.replace(";", "/")
    mhcgnomes_ver = package_version("mhcgnomes")
    try:
        identity = native_allele_identity(upstream)
    except AlleleIdentityError as error:
        parsed = parse(normalized, **PARSER_OPTIONS)
        if isinstance(parsed, Allele):
            if parsed.is_allele_group:
                target_kind = "allele_group"
            elif set(parsed.annotations).intersection("GP"):
                target_kind = "typing_group"
            else:
                target_kind = "allele"
        elif isinstance(parsed, Pair):
            target_kind = "invalid_pair"
        else:
            target_kind = type(parsed).__name__.lower()
        return MixTCRpredMHCTarget(
            upstream_mhc=upstream,
            canonical_target=parsed.to_string(),
            target_kind=target_kind,
            is_specific_allele=False,
            reason=str(error),
            mhcgnomes_version=mhcgnomes_ver,
        )

    try:
        parsed = parse(normalized, **PARSER_OPTIONS)
    except ParseError:
        target_kind = "raw"
    else:
        if isinstance(parsed, Allele) and parsed.is_allele_group:
            parsed = parse(
                normalized, **dict(PARSER_OPTIONS, use_allele_aliases=True))
        target_kind = "pair" if isinstance(parsed, Pair) else "allele"
    return MixTCRpredMHCTarget(
        upstream_mhc=upstream,
        canonical_target=identity,
        target_kind=target_kind,
        is_specific_allele=True,
        reason="",
        mhcgnomes_version=mhcgnomes_ver,
    )


def resolve_non_specific_mhc_target(raw_mhc):
    """Parse a caller-supplied non-specific MHC target query.

    Accepts haplotypes (e.g. ``H-2q`` or ``H2-q``), allele groups
    (``HLA-B*08``, ``HLA-DQA*01`` or ``HLA-DQA1*01``), and serotypes
    (``HLA-DQ2``). Never chooses a representative allele. Raises
    :class:`ValueError` when *raw_mhc* identifies a specific single allele
    (use :meth:`MixTCRpred.resolve_model` instead) and
    :class:`mhcgnomes.ParseError` when *raw_mhc* cannot be parsed.
    """
    upstream = str(raw_mhc).strip()
    normalized = upstream.replace(";", "/")
    parse(normalized, **PARSER_OPTIONS)
    classified = classify_catalog_mhc(upstream)
    if classified.is_specific_allele:
        raise ValueError(
            "%r resolves to specific allele %s; use "
            "MixTCRpred.resolve_model() for single-allele targets"
            % (raw_mhc, classified.canonical_target))
    return classified


def _manager_for_home(home):
    return (
        "mhctools"
        if (Path(home) / ".mhctools-artifact.json").is_file()
        else "user"
    )


def partition_model_catalog(mixtcrpred_path=None, data_dir=None):
    """Partition the MixTCRpred catalog into specific and non-specific models.

    Returns a :class:`MixTCRpredCatalog` whose ``models`` tuple contains only
    strictly verified single-allele / explicit-pair models and whose
    ``non_specific_models`` tuple records haplotype, allele-group, and
    serotype catalog entries with structured classification and provenance.
    """
    home = resolve_artifact_path("mixtcrpred", mixtcrpred_path, data_dir=data_dir)
    manager = _manager_for_home(home)
    catalog_path = home / "pretrained_models" / "info_models.csv"
    specific = []
    non_specific = []
    with catalog_path.open(newline="") as input_file:
        for row in csv.DictReader(input_file):
            name = row["MixTCRpred_model_name"].strip()
            training_tcrs = int(row["Number_training_abTCR"])
            # The pinned table contains one zero-training placeholder
            # (A0201_SLLMWITQC) with no checkpoint in either the repository or
            # immutable Zenodo record. Upstream documents 146 usable models,
            # so do not advertise that placeholder as fetchable inference.
            if training_tcrs <= 0:
                continue
            checkpoint = (
                home / "pretrained_models" / ("model_%s.ckpt" % name))
            upstream_class = row["MHC_class"].strip()
            mhc_class = {"MHCI": "I", "MHCII": "II"}.get(
                upstream_class, upstream_class)
            peptide = row["Peptide"].strip()
            host_species = row["Host_species"].strip()
            origin = row["Origin"].strip()
            auc_5fold = float(row["AUC_5fold"])
            # Match upstream's --download_high implementation.
            high_confidence = training_tcrs >= 50
            status = "ready" if checkpoint.is_file() else "missing"
            target = classify_catalog_mhc(row["MHC"].strip())
            if target.is_specific_allele:
                specific.append(MixTCRpredModel(
                    name=name,
                    peptide=peptide,
                    allele=target.canonical_target,
                    upstream_allele=target.upstream_mhc,
                    mhc_class=mhc_class,
                    host_species=host_species,
                    origin=origin,
                    training_tcrs=training_tcrs,
                    auc_5fold=auc_5fold,
                    high_confidence=high_confidence,
                    status=status,
                    manager=manager,
                    path=str(checkpoint),
                ))
            else:
                non_specific.append(MixTCRpredNonSpecificModel(
                    name=name,
                    peptide=peptide,
                    mhc_target=target.canonical_target,
                    mhc_target_kind=target.target_kind,
                    upstream_mhc=target.upstream_mhc,
                    mhc_class=mhc_class,
                    host_species=host_species,
                    origin=origin,
                    training_tcrs=training_tcrs,
                    auc_5fold=auc_5fold,
                    high_confidence=high_confidence,
                    status=status,
                    manager=manager,
                    path=str(checkpoint),
                    reason=target.reason,
                    upstream_revision=UPSTREAM_REVISION,
                    mhcgnomes_version=target.mhcgnomes_version,
                ))
    return MixTCRpredCatalog(
        models=tuple(specific),
        non_specific_models=tuple(non_specific),
    )


def model_catalog(mixtcrpred_path=None, data_dir=None):
    """Return all single-allele pMHC-specific MixTCRpred models and weight status.

    Catalog rows whose MHC target is a haplotype, allele group, or serotype
    are excluded with a warning and remain available via
    :func:`non_specific_model_catalog`.
    """
    catalog = partition_model_catalog(mixtcrpred_path, data_dir=data_dir)
    for entry in catalog.non_specific_models:
        logger.warning("Excluding MixTCRpred model %s: %s", entry.name, entry.reason)
    return catalog.models


def non_specific_model_catalog(mixtcrpred_path=None, data_dir=None):
    """Return MixTCRpred catalog models with non-specific MHC targets.

    Returns :class:`MixTCRpredNonSpecificModel` records for haplotype,
    allele-group, and serotype targets without choosing representative
    alleles.
    """
    return partition_model_catalog(
        mixtcrpred_path, data_dir=data_dir).non_specific_models


def select_catalog_models(
        catalog,
        models=None,
        all_models=False,
        high_confidence=False,
        include_non_specific=False,
        non_specific_catalog=()):
    """Select MixTCRpred catalog models while enforcing non-specific target opt-in.

    Parameters
    ----------
    catalog : MixTCRpredCatalog or iterable of MixTCRpredModel
        Single-allele catalog entries, or a partitioned :class:`MixTCRpredCatalog`.
    models : str or iterable of str, optional
        Explicit model names to select.
    all_models : bool, optional
        Select all catalog models (single-allele only unless
        ``include_non_specific=True``).
    high_confidence : bool, optional
        Select upstream high-confidence models (single-allele only unless
        ``include_non_specific=True``).
    include_non_specific : bool, optional
        Explicit opt-in allowing non-specific group/haplotype/serotype models
        to be selected. Default ``False``.
    non_specific_catalog : iterable of MixTCRpredNonSpecificModel, optional
        Non-specific entries when *catalog* is not a :class:`MixTCRpredCatalog`.

    Returns
    -------
    list of MixTCRpredModel or MixTCRpredNonSpecificModel
    """
    if isinstance(catalog, MixTCRpredCatalog):
        specific_entries = catalog.models
        non_specific_entries = catalog.non_specific_models
    else:
        specific_entries = tuple(catalog)
        non_specific_entries = tuple(non_specific_catalog)

    if sum(bool(value) for value in (models, all_models, high_confidence)) > 1:
        raise ValueError(
            "Choose model names, all_models, or high_confidence, not multiple")
    pool = (
        specific_entries + non_specific_entries
        if include_non_specific else specific_entries
    )
    if all_models:
        return list(pool)
    if high_confidence:
        return [model for model in pool if model.high_confidence]
    if not models:
        return []
    if isinstance(models, str):
        models = [models]
    by_name = {model.name.lower(): model for model in specific_entries}
    non_specific_by_name = {
        model.name.lower(): model for model in non_specific_entries
    }
    selected = []
    for requested in models:
        key = str(requested).lower()
        if key in by_name:
            model = by_name[key]
        elif key in non_specific_by_name:
            non_specific_entry = non_specific_by_name[key]
            if not include_non_specific:
                raise AlleleIdentityError(
                    "MixTCRpred model %s targets non-specific MHC %s (%s, "
                    "upstream %r) and is excluded from the single-allele "
                    "catalog. mhctools does not relabel non-specific models "
                    "with a representative allele. Inspect with "
                    "MixTCRpred.non_specific_catalog() or opt in explicitly "
                    "with allow_non_specific_target=True / "
                    "include_non_specific=True." % (
                        non_specific_entry.name,
                        non_specific_entry.mhc_target,
                        non_specific_entry.mhc_target_kind,
                        non_specific_entry.upstream_mhc,
                    )
                )
            model = non_specific_entry
        else:
            raise ValueError(
                "Unknown MixTCRpred model %r. Use `mhctools ls mixtcrpred "
                "--models` to inspect the catalog." % requested)
        if model not in selected:
            selected.append(model)
    return selected


def _retry_zenodo(description, operation):
    """Run an idempotent Zenodo read with bounded transient retries."""
    attempts = len(ZENODO_RETRY_DELAYS) + 1
    for attempt in range(1, attempts + 1):
        try:
            return operation()
        except OfflineError:
            raise
        except (OSError, ValueError) as error:
            retryable = (
                not isinstance(error, HTTPError)
                or error.code in _RETRYABLE_HTTP_STATUS
            )
            if not retryable or attempt == attempts:
                raise RuntimeError(
                    "Could not %s after %d attempt%s: %s" % (
                        description,
                        attempt,
                        "" if attempt == 1 else "s",
                        error,
                    )) from error
            delay = ZENODO_RETRY_DELAYS[attempt - 1]
            print(
                "MixTCRpred %s attempt %d/%d failed: %s; retrying in %ss"
                % (description, attempt, attempts, error, delay),
                file=sys.stderr,
            )
            time.sleep(delay)


def _zenodo_files():
    require_online(ZENODO_API_URL)
    if urlopen is _default_urlopen:
        try:
            raw = fetch_bytes(
                ZENODO_API_URL,
                timeout=60,
                max_retries=len(ZENODO_RETRY_DELAYS),
                retry_backoff=(
                    ZENODO_RETRY_DELAYS[0] if ZENODO_RETRY_DELAYS else 0
                ),
            )
            record = json.loads(raw.decode("utf-8"))
        except OfflineError:
            raise
        except (OSError, ValueError) as error:
            raise RuntimeError(
                "Could not read MixTCRpred model metadata from %s: %s"
                % (ZENODO_API_URL, error)) from error
    else:
        def read_metadata():
            with urlopen(ZENODO_API_URL, timeout=60) as response:
                return json.load(response)

        record = _retry_zenodo(
            "read MixTCRpred model metadata from %s" % ZENODO_API_URL,
            read_metadata,
        )
    if str(record.get("id")) != ZENODO_RECORD:
        raise RuntimeError("Zenodo returned unexpected MixTCRpred record")
    return {entry["key"]: entry for entry in record.get("files", ())}


def _download_zenodo_file(url, destination, expected_size=None):
    require_online(url)
    if urlopen is _default_urlopen:
        try:
            fetch_file(
                url,
                destination=destination,
                force=True,
                timeout=120,
                expected_size=expected_size,
                max_bytes=expected_size,
                max_retries=len(ZENODO_RETRY_DELAYS),
                retry_backoff=(
                    ZENODO_RETRY_DELAYS[0] if ZENODO_RETRY_DELAYS else 0
                ),
                raw=True,
            )
        except OfflineError:
            raise
        except (OSError, ValueError) as error:
            raise RuntimeError(
                "Could not download MixTCRpred checkpoint from %s: %s"
                % (url, error)) from error
        return

    def download():
        # Truncate on every attempt so a partial response is never appended to
        # the retry. The caller checksum-verifies before atomic installation.
        with urlopen(url, timeout=120) as response:
            with Path(destination).open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)

    _retry_zenodo("download MixTCRpred checkpoint from %s" % url, download)


def _hashes(path):
    # Zenodo publishes MD5 as its file identity; SHA-256 is recorded too.
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with Path(path).open("rb") as input_file:
        for block in iter(lambda: input_file.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def _read_model_manifest(model_dir):
    path = Path(model_dir) / _MODEL_MANIFEST
    if not path.is_file():
        return {"record": ZENODO_RECORD, "models": {}}
    try:
        with path.open() as input_file:
            value = json.load(input_file)
    except (OSError, ValueError) as error:
        raise RuntimeError("Invalid MixTCRpred model manifest: %s" % path) from error
    if str(value.get("record")) != ZENODO_RECORD:
        raise RuntimeError("Unexpected Zenodo record in %s" % path)
    value.setdefault("models", {})
    return value


def _write_model_manifest(model_dir, manifest):
    model_dir = Path(model_dir)
    destination = model_dir / _MODEL_MANIFEST
    with tempfile.NamedTemporaryFile(
            mode="w", dir=str(model_dir), prefix=".manifest-",
            delete=False) as output:
        temporary = Path(output.name)
        json.dump(manifest, output, indent=2, sort_keys=True)
        output.write("\n")
    os.replace(str(temporary), str(destination))


def _fetch_one_model(model, files, manifest):
    filename = "model_%s.ckpt" % model.name
    try:
        remote = files[filename]
    except KeyError:
        raise RuntimeError(
            "%s is in the upstream catalog but absent from Zenodo record %s"
            % (model.name, ZENODO_RECORD))
    expected_md5 = remote["checksum"].removeprefix("md5:")
    expected_size = int(remote["size"])
    destination = Path(model.path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.is_file():
        actual_md5, actual_sha256 = _hashes(destination)
        if destination.stat().st_size != expected_size or actual_md5 != expected_md5:
            raise RuntimeError(
                "Existing MixTCRpred checkpoint does not match Zenodo: %s. "
                "Move it aside before fetching again." % destination)
    else:
        url = remote["links"]["self"]
        print("Downloading %s from Zenodo" % model.name, file=sys.stderr)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(
                    mode="wb", dir=str(destination.parent),
                    prefix=".%s-" % filename, delete=False) as output:
                temporary = Path(output.name)
            _download_zenodo_file(url, temporary, expected_size=expected_size)
            actual_md5, actual_sha256 = _hashes(temporary)
            if temporary.stat().st_size != expected_size or actual_md5 != expected_md5:
                raise RuntimeError(
                    "Downloaded MixTCRpred checkpoint failed its Zenodo "
                    "size/checksum verification: %s" % model.name)
            os.replace(str(temporary), str(destination))
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()

    manifest["models"][model.name] = {
        "filename": filename,
        "md5": actual_md5,
        "sha256": actual_sha256,
        "size": expected_size,
        "url": remote["links"]["self"],
    }


def fetch_models(
        mixtcrpred_path=None,
        data_dir=None,
        models=None,
        all_models=False,
        high_confidence=False,
        refresh_metadata=False,
        include_non_specific=False):
    """Fetch weights; verify recorded downloads locally unless refreshing.

    An immutable Zenodo record and its saved MD5/SHA-256/size manifest
    identify existing downloads. Unrecorded weights still require Zenodo
    verification. ``refresh_metadata=True`` rechecks remote metadata too.
    Non-specific group/haplotype/serotype models are excluded by default and
    require ``include_non_specific=True``.
    """
    home = resolve_artifact_path("mixtcrpred", mixtcrpred_path, data_dir=data_dir)
    catalog = partition_model_catalog(mixtcrpred_path=home)
    selected = select_catalog_models(
        catalog,
        models=models,
        all_models=all_models,
        high_confidence=high_confidence,
        include_non_specific=include_non_specific,
    )
    if not selected:
        return ()
    model_dir = home / "pretrained_models"
    manifest = _read_model_manifest(model_dir)
    pending = []
    for model in selected:
        path = Path(model.path)
        recorded = manifest["models"].get(model.name)
        if refresh_metadata or recorded is None or not path.is_file():
            pending.append(model)
            continue
        md5, sha256 = _hashes(path)
        if (recorded.get("filename") != path.name or
                recorded.get("size") != path.stat().st_size or
                recorded.get("md5") != md5 or recorded.get("sha256") != sha256):
            raise RuntimeError(
                "Existing MixTCRpred checkpoint does not match its recorded "
                "manifest: %s. Move it aside before fetching again." % path)
    files = _zenodo_files() if pending else {}
    for model in pending:
        _fetch_one_model(model, files, manifest)
        _write_model_manifest(model_dir, manifest)
    refreshed_catalog = partition_model_catalog(mixtcrpred_path=home)
    refreshed = {
        model.name: model
        for model in (
            refreshed_catalog.models + refreshed_catalog.non_specific_models
        )
    }
    return tuple(refreshed[model.name] for model in selected)


def print_model_catalog(
        mixtcrpred_path=None,
        data_dir=None,
        downloaded=False,
        high_confidence=False,
        json_output=False,
        non_specific=False):
    """Print the model catalog for ``mhctools ls mixtcrpred --models``."""
    if non_specific:
        models = non_specific_model_catalog(mixtcrpred_path, data_dir=data_dir)
    else:
        models = model_catalog(mixtcrpred_path, data_dir=data_dir)
    if downloaded:
        models = tuple(model for model in models if model.status == "ready")
    if high_confidence:
        models = tuple(model for model in models if model.high_confidence)
    if json_output:
        print(json.dumps([model.to_dict() for model in models], indent=2))
        return
    if non_specific:
        headers = (
            "MODEL", "PEPTIDE", "MHC_TARGET", "TARGET_KIND", "UPSTREAM_MHC",
            "CLASS", "SPECIES", "TRAINING", "AUC", "CONFIDENCE", "STATUS",
            "MANAGER", "PATH")
        rows = [(
            model.name,
            model.peptide,
            model.mhc_target,
            model.mhc_target_kind,
            model.upstream_mhc,
            model.mhc_class,
            model.host_species,
            str(model.training_tcrs),
            "%.3f" % model.auc_5fold,
            "high" if model.high_confidence else "low",
            model.status,
            model.manager,
            model.path,
        ) for model in models]
    else:
        headers = (
            "MODEL", "PEPTIDE", "MHC", "CLASS", "SPECIES", "TRAINING",
            "AUC", "CONFIDENCE", "STATUS", "MANAGER", "PATH")
        rows = [(
            model.name,
            model.peptide,
            model.allele,
            model.mhc_class,
            model.host_species,
            str(model.training_tcrs),
            "%.3f" % model.auc_5fold,
            "high" if model.high_confidence else "low",
            model.status,
            model.manager,
            model.path,
        ) for model in models]
    from .cli.formatting import format_table
    print(format_table(headers, rows, stream=sys.stdout))


# Third-party modules reachable from the sidecar's own import chain
# (src/dataloaders, src/models, src/utils). Deliberately not every import under
# upstream's src/: natsort, for one, is only used by a V-gene helper that
# inference never loads, and requiring it would skip a runtime that works.
MIXTCRPRED_SIDECAR_MODULES = (
    "torch", "torchvision", "pytorch_lightning",
    "numpy", "pandas", "scipy", "sklearn",
)


def mixtcrpred_runtime_available(mixtcrpred_python=None, timeout=120):
    """Whether that interpreter can run MixTCRpred inference.

    A downloaded checkpoint reports ``ready`` once its files are located, which
    says nothing about the interpreter: upstream's ``src/`` imports torchvision
    and pytorch_lightning, which the interpreter running mhctools need not have.

    Probes under the environment inference actually uses. `PYTHONNOUSERSITE=1`
    hides `pip install --user` packages from the sidecar, so an ambient probe
    would see a torch the run cannot.

    A `MIXTCRPRED_PYTHON` that does not exist raises rather than reporting an
    unavailable runtime: that is a misconfiguration to fix, not a backend to
    skip over.
    """
    environment = os.environ.copy()
    environment["PYTHONNOUSERSITE"] = "1"
    return python_imports_available(
        resolve_sidecar_python(mixtcrpred_python, "MIXTCRPRED_PYTHON"), MIXTCRPRED_SIDECAR_MODULES,
        timeout=timeout, env=environment)


def _validate_tcr(tcr):
    if not isinstance(tcr, TCR):
        raise TypeError("Expected mhctools.TCR, got %r" % type(tcr))
    for name in ("cdr3a", "cdr3b"):
        sequence = getattr(tcr, name)
        if not sequence:
            raise ValueError("MixTCRpred requires %s" % name)
        if len(sequence) > 20:
            raise ValueError(
                "MixTCRpred %s is limited to 20 residues: %r" % (
                    name, sequence))
        invalid = sorted(set(sequence) - _STANDARD_AMINO_ACIDS)
        if invalid:
            raise UnsupportedResidueError(
                "MixTCRpred %s contains unsupported residues %s: %r" % (
                    name, ", ".join(invalid), sequence),
                peptides=(sequence,),
                residues=invalid,
                field=name,
                predictor_name="MixTCRpred",
            )


class MixTCRpred:
    """Score paired alpha/beta TCRs for one fixed pMHC target.

    Parameters
    ----------
    model : str
        A model name from :meth:`catalog`, for example
        ``"A0201_GILGFVFTL"``. Each model fixes both peptide and MHC.
    mixtcrpred_path : str, optional
        Licensed upstream checkout. Otherwise ``MIXTCRPRED_HOME``, common
        checkout paths, and the mhctools artifact cache are searched.
    mixtcrpred_python : str, optional
        Python interpreter containing upstream's PyTorch/Lightning runtime.
        Defaults to ``MIXTCRPRED_PYTHON`` and then the current interpreter.
    checkpoint_path : str, optional
        Explicit trusted checkpoint override. PyTorch checkpoints may execute
        code while loading; the normal path is resolved from the pinned or
        user-managed catalog instead.
    batch_size : int, optional
        Deterministic CPU inference batch size.
    allow_non_specific_target : bool, optional
        Explicit opt-in for catalog models whose upstream MHC target is a
        haplotype, allele group, or serotype rather than a single canonical
        allele. Default ``False``. Non-specific models are never advertised as
        ``single_allele`` and leave ``Prediction.allele`` empty.

    Notes
    -----
    MixTCRpred accepts J-gene columns but its released network consumes only
    CDR3alpha/beta and CDR1/2 sequences derived from the V genes. J genes are
    retained in :class:`TCR` and reported by QC, but do not affect scores.
    """

    @classmethod
    def fetch(
            cls,
            models=None,
            version=None,
            data_dir=None,
            accept_license=False,
            all_models=False,
            high_confidence=False,
            include_non_specific=False):
        """Fetch licensed upstream code and selected CC-BY model weights."""
        from .artifacts import fetch
        return fetch(
            "mixtcrpred",
            version=version,
            data_dir=data_dir,
            accept_license=accept_license,
            models=models,
            all_models=all_models,
            high_confidence=high_confidence,
            include_non_specific=include_non_specific,
        )

    @classmethod
    def catalog(cls, mixtcrpred_path=None, data_dir=None):
        """Return the available single-allele pMHC model catalog."""
        return model_catalog(mixtcrpred_path, data_dir=data_dir)

    @classmethod
    def non_specific_catalog(cls, mixtcrpred_path=None, data_dir=None):
        """Return catalog models with non-specific group/haplotype/serotype targets."""
        return non_specific_model_catalog(mixtcrpred_path, data_dir=data_dir)

    @classmethod
    def resolve_model(
            cls, peptide, allele, mixtcrpred_path=None, data_dir=None):
        """Resolve a unique single-allele catalog model from its peptide/MHC target."""
        normalized = normalize_allele_name_or_raw(allele)
        matches = [
            model for model in cls.catalog(mixtcrpred_path, data_dir=data_dir)
            if model.peptide == str(peptide).upper()
            and model.allele == normalized
        ]
        if len(matches) != 1:
            raise ValueError(
                "Expected one MixTCRpred model for %s / %s, found %d" % (
                    peptide, normalized, len(matches)))
        return matches[0]

    @classmethod
    def resolve_non_specific_model(
            cls, peptide, mhc_target, mixtcrpred_path=None, data_dir=None):
        """Resolve a unique non-specific catalog model from peptide and MHC group/haplotype/serotype."""
        target = resolve_non_specific_mhc_target(mhc_target)
        matches = [
            model
            for model in cls.non_specific_catalog(
                mixtcrpred_path, data_dir=data_dir)
            if model.peptide == str(peptide).upper()
            and model.mhc_target == target.canonical_target
            and model.mhc_target_kind == target.target_kind
        ]
        if len(matches) != 1:
            raise ValueError(
                "Expected one non-specific MixTCRpred model for %s / %s (%s), "
                "found %d" % (
                    peptide,
                    target.canonical_target,
                    target.target_kind,
                    len(matches),
                ))
        return matches[0]

    def __init__(
            self,
            model,
            mixtcrpred_path=None,
            mixtcrpred_python=None,
            checkpoint_path=None,
            batch_size=256,
            allow_non_specific_target=False):
        self.mixtcrpred_home = resolve_artifact_path("mixtcrpred", mixtcrpred_path)
        self.mixtcrpred_python = resolve_sidecar_python(mixtcrpred_python, "MIXTCRPRED_PYTHON")
        catalog = partition_model_catalog(mixtcrpred_path=self.mixtcrpred_home)
        selected = select_catalog_models(
            catalog,
            models=[model],
            include_non_specific=allow_non_specific_target,
        )
        self.model_info = selected[0]
        self.is_non_specific_target = isinstance(
            self.model_info, MixTCRpredNonSpecificModel)
        self.checkpoint_path = str(
            Path(checkpoint_path).expanduser().resolve()
            if checkpoint_path else self.model_info.path)
        if not Path(self.checkpoint_path).is_file():
            non_specific_flag = (
                " --allow-non-specific" if self.is_non_specific_target else ""
            )
            raise FileNotFoundError(
                "MixTCRpred checkpoint is missing: %s. Run `mhctools fetch "
                "mixtcrpred --model %s%s --accept-license`." % (
                    self.checkpoint_path,
                    self.model_info.name,
                    non_specific_flag,
                ))
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.batch_size = int(batch_size)
        self.last_qc = pd.DataFrame()

    def __repr__(self):
        if self.is_non_specific_target:
            return (
                "MixTCRpred(model=%s, mhc_target=%s, mhc_target_kind=%s)" % (
                    self.model_info.name,
                    self.model_info.mhc_target,
                    self.model_info.mhc_target_kind,
                )
            )
        return "MixTCRpred(model=%s)" % self.model_info.name

    def __str__(self):
        return repr(self)

    def _predictor_name(self):
        return "mixtcrpred"

    @property
    def predictor_version(self):
        return "%s+zenodo.%s:%s" % (
            UPSTREAM_REVISION, ZENODO_RECORD, self.model_info.name)

    def kind_support(self):
        if self.is_non_specific_target:
            return {
                Kind.pMHC_TCR_binding: {
                    "mhc_dependence": "mhc_group",
                    "mhc_class": self.model_info.mhc_class,
                    "mhc_target": self.model_info.mhc_target,
                    "mhc_target_kind": self.model_info.mhc_target_kind,
                }
            }
        return {
            Kind.pMHC_TCR_binding: {
                "mhc_dependence": "single_allele",
                "mhc_class": self.model_info.mhc_class,
            }
        }

    @property
    def supported_kinds(self):
        return tuple(self.kind_support())

    def _run_sidecar(self, tcrs):
        if not tcrs:
            self.last_qc = pd.DataFrame()
            return pd.DataFrame(columns=(
                "__mhctools_id", "score", "percentile_rank"))
        rows = []
        for index, tcr in enumerate(tcrs):
            _validate_tcr(tcr)
            rows.append({
                "__mhctools_id": index,
                "cdr3_TRA": tcr.cdr3a,
                "cdr3_TRB": tcr.cdr3b,
                "TRAV": tcr.trav,
                "TRAJ": tcr.traj,
                "TRBV": tcr.trbv,
                "TRBJ": tcr.trbj,
            })
        with tempfile.TemporaryDirectory(prefix="mhctools_mixtcrpred_") as tmp:
            input_path = Path(tmp) / "input.csv"
            output_path = Path(tmp) / "output.csv"
            pd.DataFrame(rows).to_csv(input_path, index=False)
            sidecar = Path(__file__).with_name("mixtcrpred_sidecar.py")
            command = [
                self.mixtcrpred_python,
                "-c",
                (
                    "import runpy,sys; sys.argv=sys.argv[1:]; "
                    "runpy.run_path(sys.argv[0],run_name='__main__')"
                ),
                str(sidecar),
                "--home", str(self.mixtcrpred_home),
                "--model", self.model_info.name,
                "--peptide", self.model_info.peptide,
                "--checkpoint", self.checkpoint_path,
                "--input", str(input_path),
                "--output", str(output_path),
                "--batch-size", str(self.batch_size),
                "--host", self.model_info.host_species,
            ]
            environment = os.environ.copy()
            # PyTorch 2.6 changed the default and otherwise rejects upstream's
            # Lightning 1.6 checkpoints before loading their state dict. A
            # checkpoint can execute code, so managed files are source-pinned
            # or checksummed and custom/user-managed paths must be trusted.
            environment["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
            environment["PYTHONNOUSERSITE"] = "1"
            process = subprocess.run(
                command,
                cwd=str(self.mixtcrpred_home),
                env=environment,
                capture_output=True,
                text=True,
            )
            if process.returncode != 0:
                raise RuntimeError(
                    "MixTCRpred inference failed (exit %d):\n%s" % (
                        process.returncode,
                        (process.stderr or process.stdout).strip()))
            output = pd.read_csv(output_path, keep_default_na=False)

        expected_ids = list(range(len(tcrs)))
        if sorted(output["__mhctools_id"].tolist()) != expected_ids:
            raise RuntimeError(
                "MixTCRpred output did not preserve every input TCR")
        output = output.sort_values("__mhctools_id").reset_index(drop=True)
        qc_columns = [
            column for column in output.columns
            if column not in ("score", "percentile_rank")
        ]
        self.last_qc = output[qc_columns].copy()
        return output

    def predict_tcrs(self, tcrs):
        """Return one :class:`Prediction` per input TCR, in input order."""
        if isinstance(tcrs, TCR):
            tcrs = [tcrs]
        tcrs = list(tcrs)
        output = self._run_sidecar(tcrs)
        if self.is_non_specific_target:
            allele = ""
            measurement_context = intern_measurement_context(MeasurementContext(
                estimate_type="ml_predicted",
                detail=(
                    "non_specific_mhc_target:%s=%s (upstream=%s)" % (
                        self.model_info.mhc_target_kind,
                        self.model_info.mhc_target,
                        self.model_info.upstream_mhc,
                    )
                ),
            ))
        else:
            allele = self.model_info.allele
            measurement_context = None
        return [
            Prediction(
                kind=Kind.pMHC_TCR_binding,
                score=float(row.score),
                peptide=self.model_info.peptide,
                allele=allele,
                tcr=tcr.identifier,
                percentile_rank=float(row.percentile_rank),
                predictor_name=self._predictor_name(),
                predictor_version=self.predictor_version,
                measurement_context=measurement_context,
            )
            for tcr, row in zip(tcrs, output.itertuples(index=False))
        ]

    def predict(self, tcrs):
        """Return one single-prediction :class:`PeptideResult` per TCR."""
        return [
            PeptideResult(preds=(prediction,))
            for prediction in self.predict_tcrs(tcrs)
        ]

    def predict_dataframe(self, tcrs, sample_name=""):
        """Return standard mhctools prediction columns for the input TCRs."""
        if isinstance(tcrs, pd.DataFrame):
            tcrs = [TCR.from_series(row) for _, row in tcrs.iterrows()]
        rows = [
            prediction.to_row(sample_name)
            for prediction in self.predict_tcrs(tcrs)
        ]
        return pd.DataFrame(rows, columns=COLUMNS)

    def annotate_dataframe(self, dataframe):
        """Preserve a TCR table and append scores, target metadata, and QC.

        This captures the two primary upstream outputs (``score`` and
        ``perc_rank``), its model-level metadata, and the corrected V/J genes,
        V-derived CDR1/2 loops, and warning status otherwise written only to
        upstream's optional logfile.
        """
        dataframe = dataframe.copy()
        tcrs = [TCR.from_series(row) for _, row in dataframe.iterrows()]
        output = self._run_sidecar(tcrs)
        dataframe["mixtcrpred_model"] = self.model_info.name
        dataframe["mixtcrpred_peptide"] = self.model_info.peptide
        if self.is_non_specific_target:
            dataframe["mixtcrpred_allele"] = ""
            dataframe["mixtcrpred_mhc_target"] = self.model_info.mhc_target
            dataframe["mixtcrpred_mhc_target_kind"] = (
                self.model_info.mhc_target_kind)
            dataframe["mixtcrpred_upstream_mhc"] = self.model_info.upstream_mhc
        else:
            dataframe["mixtcrpred_allele"] = self.model_info.allele
        dataframe["mixtcrpred_score"] = output["score"].to_numpy()
        dataframe["mixtcrpred_percentile_rank"] = output[
            "percentile_rank"].to_numpy()
        for column in (
                "trav_corrected", "traj_corrected", "trbv_corrected",
                "trbj_corrected", "cdr1a_derived", "cdr2a_derived",
                "cdr1b_derived", "cdr2b_derived", "warning"):
            dataframe["mixtcrpred_%s" % column] = output[column].to_numpy()
        return dataframe
