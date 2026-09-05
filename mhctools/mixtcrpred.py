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
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.request import urlopen

import pandas as pd

from .allele_normalization import normalize_allele_name_or_raw
from .pred import COLUMNS, Kind, PeptideResult, Prediction
from .tcr import TCR


UPSTREAM_REVISION = "acd6f57444bde675840890207c74ca3b0c7ffac2"
ZENODO_RECORD = "7930623"
ZENODO_API_URL = "https://zenodo.org/api/records/%s" % ZENODO_RECORD
_MODEL_MANIFEST = ".mhctools-mixtcrpred-models.json"
_STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYX-")
_REQUIRED_HOME_PATHS = (
    "MixTCRpred.py",
    "src/dataloaders.py",
    "src/models.py",
    "src/utils.py",
    "pretrained_models/anchors_perc_rank.pickle",
    "pretrained_models/info_models.csv",
)


@dataclass(frozen=True)
class MixTCRpredModel:
    """One pMHC-specific model in the upstream MixTCRpred catalog."""

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


def _valid_home(path):
    path = Path(path)
    return all((path / relative).is_file() for relative in _REQUIRED_HOME_PATHS)


def _resolve_home(mixtcrpred_path=None, data_dir=None):
    """Resolve a user-managed or mhctools-managed MixTCRpred checkout."""
    for source, value in (
            ("mixtcrpred_path argument", mixtcrpred_path),
            ("MIXTCRPRED_HOME", os.environ.get("MIXTCRPRED_HOME"))):
        if value:
            path = Path(value).expanduser().resolve()
            if not _valid_home(path):
                raise FileNotFoundError(
                    "MixTCRpred directory from %s is incomplete: %s" % (
                        source, path))
            return path

    for value in ("~/MixTCRpred", "~/code/MixTCRpred"):
        path = Path(value).expanduser()
        if _valid_home(path):
            return path.resolve()

    from .artifacts import artifact_status
    status = artifact_status("mixtcrpred", data_dir=data_dir)
    if status.status == "ready":
        return Path(status.path).resolve()
    raise FileNotFoundError(
        "MixTCRpred is not installed. Run `mhctools fetch mixtcrpred "
        "--accept-license`, or set MIXTCRPRED_HOME to a licensed checkout.")


def _manager_for_home(home):
    return (
        "mhctools"
        if (Path(home) / ".mhctools-artifact.json").is_file()
        else "user"
    )


def model_catalog(mixtcrpred_path=None, data_dir=None):
    """Return all pMHC-specific MixTCRpred models and local weight status."""
    home = _resolve_home(mixtcrpred_path, data_dir=data_dir)
    manager = _manager_for_home(home)
    catalog_path = home / "pretrained_models" / "info_models.csv"
    result = []
    with catalog_path.open(newline="") as input_file:
        for row in csv.DictReader(input_file):
            name = row["MixTCRpred_model_name"].strip()
            checkpoint = (
                home / "pretrained_models" / ("model_%s.ckpt" % name))
            upstream_class = row["MHC_class"].strip()
            result.append(MixTCRpredModel(
                name=name,
                peptide=row["Peptide"].strip(),
                allele=normalize_allele_name_or_raw(row["MHC"].strip()),
                upstream_allele=row["MHC"].strip(),
                mhc_class={"MHCI": "I", "MHCII": "II"}.get(
                    upstream_class, upstream_class),
                host_species=row["Host_species"].strip(),
                origin=row["Origin"].strip(),
                training_tcrs=int(row["Number_training_abTCR"]),
                auc_5fold=float(row["AUC_5fold"]),
                # Match upstream's --download_high implementation.
                high_confidence=int(row["Number_training_abTCR"]) >= 50,
                status="ready" if checkpoint.is_file() else "missing",
                manager=manager,
                path=str(checkpoint),
            ))
    return tuple(result)


def _select_catalog_models(
        catalog,
        models=None,
        all_models=False,
        high_confidence=False):
    if sum(bool(value) for value in (models, all_models, high_confidence)) > 1:
        raise ValueError(
            "Choose model names, all_models, or high_confidence, not multiple")
    if all_models:
        return list(catalog)
    if high_confidence:
        return [model for model in catalog if model.high_confidence]
    if not models:
        return []
    if isinstance(models, str):
        models = [models]
    by_name = {model.name.lower(): model for model in catalog}
    selected = []
    for requested in models:
        try:
            model = by_name[str(requested).lower()]
        except KeyError:
            raise ValueError(
                "Unknown MixTCRpred model %r. Use `mhctools ls mixtcrpred "
                "--models` to inspect the catalog." % requested)
        if model not in selected:
            selected.append(model)
    return selected


def _zenodo_files():
    try:
        with urlopen(ZENODO_API_URL, timeout=60) as response:
            record = json.load(response)
    except (OSError, ValueError) as error:
        raise RuntimeError(
            "Could not read MixTCRpred model metadata from %s" %
            ZENODO_API_URL) from error
    if str(record.get("id")) != ZENODO_RECORD:
        raise RuntimeError("Zenodo returned unexpected MixTCRpred record")
    return {entry["key"]: entry for entry in record.get("files", ())}


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
                with urlopen(url, timeout=120) as response:
                    shutil.copyfileobj(response, output, length=1024 * 1024)
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
        high_confidence=False):
    """Fetch selected MixTCRpred weights and return their catalog entries."""
    home = _resolve_home(mixtcrpred_path, data_dir=data_dir)
    catalog = model_catalog(mixtcrpred_path=home)
    selected = _select_catalog_models(
        catalog,
        models=models,
        all_models=all_models,
        high_confidence=high_confidence,
    )
    if not selected:
        return ()
    files = _zenodo_files()
    model_dir = home / "pretrained_models"
    manifest = _read_model_manifest(model_dir)
    for model in selected:
        _fetch_one_model(model, files, manifest)
        _write_model_manifest(model_dir, manifest)
    refreshed = {model.name: model for model in model_catalog(mixtcrpred_path=home)}
    return tuple(refreshed[model.name] for model in selected)


def print_model_catalog(
        mixtcrpred_path=None,
        data_dir=None,
        downloaded=False,
        high_confidence=False,
        json_output=False):
    """Print the model catalog for ``mhctools ls mixtcrpred --models``."""
    models = model_catalog(mixtcrpred_path, data_dir=data_dir)
    if downloaded:
        models = tuple(model for model in models if model.status == "ready")
    if high_confidence:
        models = tuple(model for model in models if model.high_confidence)
    if json_output:
        print(json.dumps([model.to_dict() for model in models], indent=2))
        return
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
    widths = [
        max([len(headers[index])] + [len(row[index]) for row in rows])
        for index in range(len(headers))
    ]
    print("  ".join(
        value.ljust(widths[index]) for index, value in enumerate(headers)))
    for row in rows:
        print("  ".join(
            value.ljust(widths[index]) for index, value in enumerate(row)))


def _resolve_python(mixtcrpred_python=None):
    value = mixtcrpred_python or os.environ.get("MIXTCRPRED_PYTHON")
    if value:
        path = Path(value).expanduser()
        if path.is_file():
            # Do not resolve a virtualenv's interpreter symlink: invoking the
            # target binary directly would discard the environment prefix.
            return str(path.absolute())
        executable = shutil.which(str(value))
        if executable:
            return executable
        raise FileNotFoundError("MixTCRpred Python does not exist: %s" % value)
    return sys.executable


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
            raise ValueError(
                "MixTCRpred %s contains unsupported residues %s: %r" % (
                    name, ", ".join(invalid), sequence))


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
            high_confidence=False):
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
        )

    @classmethod
    def catalog(cls, mixtcrpred_path=None, data_dir=None):
        """Return the available pMHC model catalog."""
        return model_catalog(mixtcrpred_path, data_dir=data_dir)

    @classmethod
    def resolve_model(
            cls, peptide, allele, mixtcrpred_path=None, data_dir=None):
        """Resolve a unique catalog model from its peptide/MHC target."""
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

    def __init__(
            self,
            model,
            mixtcrpred_path=None,
            mixtcrpred_python=None,
            checkpoint_path=None,
            batch_size=256):
        self.mixtcrpred_home = _resolve_home(mixtcrpred_path)
        self.mixtcrpred_python = _resolve_python(mixtcrpred_python)
        catalog = model_catalog(mixtcrpred_path=self.mixtcrpred_home)
        selected = _select_catalog_models(catalog, models=[model])
        self.model_info = selected[0]
        self.checkpoint_path = str(
            Path(checkpoint_path).expanduser().resolve()
            if checkpoint_path else self.model_info.path)
        if not Path(self.checkpoint_path).is_file():
            raise FileNotFoundError(
                "MixTCRpred checkpoint is missing: %s. Run `mhctools fetch "
                "mixtcrpred --model %s --accept-license`." % (
                    self.checkpoint_path, self.model_info.name))
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.batch_size = int(batch_size)
        self.last_qc = pd.DataFrame()

    def __repr__(self):
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
        return [
            Prediction(
                kind=Kind.pMHC_TCR_binding,
                score=float(row.score),
                peptide=self.model_info.peptide,
                allele=self.model_info.allele,
                tcr=tcr.identifier,
                percentile_rank=float(row.percentile_rank),
                predictor_name=self._predictor_name(),
                predictor_version=self.predictor_version,
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
