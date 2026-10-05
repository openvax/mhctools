# Licensed under the Apache License, Version 2.0 (the "License");
# https://www.apache.org/licenses/LICENSE-2.0

"""CleaveNet's native whole-substrate MMP scores, separate from bond cleavage."""

from dataclasses import asdict, dataclass, replace
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
from typing import Tuple

import pandas as pd

from .optional_backend import (
    BackendInventory, BackendSpec, backend_inventory, inspect_artifact,
    prediction_cache_key, run_python_sidecar,
)
from .peptide_input import PeptideInput, sequence_only_chemistry_error
from .pred import Kind

UPSTREAM_REVISION = "4dac67defc99ca35d967ddc76eca0fe8b74afdad"
PAPER = "https://doi.org/10.1038/s41467-025-67226-1"
ENZYMES = (
    "MMP1", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16",
    "MMP17", "MMP19", "MMP2", "MMP20", "MMP24", "MMP25", "MMP3", "MMP7",
    "MMP8", "MMP9",
)
ENDPOINT = "whole-substrate relative cleavage Z-score"
ASSAY = "mRNA-display substrate cleavage screen; isolated MMPs"

CLEAVENET_BACKEND_SPEC = BackendSpec(
    name="cleavenet", endpoint=ENDPOINT,
    developed_against="microsoft/cleavenet@" + UPSTREAM_REVISION,
    license="code: MIT; data: CDLA-Permissive-2.0",
    serialization="HDF5 weights; Python source; CSV vocabulary splits",
    entry_point="prediction_only",
    supported_platforms=("macos-arm64-cpu", "linux-x86_64-cpu"),
    supported_interpreters=("python3.11", "python3.12"),
)

# Exact runtime inventory at UPSTREAM_REVISION; no generator weights are loaded.
_ARTIFACTS = {'cleavenet/__init__.py': ('runtime_source',
                           '2be443e02a82bcc74a01945903529aa2c61148e607a51456338d4c7cc6f75aa4',
                           'python_source'),
 'cleavenet/analysis.py': ('runtime_source',
                           'd896991e8461571bad17e7d6b20f17ea7645f63051e4e344af15bc6b399d2011',
                           'python_source'),
 'cleavenet/data.py': ('runtime_source',
                       'b512de422da06f5a3c872750445164640c818a15265cb94f9c89aaa991e2567b',
                       'python_source'),
 'cleavenet/models.py': ('runtime_source',
                         '3c77fafc9447db34079c122002d81542f158b391adbbbc603884a978d50da36a',
                         'python_source'),
 'cleavenet/plotter.py': ('runtime_source',
                          '6db3d72a5f214caf16cad3e024cfba36e432285170d5cf4cd2083bdd45aa3c32',
                          'python_source'),
 'cleavenet/utils.py': ('runtime_source',
                        '0eb90577f6f9a23ee7f0a5a20fc72ffa604a90e9d878c0666ca5518394aaa2b1',
                        'python_source'),
 'splits/kukreja/X_all.csv': ('vocabulary_split',
                              '31585e83286eceb5784c0b0c69ac93f39c37cc581c1aa4803705532d3ab157ba',
                              'csv'),
 'splits/kukreja/y_all.csv': ('vocabulary_split',
                              'c5c772aba6574e5155d85a3573108e572a22e1254c5d82aa25b3a981909f8b09',
                              'csv'),
 'LICENSE': ('runtime_provenance',
             'c2cfccb812fe482101a8f04597dfc5a9991a6b2748266c47ac91b6a5aae15383',
             'text'),
 'data/LICENSE': ('runtime_provenance',
                  '9d242f2775d7a1d61249e49ed08bfcca3d8759782ad4bfdaa6dcfff830765054',
                  'text'),
 'requirements.txt': ('runtime_provenance',
                      'c500243da706c323d6860975c5ab90946355ec0a89bcf4bc7f584598a0db20d0',
                      'text'),
 'weights/transformer_0/model.h5': ('transformer_checkpoint',
                                    'f257f8090dbd2ffb123955539dabdb8ac8aef3b043187182b2849d9557ca4d5f',
                                    'HDF5_weights'),
 'weights/transformer_1/model.h5': ('transformer_checkpoint',
                                    '0c726c98f9064be2023d95add22bf96bfc79b999efd5e50994e992280cc983e9',
                                    'HDF5_weights'),
 'weights/transformer_2/model.h5': ('transformer_checkpoint',
                                    '9021e505fedcc5d6743b82a12840813942db43a9991345ae12291b034dd6590b',
                                    'HDF5_weights'),
 'weights/transformer_3/model.h5': ('transformer_checkpoint',
                                    '7a6843ea23145256bd1f6377f9a0ba83e112e5a4def762d5177f196730e4ca42',
                                    'HDF5_weights'),
 'weights/transformer_4/model.h5': ('transformer_checkpoint',
                                    'acd91ea0dffe860b6ebb46751873214a7b02374e5c89a305f69800d411a5ff2d',
                                    'HDF5_weights')}


@dataclass(frozen=True)
class CleaveNetScore:
    """One enzyme's native Z-score and spread across the five models.

    The standard deviation is not a calibrated interval or a probability.
    """

    enzyme: str
    z_score: float
    ensemble_sd: float

    def __post_init__(self):
        if self.enzyme not in ENZYMES:
            raise ValueError("Unknown CleaveNet enzyme: %s" % self.enzyme)
        for name in ("z_score", "ensemble_sd"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError("%s must be a finite number" % name)
        if self.ensemble_sd < 0:
            raise ValueError("ensemble_sd cannot be negative")


@dataclass(frozen=True)
class CleaveNetResult:
    """Whole-substrate evidence preserving input, padding and model identity.

    Scores do not locate a scissile bond or quantify intact-peptide survival.
    Source coordinates identify the supplied window, never a cleavage site.
    """

    peptide_input: PeptideInput
    padded_sequence: str
    scores: Tuple[CleaveNetScore, ...]
    inventory: BackendInventory
    cache_key: str
    runtime: Tuple[Tuple[str, str], ...]

    @property
    def applicability(self):
        return ("ten-residue input" if len(self.peptide_input.sequence) == 10
                else "short input with centered padding; limited validation")

    def to_dict(self):
        """Return JSON-compatible evidence without adding per-bond labels."""
        return dict(
            peptide_input=self.peptide_input.to_dict(),
            padded_sequence=self.padded_sequence,
            scores=[asdict(score) for score in self.scores],
            kind=Kind.substrate_cleavage, endpoint=ENDPOINT, score_units="dimensionless Z-score", assay=ASSAY, reference=PAPER,
            assay_conditions="model uses sequence only; no enzyme dose, incubation, pH or matrix inputs",
            training_source="Kukreja et al. 2015 mRNA-display data; pinned splits/kukreja",
            validation="official-source conformance only; no independently audited held-out benchmark",
            applicability=self.applicability,
            predictor_version=self.inventory.predictor_version,
            inventory=self.inventory.to_dict(), cache_key=self.cache_key,
            runtime=dict(self.runtime))


class CleaveNet:
    """Optional local inference for canonical sequences of 1–10 residues.

    Parameters
    ----------
    cleavenet_home : str, optional
        Pinned upstream source/weights, or ``CLEAVENET_HOME``. Otherwise use
        the managed installation from ``mhctools fetch cleavenet``.
    cleavenet_python : str, optional
        Isolated interpreter containing upstream requirements (TensorFlow
        2.18.0). Defaults to ``CLEAVENET_PYTHON`` or the current interpreter.
    subprocess_timeout : float
        Maximum inference duration in seconds.

    Notes
    -----
    This adapter returns dedicated substrate results rather than site tracks.
    The assay endpoint is not serum half-life, uptake or antigen presentation.
    Construction verifies every runtime file against the pinned checksums.
    """

    @classmethod
    def fetch(cls, version=None, data_dir=None, accept_license=False):
        """Fetch the pinned prediction assets and both license agreements."""
        from .artifacts import fetch
        return fetch("cleavenet", version=version, data_dir=data_dir,
                     accept_license=accept_license)

    def __init__(self, cleavenet_home=None, cleavenet_python=None, subprocess_timeout=300):
        from .artifacts import artifact_status
        candidate = cleavenet_home or os.environ.get("CLEAVENET_HOME")
        if not candidate:
            status = artifact_status("cleavenet")
            if status.status == "ready":
                candidate = status.path
        if not candidate:
            raise FileNotFoundError("Run `mhctools fetch cleavenet` or set CLEAVENET_HOME")
        self.home = Path(candidate).expanduser().resolve()
        interpreter = cleavenet_python or os.environ.get("CLEAVENET_PYTHON") or sys.executable
        path = Path(interpreter).expanduser()
        self.python = str(path.absolute()) if path.is_file() else shutil.which(str(interpreter))
        if not self.python:
            raise FileNotFoundError("CleaveNet interpreter not found: %s" % interpreter)
        self.subprocess_timeout = subprocess_timeout
        self.artifact_inventory = backend_inventory(
            CLEAVENET_BACKEND_SPEC,
            [inspect_artifact(name=name, role=role, path=self.home / name,
                              expected_sha256=digest, serialization=serialization)
             for name, (role, digest, serialization) in _ARTIFACTS.items()],
            settings={"architecture": "transformer", "ensemble_size": 5,
                      "padding": "centered to 10; odd remainder on right",
                      "tensorflow": "2.18.0", "ensemble_sd_ddof": 0})
        self.artifact_inventory.require_usable()

    def kind_support(self):
        """Native whole-substrate endpoint and its MHC independence."""
        return {Kind.substrate_cleavage: {"mhc_dependence": "none", "mhc_class": "none"}}

    @property
    def supported_kinds(self):
        return tuple(self.kind_support())

    def _run_sidecar(self, sequences):
        with tempfile.TemporaryDirectory(prefix="mhctools-cleavenet-") as directory:
            input_path, output_path = [Path(directory) / name for name in ("input.json", "output.json")]
            input_path.write_text(json.dumps(sequences))
            environment = dict(os.environ)
            environment.pop("TF_USE_LEGACY_KERAS", None)
            environment.update(TF_NUM_INTEROP_THREADS="1", TF_NUM_INTRAOP_THREADS="1")
            run_python_sidecar(
                "CleaveNet", self.python, Path(__file__).with_name("cleavenet_sidecar.py"),
                [self.home, input_path, output_path], timeout=self.subprocess_timeout,
                environment=environment)
            output = json.loads(output_path.read_text())
        if output.get("sequences") != sequences or output.get("enzymes") != list(ENZYMES):
            raise RuntimeError("CleaveNet input identity or enzyme order mismatch")
        for name in ("means", "ensemble_sd"):
            rows = output.get(name)
            if not isinstance(rows, list) or len(rows) != len(sequences) or any(
                    not isinstance(row, list) or len(row) != len(ENZYMES) for row in rows):
                raise RuntimeError("CleaveNet %s output shape mismatch" % name)
        if output.get("runtime", {}).get("tensorflow") != "2.18.0":
            raise RuntimeError("CleaveNet runtime version mismatch")
        return output

    def predict(self, peptides):
        """Return ordered substrate results, preserving duplicate occurrences.

        Accept strings or exact ``PeptideInput`` records. Reject noncanonical
        sequences, longer inputs and unsupported chemistry before inference.
        Short inputs are centered in ten positions with the odd gap on the
        right, following upstream's README. Empty batches return an empty list.
        """
        if isinstance(peptides, (str, PeptideInput)):
            peptides = [peptides]
        inputs = [PeptideInput(item) if isinstance(item, str) else item for item in peptides]
        for item in inputs:
            if not isinstance(item, PeptideInput):
                raise TypeError("Expected PeptideInput or canonical peptide string")
            error = sequence_only_chemistry_error(item)
            if error:
                raise ValueError(error)
            if len(item.sequence) > 10:
                raise ValueError("CleaveNet supports 1–10 residues; use predict_windows for longer inputs")
        if not inputs:
            return []
        sequences = ["-" * ((10 - len(item.sequence)) // 2) + item.sequence
                     + "-" * ((11 - len(item.sequence)) // 2) for item in inputs]
        output = self._run_sidecar(sequences)
        results = []
        for item, padded, means, deviations in zip(inputs, sequences, output["means"], output["ensemble_sd"]):
            scores = tuple(CleaveNetScore(enzyme, mean, deviation)
                           for enzyme, mean, deviation in zip(ENZYMES, means, deviations))
            results.append(CleaveNetResult(
                item, padded, scores, self.artifact_inventory.with_inference_reproduced(),
                prediction_cache_key(item, self.artifact_inventory),
                tuple(sorted(output["runtime"].items()))))
        self.artifact_inventory = self.artifact_inventory.with_inference_reproduced()
        return results

    def predict_windows(self, peptide):
        """Score every ten-residue window with zero-based source coordinates.

        These conditional sequence scores do not assert that the window is
        physically released or that an internal/central bond is cleaved.
        """
        item = PeptideInput(peptide) if isinstance(peptide, str) else peptide
        if not isinstance(item, PeptideInput):
            raise TypeError("Expected PeptideInput or canonical peptide string")
        error = sequence_only_chemistry_error(item)
        if error:
            raise ValueError(error)
        if len(item.sequence) < 10:
            raise ValueError("Window scanning requires at least 10 residues")
        return self.predict([
            replace(item, sequence=item.sequence[start:start + 10],
                    source_start=item.source_start + start)
            for start in range(len(item.sequence) - 9)])

    def predict_dataframe(self, peptides):
        """Flatten native substrate evidence to one row per input/enzyme."""
        columns = ["peptide", "source_sequence_name", "source_start", "source_end",
                   "padded_sequence", "enzyme", "z_score", "ensemble_sd", "endpoint",
                   "assay", "applicability", "predictor_version", "cache_key", "reference"]
        rows = []
        for result in self.predict(peptides):
            item = result.peptide_input
            for score in result.scores:
                rows.append(dict(
                    peptide=item.sequence, source_sequence_name=item.source_sequence_name,
                    source_start=item.source_start, source_end=item.source_start + len(item.sequence),
                    padded_sequence=result.padded_sequence, **asdict(score), endpoint=ENDPOINT,
                    assay=ASSAY, applicability=result.applicability,
                    predictor_version=result.inventory.predictor_version,
                    cache_key=result.cache_key, reference=PAPER))
        return pd.DataFrame(rows, columns=columns)
