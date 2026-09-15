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

"""Wrapper for PeptiVerse's serum half-life endpoint.

PeptiVerse is a multi-property peptide platform; this wrapper deliberately
exposes **one** of its endpoints — parent-peptide half-life in human serum, in
hours — and emits ``Kind.peptide_half_life`` with that matrix in its context.

That is not ``Kind.pMHC_stability``. NetMHCstabpan measures how long an
assembled peptide-MHC complex holds together; this measures how long the free
peptide survives in serum before proteases degrade it. Different molecule,
different assay, different matrix. It is also not a cleavage-site map: a single
number per peptide says nothing about *where* a peptide gets cut, and must never
be turned into per-bond probabilities.

Sequence input only
-------------------
Upstream also offers SMILES half-life models (``chemberta`` / ``peptideclm``
embeddings), and this wrapper does not use them, for two reasons. First,
``inference.py`` applies the ``expm1`` inverse transform only when
``col == "wt"`` and the model name contains ``log``, so the SMILES models return
an untransformed number that is *not* on the hours scale the sequence model
reports — pooling them into one field would mix units. Second, this adapter
does not consume the exact chemical form now carried by ``PeptideInput``.
Modified inputs are therefore rejected rather than reduced to sequence.

Installation
------------
Snapshot the upstream model repository and the pinned ESM2 model, then point
``PEPTIVERSE_HOME`` and ``PEPTIVERSE_ESM_HOME`` at them::

    git clone https://huggingface.co/ChatterjeeLab/PeptiVerse
    cd PeptiVerse && git checkout 8cf0b21dae356278ae96b414a088e4360357d16c
    huggingface-cli download facebook/esm2_t33_650M_UR50D \
        --revision 08e4846e537177426273712802403f7ba8261b6c \
        --local-dir /models/esm2_t33_650M_UR50D

The upstream stack (torch, transformers==4.46.0, xgboost, lightning) is kept out
of the mhctools environment: inference runs offline in a subprocess under
``PEPTIVERSE_PYTHON``. The wrapper supplies the local ESM2 snapshot and prevents
upstream from constructing its unused PeptideCLM and ChemBERTa embedders.

Provenance and limits
---------------------
Upstream: https://huggingface.co/ChatterjeeLab/PeptiVerse
Cite: https://pmc.ncbi.nlm.nih.gov/articles/PMC12773018/

The retrieved version is a January 2026 preprint. The sequence half-life model
was fit on **130 examples** and evaluated by cross-validation only, with no
external test set — a weak evidence base for ranking decisions, and untested on
the long peptides used in vaccine constructs. Upstream's model card declares
Apache-2.0 while its README declares MIT; both are permissive, but the artifacts
carry no single unambiguous grant. The separate UI Space is CC-BY-NC-ND and is
not used here.

Checkpoints load through ``torch.load(..., weights_only=False)``, which executes
pickled code. Point ``PEPTIVERSE_HOME`` only at a snapshot you trust.
"""

import os
import math
from pathlib import Path
import shutil
import sys
import tempfile

import pandas as pd

from .optional_backend import (
    BackendSpec,
    backend_inventory,
    inspect_artifact,
    prediction_cache_key,
    run_python_sidecar,
)
from .peptide_input import coerce_peptide_inputs, sequence_only_chemistry_error
from .pred import Kind, MeasurementContext, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor

#: Upstream revision this wrapper was written and verified against.
UPSTREAM_REVISION = "8cf0b21dae356278ae96b414a088e4360357d16c"
ESM2_REVISION = "08e4846e537177426273712802403f7ba8261b6c"

#: WTEmbedder allows 1022 tokens, including the two ESM special tokens.
PEPTIVERSE_MAX_PEPTIDE_LENGTH = 1020

_VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# A one-row manifest in upstream's format. Restricting it to Halflife (with the
# SMILES column blanked out with "-") makes `_load_all_best_models` load the
# single sequence half-life model rather than every property's best model.
# Name the exact variant. Upstream's generic "Transformer" alias silently falls
# back from transformer_wt_log to transformer_wt when the expected directory is
# absent, which changes the output transform and invalidates provenance.
_HALF_LIFE_MANIFEST = (
    "Properties, Best_Model_WT, Best_Model_SMILES, Type, "
    "Threshold_WT, Threshold_SMILES,\n"
    "Halflife, Transformer_WT_Log, -, Regression, -, -,\n"
)

_MODEL_DIRECTORY = Path("training_classifiers/half_life/transformer_wt_log")

_PEPTIVERSE_ARTIFACTS = {
    "inference.py": (
        "inference_and_feature_code",
        "b899369f020cc73f939a7113fd9aff471a923b109e77ae88c22aa706e66d119a",
        "python_source"),
    str(_MODEL_DIRECTORY / "best_model.pt"): (
        "half_life_model_weights",
        "d05a90b794ed45246cf2778e677f5ef312db90131c0b39367b3102fde59b09fa",
        "torch_pickle"),
    str(_MODEL_DIRECTORY / "best_params.json"): (
        "model_configuration",
        "b5fc9465638ad21e554310e8a03650c12dc7ea46febac7c75dea34251fa96b84",
        "json"),
    str(_MODEL_DIRECTORY / "mapie_calibration.joblib"): (
        "uncertainty_calibration",
        "bc801e894da7d65c906d8c5d4a14dd680d6ca1ed28eb681ed01666a8fc34b66e",
        "joblib_pickle"),
}

_ESM2_ARTIFACTS = {
    "config.json": (
        "embedding_model_configuration",
        "539095c22efc52a09d6147074ba4ca119f76a890df5901213b2b55f7d2f96b2b",
        "json"),
    "tokenizer_config.json": (
        "tokenizer_configuration",
        "7e9161ecdb548ec45a41cbc6b24aa4476fdd418461f491c4207baa99419a29ad",
        "json"),
    "special_tokens_map.json": (
        "tokenizer_configuration",
        "3aedcd4211c0d43aec4e607ff60a63255f3174ead795e997350f09a5f8cd9ee1",
        "json"),
    "vocab.txt": (
        "tokenizer_vocabulary",
        "0b82cc0a7c7cf9e567b1e5892d793285b9fbae822c964ca48696f7db44598e03",
        "text"),
}

_ESM2_WEIGHT_ARTIFACTS = {
    "model.safetensors": (
        "a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0",
        "safetensors"),
    "pytorch_model.bin": (
        "c874668852c7275a159e2c7ceb6069671d7b1ba2c7b52f59600b34ce0f721008",
        "torch_pickle"),
}

PEPTIVERSE_BACKEND_SPEC = BackendSpec(
    name="peptiverse",
    endpoint="human_serum_half_life_hours",
    developed_against="peptiverse@%s+esm2@%s" % (
        UPSTREAM_REVISION, ESM2_REVISION),
    license="PeptiVerse: Apache-2.0/MIT ambiguity; ESM2: MIT",
    serialization="torch pickle checkpoint + safetensors/torch weights",
    entry_point="prediction_only",
    # The conformance path is covered on Linux/macOS and Python 3.10-3.12, but
    # real-model inference remains an opt-in smoke test, so do not advertise a
    # validated runtime combination here yet.
    supported_platforms=(),
    supported_interpreters=(),
)


def _find_peptiverse_home(peptiverse_home=None):
    """Resolve the PeptiVerse snapshot directory.

    Checks, in order: the *peptiverse_home* argument, ``$PEPTIVERSE_HOME``,
    then ``~/PeptiVerse``.
    """
    candidate = peptiverse_home or os.environ.get("PEPTIVERSE_HOME")
    if not candidate:
        home = Path.home() / "PeptiVerse"
        if home.is_dir():
            candidate = str(home)
    if not candidate:
        raise FileNotFoundError(
            "PeptiVerse not found. Set PEPTIVERSE_HOME or pass "
            "peptiverse_home= to the constructor. Clone from "
            "https://huggingface.co/ChatterjeeLab/PeptiVerse")
    candidate = str(Path(candidate).expanduser().resolve())
    if not Path(candidate, "inference.py").is_file():
        raise FileNotFoundError(
            "inference.py not found in %r — is this a PeptiVerse snapshot?"
            % candidate)
    weights = Path(candidate) / _MODEL_DIRECTORY
    if not weights.is_dir():
        raise FileNotFoundError(
            "%s not found in %r — the exact log-scale half-life model is "
            "required; mhctools will not fall back to another variant"
            % (_MODEL_DIRECTORY, candidate))
    return candidate


def _find_esm_home(peptiverse_home, peptiverse_esm_home=None):
    """Resolve a local snapshot of the pinned ESM2 feature model."""
    candidate = peptiverse_esm_home or os.environ.get("PEPTIVERSE_ESM_HOME")
    if not candidate:
        colocated = Path(peptiverse_home) / "esm2_t33_650M_UR50D"
        if colocated.is_dir():
            candidate = colocated
    if not candidate:
        cached = (
            Path.home() / ".cache" / "huggingface" / "hub" /
            "models--facebook--esm2_t33_650M_UR50D" / "snapshots" /
            ESM2_REVISION)
        if cached.is_dir():
            candidate = cached
    if not candidate:
        raise FileNotFoundError(
            "Pinned ESM2 snapshot not found. Set PEPTIVERSE_ESM_HOME or pass "
            "peptiverse_esm_home=; inference is offline and will not download it")
    return str(Path(candidate).expanduser().resolve())


def _artifact_inventory(
        peptiverse_home, esm_home, max_peptide_length, device=None):
    artifacts = []
    root = Path(peptiverse_home)
    for relative, (role, expected_sha256, serialization) in (
            _PEPTIVERSE_ARTIFACTS.items()):
        artifacts.append(inspect_artifact(
            name="peptiverse/%s" % relative,
            role=role,
            path=root / relative,
            expected_sha256=expected_sha256,
            serialization=serialization))

    esm_root = Path(esm_home)
    for relative, (role, expected_sha256, serialization) in (
            _ESM2_ARTIFACTS.items()):
        artifacts.append(inspect_artifact(
            name="esm2/%s" % relative,
            role=role,
            path=esm_root / relative,
            expected_sha256=expected_sha256,
            serialization=serialization))

    # Transformers prefers safetensors when both formats are present. Inventory
    # only the file it will actually load, while accepting the pinned PyTorch
    # alternative when safetensors was not provisioned.
    weight_name = "model.safetensors"
    if not (esm_root / weight_name).is_file() and (
            esm_root / "pytorch_model.bin").is_file():
        weight_name = "pytorch_model.bin"
    expected_sha256, serialization = _ESM2_WEIGHT_ARTIFACTS[weight_name]
    artifacts.append(inspect_artifact(
        name="esm2/%s" % weight_name,
        role="embedding_model_weights",
        path=esm_root / weight_name,
        expected_sha256=expected_sha256,
        serialization=serialization))

    return backend_inventory(
        spec=PEPTIVERSE_BACKEND_SPEC,
        artifacts=artifacts,
        settings={
            "embedding_model": "facebook/esm2_t33_650M_UR50D",
            "device": device or "auto",
            "max_peptide_length": max_peptide_length,
            "model_variant": "transformer_wt_log",
            "output_units": "hours",
        })


def _resolve_python(peptiverse_python=None):
    value = peptiverse_python or os.environ.get("PEPTIVERSE_PYTHON")
    if not value:
        return sys.executable
    path = Path(value).expanduser()
    if path.is_file():
        # Do not resolve a virtualenv interpreter symlink: invoking the target
        # binary directly would discard the environment prefix.
        return str(path.absolute())
    executable = shutil.which(str(value))
    if executable:
        return executable
    raise FileNotFoundError("PeptiVerse Python does not exist: %s" % value)


class PeptiVerse(AlleleFreePredictor):
    """Serum half-life predictions from a local PeptiVerse snapshot.

    Parameters
    ----------
    peptiverse_home : str, optional
        Path to a PeptiVerse snapshot (holds ``inference.py`` and
        ``training_classifiers/``). Resolved from the argument, then
        ``$PEPTIVERSE_HOME``, then ``~/PeptiVerse``.
    peptiverse_python : str, optional
        Interpreter that has PeptiVerse's dependencies. Resolved from the
        argument, then ``$PEPTIVERSE_PYTHON``, then the current interpreter.
    peptiverse_esm_home : str, optional
        Local snapshot of ``facebook/esm2_t33_650M_UR50D`` at the pinned
        revision. Resolved from the argument, ``$PEPTIVERSE_ESM_HOME``, a
        directory colocated under ``PEPTIVERSE_HOME``, then the exact revision
        in the HuggingFace cache. It is never downloaded during prediction.
    device : str, optional
        Passed to upstream (``"cpu"``, ``"cuda"``, ...). Default lets upstream
        choose: CUDA when available, otherwise CPU.
    uncertainty : bool
        Ask upstream for its conformal interval alongside each prediction, kept
        in :attr:`last_qc` with upstream's own status string. Off by default and
        never surfaced in the returned ``Prediction``, for two reasons.

        In practice it does not load: the shipped
        ``half_life/transformer_wt_log/mapie_calibration.joblib`` unpickles a
        ``__main__.PassthroughRegressor``, a class defined in the training
        script rather than in any importable module, so upstream reports
        ``"unavailable (no MAPIE bundle and no seed ensemble)"``.

        And were it loadable, the bounds would not be on the hours scale the
        point estimate uses: upstream calls ``_compute_uncertainty`` with the
        already-``expm1``-ed score and adds a quantile calibrated on log-space
        residuals. Treat anything that appears here as a raw upstream
        diagnostic, not a calibrated confidence interval.
    max_peptide_length : int
        Peptides longer than this are rejected instead of being silently
        truncated by ESM2. Default 1020 (1022 tokens minus two special tokens).
    allow_unverified_assets : bool
        Permit explicitly supplied files whose checksums differ from the pinned
        inventory. False by default. The actual content identity is still
        recorded in :attr:`artifact_inventory` and every prediction.
    subprocess_timeout : float, optional
        Maximum seconds allowed for the isolated inference process.
    """

    mhc_class = "none"

    def __init__(
            self,
            peptiverse_home=None,
            peptiverse_python=None,
            peptiverse_esm_home=None,
            device=None,
            uncertainty=False,
            max_peptide_length=PEPTIVERSE_MAX_PEPTIDE_LENGTH,
            allow_unverified_assets=False,
            subprocess_timeout=3600):
        if max_peptide_length > PEPTIVERSE_MAX_PEPTIDE_LENGTH:
            raise ValueError(
                "ESM2 accepts at most %d residues; max_peptide_length cannot "
                "exceed that (got %d)"
                % (PEPTIVERSE_MAX_PEPTIDE_LENGTH, max_peptide_length))
        self.peptiverse_home = _find_peptiverse_home(peptiverse_home)
        self.peptiverse_esm_home = _find_esm_home(
            self.peptiverse_home, peptiverse_esm_home)
        self.peptiverse_python = _resolve_python(peptiverse_python)
        self.device = device
        self.uncertainty = uncertainty
        self.max_peptide_length = max_peptide_length
        self.subprocess_timeout = subprocess_timeout
        self.artifact_inventory = _artifact_inventory(
            self.peptiverse_home,
            self.peptiverse_esm_home,
            self.max_peptide_length,
            self.device)
        self.artifact_inventory.require_usable(
            allow_unverified=allow_unverified_assets)
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return "PeptiVerse(peptiverse_home=%r, device=%r)" % (
            self.peptiverse_home, self.device)

    def _default_pred_kind(self):
        return Kind.peptide_half_life

    def _predictor_name(self):
        return "peptiverse"

    @property
    def predictor_version(self):
        return self.artifact_inventory.predictor_version

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if not peptide:
                raise ValueError("Empty peptide is not allowed")
            if len(peptide) > self.max_peptide_length:
                raise ValueError(
                    "PeptiVerse supports peptides up to %d residues; got %r "
                    "(length %d)"
                    % (self.max_peptide_length, peptide, len(peptide)))
            invalid = set(peptide) - _VALID_AMINO_ACIDS
            if invalid:
                raise ValueError(
                    "Peptide %r contains characters the sequence model cannot "
                    "encode: %s. PeptiVerse's chemical-form (SMILES) models are "
                    "not wrapped, so modified peptides are rejected rather than "
                    "scored as their unmodified sequence."
                    % (peptide, "".join(sorted(invalid))))

    def _input_error(self, peptide_input):
        error = sequence_only_chemistry_error(peptide_input)
        if error:
            return error
        try:
            self._check_peptides([peptide_input.sequence])
        except ValueError as error:
            return str(error)
        return None

    def _measurement_context(self, status="available", detail=None):
        return MeasurementContext(
            estimate_type="ml_predicted",
            status=status,
            analyte="parent peptide",
            matrix="human serum",
            unit="hours" if status == "available" else None,
            transform="linear" if status == "available" else None,
            score_semantics=(
                "same linear half-life in hours as value"
                if status == "available" else None),
            detail=detail,
        )

    def predict(self, peptides, on_unsupported="raise"):
        """Predict serum half-life for a list of peptides.

        Inputs may be strings (the backward-compatible shorthand for a
        canonical unmodified L-peptide with free termini) or exact
        :class:`~mhctools.peptide_input.PeptideInput` records. Modified forms
        are never reduced to sequence. By default the first unsupported input
        raises; ``on_unsupported="record"`` returns an explicit unavailable
        prediction in that position and scores the supported remainder.

        Returns
        -------
        list of PeptideResult
            One entry per input peptide, in input order, each holding a single
            ``Kind.peptide_half_life`` prediction with an empty ``allele``.
            Its measurement context identifies human serum. Both
            ``score`` and ``value`` carry the predicted half-life in **hours**
            (higher = longer-lived); ``value`` is the units-bearing field and
            ``score`` repeats it so that rank-based consumers work without
            knowing the unit.
        """
        if on_unsupported not in ("raise", "record"):
            raise ValueError("on_unsupported must be 'raise' or 'record'")
        peptide_inputs = coerce_peptide_inputs(peptides)
        errors = [self._input_error(item) for item in peptide_inputs]
        if on_unsupported == "raise" and any(errors):
            index = next(i for i, error in enumerate(errors) if error)
            raise ValueError(f"Input {index} is unsupported: {errors[index]}")
        supported = [
            item for item, error in zip(peptide_inputs, errors) if not error]
        if not peptide_inputs:
            self.last_qc = pd.DataFrame()
            return []

        if supported:
            output = self._run_sidecar(
                [item.sequence for item in supported])
        else:
            self.last_qc = pd.DataFrame()
            output = pd.DataFrame(columns=["hours"])
        hours = iter(output["hours"])
        results = []
        for peptide_input, error in zip(peptide_inputs, errors):
            common = {
                "kind": Kind.peptide_half_life,
                "peptide": peptide_input.sequence,
                "predictor_name": self._predictor_name(),
                "predictor_version": self.predictor_version,
                "peptide_input": peptide_input,
                "cache_key": prediction_cache_key(
                    peptide_input, self.artifact_inventory),
            }
            if error:
                prediction = Prediction(
                    score=None,
                    measurement_context=self._measurement_context(
                        status="unsupported", detail=error),
                    **common)
            else:
                value = next(hours)
                prediction = Prediction(
                    score=value,
                    value=value,
                    measurement_context=self._measurement_context(),
                    **common)
            results.append(PeptideResult(preds=(prediction,)))
        return results

    def _run_sidecar(self, peptide_list):
        with tempfile.TemporaryDirectory(prefix="mhctools_peptiverse_") as tmp:
            tmp_path = Path(tmp)
            input_path = tmp_path / "input.csv"
            output_path = tmp_path / "output.csv"
            manifest_path = tmp_path / "half_life_only_manifest.txt"
            manifest_path.write_text(_HALF_LIFE_MANIFEST)
            pd.DataFrame({
                "__mhctools_id": range(len(peptide_list)),
                "peptide": peptide_list,
            }).to_csv(input_path, index=False)

            sidecar = Path(__file__).with_name("peptiverse_sidecar.py")
            arguments = [
                "--home", self.peptiverse_home,
                "--esm-home", self.peptiverse_esm_home,
                "--manifest", str(manifest_path),
                "--input", str(input_path),
                "--output", str(output_path),
                "--device", self.device or "",
            ]
            if self.uncertainty:
                arguments.append("--uncertainty")

            run_python_sidecar(
                backend_name="PeptiVerse",
                python=self.peptiverse_python,
                sidecar=sidecar,
                args=arguments,
                cwd=self.peptiverse_home,
                timeout=self.subprocess_timeout,
            )
            output = parse_peptiverse_results(output_path, peptide_list)

        self.artifact_inventory = (
            self.artifact_inventory.with_inference_reproduced())
        self.last_qc = output.drop(columns=["hours"])
        return output


def parse_peptiverse_results(filename, peptide_list):
    """Parse sidecar output into a DataFrame aligned with *peptide_list*.

    Every input peptide must come back exactly once, identified by the
    ``__mhctools_id`` the caller assigned, and with the peptide it was scored
    under unchanged. This checks the echoed input identity; the sidecar checks
    the actual tokenized residue count before inference. Duplicate input
    peptides keep separate rows.

    Returns
    -------
    pandas.DataFrame
        Sorted by ``__mhctools_id``, with ``hours`` as float.
    """
    output = pd.read_csv(filename, keep_default_na=False)
    for column in ("__mhctools_id", "peptide", "hours"):
        if column not in output.columns:
            raise ValueError(
                "PeptiVerse output missing column %r; got %s"
                % (column, list(output.columns)))
    output["__mhctools_id"] = output["__mhctools_id"].astype(int)
    if sorted(output["__mhctools_id"].tolist()) != list(range(len(peptide_list))):
        raise RuntimeError(
            "PeptiVerse output did not preserve every input peptide: expected "
            "%d rows with ids 0..%d, got ids %s"
            % (len(peptide_list), len(peptide_list) - 1,
               sorted(output["__mhctools_id"].tolist())))
    output = output.sort_values("__mhctools_id").reset_index(drop=True)
    returned = output["peptide"].tolist()
    if returned != list(peptide_list):
        mismatched = [
            (i, sent, got)
            for i, (sent, got) in enumerate(zip(peptide_list, returned))
            if sent != got
        ]
        raise RuntimeError(
            "PeptiVerse returned a different peptide than it was given "
            "(id, sent, got): %s" % mismatched[:5])
    output["hours"] = output["hours"].astype(float)
    invalid = output["hours"].map(lambda value: not math.isfinite(value) or value < 0)
    if invalid.any():
        rows = output.loc[invalid, ["__mhctools_id", "peptide", "hours"]]
        raise RuntimeError(
            "PeptiVerse returned invalid half-lives (expected finite, "
            "nonnegative hours): %s" % rows.to_dict("records"))
    return output
