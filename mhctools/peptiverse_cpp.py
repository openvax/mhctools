# Copyright (c) 2026 Mount Sinai School of Medicine
# Licensed under the Apache License, Version 2.0 (the "License");
# you may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0

"""Optional PeptiVerse canonical CPP classifier, using its released CPU SVC.

Native P(CPP) is a class score, not a fraction entering cells or a prediction
of APC delivery, intracellular localization or productive presentation.
"""

import math
from pathlib import Path
import tempfile

import pandas as pd

from .optional_backend import (
    BackendSpec, backend_inventory, inspect_artifact, prediction_cache_key,
    run_python_sidecar,
)
from .peptide_input import coerce_peptide_inputs, sequence_only_chemistry_error
from .peptiverse import (
    ESM2_REVISION, UPSTREAM_REVISION, PEPTIVERSE_MAX_PEPTIDE_LENGTH,
    _PEPTIVERSE_ARTIFACTS, _esm2_artifacts, _find_esm_home,
    _resolve_peptiverse_source, _resolve_python,
)
from .pred import Kind, MeasurementContext, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor


CPP_THRESHOLD = 0.5493
CPP_SCIKIT_VERSION = "1.7.2"
_CPP_DIRECTORY = Path("training_classifiers/permeability_penetrance/svm_gpu_wt")
_CPP_ARTIFACTS = {
    "inference.py": _PEPTIVERSE_ARTIFACTS["inference.py"],
    str(_CPP_DIRECTORY / "best_model.joblib"): (
        "cpp_classification_weights",
        "087f7c6aa5545eb6ef555c9de10cd2fd5cc84cd44c31601e2ee6ee017a9c7277",
        "joblib_pickle"),
    "best_models.txt": (
        "native_endpoint_thresholds",
        "8a515bd115bcc7a6e4ccd6ecba76ed2b66eaf767e66fa3de5ae1c1111ebb90f9",
        "text"),
}
_CPP_MANIFEST = (
    "Properties, Best_Model_WT, Best_Model_SMILES, Type, Threshold_WT, Threshold_SMILES,\n"
    "Permeability (Penetrance), SVM, -, Classifier, 0.5493, -,\n")
_CPP_SPEC = BackendSpec(
    name="peptiverse-cpp", endpoint="canonical_cpp_classification",
    developed_against=f"peptiverse@{UPSTREAM_REVISION}+esm2@{ESM2_REVISION}",
    license="PeptiVerse: Apache-2.0/MIT ambiguity; ESM2: MIT",
    serialization="joblib pickle SVC + safetensors/torch embedding weights",
    entry_point="prediction_only", supported_platforms=("linux", "macos"),
    supported_interpreters=("Python 3.11; scikit-learn 1.7.2",))
_SCOPE = (
    "Native canonical-sequence CPP classifier; score is P(CPP), not measured "
    "uptake fraction, cell-specific/APC delivery, intracellular localization, "
    "productive antigen presentation or exposure. Source metadata training "
    "lengths 3-61 residues; long-vaccine accuracy unestablished."
)


class PeptiVerseCPP(AlleleFreePredictor):
    """Offline CPP classification in a separate PeptiVerse interpreter.

    Uses the exact CPU sklearn SVC in ``svm_gpu_wt`` and native threshold
    0.5493. The selected interpreter must have scikit-learn 1.7.2. Source,
    interpreter and ESM2 resolution match :class:`PeptiVerse`. Modified
    chemistry is rejected; computational capacity is not biological validity.
    """

    def __init__(
            self, peptiverse_home=None, peptiverse_python=None,
            peptiverse_esm_home=None, device=None, uncertainty=False,
            max_peptide_length=PEPTIVERSE_MAX_PEPTIDE_LENGTH,
            allow_unverified_assets=False, subprocess_timeout=3600):
        if (isinstance(max_peptide_length, bool) or
                not isinstance(max_peptide_length, int) or
                not 1 <= max_peptide_length <= PEPTIVERSE_MAX_PEPTIDE_LENGTH):
            raise ValueError("max_peptide_length must be an integer in [1, 1020]")
        self.peptiverse_home = _resolve_peptiverse_source(peptiverse_home)
        self.peptiverse_esm_home = _find_esm_home(self.peptiverse_home, peptiverse_esm_home)
        self.peptiverse_python = _resolve_python(peptiverse_python)
        self.device, self.uncertainty = device, uncertainty
        self.max_peptide_length = max_peptide_length
        self.subprocess_timeout = subprocess_timeout
        artifacts = [inspect_artifact(
            name="peptiverse/%s" % relative, role=role,
            path=Path(self.peptiverse_home) / relative,
            expected_sha256=sha, serialization=serialization)
            for relative, (role, sha, serialization) in _CPP_ARTIFACTS.items()]
        artifacts.extend(_esm2_artifacts(self.peptiverse_esm_home))
        self.artifact_inventory = backend_inventory(
            spec=_CPP_SPEC, artifacts=artifacts, settings={
                "embedding_model": "facebook/esm2_t33_650M_UR50D",
                "device": device or "auto", "uncertainty": bool(uncertainty),
                "max_peptide_length": max_peptide_length,
                "model_variant": "svm_gpu_wt (CPU sklearn SVC)",
                "threshold": CPP_THRESHOLD, "scikit_learn": CPP_SCIKIT_VERSION,
                "output_semantics": "P(CPP) class score; no physical uptake value",
            })
        self.artifact_inventory.require_usable(allow_unverified=allow_unverified_assets)
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return "PeptiVerseCPP(peptiverse_home=%r)" % self.peptiverse_home

    def _default_pred_kind(self):
        return Kind.cpp_classification

    @property
    def predictor_version(self):
        return self.artifact_inventory.predictor_version

    def predict(self, peptides, on_unsupported="raise"):
        """Return one CPP-class score per exact input, preserving occurrences.

        ``on_unsupported='record'`` retains unsupported chemistry and
        over-capacity inputs as unavailable records. ``value`` is always None.
        Source-domain length limits are applicability flags, not hard cutoffs.
        """
        self.last_qc = pd.DataFrame()
        if on_unsupported not in ("raise", "record"):
            raise ValueError("on_unsupported must be 'raise' or 'record'")
        inputs = coerce_peptide_inputs(peptides)
        errors = [sequence_only_chemistry_error(item) or (
            f"Configured CPP sequence limit exceeded: {len(item.sequence)} > {self.max_peptide_length}"
            if len(item.sequence) > self.max_peptide_length else None) for item in inputs]
        if on_unsupported == "raise" and any(errors):
            i = next(i for i, error in enumerate(errors) if error)
            raise ValueError(f"Input {i} is unsupported: {errors[i]}")
        sequences = [item.sequence for item, error in zip(inputs, errors) if not error]
        output = self._run_sidecar(sequences) if sequences else pd.DataFrame()
        rows = iter(output.to_dict("records"))
        results, qc = [], []
        for i, (item, error) in enumerate(zip(inputs, errors)):
            row = {} if error else next(rows)
            outside = not 3 <= len(item.sequence) <= 61
            detail = f"{_SCOPE} Outside published training-length range: {outside}."
            if error:
                detail += " Unsupported: " + error
            else:
                detail += f" Native threshold={row['threshold']}; class={row['label']}."
            context = MeasurementContext(
                estimate_type="ml_predicted", status="unsupported" if error else "available",
                analyte="parent peptide", class_label=None if error else (
                    "CPP" if row["label"] else "non-CPP"),
                score_semantics=None if error else "native SVC positive-class P(CPP)",
                detail=detail)
            pred = Prediction(
                kind=Kind.cpp_classification, peptide=item.sequence,
                score=row.get("score"), measurement_context=context,
                predictor_name="peptiverse-cpp", predictor_version=self.predictor_version,
                peptide_input=item, cache_key=prediction_cache_key(item, self.artifact_inventory))
            results.append(PeptideResult(preds=(pred,)))
            qc.append(dict(
                row, input_index=i, peptide=item.sequence,
                occurrence_id=item.occurrence_id, status=context.status,
                outside_source_training_lengths=outside, reason=error))
        self.last_qc = pd.DataFrame(qc)
        return results

    def _run_sidecar(self, sequences):
        with tempfile.TemporaryDirectory(prefix="mhctools_peptiverse_cpp_") as tmp:
            root = Path(tmp)
            pd.DataFrame({"__mhctools_id": range(len(sequences)), "peptide": sequences}).to_csv(
                root / "input.csv", index=False)
            (root / "manifest.txt").write_text(_CPP_MANIFEST)
            run_python_sidecar(
                backend_name="PeptiVerse CPP", python=self.peptiverse_python,
                sidecar=Path(__file__).with_name("peptiverse_sidecar.py"),
                cwd=self.peptiverse_home, timeout=self.subprocess_timeout,
                args=["--home", self.peptiverse_home, "--esm-home", self.peptiverse_esm_home,
                      "--manifest", str(root / "manifest.txt"), "--input", str(root / "input.csv"),
                      "--output", str(root / "output.csv"), "--device", self.device or "",
                      "--endpoint", "cpp"] + (["--uncertainty"] if self.uncertainty else []))
            output = parse_cpp_results(root / "output.csv", sequences)
        self.artifact_inventory = self.artifact_inventory.with_inference_reproduced()
        return output


def parse_cpp_results(filename, sequences):
    """Verify occurrence IDs, sequences, native scores/classes and entropy."""
    frame = pd.read_csv(filename, keep_default_na=False)
    required = {"__mhctools_id", "peptide", "score", "label", "threshold", "emb_tag",
                "uncertainty", "uncertainty_type"}
    if not required <= set(frame.columns):
        raise ValueError("PeptiVerse CPP output missing columns: %s" % sorted(required - set(frame.columns)))
    ids = pd.to_numeric(frame["__mhctools_id"], errors="raise")
    if sorted(ids.tolist()) != list(range(len(sequences))):
        raise RuntimeError("PeptiVerse CPP output did not preserve every input ID")
    frame = frame.assign(__mhctools_id=ids.astype(int)).sort_values("__mhctools_id").reset_index(drop=True)
    if frame.peptide.tolist() != list(sequences):
        raise RuntimeError("PeptiVerse CPP returned a different peptide than it was given")
    if not (frame.emb_tag == "wt").all():
        raise RuntimeError("PeptiVerse CPP returned a different embedding endpoint")
    for field in ("score", "label", "threshold"):
        frame[field] = pd.to_numeric(frame[field], errors="raise")
    if (not frame.score.map(lambda p: math.isfinite(p) and 0 <= p <= 1).all() or
            not (frame.threshold == CPP_THRESHOLD).all() or
            not (frame.label == (frame.score >= CPP_THRESHOLD).astype(int)).all()):
        raise RuntimeError("Invalid PeptiVerse CPP scores, labels or native threshold")
    entropy = []
    for score, value, kind in zip(frame.score, frame.uncertainty, frame.uncertainty_type):
        if value == "":
            if kind:
                raise RuntimeError("CPP entropy type without a value")
            entropy.append(None)
            continue
        value = float(value)
        p = min(max(score, 1e-9), 1 - 1e-9)
        expected = -(p * math.log(p) + (1 - p) * math.log(1 - p))
        if (kind != "binary_predictive_entropy_single_model" or not math.isfinite(value) or
                not math.isclose(value, expected, rel_tol=1e-9, abs_tol=1e-10)):
            raise RuntimeError("Invalid CPP single-model predictive entropy")
        entropy.append(value)
    frame["predictive_entropy_nats"] = entropy
    return frame.drop(columns=["uncertainty"])
