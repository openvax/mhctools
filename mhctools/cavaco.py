# Copyright (c) 2026 Mount Sinai School of Medicine
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0

"""Cavaco 2021 published Equation 1 with explicit descriptor provenance.

This independent baseline uses the published coefficients, not the author
app's broken W/Y counting. Default pI reproduces the app's free-terminal
algorithm; it is not a claim to reproduce every training-set descriptor.
Training assays were heterogeneous. Long vaccine peptides have not been
validated, and this output is neither a cleavage map nor serum survival.
"""

from collections import Counter
from functools import lru_cache
import hashlib
from importlib import resources
import json
import math
from numbers import Real
from types import MappingProxyType

import pandas as pd

from .peptide_input import coerce_peptide_inputs, sequence_only_chemistry_error
from .pred import Kind, MeasurementContext, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor


MODEL_SHA256 = "e7a652c6a411b8278da2570ae934f27617cadcd5b98471e0d5adb671c56ab4cc"
MODEL_VERSION = "2021-equation-1+author-app-descriptors@" + MODEL_SHA256
_SCOPE = (
    "Published Equation 1 / Table S20; mixed proteolysis/stability training "
    "assays (129 peptides), matrix unspecified. Source validation: sixteen "
    "related 15-residue C-terminal carboxamides in 50% human serum at 37 C; "
    "long vaccine peptide accuracy unestablished. Sequence-only baseline; "
    "no chemical modification or combined serum-survival model."
)


@lru_cache(maxsize=1)
def _model():
    raw = resources.files("mhctools.data").joinpath(
        "cavaco_published_model.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != MODEL_SHA256:
        raise ValueError("Cavaco bundled model data checksum mismatch")
    data = json.loads(raw)
    # No mutable model state escapes the loader.
    return MappingProxyType({
        "pka": MappingProxyType({
            aa: MappingProxyType(values) for aa, values in data["pka"].items()}),
        "coefficients": MappingProxyType(data["coefficients"]),
        "nonpolar_residues": frozenset(data["nonpolar_residues"]),
    })


def _author_app_pi(sequence, counts):
    """Reproduce unmodified app calc.js pI, including JS rounding/order."""
    pka = _model()["pka"]
    n_term = 10.0 ** pka[sequence[0]]["alpha_amino"]
    c_term = 10.0 ** pka[sequence[-1]]["alpha_carboxy"]
    ph, step = 7.0, 3.5
    for _ in range(100):
        h = 10.0 ** ph
        charge = n_term / (h + n_term) - h / (h + c_term)
        # Counter preserves first-occurrence order, as the app's object does.
        for aa, count in counts.items():
            side = 10.0 ** pka[aa]["side_chain"]
            if aa in "DCEY":
                charge -= count * h / (h + side)
            elif aa in "RKH":
                charge += count * side / (h + side)
        if math.floor(charge * 10000 + 0.5) == 0:
            return ph
        ph += -step if charge < 0 else step
        step /= 2
    return ph


def _validate_pi(value):
    if value is None:
        return None
    if (isinstance(value, bool) or not isinstance(value, Real) or
            not 0 <= value <= 14 or not math.isfinite(value)):
        raise ValueError("Each isoelectric point must be None or a finite number in [0, 14]")
    return float(value)


class CavacoHalfLife(AlleleFreePredictor):
    """Offline published-equation half-life baseline with no extra runtime.

    Only canonical L-peptides with free termini are accepted. The equation
    imposes no input length cutoff; this does not establish long-peptide
    accuracy. ``score`` is ln(half-life in minutes), ``value`` is hours.
    The assay matrix is unspecified, rather than assumed to be human serum.
    """

    predictor_version = MODEL_VERSION

    def __init__(self):
        _model()
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return "CavacoHalfLife(published_equation=1, pI='author-app or provided')"

    def _default_pred_kind(self):
        return Kind.peptide_half_life

    def predict(self, peptides, isoelectric_points=None, on_unsupported="raise"):
        """Estimate half-life while preserving exact inputs and occurrences.

        Parameters
        ----------
        peptides : str, PeptideInput, or iterable
            Canonical sequences or exact peptide records.
        isoelectric_points : iterable of float or None, optional
            Per-input pI values in [0, 14]. None entries use the author-app
            descriptor. Provided values, their source choice and model data
            identity are included in prediction provenance and cache keys.
        on_unsupported : {'raise', 'record'}
            Reject modified chemistry, or retain explicit unsupported records.

        Returns
        -------
        list of PeptideResult
            One allele-free half-life prediction per input, in input order.
        """
        self.last_qc = pd.DataFrame()
        if on_unsupported not in ("raise", "record"):
            raise ValueError("on_unsupported must be 'raise' or 'record'")
        inputs = coerce_peptide_inputs(peptides)
        supplied = ([None] * len(inputs) if isoelectric_points is None
                    else list(isoelectric_points))
        if len(supplied) != len(inputs):
            raise ValueError("isoelectric_points must have one entry per input")
        supplied = [_validate_pi(value) for value in supplied]
        errors = [sequence_only_chemistry_error(item) for item in inputs]
        if on_unsupported == "raise" and any(errors):
            index = next(i for i, error in enumerate(errors) if error)
            raise ValueError(f"Input {index} is unsupported: {errors[index]}")

        data = _model()
        coefficients = data["coefficients"]
        results, qc = [], []
        for index, (item, provided, error) in enumerate(zip(inputs, supplied, errors)):
            pi_source = "author-app-free-terminal" if provided is None else "provided"
            version = self.predictor_version + ":pI=" + pi_source
            fingerprint = json.dumps({
                "input": item.inference_identity_sha256,
                "model": version,
                "provided_pi": provided,
            }, sort_keys=True, separators=(",", ":"))
            common = dict(
                kind=Kind.peptide_half_life, peptide=item.sequence,
                predictor_name="cavaco", predictor_version=version,
                peptide_input=item,
                cache_key=hashlib.sha256(fingerprint.encode()).hexdigest())
            pi, np_percent, score, hours = None, None, None, None
            counts = Counter(item.sequence)
            if not error:
                pi = (_author_app_pi(item.sequence, counts)
                      if provided is None else provided)
                np_percent = sum(
                    n for aa, n in counts.items() if aa in data["nonpolar_residues"]
                ) / len(item.sequence) * 100.0
                score = (
                    coefficients["intercept"] +
                    coefficients["nonpolar_percent"] * np_percent +
                    coefficients["w_present"] * (counts["W"] >= 1) +
                    coefficients["y_two_or_more"] * (counts["Y"] >= 2) +
                    coefficients["pi_ge_10"] * (pi >= 10))
                hours = math.exp(score) / 60.0
            detail = f"{_SCOPE} pI source={pi_source}; pI={pi}; provided pI={provided}."
            if error:
                detail += " Unsupported: " + error
            context = MeasurementContext(
                estimate_type="ml_predicted", analyte="parent peptide",
                status="unsupported" if error else "available",
                unit=None if error else "hours",
                transform=None if error else "linear",
                score_semantics=None if error else "published Equation 1 ln(half-life in minutes)",
                detail=detail)
            results.append(PeptideResult(preds=(Prediction(
                score=score, value=hours, measurement_context=context, **common),)))
            qc.append(dict(
                input_index=index, peptide=item.sequence,
                occurrence_id=item.occurrence_id, status=context.status,
                pI=pi, pI_source=pi_source, provided_pI=provided,
                nonpolar_percent=np_percent, w_count=counts["W"],
                y_count=counts["Y"], ln_minutes=score, hours=hours,
                reason=error))
        self.last_qc = pd.DataFrame(qc)
        return results
