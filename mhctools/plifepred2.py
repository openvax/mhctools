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

"""Wrapper for PlifePred2's natural-peptide half-life model.

.. warning::

   **This endpoint's semantics are not established.** PlifePred2 ships no
   publication, no training data and no target definition, so its units,
   transform, species and assay matrix are all inferred from the artifacts.
   By default this wrapper therefore reports only the model's native output as
   ``score`` and leaves ``value`` empty. Pass ``assume_log10_seconds=True`` to
   opt in to a half-life in hours, accepting the inference below.

What is actually known
----------------------
Both shipped models are ``RandomForestRegressor``. That much is verified by
loading them, and it is enough to establish one thing: the output is **not** a
class probability, so upstream's README ("Halflife ... Predicted probability")
and its CLI's ``predict_proba`` branch are both wrong — the branch is dead code
and ``model.predict`` always runs. The output is a regression target on some
monotone-in-half-life scale, which makes it usable for ranking whatever the
transform turns out to be.

What is inferred, and how strongly
----------------------------------
The evidence for ``log10(half-life in seconds)`` is circumstantial but not
weak. Inverting the extreme leaf values of the shipped forests under that
reading gives round durations — the natural model's maximum is 604800.0 s,
exactly 7.000 days, and the modified model's is 8208000 s, exactly 95.0 days —
agreeing to about seven significant figures. Under log2 or ln the same extrema
invert to a few seconds up to a couple of minutes, which no half-life dataset
would span. The minimum also inverts to 20.2 s, matching the 20-second floor
documented by the lineage paper (Mathur et al. 2018, PLoS ONE 13(6):e0196829).

That last point is the weakest of the three, and it is worth being precise
about why. The same forests hold targets out to 7 and 95 days, which violates
the 24-hour ceiling that paper also documents. So PlifePred2 was trained on a
*different* dataset, and a filter from the old one cannot independently
establish the new one's target. The round-number extrema stand on their own;
the floor agreement is corroboration, not proof.

Note also that the lineage paper states log2, not log10 — anyone carrying
assumptions across from it will be wrong by a factor of log2(10) ~ 3.32.

What is not established at all
------------------------------
The **species and assay matrix**. ``Kind.blood_half_life`` is assigned because
the lineage paper drew its data from PEPlife filtered to mammalian whole blood,
and that is the best available guide. But it is inherited from a dataset this
model demonstrably does not use. Whether PlifePred2's own data is whole blood,
plasma, serum or a mixture, from which species, and ex vivo or in vivo, is
unknown. Do not report this as a measured whole-blood property, and do not
treat it as interchangeable with :mod:`mhctools.peptiverse`'s human-serum
endpoint. Resolving this needs the authors or a model-specific publication.

Natural peptides only
---------------------
Upstream ships a second model for modified peptides. It is not wrapped. Its
modification flags (D-amino acid, terminal modifications, cyclization, PTM) are
applied uniformly to every sequence in one invocation rather than per peptide,
and this adapter does not consume the exact per-input chemical form carried by
``PeptideInput``. Modified inputs and non-standard residues are rejected rather
than reduced to the natural-peptide sequence model.

Installation
------------
Two GPLv3 installations, neither vendored:

- ``PLIFEPRED2_HOME`` — the installed ``plifepred2`` package directory, holding
  ``models/plifepred2_natural_model.sav`` (``pip install plifepred2`` into any
  environment, then point at its ``site-packages/plifepred2``).
- ``PFEATURE_HOME`` — a checkout of ``raghavagps/Pfeature``'s ``Standalone``
  directory at revision ``93636eb95bed9df2893b7a0c56b1215e648ecdbf``,
  holding ``pfeature_comp.py`` and ``Data/``.

The Linux-only ``pfeature_comp`` binary that PlifePred2 bundles is **not** used.
It is a PyInstaller freeze of Pfeature's ``pfeature_comp.py``; the pinned plain
Python source is used in the isolated environment instead.

Provenance and limits
---------------------
Upstream: https://pypi.org/project/plifepred2/ and
https://github.com/shindebpratik/plifepred2
Lineage: https://doi.org/10.1371/journal.pone.0196829

PlifePred2 itself cites no publication, so its training set is unverified beyond
what the artifacts reveal; the ceiling change suggests a larger, later dataset
than the 2018 one. In the lineage paper the composition-based natural model was
the weaker of the pair (r = 0.643, against 0.743 for chemical descriptors), and
sequences up to 90% similar were deliberately kept in the data, so reported
accuracy is optimistic for novel peptides.
"""

import math
import os
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

#: Upstream package version this wrapper was written and verified against.
UPSTREAM_VERSION = "1.0"
PFEATURE_REVISION = "93636eb95bed9df2893b7a0c56b1215e648ecdbf"

#: sha256 of the ``plifepred2-1.0-py3-none-any.whl`` these paths were read from.
UPSTREAM_WHEEL_SHA256 = (
    "e0eb32928b7970d5a113d8826ef9e451c94f0156e0b628dddb0da7103ac2c32e")

# Upstream's CLI silently drops anything outside this range; the wrapper rejects
# instead, so a caller never gets a short result set with no explanation.
PLIFEPRED2_MIN_PEPTIDE_LENGTH = 12
PLIFEPRED2_MAX_PEPTIDE_LENGTH = 100

_VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

_NATURAL_MODEL = os.path.join("models", "plifepred2_natural_model.sav")

_PLIFEPRED2_ARTIFACTS = {
    _NATURAL_MODEL: (
        "natural_peptide_random_forest",
        "1c9009319685c7ba9a0ceb1fe5e00322590039cca4589c24f89bf3995f83080b",
        "joblib_pickle"),
}

_PFEATURE_ARTIFACTS = {
    "pfeature_comp.py": (
        "qso_feature_implementation",
        "f02f39ed2aa95dd282b0f20684c0f8570522a7a886c987b6acb3cf4a18a07936",
        "python_source"),
    os.path.join("Data", "Schneider-Wrede.csv"): (
        "qso_distance_matrix",
        "5ae5c4e83e41490464ac0fc120e209de2e87a2b37b8ed3350b028b7760b89a56",
        "csv"),
    os.path.join("Data", "Grantham.csv"): (
        "qso_distance_matrix",
        "87baa46a74c6510d466011358973d340d943844aeb8ce5fb46541db97aaa883c",
        "csv"),
    os.path.join("Data", "PhysicoChemical.csv"): (
        "pfeature_startup_resource",
        "f8e94cf9464539ddd0e9e97a20da2fe0ffbdcb8a48c1a632d7117153bea11c04",
        "csv"),
    os.path.join("Data", "aaindex.csv"): (
        "pfeature_startup_resource",
        "c18b4e54a40136b56533100e69c6f8251d80eb5923a02a00db6d00b14d843b65",
        "csv"),
    os.path.join("Data", "AAIndexNames.csv"): (
        "pfeature_startup_resource",
        "608dca782677b57349c1a6a208be2411278db35fa739c26bdb522b275ebe0f43",
        "csv"),
}

PLIFEPRED2_BACKEND_SPEC = BackendSpec(
    name="plifepred2",
    endpoint="undocumented_blood_half_life_native_regression",
    developed_against="plifepred2@%s+pfeature@%s" % (
        UPSTREAM_VERSION, PFEATURE_REVISION),
    license="PlifePred2: GPL-3.0; Pfeature: GPL-3.0",
    serialization="joblib pickle random forest",
    entry_point="prediction_only",
    # Real-model smoke recorded in docs/optional-backends.md. CI continues to
    # exercise only the conformance path and never unpickles third-party data.
    supported_platforms=("macos-arm64-cpu",),
    supported_interpreters=("python3.12",),
)

_SECONDS_PER_HOUR = 3600.0


def half_life_hours(log10_seconds):
    """Convert a raw PlifePred2 prediction to hours.

    This applies the **inferred** ``log10(seconds)`` transform; see the module
    docstring for the evidence and its limits. Upstream documents no target,
    so a caller reaching for this is accepting that inference.
    """
    return (10.0 ** log10_seconds) / _SECONDS_PER_HOUR


def _find_plifepred2_home(plifepred2_home=None):
    """Resolve the installed ``plifepred2`` package directory."""
    candidate = plifepred2_home or os.environ.get("PLIFEPRED2_HOME")
    if not candidate:
        home = Path.home() / "plifepred2"
        if home.is_dir():
            candidate = str(home)
    if not candidate:
        raise FileNotFoundError(
            "PlifePred2 not found. Set PLIFEPRED2_HOME or pass "
            "plifepred2_home= to the constructor. `pip install plifepred2` "
            "into an environment and point at its site-packages/plifepred2.")
    candidate = str(Path(candidate).expanduser().resolve())
    if not Path(candidate, _NATURAL_MODEL).is_file():
        raise FileNotFoundError(
            "%s not found in %r — is this an installed plifepred2 package?"
            % (_NATURAL_MODEL, candidate))
    return candidate


def _find_pfeature_home(pfeature_home=None):
    """Resolve the Pfeature ``Standalone`` directory."""
    candidate = pfeature_home or os.environ.get("PFEATURE_HOME")
    if not candidate:
        home = Path.home() / "Pfeature" / "Standalone"
        if home.is_dir():
            candidate = str(home)
    if not candidate:
        raise FileNotFoundError(
            "Pfeature not found. Set PFEATURE_HOME or pass pfeature_home= to "
            "the constructor. Clone https://github.com/raghavagps/Pfeature and "
            "point at its Standalone directory.")
    candidate = str(Path(candidate).expanduser().resolve())
    if not Path(candidate, "pfeature_comp.py").is_file():
        raise FileNotFoundError(
            "pfeature_comp.py not found in %r — point at Pfeature's Standalone "
            "directory, not the repository root." % candidate)
    for relative in _PFEATURE_ARTIFACTS:
        if relative == "pfeature_comp.py":
            continue
        if not Path(candidate, relative).is_file():
            raise FileNotFoundError(
                "%s not found in %r — the pinned Pfeature process reads this "
                "resource before or during QSO extraction."
                % (relative, candidate))
    return candidate


def _artifact_inventory(
        plifepred2_home, pfeature_home, assume_log10_seconds):
    artifacts = []
    for relative, (role, expected_sha256, serialization) in (
            _PLIFEPRED2_ARTIFACTS.items()):
        artifacts.append(inspect_artifact(
            name="plifepred2/%s" % relative,
            role=role,
            path=Path(plifepred2_home) / relative,
            expected_sha256=expected_sha256,
            serialization=serialization))
    for relative, (role, expected_sha256, serialization) in (
            _PFEATURE_ARTIFACTS.items()):
        artifacts.append(inspect_artifact(
            name="pfeature/%s" % relative,
            role=role,
            path=Path(pfeature_home) / relative,
            expected_sha256=expected_sha256,
            serialization=serialization))
    return backend_inventory(
        spec=PLIFEPRED2_BACKEND_SPEC,
        artifacts=artifacts,
        settings={
            "assume_log10_seconds": assume_log10_seconds,
            "descriptor": "QSO",
            "model_variant": "natural",
            "qso_lag": 1,
            "qso_weight": 0.1,
            "scikit_learn": "1.4.2",
        })


def _resolve_python(plifepred2_python=None):
    value = plifepred2_python or os.environ.get("PLIFEPRED2_PYTHON")
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
    raise FileNotFoundError("PlifePred2 Python does not exist: %s" % value)


class PlifePred2(AlleleFreePredictor):
    """Whole-blood half-life predictions from local PlifePred2 + Pfeature.

    Parameters
    ----------
    plifepred2_home : str, optional
        Installed ``plifepred2`` package directory. Resolved from the argument,
        then ``$PLIFEPRED2_HOME``, then ``~/plifepred2``.
    pfeature_home : str, optional
        Pfeature ``Standalone`` directory. Resolved from the argument, then
        ``$PFEATURE_HOME``, then ``~/Pfeature/Standalone``.
    plifepred2_python : str, optional
        Interpreter with scikit-learn, joblib, pandas and tqdm. Resolved from
        the argument, then ``$PLIFEPRED2_PYTHON``, then the current interpreter.
        Upstream pins ``scikit-learn==1.4.2``; a different version loads the
        forest but warns, so give this its own environment to be exact.
    assume_log10_seconds : bool
        Opt in to reporting a half-life in hours by treating the model's output
        as ``log10(seconds)``. Off by default: that transform is inferred from
        the shipped forests, not documented by upstream (see the module
        docstring). While off, ``value`` is left empty and only the native
        output is reported, in ``score``. Turning it on is an assertion that
        you accept the inference.
    allow_unverified_assets : bool
        Permit explicitly supplied files whose checksums differ from the pinned
        model and Pfeature inventory. False by default. Their actual identity
        is still recorded in :attr:`artifact_inventory` and every prediction.
    subprocess_timeout : float, optional
        Maximum seconds allowed for the isolated inference process.
    """

    mhc_class = "none"

    def __init__(
            self,
            plifepred2_home=None,
            pfeature_home=None,
            plifepred2_python=None,
            assume_log10_seconds=False,
            allow_unverified_assets=False,
            subprocess_timeout=300):
        self.plifepred2_home = _find_plifepred2_home(plifepred2_home)
        self.pfeature_home = _find_pfeature_home(pfeature_home)
        self.plifepred2_python = _resolve_python(plifepred2_python)
        self.assume_log10_seconds = assume_log10_seconds
        self.subprocess_timeout = subprocess_timeout
        self.artifact_inventory = _artifact_inventory(
            self.plifepred2_home,
            self.pfeature_home,
            self.assume_log10_seconds)
        self.artifact_inventory.require_usable(
            allow_unverified=allow_unverified_assets)
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return ("PlifePred2(plifepred2_home=%r, pfeature_home=%r, "
                "assume_log10_seconds=%r)") % (
            self.plifepred2_home, self.pfeature_home,
            self.assume_log10_seconds)

    def _default_pred_kind(self):
        return Kind.blood_half_life

    def _predictor_name(self):
        return "plifepred2"

    @property
    def predictor_version(self):
        return self.artifact_inventory.predictor_version

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if not peptide:
                raise ValueError("Empty peptide is not allowed")
            if not (PLIFEPRED2_MIN_PEPTIDE_LENGTH
                    <= len(peptide) <= PLIFEPRED2_MAX_PEPTIDE_LENGTH):
                raise ValueError(
                    "PlifePred2 supports peptides of %d-%d residues; got %r "
                    "(length %d). Upstream drops out-of-range sequences "
                    "silently, so they are rejected here instead."
                    % (PLIFEPRED2_MIN_PEPTIDE_LENGTH,
                       PLIFEPRED2_MAX_PEPTIDE_LENGTH, peptide, len(peptide)))
            invalid = set(peptide) - _VALID_AMINO_ACIDS
            if invalid:
                raise ValueError(
                    "Peptide %r contains non-standard residues: %s. Only the "
                    "natural-peptide model is wrapped, so modified peptides "
                    "are rejected rather than scored as their unmodified "
                    "sequence." % (peptide, "".join(sorted(invalid))))

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
            analyte="free parent peptide",
            matrix=None,
            unit=(
                "hours"
                if status == "available" and self.assume_log10_seconds
                else None),
            transform=(
                "linear"
                if status == "available" and self.assume_log10_seconds
                else None),
            score_semantics=(
                "native regression output; target transform and assay "
                "context unresolved"
                if status == "available" else None),
            detail=detail,
        )

    def predict(self, peptides, on_unsupported="raise"):
        """Predict whole-blood half-life for a list of peptides.

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
            ``Kind.blood_half_life`` prediction with an empty ``allele``.

            ``score`` is always the model's **native output** — higher means
            longer-lived, which is all that is needed to rank. ``value`` is the
            half-life in hours, and is filled **only** when the predictor was
            constructed with ``assume_log10_seconds=True``; otherwise it is
            ``None``, because the transform from the native output to a
            duration is inferred rather than documented.
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
            output = pd.DataFrame(columns=["log10_seconds"])
        native_values = iter(output["log10_seconds"])
        results = []
        for peptide_input, error in zip(peptide_inputs, errors):
            common = {
                "kind": Kind.blood_half_life,
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
                native = next(native_values)
                prediction = Prediction(
                    score=native,
                    value=(half_life_hours(native)
                           if self.assume_log10_seconds else None),
                    measurement_context=self._measurement_context(),
                    **common)
            results.append(PeptideResult(preds=(prediction,)))
        return results

    def _run_sidecar(self, peptide_list):
        with tempfile.TemporaryDirectory(prefix="mhctools_plifepred2_") as tmp:
            tmp_path = Path(tmp)
            input_path = tmp_path / "input.csv"
            output_path = tmp_path / "output.csv"
            pd.DataFrame({
                "__mhctools_id": range(len(peptide_list)),
                "peptide": peptide_list,
            }).to_csv(input_path, index=False)

            sidecar = Path(__file__).with_name("plifepred2_sidecar.py")
            arguments = [
                "--pfeature-home", self.pfeature_home,
                "--model", str(Path(self.plifepred2_home, _NATURAL_MODEL)),
                "--input", str(input_path),
                "--output", str(output_path),
            ]
            if self.subprocess_timeout is not None:
                arguments.extend([
                    "--pfeature-timeout", str(self.subprocess_timeout)])
            run_python_sidecar(
                backend_name="PlifePred2",
                python=self.plifepred2_python,
                sidecar=sidecar,
                args=arguments,
                timeout=self.subprocess_timeout)
            output = parse_plifepred2_results(output_path, peptide_list)

        self.artifact_inventory = (
            self.artifact_inventory.with_inference_reproduced())
        self.last_qc = output
        return output


def parse_plifepred2_results(filename, peptide_list):
    """Parse sidecar output into a DataFrame aligned with *peptide_list*.

    Pfeature links its feature rows to inputs by position only, so every input
    peptide must come back exactly once, under the ``__mhctools_id`` it was
    given and with its sequence unchanged. Anything else means the rows have
    shifted relative to the peptides and the scores belong to the wrong ones.

    Returns
    -------
    pandas.DataFrame
        Sorted by ``__mhctools_id``, carrying the model's native output as
        ``log10_seconds`` and, as ``hours_if_log10_seconds``, what that would
        mean in hours under the inferred transform. The column is named for
        its assumption on purpose: nothing upstream documents the target.
    """
    output = pd.read_csv(filename, keep_default_na=False)
    for column in ("__mhctools_id", "peptide", "log10_seconds"):
        if column not in output.columns:
            raise ValueError(
                "PlifePred2 output missing column %r; got %s"
                % (column, list(output.columns)))
    output["__mhctools_id"] = output["__mhctools_id"].astype(int)
    if sorted(output["__mhctools_id"].tolist()) != list(range(len(peptide_list))):
        raise RuntimeError(
            "PlifePred2 output did not preserve every input peptide: expected "
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
            "PlifePred2 returned a different peptide than it was given "
            "(id, sent, got): %s" % mismatched[:5])
    output["log10_seconds"] = output["log10_seconds"].astype(float)
    if not output["log10_seconds"].map(math.isfinite).all():
        raise RuntimeError(
            "PlifePred2 returned a non-finite prediction; the feature matrix "
            "is probably misaligned with the model")
    output["hours_if_log10_seconds"] = output["log10_seconds"].map(
        half_life_hours)
    return output
