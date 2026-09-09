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
exposes **one** of its endpoints — the degradation half-life of a free peptide
in human serum, in hours — and emits ``Kind.serum_half_life``.

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
reports — pooling them into one field would mix units. Second, mhctools
identifies a peptide by its residue sequence and has nowhere to record a
chemical form, so a modified construct would silently be reported under its
unmodified sequence. Modified peptides are rejected rather than misattributed.

Installation
------------
Snapshot the upstream model repository (git-lfs; roughly 3 GB) and point
``PEPTIVERSE_HOME`` at it::

    git clone https://huggingface.co/ChatterjeeLab/PeptiVerse
    cd PeptiVerse && git checkout 8cf0b21dae356278ae96b414a088e4360357d16c

The upstream stack (torch, transformers==4.46.0, xgboost, lightning) is kept out
of the mhctools environment: inference runs in a subprocess under
``PEPTIVERSE_PYTHON``. Constructing upstream's predictor also instantiates ESM2
(``esm2_t33_650M_UR50D``), PeptideCLM and ChemBERTa, which are fetched from the
HuggingFace hub on first use and cached — so the first call needs network access
even though prediction itself does not.

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
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from .pred import Kind, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor

#: Upstream revision this wrapper was written and verified against.
UPSTREAM_REVISION = "8cf0b21dae356278ae96b414a088e4360357d16c"

#: ESM2's positional limit; upstream's WTEmbedder truncates beyond it.
PEPTIVERSE_MAX_PEPTIDE_LENGTH = 1022

_VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# A one-row manifest in upstream's format. Restricting it to Halflife (with the
# SMILES column blanked out with "-") makes `_load_all_best_models` load the
# single sequence half-life model rather than every property's best model.
# "Transformer" resolves to transformer_wt_log, whose output is expm1'd by
# upstream into hours.
_HALF_LIFE_MANIFEST = (
    "Properties, Best_Model_WT, Best_Model_SMILES, Type, "
    "Threshold_WT, Threshold_SMILES,\n"
    "Halflife, Transformer, -, Regression, -, -,\n"
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
    candidate = str(Path(candidate).expanduser())
    if not Path(candidate, "inference.py").is_file():
        raise FileNotFoundError(
            "inference.py not found in %r — is this a PeptiVerse snapshot?"
            % candidate)
    weights = Path(candidate, "training_classifiers", "half_life")
    if not weights.is_dir():
        raise FileNotFoundError(
            "training_classifiers/half_life not found in %r — the snapshot is "
            "missing its half-life weights (git-lfs not pulled?)" % candidate)
    return candidate


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
        truncated by ESM2. Default 1022.
    """

    mhc_class = "none"

    def __init__(
            self,
            peptiverse_home=None,
            peptiverse_python=None,
            device=None,
            uncertainty=False,
            max_peptide_length=PEPTIVERSE_MAX_PEPTIDE_LENGTH):
        if max_peptide_length > PEPTIVERSE_MAX_PEPTIDE_LENGTH:
            raise ValueError(
                "ESM2 accepts at most %d residues; max_peptide_length cannot "
                "exceed that (got %d)"
                % (PEPTIVERSE_MAX_PEPTIDE_LENGTH, max_peptide_length))
        self.peptiverse_home = _find_peptiverse_home(peptiverse_home)
        self.peptiverse_python = _resolve_python(peptiverse_python)
        self.device = device
        self.uncertainty = uncertainty
        self.max_peptide_length = max_peptide_length
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return "PeptiVerse(peptiverse_home=%r, device=%r)" % (
            self.peptiverse_home, self.device)

    def _default_pred_kind(self):
        return Kind.serum_half_life

    def _predictor_name(self):
        return "peptiverse"

    @property
    def predictor_version(self):
        return "%s:transformer_wt_log" % UPSTREAM_REVISION

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

    def predict(self, peptides):
        """Predict serum half-life for a list of peptides.

        Returns
        -------
        list of PeptideResult
            One entry per input peptide, in input order, each holding a single
            ``Kind.serum_half_life`` prediction with an empty ``allele``. Both
            ``score`` and ``value`` carry the predicted half-life in **hours**
            (higher = longer-lived); ``value`` is the units-bearing field and
            ``score`` repeats it so that rank-based consumers work without
            knowing the unit.
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        if not peptide_list:
            self.last_qc = pd.DataFrame()
            return []

        output = self._run_sidecar(peptide_list)
        return [
            PeptideResult(preds=(Prediction(
                kind=Kind.serum_half_life,
                score=hours,
                value=hours,
                peptide=peptide,
                predictor_name=self._predictor_name(),
                predictor_version=self.predictor_version),))
            for peptide, hours in zip(peptide_list, output["hours"])
        ]

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
            command = [
                self.peptiverse_python,
                # Run through runpy rather than passing the script path
                # directly: that would put mhctools/ on sys.path[0], where
                # mhctools/logging.py shadows the stdlib logging module that
                # torch imports.
                "-c",
                (
                    "import runpy,sys; sys.argv=sys.argv[1:]; "
                    "runpy.run_path(sys.argv[0],run_name='__main__')"
                ),
                str(sidecar),
                "--home", self.peptiverse_home,
                "--manifest", str(manifest_path),
                "--input", str(input_path),
                "--output", str(output_path),
                "--device", self.device or "",
            ]
            if self.uncertainty:
                command.append("--uncertainty")

            environment = os.environ.copy()
            environment["PYTHONNOUSERSITE"] = "1"
            process = subprocess.run(
                command,
                cwd=self.peptiverse_home,
                env=environment,
                capture_output=True,
                text=True,
            )
            if process.returncode != 0:
                raise RuntimeError(
                    "PeptiVerse inference failed (exit %d):\n%s" % (
                        process.returncode,
                        (process.stderr or process.stdout).strip()))
            output = parse_peptiverse_results(output_path, peptide_list)

        self.last_qc = output.drop(columns=["hours"])
        return output


def parse_peptiverse_results(filename, peptide_list):
    """Parse sidecar output into a DataFrame aligned with *peptide_list*.

    Every input peptide must come back exactly once, identified by the
    ``__mhctools_id`` the caller assigned, and with the peptide it was scored
    under unchanged — upstream's embedders truncate or alter unsupported input
    rather than failing, so a mismatch is an error rather than something to
    reconcile silently. Duplicate input peptides keep separate rows.

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
    return output
