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

"""Wrapper for PlifePred2's natural-peptide blood half-life model.

How long a free peptide survives incubated in **mammalian whole blood**, in
hours. That matrix matters and is why this emits ``Kind.blood_half_life`` rather
than the ``Kind.serum_half_life`` that :mod:`mhctools.peptiverse` emits: serum
is blood with the cells and clotting factors removed, and peptide stability
differs measurably between the two. Neither is ``pMHC_stability``, which is
peptide-MHC complex dissociation, and neither is a cleavage-site map.

Units
-----
Upstream's README calls the output column ``Halflife`` and describes it as a
"Predicted probability". It is neither a probability nor a raw duration. Both
shipped models are ``RandomForestRegressor``, so the ``predict_proba`` branch in
upstream's CLI is dead code, and the target is ``log10(half-life in seconds)``.

That transform is established rather than assumed. The lineage paper (Mathur et
al. 2018, PLoS ONE 13(6):e0196829) built its dataset from PEPlife filtered to
mammalian blood, discarding anything below 20 seconds. Inverting the extreme
leaf values of the shipped forests under log10-seconds reproduces that floor
(20.2 s) and yields round durations at the top (exactly 7.000 days for the
natural model, 95.0 days for the modified one); under log2 or ln the entire
training range falls below the documented 20-second floor, which is impossible.

Note that PlifePred2 departs from that paper in two ways, so its numbers are not
the paper's: the transform is log10 where the paper used log2 (a constant factor
of ~3.32), and the paper's 24-hour ceiling is gone.

Natural peptides only
---------------------
Upstream ships a second model for modified peptides. It is not wrapped. Its
modification flags (D-amino acid, terminal modifications, cyclization, PTM) are
applied uniformly to every sequence in one invocation rather than per peptide,
and mhctools identifies a peptide by its residue sequence with nowhere to record
a chemical form — so a batch of differently modified constructs would be
mislabelled. Peptides with non-standard residues are rejected instead.

Installation
------------
Two GPLv3 checkouts, neither vendored:

- ``PLIFEPRED2_HOME`` — the installed ``plifepred2`` package directory, holding
  ``models/plifepred2_natural_model.sav`` (``pip install plifepred2`` into any
  environment, then point at its ``site-packages/plifepred2``).
- ``PFEATURE_HOME`` — a checkout of ``raghavagps/Pfeature``'s ``Standalone``
  directory, holding ``pfeature_comp.py`` and ``Data/``.

The Linux-only ``pfeature_comp`` binary that PlifePred2 bundles is **not** used.
It is a PyInstaller freeze of Pfeature's ``pfeature_comp.py``, and that plain
Python source computes the same descriptor on any platform.

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
import subprocess
import sys
import tempfile

import pandas as pd

from .pred import Kind, PeptideResult, Prediction
from .wrapper_base import AlleleFreePredictor

#: Upstream package version this wrapper was written and verified against.
UPSTREAM_VERSION = "1.0"

#: sha256 of the ``plifepred2-1.0-py3-none-any.whl`` these paths were read from.
UPSTREAM_WHEEL_SHA256 = (
    "e0eb32928b7970d5a113d8826ef9e451c94f0156e0b628dddb0da7103ac2c32e")

# Upstream's CLI silently drops anything outside this range; the wrapper rejects
# instead, so a caller never gets a short result set with no explanation.
PLIFEPRED2_MIN_PEPTIDE_LENGTH = 12
PLIFEPRED2_MAX_PEPTIDE_LENGTH = 100

_VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

_NATURAL_MODEL = os.path.join("models", "plifepred2_natural_model.sav")

_SECONDS_PER_HOUR = 3600.0


def half_life_hours(log10_seconds):
    """Convert a raw PlifePred2 prediction to hours.

    The model's target is ``log10(half-life in seconds)``; see the module
    docstring for how that was established.
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
    candidate = str(Path(candidate).expanduser())
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
    candidate = str(Path(candidate).expanduser())
    if not Path(candidate, "pfeature_comp.py").is_file():
        raise FileNotFoundError(
            "pfeature_comp.py not found in %r — point at Pfeature's Standalone "
            "directory, not the repository root." % candidate)
    for data_file in ("Schneider-Wrede.csv", "Grantham.csv"):
        if not Path(candidate, "Data", data_file).is_file():
            raise FileNotFoundError(
                "Data/%s not found in %r — Pfeature's QSO descriptor reads its "
                "distance matrices from that directory." % (data_file, candidate))
    return candidate


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
    """

    mhc_class = "none"

    def __init__(
            self,
            plifepred2_home=None,
            pfeature_home=None,
            plifepred2_python=None):
        self.plifepred2_home = _find_plifepred2_home(plifepred2_home)
        self.pfeature_home = _find_pfeature_home(pfeature_home)
        self.plifepred2_python = _resolve_python(plifepred2_python)
        self.last_qc = pd.DataFrame()

    def __str__(self):
        return "PlifePred2(plifepred2_home=%r, pfeature_home=%r)" % (
            self.plifepred2_home, self.pfeature_home)

    def _default_pred_kind(self):
        return Kind.blood_half_life

    def _predictor_name(self):
        return "plifepred2"

    @property
    def predictor_version(self):
        return "%s:natural" % UPSTREAM_VERSION

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

    def predict(self, peptides):
        """Predict whole-blood half-life for a list of peptides.

        Returns
        -------
        list of PeptideResult
            One entry per input peptide, in input order, each holding a single
            ``Kind.blood_half_life`` prediction with an empty ``allele``. Both
            ``score`` and ``value`` carry the half-life in **hours**; the raw
            ``log10(seconds)`` model output is kept in :attr:`last_qc`.
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        if not peptide_list:
            self.last_qc = pd.DataFrame()
            return []

        output = self._run_sidecar(peptide_list)
        return [
            PeptideResult(preds=(Prediction(
                kind=Kind.blood_half_life,
                score=hours,
                value=hours,
                peptide=peptide,
                predictor_name=self._predictor_name(),
                predictor_version=self.predictor_version),))
            for peptide, hours in zip(peptide_list, output["hours"])
        ]

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
            command = [
                self.plifepred2_python,
                # runpy rather than the script path, so that mhctools/ does not
                # land on sys.path[0] where mhctools/logging.py would shadow
                # the stdlib logging module.
                "-c",
                ("import runpy,sys; sys.argv=sys.argv[1:]; "
                 "runpy.run_path(sys.argv[0],run_name='__main__')"),
                str(sidecar),
                "--pfeature-home", self.pfeature_home,
                "--model", str(Path(self.plifepred2_home, _NATURAL_MODEL)),
                "--input", str(input_path),
                "--output", str(output_path),
            ]
            environment = os.environ.copy()
            environment["PYTHONNOUSERSITE"] = "1"
            process = subprocess.run(
                command,
                env=environment,
                capture_output=True,
                text=True,
            )
            if process.returncode != 0:
                raise RuntimeError(
                    "PlifePred2 inference failed (exit %d):\n%s" % (
                        process.returncode,
                        (process.stderr or process.stdout).strip()))
            output = parse_plifepred2_results(output_path, peptide_list)

        self.last_qc = output.drop(columns=["hours"])
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
        Sorted by ``__mhctools_id``, with the raw ``log10_seconds`` output and
        the derived ``hours``.
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
    output["hours"] = output["log10_seconds"].map(half_life_hours)
    return output
