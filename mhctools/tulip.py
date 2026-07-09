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

"""Wrapper for TULIP-TCR (https://github.com/barthelemymp/TULIP-TCR).

TULIP predicts pMHC:TCR binding — how well a paired αβ T-cell receptor
(given by its CDR3α/CDR3β loops) recognises a peptide presented by an MHC
allele. Like :class:`~mhctools.nettcr.NetTCR` it produces the
:attr:`~mhctools.pred.Kind.pMHC_TCR_binding` kind, but unlike NetTCR it also
takes the presenting **MHC allele** as an input.

Licensing / isolation
----------------------
TULIP-TCR is distributed under the **GPLv3** and is coupled to a specific,
old transformers release (``transformers==4.32.1``); mhctools is Apache-2.0
and depends on neither torch nor transformers. To keep the two apart in both
senses, this wrapper **vendors none of TULIP** — no source, no weights, no
tokenizers. It runs the user's own TULIP checkout as a black box, in a
separate interpreter, by invoking TULIP's own ``predict.py`` via subprocess
(the same "shell out to a user-provided install" pattern used for the DTU
``netMHC*`` tools and NetTCR). Nothing here imports TULIP's GPL code.

You therefore need two things, supplied via constructor arguments or
environment variables:

* ``TULIP_HOME`` — a clone of https://github.com/barthelemymp/TULIP-TCR
  (provides ``predict.py``, ``src/``, ``aatok/``, ``mhctok/``, ``configs/``
  and the released ``model_weights/``).
* ``TULIP_PYTHON`` — a Python interpreter that can import TULIP, i.e. an
  isolated environment with ``torch`` and ``transformers==4.32.1``. Because
  ``tokenizers 0.13.x`` (which that transformers pins) has no cp312 wheel,
  build this env on **Python 3.11** so the install uses a prebuilt wheel and
  needs no Rust toolchain. See ``scripts/setup_tulip_env.sh``.

Scores are TULIP's per-pair log-likelihood score (higher = more likely to
bind), read back from the CSVs ``predict.py`` writes.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile

import pandas as pd

from .pred import COLUMNS, Kind, PeptideResult, Prediction
from .tcr import TCR


# TULIP's missing-value token; used for an absent CDR3 chain or MHC allele.
MISSING_TOKEN = "<MIS>"

_CLONE_HINT = "Clone from https://github.com/barthelemymp/TULIP-TCR"
_ENV_HINT = (
    "Build an isolated env (Python 3.11, torch + transformers==4.32.1) and "
    "point TULIP_PYTHON at its interpreter; see scripts/setup_tulip_env.sh")


def _find_tulip_home(tulip_home=None):
    """Resolve the TULIP-TCR checkout directory.

    Checks, in order: the *tulip_home* argument, ``TULIP_HOME``, then
    ``~/TULIP-TCR`` and ``~/code/TULIP-TCR``. An explicitly-provided path is
    validated up front (including that it actually contains ``predict.py``)
    so a typo fails clearly rather than deep inside a subprocess.
    """
    for source, path in (
            ("tulip_home argument", tulip_home),
            ("TULIP_HOME", os.environ.get("TULIP_HOME"))):
        if path:
            if not os.path.isdir(path):
                raise FileNotFoundError(
                    "TULIP-TCR directory from %s does not exist: %s. %s"
                    % (source, path, _CLONE_HINT))
            if not os.path.isfile(os.path.join(path, "predict.py")):
                raise FileNotFoundError(
                    "TULIP-TCR directory from %s has no predict.py: %s. %s"
                    % (source, path, _CLONE_HINT))
            return path
    home = os.path.expanduser("~")
    for candidate in (
            os.path.join(home, "TULIP-TCR"),
            os.path.join(home, "code", "TULIP-TCR")):
        if os.path.isfile(os.path.join(candidate, "predict.py")):
            return candidate
    raise FileNotFoundError(
        "TULIP-TCR not found. Set TULIP_HOME or pass tulip_home= to the "
        "constructor. %s" % _CLONE_HINT)


def _find_tulip_python(tulip_python=None):
    """Resolve the interpreter used to run TULIP.

    Checks the *tulip_python* argument then ``TULIP_PYTHON``. Falls back to
    the interpreter running mhctools only as a last resort — that works only
    if torch and transformers==4.32.1 happen to be importable there, which is
    exactly the coupling this wrapper exists to avoid, so it is not the
    recommended path.
    """
    for source, path in (
            ("tulip_python argument", tulip_python),
            ("TULIP_PYTHON", os.environ.get("TULIP_PYTHON"))):
        if path:
            if not (os.path.isfile(path) and os.access(path, os.X_OK)):
                raise FileNotFoundError(
                    "TULIP_PYTHON from %s is not an executable file: %s. %s"
                    % (source, path, _ENV_HINT))
            return path
    return sys.executable


class Tulip(object):
    """Wrapper for TULIP-TCR pMHC:TCR binding predictions.

    Parameters
    ----------
    tulip_home : str, optional
        Path to the cloned TULIP-TCR repository root. If omitted, resolved
        from ``TULIP_HOME`` or ``~/TULIP-TCR`` / ``~/code/TULIP-TCR``.
    tulip_python : str, optional
        Interpreter able to import TULIP (torch + transformers==4.32.1). If
        omitted, resolved from ``TULIP_PYTHON``, else the current interpreter.
    model_config : str, optional
        Path to the model-architecture JSON. Defaults to
        ``<tulip_home>/configs/shallow.config.json``, which matches the
        released weights.
    weights : str, optional
        Path to the pretrained ``state_dict``. Defaults to
        ``<tulip_home>/model_weights/pytorch_model.bin``.
    batch_size : int
        Passed through to ``predict.py``.

    Notes
    -----
    TULIP-TCR is GPLv3; this wrapper vendors none of it and only runs a
    user-provided installation out-of-process.
    """

    def __init__(
            self,
            tulip_home=None,
            tulip_python=None,
            model_config=None,
            weights=None,
            batch_size=512):
        self.tulip_home = _find_tulip_home(tulip_home)
        self.tulip_python = _find_tulip_python(tulip_python)
        if model_config is None:
            model_config = os.path.join(
                self.tulip_home, "configs", "shallow.config.json")
        if weights is None:
            weights = os.path.join(
                self.tulip_home, "model_weights", "pytorch_model.bin")
        for label, path in (("model_config", model_config), ("weights", weights)):
            if not os.path.isfile(path):
                raise FileNotFoundError(
                    "TULIP %s file not found: %s" % (label, path))
        self.model_config = model_config
        self.weights = weights
        self.batch_size = batch_size

    def __str__(self):
        return "Tulip(home=%s, python=%s)" % (self.tulip_home, self.tulip_python)

    def __repr__(self):
        return str(self)

    def _predictor_name(self):
        return "tulip"

    def kind_support(self):
        return {
            Kind.pMHC_TCR_binding: {
                # TULIP takes the presenting MHC allele as an input.
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            }
        }

    @property
    def supported_kinds(self):
        return tuple(self.kind_support())

    # ------------------------------------------------------------------
    # Subprocess bridge to TULIP's own predict.py
    # ------------------------------------------------------------------

    @staticmethod
    def _clean(value):
        """Normalize a sequence/allele cell: strip, uppercase-safe, and map
        empty values to TULIP's missing token."""
        if value is None:
            return MISSING_TOKEN
        text = str(value).strip()
        return text if text else MISSING_TOKEN

    def _run_predict(self, rows):
        """Run TULIP ``predict.py`` on *rows* and return a score dict.

        *rows* is a list of ``(peptide, mhc, cdr3a, cdr3b)`` tuples (already
        cleaned). Returns ``{(peptide, mhc, cdr3a, cdr3b): score}``. Because
        predict.py writes one CSV per unique peptide, in the (unshuffled)
        input order of the rows carrying that peptide, we map scores back by
        position within each peptide group — robust even to duplicate CDR3s.
        """
        if not rows:
            return {}
        with tempfile.TemporaryDirectory(prefix="mhctools_tulip_") as tmp:
            input_csv = os.path.join(tmp, "input.csv")
            # predict.py builds output paths as ``<output><peptide>.csv``, so
            # the output prefix must end in a path separator.
            output_prefix = os.path.join(tmp, "out") + os.sep
            os.makedirs(output_prefix, exist_ok=True)

            with open(input_csv, "w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["CDR3a", "CDR3b", "peptide", "MHC", "binder"])
                for peptide, mhc, cdr3a, cdr3b in rows:
                    writer.writerow([cdr3a, cdr3b, peptide, mhc, 1])

            cmd = [
                self.tulip_python,
                os.path.join(self.tulip_home, "predict.py"),
                "--test_dir", input_csv,
                "--modelconfig", self.model_config,
                "--load", self.weights,
                "--output", output_prefix,
                "--batch_size", str(self.batch_size),
            ]
            proc = subprocess.run(
                cmd,
                cwd=self.tulip_home,   # aatok/, mhctok/, src/ resolve here
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    "TULIP predict.py failed (exit %d).\nCommand: %s\n"
                    "stderr:\n%s" % (
                        proc.returncode, " ".join(cmd), proc.stderr[-4000:]))

            # Group input rows by peptide, preserving order, to align with
            # predict.py's per-peptide output ordering.
            order = {}
            for row in rows:
                order.setdefault(row[0], []).append(row)

            scores = {}
            for peptide, group in order.items():
                out_csv = os.path.join(output_prefix, peptide + ".csv")
                if not os.path.isfile(out_csv):
                    raise RuntimeError(
                        "TULIP produced no output for peptide %r (expected %s)."
                        " predict.py stdout tail:\n%s"
                        % (peptide, out_csv, proc.stdout[-2000:]))
                out = pd.read_csv(out_csv)
                if len(out) != len(group):
                    raise RuntimeError(
                        "TULIP output row count %d != input %d for peptide %r"
                        % (len(out), len(group), peptide))
                for (pep, mhc, cdr3a, cdr3b), score in zip(
                        group, out["score"].tolist()):
                    scores[(pep, mhc, cdr3a, cdr3b)] = float(score)
            return scores

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict_pairs(self, pairs):
        """Score explicit ``(peptide, TCR)`` or ``(peptide, TCR, mhc)`` items.

        Parameters
        ----------
        pairs : iterable
            Each element is ``(peptide, TCR)`` or ``(peptide, TCR, mhc)``.
            When the MHC allele is omitted it is treated as missing
            (TULIP's ``<MIS>``), i.e. MHC-agnostic scoring.

        Returns
        -------
        list of PeptideResult
            One :class:`PeptideResult` per input item, in order, each holding
            a single :class:`Prediction`.
        """
        items = []
        for pair in pairs:
            if len(pair) == 3:
                peptide, tcr, mhc = pair
            elif len(pair) == 2:
                (peptide, tcr), mhc = pair, None
            else:
                raise ValueError(
                    "Expected (peptide, TCR) or (peptide, TCR, mhc), got %r"
                    % (pair,))
            if not isinstance(tcr, TCR):
                raise TypeError(
                    "Expected mhctools.TCR instances, got %r" % type(tcr))
            items.append((peptide, tcr, mhc))

        # Build cleaned rows and remember each item's lookup key.
        rows = []
        keys = []
        for peptide, tcr, mhc in items:
            key = (
                self._clean(peptide),
                self._clean(mhc),
                self._clean(tcr.cdr3a),
                self._clean(tcr.cdr3b))
            keys.append(key)
            rows.append(key)

        scores = self._run_predict(list(dict.fromkeys(rows)))

        name = self._predictor_name()
        results = []
        for (peptide, tcr, mhc), key in zip(items, keys):
            allele = "" if key[1] == MISSING_TOKEN else key[1]
            results.append(PeptideResult(preds=(Prediction(
                kind=Kind.pMHC_TCR_binding,
                score=scores[key],
                peptide=peptide,
                allele=allele,
                tcr=tcr.identifier,
                predictor_name=name,
            ),)))
        return results

    def predict(self, peptides, tcrs, mhc=None):
        """Score every ``peptide × TCR`` combination at a given MHC allele.

        Parameters
        ----------
        peptides : str or list of str
        tcrs : TCR or list of TCR
        mhc : str or list of str, optional
            Presenting MHC allele(s). A single string applies to all
            peptides; a list must be parallel to *peptides*. Omit for
            MHC-agnostic scoring.

        Returns
        -------
        list of PeptideResult
            One :class:`PeptideResult` per peptide; each holds one
            :class:`Prediction` per TCR.
        """
        if isinstance(peptides, str):
            peptides = [peptides]
        peptides = list(peptides)
        if isinstance(tcrs, TCR):
            tcrs = [tcrs]
        tcrs = list(tcrs)

        if mhc is None or isinstance(mhc, str):
            mhc_list = [mhc] * len(peptides)
        else:
            mhc_list = list(mhc)
            if len(mhc_list) != len(peptides):
                raise ValueError(
                    "mhc list length %d != peptides length %d"
                    % (len(mhc_list), len(peptides)))

        flat = []
        for peptide, allele in zip(peptides, mhc_list):
            for tcr in tcrs:
                flat.append((peptide, tcr, allele))

        flat_results = self.predict_pairs(flat)

        results = []
        idx = 0
        for _ in peptides:
            preds = []
            for _ in tcrs:
                preds.extend(flat_results[idx].preds)
                idx += 1
            results.append(PeptideResult(preds=tuple(preds)))
        return results

    def predict_dataframe(self, peptides, tcrs, mhc=None, sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        dfs = [pp.to_dataframe(sample_name)
               for pp in self.predict(peptides, tcrs, mhc=mhc)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)
