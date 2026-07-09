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

"""Wrapper for NetCleave (https://github.com/BSC-CNS-EAPM/NetCleave).

NetCleave predicts **C-terminal antigen processing** for both the MHC-I
(proteasomal) and MHC-II (endolysosomal) pathways. It scores the cleavage
site (4 residues of the peptide C-terminus + 3 downstream residues) using a
neural network over QSAR descriptors, and emits **one score per peptide**
(unlike NetChop/pepsickle, which emit a per-position profile).

Because the output granularity is per-peptide rather than per-position, this
wrapper does not subclass :class:`ProcessingPredictor`; it shells out to
NetCleave's ``NetCleave.py --predict`` CLI (the NetChop-style subprocess
model) and parses the resulting CSV. The class-I and class-II models emit
different kinds:

* class I -> :attr:`~mhctools.pred.Kind.proteasome_cleavage`
* class II -> :attr:`~mhctools.pred.Kind.endolysosomal_cleavage`

NetCleave ships its pretrained weights (small Keras ``.h5`` files) directly
in its git repository, so only a ``git clone`` is required — no separate
download. Its ``--predict`` path needs Python (tensorflow/keras, scikit-learn,
biopython, pandas, numpy); the R dependency in NetCleave's README is only for
its training-data generation and is **not** used for prediction.

.. warning::
    NetCleave's own paper reports class-II C-terminal cleavage is a much
    weaker signal than class I (AUC ~0.66 vs ~0.91): the class-II peptide
    C-terminus lies outside the binding groove and cathepsin specificity is
    diffuse. Treat ``endolysosomal_cleavage`` scores accordingly.
"""

import logging
import os
import subprocess
import sys
import tempfile
import uuid
from collections import defaultdict

import pandas as pd

from .base_predictor import (
    _check_flank_inputs,
    _normalize_sequence_dict,
    _peptide_contexts,
)
from .pred import COLUMNS, Kind, PeptideResult, Prediction

logger = logging.getLogger(__name__)

NETCLEAVE_TIMEOUT_SECONDS = 600

# NetCleave needs 3 residues downstream of the peptide C-terminus to build
# the 4+3 cleavage site; without them it cannot score the site.
_C_FLANK_REQUIRED = 3


def _find_netcleave_dir(netcleave_path=None):
    """Resolve the NetCleave installation directory.

    Checks, in order: the *netcleave_path* argument, the ``NETCLEAVE_DIR``
    environment variable, then ``~/NetCleave`` and ``~/code/NetCleave``. An
    explicitly-provided path is validated up front.
    """
    clone_hint = "Clone from https://github.com/BSC-CNS-EAPM/NetCleave"
    for source, path in (
            ("netcleave_path argument", netcleave_path),
            ("NETCLEAVE_DIR", os.environ.get("NETCLEAVE_DIR"))):
        if path:
            if not os.path.isdir(path):
                raise FileNotFoundError(
                    "NetCleave directory from %s does not exist: %s. %s"
                    % (source, path, clone_hint))
            return path
    home = os.path.expanduser("~")
    for candidate in (
            os.path.join(home, "NetCleave"),
            os.path.join(home, "code", "NetCleave")):
        if os.path.isdir(candidate):
            return candidate
    raise FileNotFoundError(
        "NetCleave not found. Set NETCLEAVE_DIR or pass netcleave_path= to "
        "the constructor. %s" % clone_hint)


class NetCleave(object):
    """Wrapper for NetCleave C-terminal cleavage predictions.

    Parameters
    ----------
    mhc_class : str
        ``"I"`` (proteasomal) or ``"II"`` (endolysosomal). Selects the
        default model and the emitted :class:`~mhctools.pred.Kind`.
    mhc_allele : str
        Allele label used to pick the pretrained model directory
        (``data/models/{class}_mass-spectrometry_{allele}``). Default
        ``"HLA"`` (the pan-allele model for the class).
    netcleave_path : str, optional
        Path to the cloned NetCleave repository. Resolved from
        ``NETCLEAVE_DIR`` / ``~/NetCleave`` when omitted.
    model_path : str, optional
        Full path to a specific model directory, overriding
        ``mhc_class`` / ``mhc_allele`` selection.
    python_executable : str, optional
        Interpreter used to run ``NetCleave.py``. Defaults to the current
        interpreter; override if NetCleave's dependencies live in a separate
        environment.
    subprocess_timeout : int
        Timeout (seconds) for a single NetCleave invocation.

    Notes
    -----
    NetCleave is distributed under GPL-v2; this wrapper only *runs* a
    user-provided installation and vendors none of it.

    A ``NetCleave`` instance is not safe to call from multiple threads
    concurrently (each call shells out via a per-call temp file); use one
    instance per thread.
    """

    VALID_CLASSES = ("I", "II")

    def __init__(self, mhc_class="I", mhc_allele="HLA", netcleave_path=None,
                 model_path=None, python_executable=None,
                 subprocess_timeout=NETCLEAVE_TIMEOUT_SECONDS):
        if mhc_class not in self.VALID_CLASSES:
            raise ValueError(
                "mhc_class must be one of %s, got %r"
                % (self.VALID_CLASSES, mhc_class))
        self.mhc_class = mhc_class
        self.mhc_allele = mhc_allele
        self.netcleave_dir = _find_netcleave_dir(netcleave_path)
        self.python_executable = python_executable or sys.executable
        self.subprocess_timeout = subprocess_timeout

        self._script = os.path.join(self.netcleave_dir, "NetCleave.py")
        if not os.path.isfile(self._script):
            raise FileNotFoundError(
                "NetCleave.py not found in %s" % self.netcleave_dir)

        if model_path is None:
            model_path = os.path.join(
                self.netcleave_dir, "data", "models",
                "%s_mass-spectrometry_%s" % (mhc_class, mhc_allele))
        self.model_path = model_path
        weights = os.path.join(
            model_path, "%s_model.h5" % os.path.basename(model_path))
        if not os.path.isfile(weights):
            raise FileNotFoundError(
                "NetCleave model weights not found: %s. Use --mhc_options to "
                "list available models, or pass model_path=." % weights)

        self._call_counter = 0

    def __str__(self):
        return "NetCleave(mhc_class=%s, model=%s)" % (
            self.mhc_class, os.path.basename(self.model_path))

    def __repr__(self):
        return str(self)

    def _predictor_name(self):
        return "netcleave"

    def _pred_kind(self):
        return (Kind.proteasome_cleavage if self.mhc_class == "I"
                else Kind.endolysosomal_cleavage)

    def kind_support(self):
        return {
            self._pred_kind(): {
                # Cleavage is MHC-independent (no allele in the prediction).
                "mhc_dependence": "none",
                "mhc_class": self.mhc_class,
            }
        }

    @property
    def supported_kinds(self):
        return tuple(self.kind_support())

    # ------------------------------------------------------------------
    # Subprocess
    # ------------------------------------------------------------------

    def _run_netcleave(self, requests):
        """Run NetCleave (pred_input 3) on a list of scoring requests.

        Each request is ``(epitope, protein_seq, expected_site)``: the
        peptide, the sequence to search it in, and the expected 7-residue
        cleavage site (the peptide's C-terminal 4 residues + 3 downstream).

        NetCleave emits **one output row per regex occurrence** of the
        epitope in ``protein_seq`` (and drops occurrences without 3
        downstream residues), so results are matched back to requests by
        input-row id and cleavage site — never by row position.

        Returns a list of scores (float or None) aligned to *requests*;
        None means no matching valid cleavage site was produced.
        """
        if not requests:
            return []
        self._call_counter += 1
        # Basename must be unique (a shared NetCleave output/ dir is used by
        # every instance) and contain no '.' before the extension — NetCleave
        # derives the output filename via ``basename.split('.')[0]``.
        basename = "mhctools_netcleave_%d_%s" % (os.getpid(), uuid.uuid4().hex)
        tmp_dir = tempfile.mkdtemp(prefix="mhctools_netcleave_")
        input_csv = os.path.join(tmp_dir, basename + ".csv")
        output_csv = os.path.join(
            self.netcleave_dir, "output", basename + "_NetCleave.csv")

        pd.DataFrame({
            "epitope": [r[0].upper() for r in requests],
            "protein_seq": [r[1].upper() for r in requests],
            "protein_name": [str(i) for i in range(len(requests))],
        }).to_csv(input_csv, index=False)

        try:
            try:
                result = subprocess.run(
                    [self.python_executable, self._script,
                     "--predict", input_csv, "--pred_input", "3",
                     "--model_path", self.model_path,
                     "--mhc_class", self.mhc_class],
                    cwd=self.netcleave_dir,
                    capture_output=True,
                    timeout=self.subprocess_timeout,
                )
            except subprocess.TimeoutExpired as e:
                raise RuntimeError(
                    "NetCleave timed out after %d seconds on %d epitopes"
                    % (self.subprocess_timeout, len(requests))) from e
            except OSError as e:
                raise RuntimeError(
                    "Could not run NetCleave with %r: %s"
                    % (self.python_executable, e)) from e

            stderr_text = result.stderr.decode("utf-8", errors="replace").strip()
            if result.returncode != 0:
                raise RuntimeError(
                    "NetCleave exited with code %d.\nstdout: %s\nstderr: %s"
                    % (result.returncode,
                       result.stdout.decode("utf-8", errors="replace").strip(),
                       stderr_text))
            if not os.path.isfile(output_csv):
                raise RuntimeError(
                    "NetCleave produced no output file at %s.\nstderr: %s"
                    % (output_csv, stderr_text))

            out = pd.read_csv(output_csv)
            # Map input-row id -> {cleavage_site (upper): score}. NetCleave
            # carries our per-row protein_name through as ``uniprot_id`` and
            # the 4+3 site as ``cleavage_site``; invalid sites read as NaN.
            by_row = defaultdict(dict)
            row_epitope = {}
            for rid, epi, site, value in zip(
                    out["uniprot_id"], out["epitope"],
                    out["cleavage_site"], out["prediction"]):
                try:
                    score = float(value)
                except (TypeError, ValueError):
                    score = None
                if score is not None and score != score:  # NaN
                    score = None
                by_row[str(rid)][str(site).upper()] = score
                row_epitope.setdefault(str(rid), str(epi).upper())

            scores = []
            for i, (epitope, _protein_seq, expected_site) in enumerate(requests):
                key = str(i)
                # Guard against a NetCleave contract change silently
                # misaligning rows: rows for this id must be our epitope.
                seen = row_epitope.get(key)
                if seen is not None and seen != epitope.upper():
                    raise RuntimeError(
                        "NetCleave row %s epitope %r does not match the "
                        "requested peptide %r" % (key, seen, epitope))
                scores.append(by_row.get(key, {}).get(expected_site.upper()))
            return scores
        finally:
            for path in (input_csv, output_csv):
                try:
                    os.remove(path)
                except OSError:
                    pass
            try:
                os.rmdir(tmp_dir)
            except OSError:
                pass

    @staticmethod
    def _expected_site(peptide, c_flank):
        """The 7-residue cleavage site NetCleave scores: peptide C-terminal
        4 residues + 3 downstream. Empty if there isn't enough context."""
        if len(peptide) < 4 or len(c_flank) < _C_FLANK_REQUIRED:
            return ""
        return (peptide[-4:] + c_flank[:_C_FLANK_REQUIRED]).upper()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict C-terminal cleavage scores for peptides.

        NetCleave needs the residues **downstream** of the peptide to build
        the cleavage site, so ``c_flanks`` (>= 3 residues each) is required.
        ``n_flanks`` is accepted for interface parity and recorded on each
        :class:`Prediction`, but does not affect the C-terminal score.

        Parameters
        ----------
        peptides : list of str
        n_flanks : list of str, optional
        c_flanks : list of str
            C-terminal flanks (>= 3 residues) — one per peptide.

        Returns
        -------
        list of PeptideResult
            One per peptide; empty (no preds) for peptides whose cleavage
            site could not be built.
        """
        peptide_list, n_flank_list, c_flank_list = _check_flank_inputs(
            peptides, n_flanks, c_flanks)
        if not peptide_list:
            return []
        if c_flank_list is None:
            raise ValueError(
                "NetCleave.predict requires c_flanks (>= %d residues per "
                "peptide) to build the C-terminal cleavage site. To score "
                "peptides in a protein context, use predict_proteins()."
                % _C_FLANK_REQUIRED)

        requests, meta = [], []
        for i, peptide in enumerate(peptide_list):
            n_flank = n_flank_list[i] if n_flank_list is not None else ""
            c_flank = c_flank_list[i]
            expected_site = self._expected_site(peptide, c_flank)
            if not expected_site:
                logger.warning(
                    "NetCleave: peptide %r has c_flank %r shorter than %d "
                    "residues; cannot score its cleavage site",
                    peptide, c_flank, _C_FLANK_REQUIRED)
                meta.append((peptide, n_flank, c_flank, None))
                continue
            # Score each peptide in its own surrogate protein so its cleavage
            # site is well-defined regardless of repeats elsewhere.
            meta.append((peptide, n_flank, c_flank, len(requests)))
            requests.append((peptide, peptide + c_flank, expected_site))

        scores = self._run_netcleave(requests)

        results = []
        for peptide, n_flank, c_flank, req_idx in meta:
            score = None if req_idx is None else scores[req_idx]
            results.append(self._make_result(
                peptide, score, n_flank=n_flank, c_flank=c_flank))
        return results

    def predict_proteins(self, sequence_dict, peptide_lengths=None,
                         flank_length=3):
        """Score peptides scanned from full protein sequences.

        Each peptide is scored in its real protein context, so NetCleave sees
        the true downstream residues. This is the preferred entry-point.

        Parameters
        ----------
        sequence_dict : dict or str
        peptide_lengths : list of int, optional
            Default: ``[9]`` for class I, ``[15]`` for class II.
        flank_length : int
            Residues of flank recorded on each Prediction (default 3).

        Returns
        -------
        dict mapping sequence_name -> list of PeptideResult
        """
        sequence_dict = _normalize_sequence_dict(sequence_dict)
        if peptide_lengths is None:
            peptide_lengths = [9] if self.mhc_class == "I" else [15]
        if isinstance(peptide_lengths, int):
            peptide_lengths = [peptide_lengths]

        # Ensure at least 3 downstream residues are captured for the cleavage
        # site, even if the caller asks for a shorter recorded flank.
        contexts = _peptide_contexts(
            sequence_dict, peptide_lengths, flank_length,
            n_flank_length=flank_length,
            c_flank_length=max(flank_length, _C_FLANK_REQUIRED))

        # Score each peptide in a surrogate protein (peptide + its downstream
        # flank). The cleavage site depends only on those residues, so this
        # matches the full-protein score while keeping each peptide unique in
        # its row (NetCleave emits one row per regex occurrence).
        requests, meta = [], []
        for context in contexts:
            expected_site = self._expected_site(context.peptide, context.c_flank)
            if not expected_site:
                meta.append((context, None))
                continue
            meta.append((context, len(requests)))
            requests.append((
                context.peptide,
                context.peptide + context.c_flank,
                expected_site))

        scores = self._run_netcleave(requests)

        results = defaultdict(list)
        for context, req_idx in meta:
            score = None if req_idx is None else scores[req_idx]
            results[context.source_sequence_name].append(self._make_result(
                context.peptide, score,
                n_flank=context.n_flank, c_flank=context.c_flank,
                source_sequence_name=context.source_sequence_name,
                offset=context.offset))
        return dict(results)

    def _make_result(self, peptide, score, n_flank="", c_flank="",
                     source_sequence_name=None, offset=0):
        if score is None:
            return PeptideResult(preds=())
        return PeptideResult(preds=(Prediction(
            kind=self._pred_kind(),
            score=score,
            peptide=peptide,
            n_flank=n_flank,
            c_flank=c_flank,
            source_sequence_name=source_sequence_name,
            offset=offset,
            predictor_name=self._predictor_name(),
        ),))

    def predict_dataframe(self, peptides, n_flanks=None, c_flanks=None,
                          sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        dfs = [pp.to_dataframe(sample_name)
               for pp in self.predict(peptides, n_flanks, c_flanks)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict_proteins_dataframe(self, sequence_dict, peptide_lengths=None,
                                   flank_length=3, sample_name=""):
        """``predict_proteins()`` flattened to a DataFrame."""
        dfs = []
        for pp_list in self.predict_proteins(
                sequence_dict, peptide_lengths, flank_length).values():
            for pp in pp_list:
                dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)


class NetCleave_I(NetCleave):
    """NetCleave for the MHC-I (proteasomal) pathway."""
    def __init__(self, mhc_allele="HLA", **kwargs):
        NetCleave.__init__(self, mhc_class="I", mhc_allele=mhc_allele, **kwargs)


class NetCleave_II(NetCleave):
    """NetCleave for the MHC-II (endolysosomal) pathway."""
    def __init__(self, mhc_allele="HLA", **kwargs):
        NetCleave.__init__(self, mhc_class="II", mhc_allele=mhc_allele, **kwargs)
