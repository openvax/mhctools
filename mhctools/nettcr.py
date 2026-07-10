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

"""Wrapper for NetTCR-2.2 (https://github.com/mnielLab/NetTCR-2.2).

NetTCR-2.2 predicts pMHC:TCR binding — whether a paired αβ T-cell receptor
(given by its six CDR loops) recognises a given (class-I) peptide. This is a
different modality from the MHC-ligand predictors in the rest of mhctools:
the input is ``(peptide, TCR)`` rather than ``(peptide, allele)``, and the
output kind is :attr:`~mhctools.pred.Kind.pMHC_TCR_binding`.

Unlike the DTU ``netMHC*`` tools, NetTCR ships its pretrained weights
directly in its git repository as small TFLite models (the pan ensemble is
20 models of ~0.4 MB each). This wrapper runs them **in-process** through a
TFLite interpreter — it does *not* need NetTCR's conda environment. We only
require the cloned repository (for the weights) and a TFLite runtime
(``ai-edge-litert``, ``tflite-runtime``, or ``tensorflow``).

The BLOSUM50 encoding, per-feature padding lengths, and name-based input
tensor assignment here reproduce ``src/predict.py`` from NetTCR-2.2. Scores are
the mean over the pan cross-validation ensemble, matching the published usage.
"""

from __future__ import annotations

import contextlib
import glob
import os
import sys

import numpy as np
import pandas as pd

from .pred import COLUMNS, Kind, PeptideResult, Prediction
from .tcr import TCR


@contextlib.contextmanager
def _suppress_native_stderr():
    """Silence native (C-level) writes to fd 2 for the duration of the block.

    The TFLite runtimes print ``INFO: Created TensorFlow Lite XNNPACK delegate
    for CPU.`` straight to fd 2 during interpreter setup, bypassing Python's
    ``logging``. Redirecting fd 2 is process-global, so this wraps only
    construction / allocation, never the prediction loop; Python exceptions
    still propagate.
    """
    sys.stderr.flush()
    # Acquire both descriptors inside the try so a failure partway through
    # (e.g. fd exhaustion at os.open) can't leak the one already taken.
    saved_fd = None
    devnull_fd = None
    try:
        saved_fd = os.dup(2)
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_fd, 2)
        yield
    finally:
        if saved_fd is not None:
            os.dup2(saved_fd, 2)
            os.close(saved_fd)
        if devnull_fd is not None:
            os.close(devnull_fd)


# BLOSUM50, restricted to the 20 standard amino acids, in NetTCR's column
# order (A R N D C Q E G H I L K M F P S T W Y V). Copied verbatim from
# NetTCR-2.2 `src/keras_utils.py` (`blosum50_20aa`).
_BLOSUM50_20AA = {
    'A': (5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0),
    'R': (-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3),
    'N': (-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3),
    'D': (-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4),
    'C': (-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1),
    'Q': (-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3),
    'E': (-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3),
    'G': (0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4),
    'H': (-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4),
    'I': (-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4),
    'L': (-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1),
    'K': (-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3),
    'M': (-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1),
    'F': (-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1),
    'P': (-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3),
    'S': (1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2),
    'T': (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0),
    'W': (-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3),
    'Y': (-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1),
    'V': (0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5),
}
_BLOSUM50 = {aa: np.array(vec, dtype=np.float32) for aa, vec in _BLOSUM50_20AA.items()}

# NetTCR encodes residues as BLOSUM50 rows, pads missing positions with a
# constant, then divides the whole array by this factor. Residues become
# BLOSUM50/5 and padded positions become -1 (= -5 / 5).
_PAD_VALUE = -5.0
_NORM_FACTOR = 5.0

# Per-feature maximum sequence lengths (from NetTCR-2.2 `src/predict.py`).
_FEATURE_MAX_LEN = {
    "pep": 12,
    "a1": 7,
    "a2": 8,
    "a3": 22,
    "b1": 6,
    "b2": 7,
    "b3": 23,
}


def _load_interpreter(model_path):
    """Load a TFLite interpreter from whichever runtime is installed.

    Prefers the lightweight LiteRT/tflite runtimes, falling back to the
    ``tensorflow.lite`` interpreter. All expose the same interface.

    Construction is wrapped in :func:`_suppress_native_stderr` as well as
    ``allocate_tensors`` because different runtimes apply (and log) the XNNPACK
    delegate at different points — ``tensorflow`` does it in
    ``allocate_tensors`` (verified), while the LiteRT runtimes may do it at
    construction. The import itself stays outside the suppression so a genuine
    ImportError still surfaces normally.
    """
    try:
        from ai_edge_litert.interpreter import Interpreter
    except ImportError:
        pass
    else:
        with _suppress_native_stderr():
            return Interpreter(model_path=model_path)
    try:
        from tflite_runtime.interpreter import Interpreter
    except ImportError:
        pass
    else:
        with _suppress_native_stderr():
            return Interpreter(model_path=model_path)
    try:
        import tensorflow as tf
    except ImportError as e:
        raise ImportError(
            "NetTCR needs a TFLite runtime. Install one of: "
            "`ai-edge-litert` (recommended, lightweight), `tflite-runtime`, "
            "or `tensorflow`.") from e
    with _suppress_native_stderr():
        return tf.lite.Interpreter(model_path=model_path)


def _find_nettcr_dir(nettcr_path=None):
    """Resolve the NetTCR-2.2 installation directory.

    Checks, in order:
    1. The *nettcr_path* argument
    2. The ``NETTCR_DIR`` environment variable
    3. ``~/NetTCR-2.2`` and ``~/code/NetTCR-2.2``

    An explicitly-provided path (argument or ``NETTCR_DIR``) is validated up
    front so a typo fails with a clear message rather than later when no
    models are found.
    """
    clone_hint = "Clone from https://github.com/mnielLab/NetTCR-2.2"
    for source, path in (
            ("nettcr_path argument", nettcr_path),
            ("NETTCR_DIR", os.environ.get("NETTCR_DIR"))):
        if path:
            if not os.path.isdir(path):
                raise FileNotFoundError(
                    "NetTCR-2.2 directory from %s does not exist: %s. %s"
                    % (source, path, clone_hint))
            return path
    home = os.path.expanduser("~")
    for candidate in (
            os.path.join(home, "NetTCR-2.2"),
            os.path.join(home, "code", "NetTCR-2.2")):
        if os.path.isdir(candidate):
            return candidate
    raise FileNotFoundError(
        "NetTCR-2.2 not found. Set NETTCR_DIR or pass nettcr_path= to the "
        "constructor. %s" % clone_hint)


def _encode_feature(sequences, feature):
    """BLOSUM50-encode a list of sequences for one NetTCR feature.

    Returns a ``float32`` array of shape ``(n, max_len, 20)`` where residues
    are BLOSUM50/5 and padded positions are -1, matching NetTCR's
    ``enc_list_bl_max_len(...) / 5``.
    """
    max_len = _FEATURE_MAX_LEN[feature]
    n = len(sequences)
    arr = _PAD_VALUE * np.ones((n, max_len, 20), dtype=np.float32)
    for i, seq in enumerate(sequences):
        seq = seq.upper()
        if len(seq) > max_len:
            raise ValueError(
                "NetTCR feature %r max length is %d, got %d-mer %r"
                % (feature, max_len, len(seq), seq))
        for j, aa in enumerate(seq):
            try:
                arr[i, j] = _BLOSUM50[aa]
            except KeyError:
                raise ValueError(
                    "Unknown amino acid %r in NetTCR feature %r sequence %r"
                    % (aa, feature, seq))
    arr /= _NORM_FACTOR
    return arr


class NetTCR(object):
    """Wrapper for NetTCR-2.2 pMHC:TCR binding predictions.

    Runs NetTCR-2.2's pan cross-validation ensemble in-process via a TFLite
    interpreter; the reported score is the mean over the ensemble. Models are
    loaded lazily on the first :meth:`predict` call and kept in memory.

    Parameters
    ----------
    nettcr_path : str, optional
        Path to the cloned NetTCR-2.2 repository root. If omitted, resolved
        from ``NETTCR_DIR`` or ``~/NetTCR-2.2`` / ``~/code/NetTCR-2.2``.
    checkpoint_dir : str, optional
        Directory of ``*.tflite`` ensemble models. Defaults to the pan model
        (``models/nettcr_2_2_pan/checkpoint``). Only the pan model
        generalises across arbitrary peptides; the peptide-specific and
        pretrained models are out of scope for this wrapper.

    Notes
    -----
    NetTCR-2.2 is distributed under an academic software license; this
    wrapper only *runs* a user-provided installation and vendors none of it.

    The cached ensemble interpreters are stateful, so a single ``NetTCR``
    instance is not safe to call from multiple threads concurrently; use one
    instance per thread.
    """

    def __init__(self, nettcr_path=None, checkpoint_dir=None):
        self.nettcr_dir = _find_nettcr_dir(nettcr_path)
        if checkpoint_dir is None:
            checkpoint_dir = os.path.join(
                self.nettcr_dir, "models", "nettcr_2_2_pan", "checkpoint")
        self.checkpoint_dir = checkpoint_dir
        self._model_paths = sorted(glob.glob(
            os.path.join(self.checkpoint_dir, "*.tflite")))
        if not self._model_paths:
            raise FileNotFoundError(
                "No NetTCR *.tflite models found in %s" % self.checkpoint_dir)
        # Lazy-loaded ensemble of interpreters, and the batch size the
        # interpreters are currently allocated for (so repeated same-size
        # calls skip re-resizing/re-allocating tensors).
        self._interpreters = None
        self._allocated_n = None

    def __str__(self):
        loaded = "loaded" if self._interpreters is not None else "not loaded"
        return "NetTCR(models=%d, %s)" % (len(self._model_paths), loaded)

    def __repr__(self):
        return str(self)

    def _predictor_name(self):
        return "nettcr"

    def kind_support(self):
        return {
            Kind.pMHC_TCR_binding: {
                # NetTCR takes no MHC allele as input; the peptide is
                # implicitly class-I-restricted.
                "mhc_dependence": "none",
                "mhc_class": "I",
            }
        }

    @property
    def supported_kinds(self):
        return tuple(self.kind_support())

    # ------------------------------------------------------------------
    # Model loading (lazy, cached)
    # ------------------------------------------------------------------

    def _ensure_loaded(self):
        if self._interpreters is None:
            self._interpreters = [
                _load_interpreter(path) for path in self._model_paths]

    # ------------------------------------------------------------------
    # Core prediction (in-process, batched over the ensemble)
    # ------------------------------------------------------------------

    def _predict_raw(self, peptides, tcrs):
        """Score parallel lists of peptides and :class:`TCR` objects.

        Returns a 1-D numpy array of ensemble-mean scores, one per pair.
        """
        n = len(peptides)
        if n == 0:
            return np.zeros(0, dtype=np.float32)
        self._ensure_loaded()

        encoded = {"pep": _encode_feature(peptides, "pep")}
        for key in ("a1", "a2", "a3", "b1", "b2", "b3"):
            encoded[key] = _encode_feature(
                [t.cdr_dict()[key] for t in tcrs], key)

        # Resize + allocate only when the batch size changed since last call.
        reallocate = self._allocated_n != n

        total = np.zeros(n, dtype=np.float64)
        for interpreter in self._interpreters:
            inputs = interpreter.get_input_details()
            output = interpreter.get_output_details()[0]
            if reallocate:
                for det in inputs:
                    interpreter.resize_tensor_input(
                        det["index"], [n, det["shape"][1], det["shape"][2]])
                interpreter.resize_tensor_input(
                    output["index"], [n, output["shape"][1]])
                # allocate_tensors() is where the TFLite runtime prints its
                # one-time "Created ... XNNPACK delegate for CPU." INFO line to
                # native stderr; keep that out of callers' output.
                with _suppress_native_stderr():
                    interpreter.allocate_tensors()
            for det in inputs:
                # NetTCR names inputs "serving_default_<feature>:0"; match by
                # the trailing feature token rather than by tensor order.
                key = det["name"].split(":")[0].split("_")[-1]
                try:
                    tensor = encoded[key]
                except KeyError:
                    raise ValueError(
                        "NetTCR model has an unexpected input tensor %r "
                        "(parsed feature %r); expected one of %s"
                        % (det["name"], key, sorted(encoded)))
                interpreter.set_tensor(det["index"], tensor)
            interpreter.invoke()
            total += interpreter.get_tensor(output["index"]).reshape(n)
        self._allocated_n = n
        return (total / len(self._interpreters)).astype(np.float32)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @staticmethod
    def _tcrs_from_dataframe(df, cdr_cols=None):
        default_cols = {
            "cdr1a": "cdr1a",
            "cdr2a": "cdr2a",
            "cdr3a": "cdr3a",
            "cdr1b": "cdr1b",
            "cdr2b": "cdr2b",
            "cdr3b": "cdr3b",
        }
        if cdr_cols is None:
            cdr_cols = default_cols
        else:
            unknown = sorted(set(cdr_cols) - set(default_cols))
            if unknown:
                raise ValueError(
                    "Unknown NetTCR CDR field(s): %s"
                    % ", ".join(unknown))
            cdr_cols = {**default_cols, **cdr_cols}

        missing = [
            column for column in cdr_cols.values()
            if column not in df.columns
        ]
        if missing:
            raise ValueError(
                "NetTCR dataframe missing CDR column(s): %s"
                % ", ".join(sorted(missing)))

        tcrs = []
        for row in df.itertuples(index=False):
            row_values = row._asdict()
            tcrs.append(TCR.from_dict({
                cdr: row_values[column]
                for cdr, column in cdr_cols.items()
            }))
        return tcrs

    def predict_pairs(self, pairs):
        """Score explicit ``(peptide, TCR)`` pairs.

        Parameters
        ----------
        pairs : iterable of (str, TCR)

        Returns
        -------
        list of PeptideResult
            One :class:`PeptideResult` per input pair (in order), each
            holding a single :class:`Prediction`.
        """
        pairs = list(pairs)
        peptides = [pep for pep, _ in pairs]
        tcrs = [tcr for _, tcr in pairs]
        for tcr in tcrs:
            if not isinstance(tcr, TCR):
                raise TypeError(
                    "Expected mhctools.TCR instances, got %r" % type(tcr))
        scores = self._predict_raw(peptides, tcrs)
        name = self._predictor_name()
        results = []
        for pep, tcr, score in zip(peptides, tcrs, scores):
            results.append(PeptideResult(preds=(Prediction(
                kind=Kind.pMHC_TCR_binding,
                score=float(score),
                peptide=pep,
                tcr=tcr.identifier,
                predictor_name=name,
            ),)))
        return results

    def predict(self, peptides, tcrs):
        """Score every ``peptide × TCR`` combination.

        Parameters
        ----------
        peptides : str or list of str
        tcrs : TCR or list of TCR

        Returns
        -------
        list of PeptideResult
            One :class:`PeptideResult` per peptide; each holds one
            :class:`Prediction` per TCR.
        """
        if isinstance(peptides, str):
            peptides = [peptides]
        if isinstance(tcrs, TCR):
            tcrs = [tcrs]
        tcrs = list(tcrs)

        flat_peptides = []
        flat_tcrs = []
        for pep in peptides:
            for tcr in tcrs:
                flat_peptides.append(pep)
                flat_tcrs.append(tcr)

        flat_results = self.predict_pairs(zip(flat_peptides, flat_tcrs))

        results = []
        idx = 0
        for _ in peptides:
            preds = []
            for _ in tcrs:
                preds.extend(flat_results[idx].preds)
                idx += 1
            results.append(PeptideResult(preds=tuple(preds)))
        return results

    def predict_dataframe(
            self,
            peptides,
            tcrs=None,
            sample_name="",
            peptide_col="peptide",
            cdr_cols=None):
        """Predict flattened to a DataFrame.

        Existing peptide-list calls still use ``predict(peptides, tcrs)``.
        When *peptides* is a DataFrame, each row is treated as one
        ``(peptide, TCR)`` pair and TCRs are built from CDR columns.
        """
        if isinstance(peptides, pd.DataFrame):
            if tcrs is not None:
                raise ValueError(
                    "tcrs must be omitted when peptides is a DataFrame")
            if peptide_col not in peptides.columns:
                raise ValueError(
                    "NetTCR dataframe missing peptide column %r"
                    % peptide_col)
            peptide_list = peptides[peptide_col].tolist()
            row_tcrs = self._tcrs_from_dataframe(peptides, cdr_cols=cdr_cols)
            dfs = [pp.to_dataframe(sample_name)
                   for pp in self.predict_pairs(zip(peptide_list, row_tcrs))]
        else:
            if tcrs is None:
                raise ValueError(
                    "tcrs is required unless peptides is a DataFrame")
            dfs = [pp.to_dataframe(sample_name)
                   for pp in self.predict(peptides, tcrs)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)
