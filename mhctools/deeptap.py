# Copyright (c) 2016. Mount Sinai School of Medicine
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

"""Wrapper for DeepTAP — a BiGRU predictor of TAP transport / binding.

TAP (transporter associated with antigen processing) shuttles cytosolic
peptides into the ER for MHC-I loading, and is a distinct step of the
processing pathway from proteasomal cleavage (NetChop / NetCleave / Pepsickle).
mhctools otherwise has no TAP predictor; DeepTAP fills that gap and emits
``Kind.tap_transport``.

Like the cleavage predictors, DeepTAP is allele-independent: it scores each
peptide once, so ``predict()`` returns one prediction per peptide (with an empty
``allele``).

DeepTAP ships its pretrained weights in-repo (Apache-2.0) but pins an old stack
(pytorch-lightning==1.9.2, torch). To keep that out of the mhctools environment,
this wrapper shells out to DeepTAP's own ``deeptap.py`` CLI in a user-provided
checkout (``DEEPTAP_HOME``), run by a user-provided interpreter
(``DEEPTAP_PYTHON``, default the current one); the checkpoints also load under
modern Lightning.

Upstream: https://github.com/zjupgx/DeepTAP
Cite: Chen et al., Comput. Biol. Med. 2023 — "DeepTAP: An RNN-based method of
TAP-binding peptide prediction in the selection of tumor neoantigens".

DeepTAP's evaluation is self-reported and no independent TAP benchmark exists —
a useful pathway signal, not a validated oracle.
"""

import os
import sys
from os.path import exists, isdir, isfile, join
from tempfile import mkdtemp

import pandas as pd

from .cleanup_context import CleanupFiles
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command
from .wrapper_base import AlleleFreePredictor

# DeepTAP one-hot encodes the 20 standard amino acids plus "X" (unknown /
# padding) and pads every peptide to a fixed width of 17; longer peptides
# overflow its input tensor, so 17 is a hard model limit.
DEEPTAP_MAX_PEPTIDE_LENGTH = 17
_VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYX")

_TASK_TYPES = ("cla", "reg")


def _find_deeptap_home(deeptap_home=None):
    """Resolve the DeepTAP checkout directory (holds ``deeptap.py`` + ``model/``).

    Checks, in order: the *deeptap_home* argument, ``$DEEPTAP_HOME``, then
    ``~/DeepTAP``.
    """
    candidate = deeptap_home or os.environ.get("DEEPTAP_HOME")
    if not candidate:
        home = join(os.path.expanduser("~"), "DeepTAP")
        if isdir(home):
            candidate = home
    if not candidate:
        raise FileNotFoundError(
            "DeepTAP not found. Set DEEPTAP_HOME or pass deeptap_home= to the "
            "constructor. Clone from https://github.com/zjupgx/DeepTAP")
    if not isfile(join(candidate, "deeptap.py")):
        raise FileNotFoundError(
            "deeptap.py not found in %r — is this a DeepTAP checkout?"
            % candidate)
    return candidate


class DeepTAP(AlleleFreePredictor):
    """Wrapper for the DeepTAP TAP-transport predictor.

    Parameters
    ----------
    task_type : str
        ``"cla"`` (classification; ``score`` in 0-1, higher = stronger TAP
        binding) or ``"reg"`` (regression; also fills ``value`` with the
        predicted affinity in nM, lower = stronger). Default ``"cla"``.
    deeptap_home : str, optional
        Path to a DeepTAP checkout. Resolved from the argument, then
        ``$DEEPTAP_HOME``, then ``~/DeepTAP``.
    deeptap_python : str, optional
        Interpreter that can run DeepTAP (has torch + pytorch-lightning).
        Resolved from the argument, then ``$DEEPTAP_PYTHON``, then the current
        interpreter (``sys.executable``).
    max_peptide_length : int
        Peptides longer than this are rejected. Capped at DeepTAP's hard limit
        of 17. Default 17.
    """

    def __init__(
            self,
            task_type="cla",
            deeptap_home=None,
            deeptap_python=None,
            max_peptide_length=DEEPTAP_MAX_PEPTIDE_LENGTH):
        if task_type not in _TASK_TYPES:
            raise ValueError(
                "task_type must be one of %s, got %r"
                % (_TASK_TYPES, task_type))
        if max_peptide_length > DEEPTAP_MAX_PEPTIDE_LENGTH:
            raise ValueError(
                "DeepTAP pads peptides to %d residues; max_peptide_length "
                "cannot exceed that (got %d)"
                % (DEEPTAP_MAX_PEPTIDE_LENGTH, max_peptide_length))
        self.task_type = task_type
        self.deeptap_home = _find_deeptap_home(deeptap_home)
        self.deeptap_python = (
            deeptap_python
            or os.environ.get("DEEPTAP_PYTHON")
            or sys.executable)
        self.max_peptide_length = max_peptide_length

    def __str__(self):
        return "DeepTAP(task_type=%r, deeptap_home=%r)" % (
            self.task_type, self.deeptap_home)

    def _default_pred_kind(self):
        return Kind.tap_transport

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if not peptide:
                raise ValueError("Empty peptide is not allowed")
            if len(peptide) > self.max_peptide_length:
                raise ValueError(
                    "DeepTAP supports peptides up to %d residues; got %r "
                    "(length %d)"
                    % (self.max_peptide_length, peptide, len(peptide)))
            invalid = set(peptide) - _VALID_AMINO_ACIDS
            if invalid:
                raise ValueError(
                    "Peptide %r contains characters DeepTAP cannot encode: %s"
                    % (peptide, "".join(sorted(invalid))))

    def predict(self, peptides):
        """Predict TAP transport for a list of peptides.

        Returns
        -------
        list of PeptideResult
            One entry per input peptide, each holding a single
            ``Kind.tap_transport`` prediction with an empty ``allele``. In
            ``"reg"`` mode the prediction's ``value`` is the predicted affinity
            in nM.
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        if not peptide_list:
            return []

        temp_dir = mkdtemp(prefix="mhctools", suffix="deeptap")
        # DeepTAP derives the output basename from the input filename up to the
        # first ".", so keep the stem dot-free.
        input_file_path = join(temp_dir, "deeptap_input.csv")
        output_file_path = join(
            temp_dir, "deeptap_input_DeepTAP_%s_predresult.csv" % self.task_type)
        rank_file_path = join(
            temp_dir,
            "deeptap_input_DeepTAP_%s_predresult_rank.csv" % self.task_type)
        pd.DataFrame({"peptide": peptide_list}).to_csv(
            input_file_path, index=False)

        args = [
            self.deeptap_python,
            join(self.deeptap_home, "deeptap.py"),
            "-t", self.task_type,
            "-f", input_file_path,
            "-o", temp_dir,
        ]

        with CleanupFiles(
                filenames=[input_file_path, output_file_path, rank_file_path],
                directories=[temp_dir]):
            run_command(args, suppress_stderr=False)
            if not exists(output_file_path):
                raise ValueError(
                    "DeepTAP produced no output file %r" % output_file_path)
            preds_by_peptide = parse_deeptap_results(
                output_file_path, self.task_type)

        return [
            PeptideResult(preds=tuple(preds_by_peptide.get(peptide, ())))
            for peptide in peptide_list
        ]


def parse_deeptap_results(filename, task_type="cla"):
    """Parse a DeepTAP prediction CSV into ``{peptide: [Prediction]}``.

    DeepTAP's ``cla`` output has columns ``peptide,pred_score,pred_label`` and
    its ``reg`` output ``peptide,pred_score,pred_affinity,pred_label``. In both,
    ``pred_score`` is in 0-1 with higher = stronger TAP binding; in ``reg`` mode
    the predicted affinity (nM, lower = stronger) is also surfaced as ``value``.

    Parameters
    ----------
    filename : str
    task_type : str
        ``"cla"`` or ``"reg"`` — selects whether ``pred_affinity`` is expected.

    Returns
    -------
    dict mapping peptide -> list with a single Prediction
    """
    df = pd.read_csv(filename)
    columns = list(df.columns)
    required = ["peptide", "pred_score"]
    if task_type == "reg":
        required.append("pred_affinity")
    missing = [c for c in required if c not in columns]
    if missing:
        raise ValueError(
            "DeepTAP output missing column(s) %s; got %s"
            % (missing, columns))

    predictor_name = "deeptap_%s" % task_type
    results = {}
    for row in df.itertuples(index=False):
        peptide = str(getattr(row, "peptide")).strip().upper()
        value = None
        if task_type == "reg":
            value = float(getattr(row, "pred_affinity"))
        results[peptide] = [Prediction(
            kind=Kind.tap_transport,
            score=float(getattr(row, "pred_score")),
            peptide=peptide,
            value=value,
            predictor_name=predictor_name)]
    return results
