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

"""Wrapper for DeepImmuno — a CNN predictor of class-I CD8 immunogenicity.

DeepImmuno scores whether a presented class-I peptide is immunogenic (elicits a
CD8+ T-cell response), from the peptide sequence and its HLA-A/B/C allele. It
joins the other immunogenicity predictors mhctools wraps (Calis, PRIME,
BigMHC_IM), emitting ``Kind.immunogenicity`` per (peptide, allele).

DeepImmuno ships its trained CNN weights in-repo (MIT-licensed) but its
``deepimmuno-cnn.py`` script rebuilds the network in code and loads those
weights with an old Keras 2 / TensorFlow stack. To keep that dependency out of
the mhctools environment, this wrapper shells out to that script in a
user-provided checkout (``DEEPIMMUNO_HOME``), run by a user-provided interpreter
(``DEEPIMMUNO_PYTHON``, default the current one). On newer TensorFlow the
interpreter only needs the ``tf-keras`` shim installed — this wrapper sets
``TF_USE_LEGACY_KERAS=1`` for the subprocess so DeepImmuno's Keras-2 model
loads.

Upstream: https://github.com/frankligy/DeepImmuno
Cite: Li et al., Briefings in Bioinformatics 2021 — "DeepImmuno: deep
learning-empowered prediction and generation of immunogenic peptides for T-cell
immunity".

Like every current CD8 immunogenicity predictor, DeepImmuno ranks better than
it generalizes to novel neoepitopes (independent benchmarks put the field near
AUC 0.5-0.65) — a prioritization aid, not ground truth.
"""

import os
import sys
from os.path import exists, isdir, isfile, join
from tempfile import mkdtemp

import pandas as pd

from .allele_normalization import normalize_allele_name
from .cleanup_context import CleanupFiles
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command
from .wrapper_base import NewModelPredictorMixin

# DeepImmuno's peptide encoder only handles 9- and 10-mers: a 9mer is padded to
# 10 with a gap inserted after position 5, and any other length leaves the
# encoding undefined. So 9 and 10 are the only lengths it can score.
DEEPIMMUNO_PEPTIDE_LENGTHS = (9, 10)


def _find_deepimmuno_home(deepimmuno_home=None):
    """Resolve the DeepImmuno checkout directory (holds ``deepimmuno-cnn.py``).

    Checks, in order: the *deepimmuno_home* argument, ``$DEEPIMMUNO_HOME``, then
    ``~/DeepImmuno``.
    """
    candidate = deepimmuno_home or os.environ.get("DEEPIMMUNO_HOME")
    if not candidate:
        home = join(os.path.expanduser("~"), "DeepImmuno")
        if isdir(home):
            candidate = home
    if not candidate:
        raise FileNotFoundError(
            "DeepImmuno not found. Set DEEPIMMUNO_HOME or pass deepimmuno_home= "
            "to the constructor. Clone from "
            "https://github.com/frankligy/DeepImmuno")
    if not isfile(join(candidate, "deepimmuno-cnn.py")):
        raise FileNotFoundError(
            "deepimmuno-cnn.py not found in %r — is this a DeepImmuno checkout?"
            % candidate)
    return candidate


def _deepimmuno_allele(allele):
    """Format an allele the way DeepImmuno keys its paratope table (HLA-A*0201).

    mhctools canonicalizes to ``HLA-A*02:01``; DeepImmuno's keys drop the colon.
    DeepImmuno itself snaps an unknown allele to the nearest one it knows
    (``rescue_unknown_hla``), so only the punctuation has to match.
    """
    return normalize_allele_name(allele).replace(":", "")


def _subprocess_env():
    """Environment for the DeepImmuno subprocess.

    ``TF_USE_LEGACY_KERAS=1`` routes ``tensorflow.keras`` to the ``tf-keras``
    (Keras 2) shim on TensorFlow >= 2.16 so the committed checkpoint loads; it
    is ignored by older TensorFlow. ``TF_CPP_MIN_LOG_LEVEL=3`` quiets TF's C++
    logging. Both only take effect if not already set by the caller.
    """
    env = dict(os.environ)
    env.setdefault("TF_USE_LEGACY_KERAS", "1")
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    return env


class DeepImmuno(NewModelPredictorMixin):
    """Wrapper for the DeepImmuno-CNN immunogenicity predictor.

    Parameters
    ----------
    alleles : list of str
        Class-I HLA-A/B/C alleles (e.g. ``["HLA-A*02:01", "HLA-B*07:02"]``).
        DeepImmuno supports a fixed set of ~62 alleles and snaps anything else
        to the nearest one it knows.
    deepimmuno_home : str, optional
        Path to a DeepImmuno checkout. Resolved from the argument, then
        ``$DEEPIMMUNO_HOME``, then ``~/DeepImmuno``.
    deepimmuno_python : str, optional
        Interpreter that can run DeepImmuno (TensorFlow with Keras 2, or newer
        TensorFlow plus the ``tf-keras`` shim). Resolved from the argument, then
        ``$DEEPIMMUNO_PYTHON``, then the current interpreter
        (``sys.executable``).
    """

    def __init__(self, alleles, deepimmuno_home=None, deepimmuno_python=None):
        if isinstance(alleles, str):
            alleles = [alleles]
        if not alleles:
            raise ValueError("DeepImmuno requires at least one allele")
        self.alleles = [normalize_allele_name(a) for a in alleles]
        self.deepimmuno_home = _find_deepimmuno_home(deepimmuno_home)
        self.deepimmuno_python = (
            deepimmuno_python
            or os.environ.get("DEEPIMMUNO_PYTHON")
            or sys.executable)

    def __str__(self):
        return "DeepImmuno(alleles=%s, deepimmuno_home=%r)" % (
            self.alleles, self.deepimmuno_home)

    def _default_pred_kind(self):
        return Kind.immunogenicity

    def kind_support(self):
        return {
            Kind.immunogenicity: {
                "mhc_dependence": "single_allele",
                "mhc_class": "I",
            },
        }

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if len(peptide) not in DEEPIMMUNO_PEPTIDE_LENGTHS:
                raise ValueError(
                    "DeepImmuno only scores 9- and 10-mers; got %r (length %d)"
                    % (peptide, len(peptide)))

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict immunogenicity for a list of peptides.

        Flanks are accepted for a uniform API but ignored (DeepImmuno does not
        use flanking context).

        Returns
        -------
        list of PeptideResult
            One entry per input peptide; each holds one ``Kind.immunogenicity``
            prediction per allele (``score`` in 0-1, higher = more immunogenic).
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        if not peptide_list:
            return []

        # (peptide, allele) grid, alleles in DeepImmuno's HLA-A*0201 key format.
        rows = [
            (peptide, _deepimmuno_allele(allele))
            for peptide in peptide_list
            for allele in self.alleles
        ]

        temp_dir = mkdtemp(prefix="mhctools", suffix="deepimmuno")
        input_file_path = join(temp_dir, "deepimmuno_input.csv")
        # deepimmuno-cnn.py writes a fixed basename into --outdir.
        output_file_path = join(temp_dir, "deepimmuno-cnn-result.txt")
        # DeepImmuno reads a header-less CSV of peptide,HLA.
        pd.DataFrame(rows).to_csv(input_file_path, index=False, header=False)

        args = [
            self.deepimmuno_python,
            join(self.deepimmuno_home, "deepimmuno-cnn.py"),
            "--mode", "multiple",
            "--intdir", input_file_path,
            "--outdir", temp_dir,
        ]

        with CleanupFiles(
                filenames=[input_file_path, output_file_path],
                directories=[temp_dir]):
            # DeepImmuno hardcodes ./data and ./models, so run in its own dir.
            run_command(
                args,
                suppress_stderr=False,
                cwd=self.deepimmuno_home,
                env=_subprocess_env())
            if not exists(output_file_path):
                raise ValueError(
                    "DeepImmuno produced no output file %r" % output_file_path)
            scores = parse_deepimmuno_results(output_file_path)

        if len(scores) != len(rows):
            raise ValueError(
                "DeepImmuno returned %d rows for %d (peptide, allele) pairs"
                % (len(scores), len(rows)))

        # DeepImmuno preserves input row order, so consume the scores in the
        # same nested (peptide, allele) order we wrote them.
        idx = 0
        results = []
        for peptide in peptide_list:
            preds = []
            for allele in self.alleles:
                preds.append(Prediction(
                    kind=Kind.immunogenicity,
                    score=scores[idx],
                    peptide=peptide,
                    allele=allele,
                    predictor_name="deepimmuno"))
                idx += 1
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def parse_deepimmuno_results(filename):
    """Parse a DeepImmuno result file into an ordered list of float scores.

    ``deepimmuno-cnn.py`` writes a tab-separated file with columns
    ``peptide``, ``HLA``, ``immunogenicity`` (0-1, higher = more immunogenic),
    one row per input pair in input order.

    Returns
    -------
    list of float
        Immunogenicity scores, in file order.
    """
    df = pd.read_csv(filename, sep="\t")
    if "immunogenicity" not in df.columns:
        raise ValueError(
            "DeepImmuno output missing 'immunogenicity' column; got %s"
            % list(df.columns))
    return [float(v) for v in df["immunogenicity"]]
