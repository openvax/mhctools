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

"""Wrapper for TLimmuno2 — MHC class-II (CD4) immunogenicity.

TLimmuno2 predicts whether a class-II peptide is immunogenic to CD4+ T cells,
from the peptide and its class-II allele. It is the first **class-II**
immunogenicity predictor in mhctools — the other three (Calis, PRIME,
BigMHC_IM) are all class-I / CD8 — so it fills a real gap in the taxonomy,
emitting one ``Kind.immunogenicity`` prediction per (peptide, allele).

TLimmuno2 ships its trained weights in-repo but loads them with a Keras 2 /
TensorFlow stack, and its upstream license is ambiguous (an Apache-2.0 README
badge with no LICENSE file), so mhctools does **not** vendor it: this wrapper
shells out to TLimmuno2's own ``Python/TLimmuno2.py`` CLI in a user-provided
checkout (``TLIMMUNO2_HOME``), run by a user-provided interpreter
(``TLIMMUNO2_PYTHON``, default the current one). On newer TensorFlow the
interpreter only needs the ``tf-keras`` shim — this wrapper sets
``TF_USE_LEGACY_KERAS=1`` for the subprocess so TLimmuno2's Keras-2 SavedModels
load.

TLimmuno2's %Rank is computed by scoring the query against ~90,000 background
peptides *per distinct allele*, so a call takes on the order of a minute per
allele regardless of how many peptides are queried — batch peptides per allele.

Upstream: https://github.com/XSLiuLab/TLimmuno2
Cite: Wang et al., Briefings in Bioinformatics 2023 — "TLimmuno2: predicting
MHC class II antigen immunogenicity through transfer learning".

Like every current immunogenicity predictor, TLimmuno2 ranks better than it
generalizes to novel epitopes (and class-II immunogenicity is noisier still) —
a prioritization aid, not ground truth.
"""

import logging
import os
import re
import sys
from os.path import exists, isdir, isfile, join
from tempfile import mkdtemp

import pandas as pd

from .cleanup_context import CleanupFiles
from .pred import Kind, PeptideResult, Prediction
from .process_helpers import run_command
from .wrapper_base import NewModelPredictorMixin

logger = logging.getLogger(__name__)

# TLimmuno2 blosum-encodes each peptide into a fixed (21, 21) tensor, so a
# peptide longer than 21 residues overflows the input — a hard model limit.
TLIMMUNO2_MAX_PEPTIDE_LENGTH = 21

_PSEUDOSEQ_RELPATH = join("Python", "data", "pseudosequence.2016.all.X.dat")


def _find_tlimmuno2_home(tlimmuno2_home=None):
    """Resolve the TLimmuno2 checkout directory (holds ``Python/TLimmuno2.py``).

    Checks, in order: the *tlimmuno2_home* argument, ``$TLIMMUNO2_HOME``, then
    ``~/TLimmuno2``.
    """
    candidate = tlimmuno2_home or os.environ.get("TLIMMUNO2_HOME")
    if not candidate:
        home = join(os.path.expanduser("~"), "TLimmuno2")
        if isdir(home):
            candidate = home
    if not candidate:
        raise FileNotFoundError(
            "TLimmuno2 not found. Set TLIMMUNO2_HOME or pass tlimmuno2_home= to "
            "the constructor. Clone from https://github.com/XSLiuLab/TLimmuno2")
    if not isfile(join(candidate, "Python", "TLimmuno2.py")):
        raise FileNotFoundError(
            "Python/TLimmuno2.py not found in %r — is this a TLimmuno2 checkout?"
            % candidate)
    return candidate


def _load_valid_alleles(tlimmuno2_home):
    """The set of allele keys TLimmuno2 knows (from its pseudosequence file)."""
    path = join(tlimmuno2_home, _PSEUDOSEQ_RELPATH)
    df = pd.read_table(path, header=None, names=("HLA", "sequence"))
    return set(df["HLA"].astype(str))


def _allele_candidates(allele):
    """Yield candidate TLimmuno2 keys for an allele, most-specific first.

    TLimmuno2 keys are NetMHCIIpan-style and inconsistent across loci: DR is a
    single chain (``DRB1_0803``) while DP/DQ are heterodimers
    (``HLA-DPA10103-DPB10101``). We try the string as-is, then a few
    punctuation transforms, then pull the DRB chain out of a DRA/DRB
    heterodimer name.
    """
    a = allele.strip()
    yield a
    core = a[4:] if a.startswith("HLA-") else a
    yield core
    yield core.replace("*", "_").replace(":", "")   # DRB1*08:03 -> DRB1_0803
    yield a.replace("*", "").replace(":", "")        # DPA1*..-DPB1*.. -> DPA10..-DPB10..
    # mhctools canonicalizes DR to an alpha/beta pair (HLA-DRA1*01:01-DRB1*08:03);
    # TLimmuno2 keys it by the beta chain only.
    match = re.search(r"DRB[1-5]\*?\d+:?\d+", a)
    if match:
        yield match.group(0).replace("*", "_").replace(":", "")


def _tlimmuno2_allele(allele, valid_alleles):
    """Map an allele onto a TLimmuno2 key, or raise if it is unsupported.

    TLimmuno2's file mode inner-joins on the allele and silently drops unknown
    ones, so we validate against its own allele list and fail loudly instead.
    """
    for candidate in _allele_candidates(allele):
        if candidate in valid_alleles:
            return candidate
    raise ValueError(
        "TLimmuno2 does not recognize allele %r. Pass a NetMHCIIpan-style name "
        "present in Python/data/pseudosequence.2016.all.X.dat, e.g. "
        "'DRB1_0803' or 'HLA-DPA10103-DPB10101'." % allele)


def _subprocess_env():
    """Environment for the TLimmuno2 subprocess (legacy-Keras + quiet TF)."""
    env = dict(os.environ)
    env.setdefault("TF_USE_LEGACY_KERAS", "1")
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    return env


class TLimmuno2(NewModelPredictorMixin):
    """Wrapper for the TLimmuno2 class-II immunogenicity predictor.

    Parameters
    ----------
    alleles : list of str
        Class-II alleles. Native NetMHCIIpan-style keys (``DRB1_0803``,
        ``HLA-DPA10103-DPB10101``, ``H-2-IAb``) pass through; common DR forms
        (``HLA-DRB1*08:03``, ``DRB1*08:03``) are converted. Anything TLimmuno2
        does not know raises a ``ValueError``.
    tlimmuno2_home : str, optional
        Path to a TLimmuno2 checkout. Resolved from the argument, then
        ``$TLIMMUNO2_HOME``, then ``~/TLimmuno2``.
    tlimmuno2_python : str, optional
        Interpreter that can run TLimmuno2 (TensorFlow with Keras 2, or newer
        TensorFlow plus the ``tf-keras`` shim). Resolved from the argument,
        then ``$TLIMMUNO2_PYTHON``, then the current interpreter
        (``sys.executable``).
    """

    def __init__(self, alleles, tlimmuno2_home=None, tlimmuno2_python=None):
        if isinstance(alleles, str):
            alleles = [alleles]
        if not alleles:
            raise ValueError("TLimmuno2 requires at least one allele")
        self.alleles = list(alleles)
        self.tlimmuno2_home = _find_tlimmuno2_home(tlimmuno2_home)
        self.tlimmuno2_python = (
            tlimmuno2_python
            or os.environ.get("TLIMMUNO2_PYTHON")
            or sys.executable)
        self._valid_alleles = None

    def __str__(self):
        return "TLimmuno2(alleles=%s, tlimmuno2_home=%r)" % (
            self.alleles, self.tlimmuno2_home)

    def _default_pred_kind(self):
        return Kind.immunogenicity

    def kind_support(self):
        return {
            Kind.immunogenicity: {
                "mhc_dependence": "single_allele",
                "mhc_class": "II",
            },
        }

    def valid_alleles(self):
        """The set of allele keys this TLimmuno2 checkout supports (cached)."""
        if self._valid_alleles is None:
            self._valid_alleles = _load_valid_alleles(self.tlimmuno2_home)
        return self._valid_alleles

    def _check_peptides(self, peptides):
        for peptide in peptides:
            if not peptide:
                raise ValueError("Empty peptide is not allowed")
            if len(peptide) > TLIMMUNO2_MAX_PEPTIDE_LENGTH:
                raise ValueError(
                    "TLimmuno2 supports peptides up to %d residues; got %r "
                    "(length %d)"
                    % (TLIMMUNO2_MAX_PEPTIDE_LENGTH, peptide, len(peptide)))

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """Predict class-II immunogenicity for a list of peptides.

        Flanks are accepted for a uniform API but ignored. Returns one
        ``PeptideResult`` per input peptide; each holds one
        ``Kind.immunogenicity`` prediction per allele (``score`` in 0-1, higher
        = more immunogenic; ``percentile_rank`` is TLimmuno2's %Rank rescaled to
        0-100, lower = more immunogenic).
        """
        peptide_list = self._normalize_peptides(peptides)
        self._check_peptides(peptide_list)
        if not peptide_list:
            return []

        valid = self.valid_alleles()
        # Map each requested allele to its TLimmuno2 key (fail loud on unknowns).
        allele_keys = [_tlimmuno2_allele(a, valid) for a in self.alleles]

        # Deduplicate the (peptide, key) pairs we send: distinct requested
        # alleles can share a key, and TLimmuno2 re-groups output rows by
        # allele, so we map results back by content rather than by position.
        pairs = sorted({
            (peptide, key)
            for peptide in peptide_list
            for key in allele_keys
        })

        if len(set(allele_keys)) > 1:
            logger.info(
                "TLimmuno2 scores ~90k background peptides per allele for its "
                "%%Rank; %d distinct alleles will take a while.",
                len(set(allele_keys)))

        temp_dir = mkdtemp(prefix="mhctools", suffix="tlimmuno2")
        input_file_path = join(temp_dir, "tlimmuno2_input.csv")
        output_file_path = join(temp_dir, "result.csv")
        # File mode reads a header-less CSV of pep,HLA.
        pd.DataFrame(pairs).to_csv(input_file_path, index=False, header=False)

        args = [
            self.tlimmuno2_python,
            join(self.tlimmuno2_home, "Python", "TLimmuno2.py"),
            "--mode", "file",
            "--intdir", input_file_path,
            "--outdir", temp_dir,
            "--gpu", "False",
        ]

        with CleanupFiles(
                filenames=[input_file_path, output_file_path],
                directories=[temp_dir]):
            # TLimmuno2 hardcodes ./Python/... paths, so run in its own dir.
            run_command(
                args,
                suppress_stderr=False,
                cwd=self.tlimmuno2_home,
                env=_subprocess_env())
            if not exists(output_file_path):
                raise ValueError(
                    "TLimmuno2 produced no output file %r" % output_file_path)
            scores_by_pair = parse_tlimmuno2_results(output_file_path)

        results = []
        for peptide in peptide_list:
            preds = []
            for allele, key in zip(self.alleles, allele_keys):
                entry = scores_by_pair.get((peptide, key))
                if entry is None:
                    continue
                score, percentile_rank = entry
                preds.append(Prediction(
                    kind=Kind.immunogenicity,
                    score=score,
                    peptide=peptide,
                    allele=allele,
                    percentile_rank=percentile_rank,
                    predictor_name="tlimmuno2"))
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def parse_tlimmuno2_results(filename):
    """Parse a TLimmuno2 ``result.csv`` into ``{(peptide, allele): (score, rank)}``.

    The file has columns ``pep``, ``HLA``, ``sequence``, ``prediction``,
    ``Rank`` (plus an unnamed index). ``prediction`` is the immunogenicity score
    (0-1, higher = more immunogenic); ``Rank`` is a 0-1 fraction (lower = more
    immunogenic) that we rescale to a 0-100 ``percentile_rank`` for consistency
    with the other predictors.

    Returns
    -------
    dict mapping (peptide, allele_key) -> (score, percentile_rank)
    """
    df = pd.read_csv(filename)
    for column in ("pep", "HLA", "prediction", "Rank"):
        if column not in df.columns:
            raise ValueError(
                "TLimmuno2 output missing %r column; got %s"
                % (column, list(df.columns)))
    results = {}
    for row in df.itertuples(index=False):
        peptide = str(getattr(row, "pep")).strip().upper()
        allele = str(getattr(row, "HLA"))
        results[(peptide, allele)] = (
            float(getattr(row, "prediction")),
            float(getattr(row, "Rank")) * 100.0)
    return results
