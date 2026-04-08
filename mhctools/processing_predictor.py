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

"""
Base class for antigen-processing predictors and built-in scoring functions.

Scoring functions
-----------------
A scoring function combines three components extracted from a per-position
cleavage profile into a single peptide-level score::

    score = scoring(c_term, n_term, internal_probs)

*c_term* (float)
    Cleavage probability after the peptide's last residue.
*n_term* (float or None)
    Cleavage probability after the residue immediately upstream of the
    peptide. ``None`` when the peptide starts at position 0 (no upstream
    residue).
*internal_probs* (list of float)
    Cleavage probabilities at positions *within* the peptide (excluding
    the C-terminal position). Cleavage at any of these would destroy the
    peptide.
"""

from collections import defaultdict

from .pred import Pred, Kind, PeptidePreds


# ------------------------------------------------------------------
# Built-in scoring functions
# ------------------------------------------------------------------

def _geomean(*values):
    product = 1.0
    for v in values:
        product *= v
    return product ** (1.0 / len(values))


def score_cterm(c_term, n_term, internal):
    """C-terminal cleavage probability only."""
    return c_term


def score_nterm_cterm(c_term, n_term, internal):
    """Geometric mean of N- and C-terminal cleavage.

    Falls back to C-terminal only when *n_term* is ``None``.
    """
    if n_term is None:
        return c_term
    return _geomean(c_term, n_term)


def score_cterm_anti_max_internal(c_term, n_term, internal):
    """``c_term * (1 - max(internal))``."""
    max_i = max(internal) if internal else 0.0
    return c_term * (1.0 - max_i)


def score_cterm_anti_mean_internal(c_term, n_term, internal):
    """``c_term * (1 - mean(internal))``."""
    mean_i = sum(internal) / len(internal) if internal else 0.0
    return c_term * (1.0 - mean_i)


def score_nterm_cterm_anti_max_internal(c_term, n_term, internal):
    """Geometric mean of C-term, N-term, and ``1 - max(internal)``.

    N-term is omitted from the mean when ``None``.
    """
    max_i = max(internal) if internal else 0.0
    components = [c_term, 1.0 - max_i]
    if n_term is not None:
        components.append(n_term)
    return _geomean(*components)


def score_nterm_cterm_anti_mean_internal(c_term, n_term, internal):
    """Geometric mean of C-term, N-term, and ``1 - mean(internal)``.

    N-term is omitted from the mean when ``None``.
    """
    mean_i = sum(internal) / len(internal) if internal else 0.0
    components = [c_term, 1.0 - mean_i]
    if n_term is not None:
        components.append(n_term)
    return _geomean(*components)


SCORING_MODES = {
    "cterm": score_cterm,
    "nterm_cterm": score_nterm_cterm,
    "cterm_max_internal": score_cterm_anti_max_internal,
    "cterm_mean_internal": score_cterm_anti_mean_internal,
    "nterm_cterm_max_internal": score_nterm_cterm_anti_max_internal,
    "nterm_cterm_mean_internal": score_nterm_cterm_anti_mean_internal,
}


def resolve_scoring(scoring):
    """Resolve a scoring argument to a callable.

    Accepts a string mode name (see :data:`SCORING_MODES`), a callable,
    or ``None`` (returns ``None`` so the caller can apply its default).
    """
    if scoring is None or callable(scoring):
        return scoring
    if isinstance(scoring, str):
        if scoring not in SCORING_MODES:
            raise ValueError(
                "Unknown scoring mode %r, choose from %s" % (
                    scoring, list(SCORING_MODES)))
        return SCORING_MODES[scoring]
    raise TypeError(
        "scoring must be a string mode, callable, or None — got %r" % type(scoring))


# ------------------------------------------------------------------
# ProcessingPredictor
# ------------------------------------------------------------------

class ProcessingPredictor:
    """
    Base class for antigen-processing predictors.

    Subclasses implement :meth:`cleavage_probs` (per-position cleavage
    probabilities for a sequence).  This class provides:

    * **Component helpers** that extract C-terminal, N-terminal, and
      internal cleavage components from a probability vector.
    * **Configurable scoring** via a callable that combines those
      components into a single peptide-level score.
    * **Flanking-sequence support** so cleavage models see context around
      the peptide (optional for ``predict``, automatic for
      ``predict_proteins``).

    Unlike :class:`BasePredictor`, this class has **no allele parameter**
    because antigen processing is allele-independent.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : str or callable, optional
        Either a mode name (``"cterm"``, ``"nterm_cterm"``,
        ``"cterm_max_internal"``, ``"cterm_mean_internal"``,
        ``"nterm_cterm_max_internal"``, ``"nterm_cterm_mean_internal"``)
        or a callable ``(c_term, n_term, internal_probs) -> float``.
        Default: ``"nterm_cterm_max_internal"``.
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None):
        if default_peptide_lengths is None:
            default_peptide_lengths = [9]
        if isinstance(default_peptide_lengths, int):
            default_peptide_lengths = [default_peptide_lengths]
        scoring = resolve_scoring(scoring)
        if scoring is None:
            scoring = score_nterm_cterm_anti_max_internal
        self.default_peptide_lengths = default_peptide_lengths
        self.scoring = scoring

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(scoring=%s)" % (
            self.__class__.__name__, getattr(self.scoring, "__name__", repr(self.scoring)))

    # ------------------------------------------------------------------
    # Abstract
    # ------------------------------------------------------------------

    def cleavage_probs(self, sequence):
        """
        Return per-position cleavage probabilities for *sequence*.

        Position *i* represents the probability of cleavage **after**
        residue *i* (0-based).

        Must be implemented by subclasses.
        """
        raise NotImplementedError(
            "%s must implement cleavage_probs" % self.__class__.__name__)

    # ------------------------------------------------------------------
    # Component helpers
    # ------------------------------------------------------------------

    @staticmethod
    def c_term_prob(probs, offset, length):
        """Cleavage probability after the peptide's last residue."""
        return probs[offset + length - 1]

    @staticmethod
    def n_term_prob(probs, offset, length):
        """Cleavage probability upstream of the peptide, or ``None``."""
        return probs[offset - 1] if offset > 0 else None

    @staticmethod
    def internal_probs(probs, offset, length):
        """Cleavage probabilities at positions within the peptide
        (excluding C-terminal)."""
        return probs[offset:offset + length - 1]

    @staticmethod
    def max_internal_prob(probs, offset, length):
        """Maximum internal cleavage probability."""
        internal = probs[offset:offset + length - 1]
        return max(internal) if internal else 0.0

    @staticmethod
    def mean_internal_prob(probs, offset, length):
        """Mean internal cleavage probability."""
        internal = probs[offset:offset + length - 1]
        return sum(internal) / len(internal) if internal else 0.0

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def _peptide_score(self, probs, offset, length):
        """Extract components and delegate to ``self.scoring``."""
        c = self.c_term_prob(probs, offset, length)
        n = self.n_term_prob(probs, offset, length)
        internal = self.internal_probs(probs, offset, length)
        return self.scoring(c, n, internal)

    def _pred_kind(self):
        """Kind value for Pred objects.  Override in subclasses."""
        return Kind.antigen_processing

    def _predictor_name(self):
        return self.__class__.__name__.lower()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """
        Predict processing scores for peptides.

        Parameters
        ----------
        peptides : list of str

        n_flanks : list of str, optional
            N-terminal flanking sequences, one per peptide.  When
            provided the cleavage model sees ``n_flank + peptide +
            c_flank`` so that edge positions get proper context.

        c_flanks : list of str, optional
            C-terminal flanking sequences, one per peptide.

        Returns
        -------
        list of PeptidePreds
        """
        results = []
        for i, peptide in enumerate(peptides):
            n_flank = n_flanks[i] if n_flanks else ""
            c_flank = c_flanks[i] if c_flanks else ""

            full_seq = n_flank + peptide + c_flank
            probs = self.cleavage_probs(full_seq)

            offset = len(n_flank)
            score = self._peptide_score(probs, offset, len(peptide))
            pred = Pred(
                kind=self._pred_kind(),
                score=score,
                peptide=peptide,
                n_flank=n_flank,
                c_flank=c_flank,
                predictor_name=self._predictor_name(),
            )
            results.append(PeptidePreds(preds=(pred,)))
        return results

    def predict_dataframe(self, peptides, n_flanks=None, c_flanks=None,
                          sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS
        dfs = [pp.to_dataframe(sample_name)
               for pp in self.predict(peptides, n_flanks, c_flanks)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict_proteins(self, sequence_dict, peptide_lengths=None,
                         flank_length=0):
        """
        Run cleavage prediction on full protein sequences and aggregate
        into peptide-level processing scores.

        This is the preferred entry-point because the full protein gives
        the cleavage model proper flanking context at every position.

        Parameters
        ----------
        sequence_dict : dict or str

        peptide_lengths : list of int, optional

        flank_length : int
            How many residues of flanking context to record in each
            :class:`Pred`'s ``n_flank`` / ``c_flank`` fields (default 0).
            This does **not** affect the cleavage model — it always sees
            the full protein.

        Returns
        -------
        dict mapping sequence_name -> list of PeptidePreds
        """
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}
        elif isinstance(sequence_dict, (list, tuple)):
            sequence_dict = {seq: seq for seq in sequence_dict}

        peptide_lengths = self._resolve_peptide_lengths(peptide_lengths)

        results = defaultdict(list)
        for name, sequence in sequence_dict.items():
            probs = self.cleavage_probs(sequence)
            for plen in peptide_lengths:
                for i in range(len(sequence) - plen + 1):
                    peptide = sequence[i:i + plen]
                    score = self._peptide_score(probs, offset=i, length=plen)
                    if flank_length:
                        n_flank = sequence[max(0, i - flank_length):i]
                        c_flank = sequence[i + plen:i + plen + flank_length]
                    else:
                        n_flank = ""
                        c_flank = ""
                    pred = Pred(
                        kind=self._pred_kind(),
                        score=score,
                        peptide=peptide,
                        n_flank=n_flank,
                        c_flank=c_flank,
                        source_sequence_name=name,
                        offset=i,
                        predictor_name=self._predictor_name(),
                    )
                    results[name].append(PeptidePreds(preds=(pred,)))
        return dict(results)

    def predict_proteins_dataframe(
            self, sequence_dict, peptide_lengths=None,
            flank_length=0, sample_name=""):
        """``predict_proteins()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS
        dfs = []
        for name, pp_list in self.predict_proteins(
                sequence_dict, peptide_lengths, flank_length).items():
            for pp in pp_list:
                dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict_cleavage_sites(self, sequence_dict):
        """
        Return raw per-position cleavage probabilities.

        Returns
        -------
        dict mapping sequence_name -> list of float
        """
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}
        return {
            name: self.cleavage_probs(seq)
            for name, seq in sequence_dict.items()
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _resolve_peptide_lengths(self, peptide_lengths):
        if peptide_lengths is None:
            peptide_lengths = self.default_peptide_lengths
        if isinstance(peptide_lengths, int):
            peptide_lengths = [peptide_lengths]
        if not peptide_lengths:
            raise ValueError("peptide_lengths must be a non-empty list of ints")
        return peptide_lengths
