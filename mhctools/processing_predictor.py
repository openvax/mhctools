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
Base class for antigen-processing predictors that produce per-position
cleavage probabilities and aggregate them into peptide-level scores.
"""

from collections import defaultdict

from .pred import Pred, Kind, PeptidePreds


class ProcessingPredictor:
    """
    Base class for antigen-processing predictors.

    Subclasses must implement :meth:`cleavage_probs`, which returns a
    list of per-position cleavage probabilities for a protein sequence.
    This class provides the aggregation logic that converts those
    per-position probabilities into peptide-level processing scores.

    Unlike :class:`BasePredictor`, this class has **no allele parameter**
    because antigen processing is allele-independent.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : str
        How to aggregate per-position cleavage probabilities into a
        single peptide processing score (all use geometric mean):

        * ``"cterm"`` -- C-terminal cleavage probability only.
        * ``"nterm_cterm"`` -- N- and C-terminal cleavage.
        * ``"cterm_max_internal"`` -- C-terminal and (1 - max internal).
        * ``"cterm_mean_internal"`` -- C-terminal and (1 - mean internal).
        * ``"nterm_cterm_max_internal"`` -- N-term, C-term, and
          (1 - max internal).
        * ``"nterm_cterm_mean_internal"`` -- N-term, C-term, and
          (1 - mean internal).

        When the N-terminal cleavage site falls outside the sequence
        (``offset == 0`` or isolated peptides), the N-term component is
        omitted from the geometric mean.
    """

    SCORING_METHODS = (
        "cterm",
        "nterm_cterm",
        "cterm_max_internal",
        "cterm_mean_internal",
        "nterm_cterm_max_internal",
        "nterm_cterm_mean_internal",
    )

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring="nterm_cterm_max_internal"):
        if default_peptide_lengths is None:
            default_peptide_lengths = [9]
        if isinstance(default_peptide_lengths, int):
            default_peptide_lengths = [default_peptide_lengths]
        if scoring not in self.SCORING_METHODS:
            raise ValueError(
                "scoring must be one of %s, got %r" % (
                    self.SCORING_METHODS, scoring))
        self.default_peptide_lengths = default_peptide_lengths
        self.scoring = scoring

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(scoring=%r)" % (self.__class__.__name__, self.scoring)

    # ------------------------------------------------------------------
    # Abstract
    # ------------------------------------------------------------------

    def cleavage_probs(self, sequence):
        """
        Return per-position cleavage probabilities for *sequence*.

        Position *i* represents the probability of cleavage **after**
        residue *i* (0-based).

        Parameters
        ----------
        sequence : str
            Amino acid sequence.

        Returns
        -------
        list of float
            One probability per residue, in ``[0, 1]``.
        """
        raise NotImplementedError(
            "%s must implement cleavage_probs" % self.__class__.__name__)

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def _peptide_score(self, probs, offset, length):
        """
        Aggregate per-position cleavage probabilities into a single
        peptide processing score.

        Components (depending on ``self.scoring``):

        * **C-terminal** -- ``probs[offset + length - 1]``: cleavage
          after the last residue of the peptide.
        * **N-terminal** -- ``probs[offset - 1]`` (when ``offset > 0``):
          cleavage after the residue immediately upstream.
        * **Internal** -- ``probs[offset : offset + length - 1]``:
          cleavage within the peptide (destroying it).  Enters the
          formula as ``1 - summary(internal)``.
        """
        c_term = probs[offset + length - 1]

        n_term = probs[offset - 1] if offset > 0 else None

        internal = probs[offset:offset + length - 1]

        components = [c_term]

        if "nterm" in self.scoring and n_term is not None:
            components.append(n_term)

        if "internal" in self.scoring and internal:
            if "max" in self.scoring:
                summary = max(internal)
            else:
                summary = sum(internal) / len(internal)
            components.append(1.0 - summary)

        product = 1.0
        for c in components:
            product *= c
        return product ** (1.0 / len(components))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def predict(self, peptides):
        """
        Predict antigen-processing scores for isolated peptides.

        .. note::

           For short peptides the C-terminal position may lack downstream
           flanking context, producing unreliable scores.  Prefer
           :meth:`predict_proteins` when the source protein is available.

        Parameters
        ----------
        peptides : list of str

        Returns
        -------
        list of PeptidePreds
        """
        results = []
        for peptide in peptides:
            probs = self.cleavage_probs(peptide)
            score = self._peptide_score(probs, offset=0, length=len(peptide))
            pred = Pred(
                kind=Kind.antigen_processing,
                score=score,
                peptide=peptide,
                predictor_name=self._predictor_name(),
            )
            results.append(PeptidePreds(preds=(pred,)))
        return results

    def predict_dataframe(self, peptides, sample_name=""):
        """``predict()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS
        dfs = [pp.to_dataframe(sample_name) for pp in self.predict(peptides)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict_proteins(self, sequence_dict, peptide_lengths=None):
        """
        Run cleavage prediction on full protein sequences and aggregate
        into peptide-level processing scores.

        This is the preferred entry-point because flanking-sequence
        context produces more accurate cleavage probabilities.

        Parameters
        ----------
        sequence_dict : dict or str
            Mapping of sequence names to amino acid strings, or a single
            sequence string.

        peptide_lengths : list of int, optional

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
                    pred = Pred(
                        kind=Kind.antigen_processing,
                        score=score,
                        peptide=peptide,
                        source_sequence_name=name,
                        offset=i,
                        predictor_name=self._predictor_name(),
                    )
                    results[name].append(PeptidePreds(preds=(pred,)))
        return dict(results)

    def predict_proteins_dataframe(
            self, sequence_dict, peptide_lengths=None, sample_name=""):
        """``predict_proteins()`` flattened to a DataFrame."""
        import pandas as pd
        from .pred import COLUMNS
        dfs = []
        for name, pp_list in self.predict_proteins(
                sequence_dict, peptide_lengths).items():
            for pp in pp_list:
                dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def predict_cleavage_sites(self, sequence_dict):
        """
        Return raw per-position cleavage probabilities.

        Parameters
        ----------
        sequence_dict : dict or str

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

    def _predictor_name(self):
        return self.__class__.__name__.lower()

    def _resolve_peptide_lengths(self, peptide_lengths):
        if peptide_lengths is None:
            peptide_lengths = self.default_peptide_lengths
        if isinstance(peptide_lengths, int):
            peptide_lengths = [peptide_lengths]
        if not peptide_lengths:
            raise ValueError("peptide_lengths must be a non-empty list of ints")
        return peptide_lengths
