# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import warnings
from collections import defaultdict, namedtuple

from typechecks import require_iterable_of
from .allele_normalization import (
    normalize_allele_name,
    normalize_allele_name_or_raw,
)

from .unsupported_allele import UnsupportedAllele
from .binding_prediction_collection import BindingPredictionCollection
from .pred import (
    Kind,
    PeptideResult,
    Prediction,
)

logger = logging.getLogger(__name__)

PeptideContext = namedtuple(
    "PeptideContext",
    ["source_sequence_name", "offset", "peptide", "n_flank", "c_flank"])


def _normalize_sequence_dict(sequence_dict):
    if isinstance(sequence_dict, str):
        return {"seq": sequence_dict}
    if isinstance(sequence_dict, (list, tuple)):
        return {seq: seq for seq in sequence_dict}
    return sequence_dict


def _check_flank_inputs(peptides, n_flanks=None, c_flanks=None):
    peptide_list = list(peptides)

    def _check_flanks(name, flanks):
        if flanks is None:
            return None
        if isinstance(flanks, str):
            raise TypeError("%s must be a sequence of strings, not a string" % name)
        flank_list = list(flanks)
        require_iterable_of(flank_list, str, name)
        if len(flank_list) != len(peptide_list):
            raise ValueError(
                "%s must have one entry per peptide, got %d flank(s) for "
                "%d peptide(s)" % (name, len(flank_list), len(peptide_list)))
        return flank_list

    return (
        peptide_list,
        _check_flanks("n_flanks", n_flanks),
        _check_flanks("c_flanks", c_flanks),
    )


def _peptide_contexts(
        sequence_dict,
        peptide_lengths,
        flank_length=0,
        n_flank_length=None,
        c_flank_length=None):
    if n_flank_length is None:
        n_flank_length = flank_length
    if c_flank_length is None:
        c_flank_length = flank_length

    contexts = []
    for name, sequence in sequence_dict.items():
        for peptide_length in peptide_lengths:
            for i in range(len(sequence) - peptide_length + 1):
                peptide = sequence[i:i + peptide_length]
                if n_flank_length or c_flank_length:
                    n_flank = sequence[max(0, i - n_flank_length):i]
                    c_flank = sequence[
                        i + peptide_length:
                        i + peptide_length + c_flank_length]
                else:
                    n_flank = ""
                    c_flank = ""
                contexts.append(PeptideContext(
                    source_sequence_name=name,
                    offset=i,
                    peptide=peptide,
                    n_flank=n_flank,
                    c_flank=c_flank,
                ))
    return contexts


class BasePredictor(object):
    """
    Base class for all MHC binding predictors.
    """
    uses_flanking_sequences = False
    flank_length = 15
    n_flank_length = None
    c_flank_length = None
    mhc_class = "I"

    def __init__(
            self,
            alleles,
            valid_alleles=None,
            default_peptide_lengths=None,
            min_peptide_length=8,
            max_peptide_length=None,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False,
            keep_unparseable_alleles=False):
        """
        Parameters
        ----------
        alleles : list
            List of strings containing names of HLA alleles we're
            making predictions for. Example:
                ["HLA-A*02:01", "HLA-B*07:02"]

        valid_alleles : list, optional
            If given, constrain HLA alleles to be contained within
            this set.

        default_peptide_lengths : list of int, optional
            When making predictions across subsequences of protein sequences,
            what peptide lengths to predict for.

        min_peptide_length : int
            Shortest peptide this predictor can handle

        max_peptide_length : int
            Longest peptide this predictor can handle

        allow_X_in_peptides : bool
            Allow unknown amino acids in peptide sequences

        allow_lowercase_in_peptides : bool
            Allow lowercase letters in peptide sequences

        keep_unparseable_alleles : bool
            If True, keep allele names our normalizer can't parse verbatim
            instead of raising (used by command-line predictors that validate
            against the tool's own supported-allele list, so that non-human
            alleles like H-2-Qa1 or BoLA-amani.1 can still be requested).
        """
        # I find myself often constructing a predictor with just one allele
        # so as a convenience, allow user to not wrap that allele as a list
        if type(alleles) is str:
            alleles = alleles.split(',')
        self.alleles = self._check_hla_alleles(
            alleles, valid_alleles, keep_unparseable=keep_unparseable_alleles)

        if type(default_peptide_lengths) is int:
            default_peptide_lengths = [default_peptide_lengths]
        require_iterable_of(default_peptide_lengths, int)
        self.default_peptide_lengths = default_peptide_lengths
        self.min_peptide_length = min_peptide_length
        self.max_peptide_length = max_peptide_length
        self.allow_X_in_peptides = allow_X_in_peptides
        self.allow_lowercase_in_peptides = allow_lowercase_in_peptides

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(alleles=%s, default_peptide_lengths=%s)" % (
            self.__class__.__name__,
            self.alleles,
            self.default_peptide_lengths)

    # --- new API ---

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        """
        Predict for a list of peptide sequences.

        n_flanks and c_flanks are accepted for a uniform flank-aware API.
        The default BindingPrediction-based implementation validates but
        ignores them; subclasses that use flanking context should override
        this method and set ``uses_flanking_sequences = True``.

        Returns
        -------
        list of PeptideResult
        """
        peptide_list, _, _ = _check_flank_inputs(
            peptides, n_flanks, c_flanks)
        collection = self.predict_peptides(peptide_list)
        return collection.to_peptide_preds(kind=self._default_pred_kind())

    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        """
        Optional flank-aware prediction path.

        The default implementation validates aligned flank lists and falls
        back to ``predict(peptides)`` so predictors that do not use flanking
        context retain their existing behavior. Subclasses whose models use
        flanks should override this method and forward the flanks to the
        underlying predictor.
        """
        peptide_list, _, _ = _check_flank_inputs(peptides, n_flanks, c_flanks)
        return self.predict(peptide_list)

    def predict_dataframe(
            self, peptides, sample_name="", n_flanks=None, c_flanks=None):
        """predict() flattened to a DataFrame."""
        import pandas as pd
        dfs = [
            pp.to_dataframe(sample_name)
            for pp in self.predict(
                peptides, n_flanks=n_flanks, c_flanks=c_flanks)
        ]
        if not dfs:
            from .pred import COLUMNS
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def _predict_protein_flank_length(self):
        return max(self._predict_protein_flank_lengths())

    def _predict_protein_flank_lengths(self):
        if not self.uses_flanking_sequences:
            return (0, 0)
        n_flank_length = (
            self.flank_length
            if self.n_flank_length is None else self.n_flank_length)
        c_flank_length = (
            self.flank_length
            if self.c_flank_length is None else self.c_flank_length)
        return (n_flank_length, c_flank_length)

    def _check_flank_inputs(self, peptides, n_flanks=None, c_flanks=None):
        return _check_flank_inputs(peptides, n_flanks, c_flanks)

    @staticmethod
    def _with_protein_location(pred, context):
        return Prediction(
            kind=pred.kind,
            score=pred.score,
            peptide=pred.peptide,
            allele=pred.allele,
            n_flank=context.n_flank,
            c_flank=context.c_flank,
            value=pred.value,
            percentile_rank=pred.percentile_rank,
            source_sequence_name=context.source_sequence_name,
            offset=context.offset,
            predictor_name=pred.predictor_name,
            predictor_version=pred.predictor_version,
        )

    def predict_proteins(self, sequence_dict, peptide_lengths=None):
        """
        Scan protein sequences and predict for all subsequences.

        Parameters
        ----------
        sequence_dict : dict or str
            Mapping of sequence names to amino acid strings.
            If a string, treated as {"seq": string}.

        peptide_lengths : list of int, optional

        Returns
        -------
        dict mapping sequence_name -> list of PeptideResult
        """
        sequence_dict = _normalize_sequence_dict(sequence_dict)

        peptide_lengths = self._check_peptide_lengths(peptide_lengths)

        n_flank_length, c_flank_length = self._predict_protein_flank_lengths()
        contexts = _peptide_contexts(
            sequence_dict,
            peptide_lengths,
            n_flank_length=n_flank_length,
            c_flank_length=c_flank_length)

        if self.uses_flanking_sequences:
            peptide_list = [context.peptide for context in contexts]
            n_flanks = [context.n_flank for context in contexts]
            c_flanks = [context.c_flank for context in contexts]
            flat_preds = self.predict_with_flanks(
                peptide_list,
                n_flanks=n_flanks,
                c_flanks=c_flanks)
            if len(flat_preds) != len(contexts):
                raise ValueError(
                    "%s.predict returned %d result(s) for %d flanked "
                    "peptide occurrence(s)" % (
                        self.__class__.__name__, len(flat_preds),
                        len(contexts)))
            results = defaultdict(list)
            for context, pp in zip(contexts, flat_preds):
                relocated = PeptideResult(preds=tuple(
                    self._with_protein_location(p, context)
                    for p in pp.preds
                ))
                results[context.source_sequence_name].append(relocated)
            return dict(results)

        peptide_set = {context.peptide for context in contexts}
        peptide_to_contexts = defaultdict(list)
        for context in contexts:
            peptide_to_contexts[context.peptide].append(context)
        peptide_list = sorted(peptide_set)
        flat_preds = self.predict(peptide_list)

        # Expand: each PeptideResult may map to multiple (name, offset) positions
        results = defaultdict(list)
        for pp in flat_preds:
            if not pp.preds:
                continue
            peptide = pp.preds[0].peptide
            for context in peptide_to_contexts.get(peptide, []):
                relocated = PeptideResult(preds=tuple(
                    self._with_protein_location(p, context)
                    for p in pp.preds
                ))
                results[context.source_sequence_name].append(relocated)
        return dict(results)

    def predict_proteins_dataframe(self, sequence_dict, peptide_lengths=None, sample_name=""):
        """predict_proteins() flattened to a DataFrame."""
        import pandas as pd
        dfs = []
        for name, pp_list in self.predict_proteins(sequence_dict, peptide_lengths).items():
            for pp in pp_list:
                dfs.append(pp.to_dataframe(sample_name))
        if not dfs:
            from .pred import COLUMNS
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

    def _default_pred_kind(self):
        """Override in subclasses to set the Kind for compat conversion."""
        return Kind.pMHC_affinity

    def kind_support(self):
        """Predictor-specific MHC context for supported prediction kinds."""
        return {
            self._default_pred_kind(): {
                "mhc_dependence": "single_allele",
                "mhc_class": self.mhc_class,
            }
        }

    @property
    def supported_kinds(self):
        """Prediction kind strings this predictor can emit."""
        return tuple(self.kind_support())

    # --- deprecated API (still works) ---

    def predict_peptides(self, peptides):
        """
        Deprecated: use predict() instead.

        Given a list of peptide sequences, returns a BindingPredictionCollection.
        """
        raise NotImplementedError(
            "%s must implement predict_peptides" % (self.__class__.__name__,))

    def predict_peptides_dataframe(self, peptides):
        """Deprecated: use predict_dataframe() instead.

        Emits the canonical prediction schema (see ``mhctools.pred.COLUMNS``)
        for parity with ``predict_proteins_dataframe`` / ``predict_dataframe``.
        Previously this returned the legacy BindingPrediction schema
        (``source_sequence_name, offset, peptide, allele, score, affinity,
        percentile_rank, prediction_method_name, length``). That schema lacked
        ``predictor_version``, ``kind``, ``value`` and used
        ``prediction_method_name`` instead of ``predictor_name`` — see #193.
        """
        warnings.warn(
            "predict_peptides_dataframe is deprecated; use predict_dataframe()",
            DeprecationWarning,
            stacklevel=2)
        return self.predict_dataframe(peptides)

    def _check_peptide_lengths(self, peptide_lengths=None):
        """
        If peptide lengths not specified, then try using the default
        lengths associated with this predictor object. If those aren't
        a valid non-empty sequence of integers, then raise an exception.
        Otherwise return the peptide lengths.
        """
        if not peptide_lengths:
            peptide_lengths = self.default_peptide_lengths

        if not peptide_lengths:
            raise ValueError(
                ("Must either provide 'peptide_lengths' argument "
                "or set 'default_peptide_lengths"))
        if isinstance(peptide_lengths, int):
            peptide_lengths = [peptide_lengths]
        require_iterable_of(peptide_lengths, int)
        for peptide_length in peptide_lengths:
            if (self.min_peptide_length is not None and
                    peptide_length < self.min_peptide_length):
                raise ValueError(
                    "Invalid peptide length %d, shorter than min %d" % (
                        peptide_length,
                        self.min_peptide_length))
            elif (self.max_peptide_length is not None and
                    peptide_length > self.max_peptide_length):
                raise ValueError(
                    "Invalid peptide length %d, longer than max %d" % (
                        peptide_length,
                        self.max_peptide_length))
        return peptide_lengths

    def _check_results(self, binding_predictions, peptides, alleles):
        # Avoid materializing the Cartesian product: scan observed pairs once,
        # checking membership against the allele/peptide sets, then compare
        # unique-pair count to the expected total.
        alleles_set = set(alleles)
        peptides_set = set(peptides)
        n_expected = len(alleles_set) * len(peptides_set)
        observed = set()
        for bp in binding_predictions:
            pair = (bp.allele, bp.peptide)
            if bp.allele not in alleles_set or bp.peptide not in peptides_set:
                raise ValueError(
                    "Unexpected binding prediction peptide='%s' allele='%s'" % (
                        bp.peptide, bp.allele))
            observed.add(pair)
        # Every observed pair is in alleles_set x peptides_set, so
        # len(observed) <= n_expected; a strict inequality means at least
        # one expected pair is missing.
        if len(observed) < n_expected:
            # Iterate the original lists for a deterministic example pair.
            for a in alleles:
                for p in peptides:
                    if (a, p) not in observed:
                        raise ValueError(
                            "Missing predictions, example peptide='%s' allele='%s'. "
                            "Expected %d unique pairs, got %d." % (
                                p, a, n_expected, len(observed)))

    def _check_peptide_inputs(self, peptides):
        """
        Check peptide sequences to make sure they are valid for this predictor.
        """
        require_iterable_of(peptides, str)
        check_X = not self.allow_X_in_peptides
        check_lower = not self.allow_lowercase_in_peptides
        check_min_length = self.min_peptide_length is not None
        min_length = self.min_peptide_length
        check_max_length = self.max_peptide_length is not None
        max_length = self.max_peptide_length
        for p in peptides:
            if not p.isalpha():
                raise ValueError("Invalid characters in peptide '%s'" % p)
            elif check_X and "X" in p:
                raise ValueError("Invalid character 'X' in peptide '%s'" % p)
            elif check_lower and not p.isupper():
                raise ValueError("Invalid lowercase letters in peptide '%s'" % p)
            elif check_min_length and len(p) < min_length:
                raise ValueError(
                    "Peptide '%s' too short (%d chars), must be at least %d" % (
                        p, len(p), min_length))
            elif check_max_length and len(p) > max_length:
                raise ValueError(
                    "Peptide '%s' too long (%d chars), must be at least %d" % (
                        p, len(p), max_length))

    def predict_subsequences(
            self,
            sequence_dict,
            peptide_lengths=None):
        """
        Deprecated: use predict_proteins() instead.

        Given a dictionary mapping sequence names to amino acid strings,
        and an optional list of peptide lengths, returns a
        BindingPredictionCollection.
        """
        sequence_dict = _normalize_sequence_dict(sequence_dict)

        peptide_lengths = self._check_peptide_lengths(peptide_lengths)

        contexts = _peptide_contexts(sequence_dict, peptide_lengths)
        peptide_set = {context.peptide for context in contexts}
        peptide_to_contexts = defaultdict(list)
        for context in contexts:
            peptide_to_contexts[context.peptide].append(context)
        peptide_list = sorted(peptide_set)

        binding_predictions = self.predict_peptides(peptide_list)

        # create BindingPrediction objects with sequence name and offset
        results = []
        for binding_prediction in binding_predictions:
            peptide = binding_prediction.peptide
            for context in peptide_to_contexts[peptide]:
                results.append(binding_prediction.clone_with_updates(
                    source_sequence_name=context.source_sequence_name,
                    offset=context.offset))
        self._check_results(
            results,
            peptides=peptide_set,
            alleles=self.alleles)
        return BindingPredictionCollection(results)

    def predict_subsequences_dataframe(
            self,
            sequence_dict,
            peptide_lengths=None):
        """Deprecated: use predict_proteins_dataframe() instead."""
        return self.predict_subsequences(
                sequence_dict=sequence_dict,
                peptide_lengths=peptide_lengths).to_dataframe()

    @staticmethod
    def _check_hla_alleles(
            alleles,
            valid_alleles=None,
            keep_unparseable=False):
        """
        Given a list of HLA alleles and an optional list of valid
        HLA alleles, return a set of alleles that we will pass into
        the MHC binding predictor.

        When keep_unparseable is True, allele names the normalizer can't
        parse are kept verbatim (original case) rather than raising, so a
        caller that validates against the tool's own supported list can still
        accept non-human alleles like H-2-Qa1 or BoLA-amani.1.
        """
        require_iterable_of(alleles, str, "HLA alleles")

        # Keep only unique alleles (don't run the predictor twice for a
        # homozygous genotype). When keep_unparseable is set (command-line
        # predictors, which validate against the tool's own -listMHC list),
        # normalize_allele_name_or_raw retains names mhcgnomes can't parse
        # using the same canonical fallback the output parser applies, so a
        # requested allele and its echoed form share one identity. Otherwise
        # normalize_allele_name raises on an unparseable name, as before.
        if keep_unparseable:
            alleles = {normalize_allele_name_or_raw(a) for a in alleles}
        else:
            alleles = {normalize_allele_name(a.strip().upper()) for a in alleles}
        if valid_alleles:
            # For some reason netMHCpan drops the '*' in names, so
            # 'HLA-A*03:01' becomes 'HLA-A03:01'
            missing_alleles = [
                allele
                for allele in alleles
                if allele not in valid_alleles
            ]
            if len(missing_alleles) > 0:
                raise UnsupportedAllele(
                    "Unsupported HLA alleles: %s" % missing_alleles)

        return list(alleles)
