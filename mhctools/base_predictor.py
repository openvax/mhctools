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
from collections import defaultdict

from typechecks import require_iterable_of
from .allele_normalization import normalize_allele_name

from .unsupported_allele import UnsupportedAllele
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Prediction, Kind, PeptideResult

logger = logging.getLogger(__name__)

class BasePredictor(object):
    """
    Base class for all MHC binding predictors.
    """
    def __init__(
            self,
            alleles,
            valid_alleles=None,
            default_peptide_lengths=None,
            min_peptide_length=8,
            max_peptide_length=None,
            allow_X_in_peptides=False,
            allow_lowercase_in_peptides=False):
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
        """
        # I find myself often constructing a predictor with just one allele
        # so as a convenience, allow user to not wrap that allele as a list
        if type(alleles) is str:
            alleles = alleles.split(',')
        self.alleles = self._check_hla_alleles(alleles, valid_alleles)

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

    def predict(self, peptides):
        """
        Predict for a list of peptide sequences.

        Returns
        -------
        list of PeptideResult
        """
        collection = self.predict_peptides(peptides)
        return collection.to_peptide_preds(kind=self._default_pred_kind())

    def predict_dataframe(self, peptides, sample_name=""):
        """predict() flattened to a DataFrame."""
        import pandas as pd
        dfs = [pp.to_dataframe(sample_name) for pp in self.predict(peptides)]
        if not dfs:
            from .pred import COLUMNS
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)

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
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}
        elif isinstance(sequence_dict, (list, tuple)):
            sequence_dict = {seq: seq for seq in sequence_dict}

        peptide_lengths = self._check_peptide_lengths(peptide_lengths)

        peptide_set = set()
        peptide_to_name_offset_pairs = defaultdict(list)

        for name, sequence in sequence_dict.items():
            for peptide_length in peptide_lengths:
                for i in range(len(sequence) - peptide_length + 1):
                    peptide = sequence[i:i + peptide_length]
                    peptide_set.add(peptide)
                    peptide_to_name_offset_pairs[peptide].append((name, i))

        peptide_list = sorted(peptide_set)
        flat_preds = self.predict(peptide_list)

        # Expand: each PeptideResult may map to multiple (name, offset) positions
        results = defaultdict(list)
        for pp in flat_preds:
            if not pp.preds:
                continue
            peptide = pp.preds[0].peptide
            for name, offset in peptide_to_name_offset_pairs.get(peptide, []):
                relocated = PeptideResult(preds=tuple(
                    Prediction(
                        kind=p.kind,
                        score=p.score,
                        peptide=p.peptide,
                        allele=p.allele,
                        value=p.value,
                        percentile_rank=p.percentile_rank,
                        source_sequence_name=name,
                        offset=offset,
                        predictor_name=p.predictor_name,
                        predictor_version=p.predictor_version,
                    ) for p in pp.preds
                ))
                results[name].append(relocated)
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

    # --- deprecated API (still works) ---

    def predict_peptides(self, peptides):
        """
        Deprecated: use predict() instead.

        Given a list of peptide sequences, returns a BindingPredictionCollection.
        """
        raise NotImplementedError(
            "%s must implement predict_peptides" % (self.__class__.__name__,))

    def predict_peptides_dataframe(self, peptides):
        """Deprecated: use predict_dataframe() instead."""
        return self.predict_peptides(peptides).to_dataframe()

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
        if len(observed) != n_expected:
            for a in alleles_set:
                for p in peptides_set:
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
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}
        elif isinstance(sequence_dict, (list, tuple)):
            sequence_dict = {seq: seq for seq in sequence_dict}

        peptide_lengths = self._check_peptide_lengths(peptide_lengths)

        # convert long protein sequences to set of peptides and
        # associated sequence name / offsets that each peptide may have come
        # from
        peptide_set = set([])
        peptide_to_name_offset_pairs = defaultdict(list)

        for name, sequence in sequence_dict.items():
            for peptide_length in peptide_lengths:
                for i in range(len(sequence) - peptide_length + 1):
                    peptide = sequence[i:i + peptide_length]
                    peptide_set.add(peptide)
                    peptide_to_name_offset_pairs[peptide].append((name, i))
        peptide_list = sorted(peptide_set)

        binding_predictions = self.predict_peptides(peptide_list)

        # create BindingPrediction objects with sequence name and offset
        results = []
        for binding_prediction in binding_predictions:
            peptide = binding_prediction.peptide
            for name, offset in peptide_to_name_offset_pairs[peptide]:
                results.append(binding_prediction.clone_with_updates(
                    source_sequence_name=name,
                    offset=offset))
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
            valid_alleles=None):
        """
        Given a list of HLA alleles and an optional list of valid
        HLA alleles, return a set of alleles that we will pass into
        the MHC binding predictor.
        """
        require_iterable_of(alleles, str, "HLA alleles")

        # Don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        alleles = {
            normalize_allele_name(allele.strip().upper())
            for allele in alleles
        }
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
