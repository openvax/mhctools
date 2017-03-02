# Copyright (c) 2014-2017. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
import logging
from collections import defaultdict

from six import string_types
from typechecks import require_iterable_of
from mhcnames import normalize_allele_name

from .unsupported_allele import UnsupportedAllele
from .binding_prediction_collection import BindingPredictionCollection

logger = logging.getLogger(__name__)

class BasePredictor(object):
    """
    Base class for all MHC binding predictors.
    """
    def __init__(
            self,
            alleles,
            valid_alleles=None,
            default_peptide_lengths=None):
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
        """
        # I find myself often constructing a predictor with just one allele
        # so as a convenience, allow user to not wrap that allele as a list
        if isinstance(alleles, string_types):
            alleles = alleles.split(',')
        self.alleles = self._check_hla_alleles(alleles, valid_alleles)

        if isinstance(default_peptide_lengths, int):
            default_peptide_lengths = [default_peptide_lengths]
        require_iterable_of(default_peptide_lengths, int)
        self.default_peptide_lengths = default_peptide_lengths

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(alleles=%s, default_peptide_lengths=%s)" % (
            self.__class__.__name__,
            self.alleles,
            self.default_peptide_lengths)

    def predict_peptides(self, peptides):
        """
        Given a list of peptide sequences, returns a BindingPredictionCollection
        """
        raise NotImplementedError("%s must implement predict_peptides" % (
            self.__class__.__name__))

    def predict_peptides_dataframe(self, peptides):
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
        return peptide_lengths

    def predict_subsequences(
            self,
            sequence_dict,
            peptide_lengths=None):
        """
        Given a dictionary mapping sequence names to amino acid strings,
        and an optional list of peptide lengths, returns a
        BindingPredictionCollection.
        """
        if isinstance(sequence_dict, string_types):
            sequence_dict = {"seq": sequence_dict}

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
        binding_predictions = self.predict_peptides(sorted(peptide_set))
        n_expected = len(peptide_set) * len(self.alleles)
        if len(binding_predictions) != n_expected:
            raise ValueError("Expected %d peptide predictions but got %d" % (
                n_expected, len(binding_predictions)))

        # create BindingPrediction objects with sequence name and offset
        results = []
        for binding_prediction in binding_predictions:
            for name, offset in peptide_to_name_offset_pairs[
                    binding_prediction.peptide]:
                results.append(binding_prediction.clone_with_updates(
                    source_sequence_name=name,
                    offset=offset))
        assert len(results) >= len(binding_predictions)
        return BindingPredictionCollection(results)

    def predict(self, sequence_dict, peptide_lengths=None):
        logger.warn("Deprecated method 'predict', use 'predict_subsequences")
        return self.predict_subsequences(sequence_dict, peptide_lengths=None)

    def predict_subsequences_dataframe(
            self,
            sequence_dict,
            peptide_lengths=None):
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
        require_iterable_of(alleles, string_types, "HLA alleles")

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
