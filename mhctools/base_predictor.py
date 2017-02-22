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

from typechecks import require_iterable_of
from mhcnames import normalize_allele_name

from .unsupported_allele import UnsupportedAllele

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
        if isinstance(alleles, str):
            alleles = alleles.split(',')
        self.alleles = self._check_hla_alleles(alleles, valid_alleles)

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
        peptide_lengths = self._check_peptide_lengths(peptide_lengths)
        peptides = []
        source_sequence_names = []
        offsets = []
        for name, sequence in sequence_dict.items():
            for peptide_length in peptide_lengths:
                for i in range(len(sequence) - peptide_length + 1):
                    peptides.append(sequence[i:i + peptide_length])
                    offsets.append(i)
                    source_sequence_names.append(name)
        binding_predictions = self.predict_peptides(peptides)
        return binding_predictions.update_fields(
            offset=offsets,
            source_sequence_name=source_sequence_names)

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
        alleles = [
            normalize_allele_name(allele.strip().upper())
            for allele in alleles
        ]
        if valid_alleles:
            # For some reason netMHCpan drops the '*' in names, so
            # 'HLA-A*03:01' becomes 'HLA-A03:01'
            missing_alleles = [
                allele
                for allele in alleles
                if allele not in valid_alleles
            ]
            if len(missing_alleles) > 0:
                raise UnsupportedAllele("Unsupported HLA alleles: %s" % missing_alleles)

        # Don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        alleles = list(set(alleles))
        return alleles
