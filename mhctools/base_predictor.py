# Copyright (c) 2014. Mount Sinai School of Medicine
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

from .alleles import normalize_allele_name


class UnsupportedAllele(ValueError):
    pass


class BasePredictor(object):
    """
    Base class for all MHC binding predictors.
    """
    def __init__(
            self,
            alleles,
            epitope_lengths,
            valid_alleles=None):
        """
        Parameters
        ----------
        alleles : list
            List of strings containing names of HLA alleles we're
            making predictions for. Example:
                ["HLA-A*02:01", "HLA-B*07:02"]

        epitope_lengths : list or int
            List of epitope lengths to make predictions for, or
            a single integer length.

        valid_alleles : list, optional
            If given, constrain HLA alleles to be contained within
            this set.
        """
        # I find myself often constructing a predictor with just one allele
        # so as a convenience, allow user to not wrap that allele as a list
        if isinstance(alleles, str):
            alleles = alleles.split(',')
        self.alleles = self._check_hla_alleles(alleles, valid_alleles)

        if isinstance(epitope_lengths, int):
            epitope_lengths = [epitope_lengths]

        if not isinstance(epitope_lengths, list):
            raise TypeError(
                'Expected epitope_lengths : list, got %s : %s' % (
                    epitope_lengths,
                    type(epitope_lengths)))
        for length in epitope_lengths:
            if not isinstance(length, int):
                raise TypeError(
                    ('Element of epitope_lengths must be int, got '
                     '%s : %s') % (length, type(length)))

        self.epitope_lengths = epitope_lengths

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "%s(alleles=%s, epitope_lengths=%s)" % (
            self.__class__.__name__,
            self.alleles,
            self.epitope_lengths)

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
