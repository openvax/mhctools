import logging

from .common import normalize_hla_allele_name


class BasePredictor(object):
    """
    Base class for all MHC binding predictors.
    """
    def __init__(self, hla_alleles, epitope_lengths,
                 valid_alleles=None):
        """
        Parameters
        ----------
        hla_alleles : list
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
        self.alleles = self._process_alleles(hla_alleles,
                                             valid_alleles)

        if isinstance(epitope_lengths, int):
            epitope_lengths = [epitope_lengths]

        if not isinstance(epitope_lengths, list):
            raise TypeError(
                'Expected epitope_lengths : list, got %s : %s' % (
                epitope_lengths, type(epitope_lengths)))
        for length in epitope_lengths:
            if not isinstance(length, int):
                raise TypeError(
                    ('Element of epitope_lengths must be int, got '
                     '%s : %s') % (length, type(length)))

        self.epitope_lengths = epitope_lengths

    @staticmethod
    def _process_alleles(hla_alleles, valid_alleles=None):
        """
        Given a list of HLA alleles and an optional list of valid
        HLA alleles, return a set of alleles that we will pass into
        the MHC binding predictor.
        """
        alleles = []
        for allele in hla_alleles:
            if not isinstance(allele, str):
                raise TypeError(
                    'Expected allele to be string, got %s : %s' % (
                        allele, type(allele)))
            allele = normalize_hla_allele_name(
                allele.strip().upper())

            # For some reason netMHCpan drops the '*' in names, so
            # 'HLA-A*03:01' becomes 'HLA-A03:01'
            if (valid_alleles and
                allele.replace('*', '') not in valid_alleles):
                logging.warning('Skipping %s (not available)'
                                % allele)
            else:
                alleles.append(allele)

        # Don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        alleles = set(alleles)
        return alleles
