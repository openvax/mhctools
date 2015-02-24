import logging

from mhc_common import normalize_hla_allele_name


class BasePredictor(object):
    """
    Base class for all MHC predictors.
    """
    def __init__(self, hla_alleles, epitope_lengths, valid_alleles=None):
        """
        Parameters
        ----------
        hla_alleles : list
            List of strings containing names of HLA alleles we're making
            predictions for. Example: ["HLA-A*02:01", "HLA-B*07:02"]

        epitope_lengths : list
            List of epitope lengths to make predictions for.

        valid_alleles : list, optional.
            If given, constrain HLA alleles to be contained within this set.
        """

        self.alleles = []
        for allele in hla_alleles:
            if not isinstance(allele, str):
                raise TypeError(
                    "Expected allele to be string, got %s : %s" % (
                        allele, type(allele)))
            allele = normalize_hla_allele_name(allele.strip().upper())
            # for some reason netMHCpan drop the "*" in names
            # such as "HLA-A*03:01" becomes "HLA-A03:01"
            if valid_alleles and allele.replace("*", "") not in valid_alleles:
                logging.warning("Skipping %s (not available)" % allele)
            else:
                self.alleles.append(allele)

        # don't run the MHC predictor twice for homozygous alleles,
        # only run it for unique alleles
        self.alleles = set(self.alleles)

        if isinstance(epitope_lengths, int):
            epitope_lengths = [epitope_lengths]

        if not isinstance(epitope_lengths, list):
            raise TypeError("Expected epitope_lengths : list, got %s : %s" % (
                epitope_lengths, type(epitope_lengths)))
        for length in epitope_lengths:
            if not isisntance(length, int):
                raise TypeError(
                    "Element of epitope_lengths must be int, got %s : %s" % (
                        length, type(length)))

        self.epitope_lengths = epitope_lengths

