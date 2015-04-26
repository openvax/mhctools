import logging
from subprocess import check_output

from .alleles import normalize_allele_name
from .base_predictor import BasePredictor
from .process_helpers import run_command


class BaseCommandlinePredictor(BasePredictor):
    """
    Base class for MHC binding predictors that run a local external
    program and write their output to a local file.
    """
    def __init__(
            self,
            name,
            command,
            alleles,
            epitope_lengths,
            supported_allele_flag='-listMHC'):
        self.name = name
        self.command = command

        if not isinstance(command, str):
            raise TypeError(
                "Expected %s command to be string, got %s : %s" % (
                    name, command, type(command)))

        try:
            run_command([command])
        except:
            raise SystemError("Failed to run %s" % command)

        valid_alleles = self._determine_valid_alleles(
            command, supported_allele_flag)

        BasePredictor.__init__(
            self,
            alleles,
            epitope_lengths,
            valid_alleles=valid_alleles)

    @staticmethod
    def _determine_valid_alleles(command, supported_allele_flag):
        """
        Try asking the commandline predictor (e.g. netMHCpan)
        which alleles it supports.
        """
        valid_alleles = None
        if supported_allele_flag:
            try:
                # convert to str since Python3 returns a `bytes` object
                valid_alleles = check_output([
                    command, supported_allele_flag
                ])
                valid_alleles_str = valid_alleles.decode("ascii", "ignore")
                assert len(valid_alleles_str) > 0, \
                    '%s returned empty allele list' % command
                valid_alleles = set([])
                for line in valid_alleles_str.split("\n"):
                    line = line.strip()
                    if not line.startswith('#') and len(line) > 0:
                        try:
                            allele = normalize_allele_name(line)
                            valid_alleles.add(allele)
                        except ValueError as error:
                            logging.info("Skipping allele %s: %s" % (
                                line, error))
                            continue
            except:
                logging.warning(
                    "Failed to run %s %s",
                    command,
                    supported_allele_flag)
                raise
        return valid_alleles
