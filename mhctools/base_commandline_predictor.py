from base_predictor import BasePredictor
from process_helpers import run_command

class BaseCommandlinePredictor(BasePredictor):
    """
    Base class for MHC predictors which run a local external program
    and write their output to a local file.
    """

    def __init__(
            self,
            name,
            command,
            hla_alleles,
            epitope_lengths,
            supported_allele_flag="-listMHC"):

        self.name = name

        if not isinstance(command, str):
            raise TypeError(
                "Expected %s command to be string, got %s : %s" % (
                    name, command, type(command)))

        try:
            run_command([command])
        except:
            raise SystemError("Failed to run %s" % command)

        self.command = command


        # try asking the commandline predictor (e.g. netMHCpan)
        # which alleles it supports
        valid_alleles = None
        if supported_allele_flag:
            try:
                valid_alleles_str = check_output([self.netmhc_command, "-listMHC"])

                assert len(valid_alleles_str) > 0, \
                    "%s returned empty allele list" % self.self.netmhc_command
                valid_alleles = set([])
                for line in valid_alleles_str.split("\n"):
                    if not line.startswith("#"):
                        valid_alleles.add(line)
            except:
                logging.warning("Failed to run %s -listMHC", self.netmhc_command)


        MHCBasePredictor.__init__(
            hla_alleles,
            epitope_lengths,
            valid_alleles=valid_alleles)