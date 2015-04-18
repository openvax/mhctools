import tempfile
import logging
import pandas as pd

from .cleanup_context import CleanupFiles
from .file_formats import create_input_fasta_file, parse_xls_file
from .base_commandline_predictor import BaseCommandlinePredictor
from .process_helpers import run_command

class NetMHCpan(BaseCommandlinePredictor):

    def __init__(
            self,
            alleles,
            netmhc_command="netMHCpan",
            epitope_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHCpan",
            command=netmhc_command,
            alleles=alleles,
            epitope_lengths=epitope_lengths)

    def predict(
            self,
            fasta_dictionary,
            raise_on_error=True):
        input_filename = create_input_fasta_file(fasta_dictionary)

        alleles_str = \
            ",".join(allele.replace("*", "") for allele in self.alleles)
        output_file = tempfile.NamedTemporaryFile(
                "r+",
                prefix="netMHCpan_output",
                delete=False)
        command = [
            self.command,
            "-xls",
            "-xlsfile", output_file.name,
            "-l", "9",
            "-f", input_filename,
            "-a", alleles_str
        ]
        logging.info(" ".join(command))

        with CleanupFiles(
                filenames=[input_filename],
                files=[output_file]):
            run_command(command)
            file_contents = output_file.read()
            epitope_collection = parse_xls_file(
                file_contents,
                fasta_dictionary=fasta_dictionary,
                prediction_method_name="netmhcpan")

        # TODO(tavi) Unwise to just return an empty DataFrame.
        if len(epitope_collection) == 0:
            if raise_on_error:
                raise ValueError("No epitopes from netMHCpan")
            return pd.DataFrame()
        return epitope_collection
