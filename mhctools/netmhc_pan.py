import tempfile
import logging

from .base_commandline_predictor import BaseCommandlinePredictor
from .cleanup_context import CleanupFiles
from .file_formats import create_input_fasta_file, parse_netmhc_stdout
from .process_helpers import AsyncProcess

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

    def predict(self, fasta_dictionary):
        input_filename, sequence_key_mapping = create_input_fasta_file(
            fasta_dictionary)

        alleles_str = \
            ",".join(allele.replace("*", "") for allele in self.alleles)
        output_file = tempfile.NamedTemporaryFile(
                "w",
                prefix="netMHCpan_output",
                delete=False)
        args = [
            self.command,
            "-l", "9",
            "-f", input_filename,
            "-a", alleles_str
        ]
        logging.info(" ".join(args))

        with CleanupFiles(
                filenames=[input_filename],
                files=[output_file]):
            process = AsyncProcess(
                args=args,
                redirect_stdout_file=output_file)
            process.wait()
            # need to flush written output and re-open for read
            output_file.close()
            with open(output_file.name, 'r') as f:
                file_contents = f.read()
                epitope_collection = parse_netmhc_stdout(
                    file_contents,
                    sequence_key_mapping=sequence_key_mapping,
                    fasta_dictionary=fasta_dictionary,
                    prediction_method_name="netmhcpan")

        if len(epitope_collection) == 0:
            logging.warn(file_contents)
            raise ValueError("No epitopes from netMHCpan")
        return epitope_collection
