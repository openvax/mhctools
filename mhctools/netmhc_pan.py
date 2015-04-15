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
            hla_alleles,
            netmhc_command="netMHCpan",
            epitope_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHCpan",
            command=netmhc_command,
            hla_alleles=hla_alleles,
            epitope_lengths=epitope_lengths)

    def predict(self, df, mutation_window_size=None):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence
        through NetMHCpan.
        If mutation_window_size is not None then only make predictions for that
        number residues away from mutations.

        Expects the input DataFrame to have the following fields:
            - SourceSequence
            - MutationStart
            - MutationEnd
            - GeneInfo
            - Gene
            - GeneMutationInfo
            - PeptideMutationInfo
            - TranscriptId
        """

        input_filename, peptide_entries = create_input_fasta_file(
            df,
            mutation_window_size=mutation_window_size
        )

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
            results = parse_xls_file(
                output_file.read(),
                peptide_entries,
                mutation_window_size=mutation_window_size)

        if len(results) == 0:
            raise ValueError("No epitopes from netMHCpan")
        return pd.DataFrame.from_records(results)
