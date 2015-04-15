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

import logging
import tempfile
import pandas as pd

from .process_helpers import run_multiple_commands_redirect_stdout
from .cleanup_context import CleanupFiles
from .file_formats import create_input_fasta_file, parse_netmhc_stdout
from .base_commandline_predictor import BaseCommandlinePredictor

class NetMHCcons(BaseCommandlinePredictor):
    def __init__(
            self,
            hla_alleles,
            netmhc_command="netMHCcons",
            epitope_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHCcons",
            command=netmhc_command,
            hla_alleles=hla_alleles,
            epitope_lengths=epitope_lengths)

    def predict(self, df, mutation_window_size=None):
        """
        Given a dataframe of mutated amino acid sequences, run each sequence
        through NetMHCcons.
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
            mutation_window_size=mutation_window_size)

        output_files = {}
        commands = {}
        dirs = []
        for i, allele in enumerate(self.alleles):
            for length in self.epitope_lengths:
                temp_dirname = tempfile.mkdtemp(
                    prefix="tmp_netmhccons_length_%d" % length)
                logging.info(
                    "Created temporary directory %s for allele %s, length %d",
                    temp_dirname,
                    allele,
                    length)
                dirs.append(temp_dirname)
                output_file = tempfile.NamedTemporaryFile(
                        "w",
                        prefix="netMHCcons_output_%d" % i,
                        delete=False)
                command = [
                    self.command,
                    "-length", str(length),
                    "-f", input_filename,
                    "-a", allele,
                    "-tdir", temp_dirname
                ]
                commands[output_file] = command

        results = []

        # Cleanup either when finished or if an exception gets raised by
        # deleting the input and output files
        filenames_to_delete = [input_filename]
        for f in output_files.keys():
            filenames_to_delete.append(f.name)

        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=dirs):
            run_multiple_commands_redirect_stdout(
                commands, print_commands=True)
            for output_file, command in commands.iteritems():
                # closing/opening looks insane
                # but I was getting empty files otherwise
                output_file.close()
                with open(output_file.name, 'r') as f:
                    rows = parse_netmhc_stdout(
                        f.read(),
                        peptide_entries,
                        mutation_window_size=mutation_window_size)
                results.extend(rows)

        assert len(results) > 0, "No epitopes from netMHCcons"

        df = pd.DataFrame.from_records(results)
        unique_alleles = set(df.Allele)
        assert len(unique_alleles) == len(self.alleles), \
            "Expected %d alleles (%s) but got %d (%s)" % (
                len(self.alleles), self.alleles,
                len(unique_alleles), unique_alleles)
        return df
