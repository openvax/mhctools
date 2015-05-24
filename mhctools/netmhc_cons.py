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
import logging
import tempfile

from .base_commandline_predictor import BaseCommandlinePredictor
from .cleanup_context import CleanupFiles
from .common import check_sequence_dictionary
from .epitope_collection import EpitopeCollection
from .file_formats import create_input_fasta_files, parse_netmhc_stdout
from .process_helpers import run_multiple_commands_redirect_stdout

class NetMHCcons(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            netmhc_command="netMHCcons",
            epitope_lengths=[9],
            max_file_records=None,
            process_limit=0):
        self.max_file_records = max_file_records
        self.process_limit = process_limit
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHCcons",
            command=netmhc_command,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            # netMHCcons does not have a supported allele flag
            supported_allele_flag=None)

    def predict(self, fasta_dictionary):
        """
        Given a dictionary mapping sequence identifiers to amino acid sequences,
        return an EpitopeCollection of binding predictions.
        """
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        input_filenames, sequence_key_mapping = create_input_fasta_files(
            fasta_dictionary, max_file_records=self.max_file_records)
        output_files = {}
        commands = {}
        dirs = []
        alleles = [allele.replace("*", "") for allele in self.alleles]
        for i, input_filename in enumerate(input_filenames):
            for j, allele in enumerate(alleles):
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
                        prefix="netMHCcons_output_%d_%d" % (i, j),
                        delete=False)
                    command = [
                        self.command,
                        "-length", str(length),
                        "-f", input_filename,
                        "-a", allele,
                        "-tdir", temp_dirname
                    ]
                    commands[output_file] = command

        epitope_collections = []

        # Cleanup either when finished or if an exception gets raised by
        # deleting the input and output files
        filenames_to_delete = input_filenames
        for f in output_files.keys():
            filenames_to_delete.append(f.name)

        with CleanupFiles(
                filenames=filenames_to_delete,
                directories=dirs):
            run_multiple_commands_redirect_stdout(
                commands, print_commands=True,
                process_limit=self.process_limit)
            for output_file, command in commands.iteritems():
                # closing/opening looks insane
                # but I was getting empty files otherwise
                output_file.close()
                with open(output_file.name, 'r') as f:
                    epitope_collection = parse_netmhc_stdout(
                        netmhc_output=f.read(),
                        fasta_dictionary=fasta_dictionary,
                        sequence_key_mapping=sequence_key_mapping,
                        prediction_method_name="netmhccons")
                    epitope_collections.append(epitope_collection)

        if len(epitope_collections) != len(commands):
            raise ValueError("Expected an epitope collection for each "
                             "command (%d), but instead there are %d" %
                             (len(commands), len(epitope_collections)))

        if len(epitope_collections) == 0:
            raise ValueError("No epitopes from netMHCcons")

        # flatten all epitope collections into a single object
        return EpitopeCollection([
            binding_prediction
            for sublist in epitope_collections
            for binding_prediction in sublist])
