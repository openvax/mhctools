# Copyright (c) 2015. Mount Sinai School of Medicine
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
import tempfile
import logging
from os import remove

from .base_commandline_predictor import BaseCommandlinePredictor
from .cleanup_context import CleanupFiles
from .common import check_sequence_dictionary
from .epitope_collection import EpitopeCollection
from .file_formats import create_input_fasta_files, parse_netmhc_stdout
from .process_helpers import AsyncProcess

class NetMHC(BaseCommandlinePredictor):

    def __init__(
            self,
            alleles,
            netmhc_command="netMHC",
            epitope_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHC",
            command=netmhc_command,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            supported_allele_flag="-A")

    def predict(self, fasta_dictionary):
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        input_filenames, sequence_key_mapping = create_input_fasta_files(
            fasta_dictionary)

        assert len(input_filenames) == 1, \
            "Unexpected number of input files: %s" % input_filenames

        # TODO: We are not currently using the file chunking
        # functionality here. See NetMHCcons.
        input_filename = input_filenames[0]

        alleles_str = \
            ",".join(allele.replace("*", "") for allele in self.alleles)

        epitope_collections = []
        # dictionary from output file paths to commands
        commands = {}
        for peptide_length in self.epitope_lengths:
            output_file = tempfile.NamedTemporaryFile(
                    "w",
                    prefix="netMHC_output_peplen%d" % peptide_length,
                    delete=False)

            command = [
                self.command,
                input_filename,
                "--peplen", str(peptide_length),
                "--nodirect",  # approximate 10mer predictions using 9mers
                "--mhc", alleles_str
            ]
            commands[output_file] = command

        for output_file, command in commands.items():
            with CleanupFiles(files=[output_file]):
                process = AsyncProcess(
                    args=command,
                    redirect_stdout_file=output_file)
                print(" ".join(command))
                process.wait()
                # need to flush written output and re-open for read
                output_file.close()
                with open(output_file.name, 'r') as f:
                    file_contents = f.read()
                    epitope_collection = parse_netmhc_stdout(
                        file_contents,
                        sequence_key_mapping=sequence_key_mapping,
                        fasta_dictionary=fasta_dictionary,
                        prediction_method_name="netmhc")

                    if len(epitope_collection) == 0:
                        logging.warn(file_contents)
                        raise ValueError("No epitopes from netMHC")
                    epitope_collections.append(epitope_collection)
        remove(input_filename)
        # flatten all epitope collections into a single object
        return EpitopeCollection([
            binding_prediction
            for sublist in epitope_collections
            for binding_prediction in sublist])
