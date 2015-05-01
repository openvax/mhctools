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
import tempfile
import logging

from .base_commandline_predictor import BaseCommandlinePredictor
from .cleanup_context import CleanupFiles
from .common import check_sequence_dictionary, seq_to_str
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
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
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
            "-l", seq_to_str(self.epitope_lengths),
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
