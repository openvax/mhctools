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

from .alleles import parse_classi_or_classii_allele_name
from .base_commandline_predictor import BaseCommandlinePredictor
from .cleanup_context import CleanupFiles
from .common import check_sequence_dictionary, seq_to_str
from .file_formats import create_input_fasta_files, parse_netmhcpan_stdout
from .process_helpers import AsyncProcess

class NetMHCIIpan(BaseCommandlinePredictor):

    def __init__(
            self,
            alleles,
            netmhc_command="netMHCIIpan",
            epitope_lengths=[15, 16, 17, 18, 19, 20]):
        BaseCommandlinePredictor.__init__(
            self,
            name="NetMHCIIpan",
            command=netmhc_command,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            supported_allele_flag="-list")

    def normalize_allele(self, allele_name):
        """
        netMHCIIpan has some unique requirements for allele formats,
        expecting the following forms:
         - DRB1_0101 (for non-alpha/beta pairs)
         - HLA-DQA10501-DQB10636 (for alpha and beta pairs)

        Other than human class II alleles, the only other alleles that
        netMHCIIpan accepts are the following mouse alleles:
         - H-2-IAb
         - H-2-IAd
        """
        parsed_alleles = parse_classi_or_classii_allele_name(allele_name)
        if len(parsed_alleles) == 1:
            # Handle mice differently, as netMHCIIpan has a different format
            # for them.
            if parsed_alleles[0].species == "H-2":
                return "%s-%s%s" % (parsed_alleles[0].species,
                                    parsed_alleles[0].gene,
                                    parsed_alleles[0].allele_code)
            return "%s_%s%s" % (parsed_alleles[0].gene,
                                parsed_alleles[0].allele_family,
                                parsed_alleles[0].allele_code)
        else:
            return "HLA-%s%s%s-%s%s%s" % (
                parsed_alleles[0].gene,
                parsed_alleles[0].allele_family,
                parsed_alleles[0].allele_code,
                parsed_alleles[1].gene,
                parsed_alleles[1].allele_family,
                parsed_alleles[1].allele_code)

    def predict(self, fasta_dictionary):
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        input_filenames, sequence_key_mapping = create_input_fasta_files(
            fasta_dictionary)
        # TODO: We are not currently using the file chunking
        # functionality here. See NetMHCcons.
        input_filename = input_filenames[0]

        alleles_str = \
            ",".join(self.normalize_allele(allele) for allele in self.alleles)
        output_file = tempfile.NamedTemporaryFile(
                "w",
                prefix="netMHCIIpan_output",
                delete=False)
        args = [
            self.command,
            "-length", seq_to_str(self.epitope_lengths, sep=" "),
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
                epitope_collection = parse_netmhcpan_stdout(
                    file_contents,
                    sequence_key_mapping=sequence_key_mapping,
                    fasta_dictionary=fasta_dictionary,
                    prediction_method_name="netmhciipan",
                    is_netmhcpanii=True)

        if len(epitope_collection) == 0:
            logging.warn(file_contents)
            raise ValueError("No epitopes from netMHCIIpan")
        return epitope_collection
