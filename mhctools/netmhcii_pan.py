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

from .alleles import parse_classi_or_classii_allele_name
from .base_commandline_predictor import BaseCommandlinePredictor
from .file_formats import parse_netmhciipan_stdout


class NetMHCIIpan(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            epitope_lengths=[15, 16, 17, 18, 19, 20],
            program_name="netMHCIIpan",
            max_file_records=None,
            process_limit=0):
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            parse_output_fn=parse_netmhciipan_stdout,
            supported_alleles_flag="-list",
            input_fasta_flag="-f",
            allele_flag="-a",
            length_flag="-length",
            tempdir_flag="-tdir",
            max_file_records=max_file_records,
            process_limit=process_limit)

    def prepare_allele_name(self, allele_name):
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
            allele = parsed_alleles[0]
            if allele.species == "H-2":
                return "%s-%s%s" % (
                    allele.species,
                    allele.gene,
                    allele.allele_code)
            else:
                # assume that we're dealing with a human DRB allele
                # which NetMHCIIpan treats differently because there is
                # supposedly no population diversity in the DR-alpha gene
                return "%s_%s%s" % (
                    allele.gene,
                    allele.allele_family,
                    allele.allele_code)
        else:
            alpha, beta = parsed_alleles
            return "HLA-%s%s%s-%s%s%s" % (
                alpha.gene,
                alpha.allele_family,
                alpha.allele_code,
                beta.gene,
                beta.allele_family,
                beta.allele_code)
