# Copyright (c) 2015-2017. Mount Sinai School of Medicine
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

from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhc4_stdout

class NetMHC4(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            program_name="netMHC",
            process_limit=0,
            default_peptide_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            parse_output_fn=parse_netmhc4_stdout,
            input_file_flag="-f",
            tempdir_flag="-tdir",
            length_flag="-l",
            allele_flag="-a",
            supported_alleles_flag="-listMHC",
            process_limit=process_limit,
            default_peptide_lengths=default_peptide_lengths)

    def prepare_allele_name(self, allele_name):
        allele_name = super(NetMHC4, self).prepare_allele_name(allele_name)
        return allele_name.replace(":", "")
