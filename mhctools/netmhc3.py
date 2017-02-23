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
from .parsing import parse_netmhc3_stdout

class NetMHC3(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            program_name="netMHC",
            default_peptide_lengths=[9]):
        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            parse_output_fn=parse_netmhc3_stdout,
            # NetMHC just expects the first arg to be an input FASTA
            input_file_flag="",
            # NetMHC doesn't have the ability to use a custom
            # temporary directory
            tempdir_flag=None,
            length_flag="--peplen",
            allele_flag="--mhc",
            extra_flags=["--nodirect"],
            supported_alleles_flag="-A",
            # because we don't have a tempdir flag, can't run more than
            # one predictor at a time
            process_limit=1,
            default_peptide_lengths=default_peptide_lengths)
