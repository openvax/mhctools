# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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
from .parsing import parse_netmhcpan4_stdout
from functools import partial

class NetMHCpan4(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            mode="binding_affinity",
            extra_flags=[]):
        """
        Wrapper for NetMHCpan4.

        The mode argument should be one of "binding_affinity" (default) or
        "elution_score".
        """

        # The -BA flag is required to predict binding affinity
        if mode == "binding_affinity":
            flags = ["-BA"]
        elif mode == "elution_score":
            flags = []
        else:
            raise ValueError("Unsupported mode", mode)

        BaseCommandlinePredictor.__init__(
            self,
            program_name=program_name,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            parse_output_fn=partial(parse_netmhcpan4_stdout, mode=mode),
            supported_alleles_flag="-listMHC",
            input_file_flag="-f",
            length_flag="-l",
            allele_flag="-a",
            extra_flags=flags + extra_flags,
            process_limit=process_limit)

        self.mode = mode

class NetMHCpan4_EL(NetMHCpan4):
    """
    Wrapper for NetMHCpan4 when the preferred mode is elution score
    """
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCpan4.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="elution_score",
            extra_flags=extra_flags)

class NetMHCpan4_BA(NetMHCpan4):
    """
    Wrapper for NetMHCpan4 when the preferred mode is binding affinity
    """
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCpan4.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="binding_affinity",
            extra_flags=extra_flags)