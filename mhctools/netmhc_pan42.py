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

from .base_commandline_predictor import BaseCommandlinePredictor
from .parsing import parse_netmhcpan41_stdout, parse_netmhcpan_to_preds
from functools import partial


class NetMHCpan42(BaseCommandlinePredictor):
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            mode="binding_affinity",
            extra_flags=[]):
        """
        Wrapper for NetMHCpan 4.2.

        Output format is identical to 4.1 (Score_EL, %Rank_EL, Score_BA,
        %Rank_BA, Aff(nM) columns).

        The mode argument should be one of "binding_affinity" (default) or
        "elution_score".
        """
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
            parse_output_fn=partial(parse_netmhcpan41_stdout, mode=mode),
            parse_to_preds_fn=parse_netmhcpan_to_preds,
            supported_alleles_flag="-listMHC",
            input_file_flag="-f",
            length_flag="-l",
            allele_flag="-a",
            extra_flags=flags + extra_flags,
            process_limit=process_limit)


class NetMHCpan42_EL(NetMHCpan42):
    """NetMHCpan 4.2 in elution score mode."""
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCpan42.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="elution_score",
            extra_flags=extra_flags)


class NetMHCpan42_BA(NetMHCpan42):
    """NetMHCpan 4.2 in binding affinity mode."""
    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            program_name="netMHCpan",
            process_limit=-1,
            extra_flags=[]):
        NetMHCpan42.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode="binding_affinity",
            extra_flags=extra_flags)
