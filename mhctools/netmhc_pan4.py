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

from .netmhc_pan4_base import NetMHCpan4Base
from .parsing import parse_netmhcpan4_stdout
from .pred import Kind


class NetMHCpan4(NetMHCpan4Base):
    parse_stdout_fn = staticmethod(parse_netmhcpan4_stdout)
    affinity_mode_kinds = (Kind.pMHC_affinity,)

    def __init__(
            self,
            alleles=None,
            default_peptide_lengths=[9],
            program_name="netMHCpan-4.0",
            process_limit=-1,
            mode="binding_affinity",
            extra_flags=[],
            max_peptides_per_file=10 ** 4,
            max_alleles_per_command="auto"):
        """
        Wrapper for NetMHCpan4.

        The mode argument should be one of "binding_affinity" (default) or
        "elution_score".
        """
        NetMHCpan4Base.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            program_name=program_name,
            process_limit=process_limit,
            mode=mode,
            extra_flags=extra_flags,
            max_peptides_per_file=max_peptides_per_file,
            max_alleles_per_command=max_alleles_per_command)


class NetMHCpan4_EL(NetMHCpan4):
    """
    Wrapper for NetMHCpan4 when the preferred mode is elution score
    """
    def __init__(
            self,
            alleles=None,
            default_peptide_lengths=[9],
            program_name="netMHCpan-4.0",
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
            alleles=None,
            default_peptide_lengths=[9],
            program_name="netMHCpan-4.0",
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
