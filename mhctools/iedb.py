# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Historical IEDB interfaces backed by local tools; no HTTP fallback.

NetMHCpan uses 4.1 BA; NetMHCIIpan uses 4.3 BA. HTTP-only arguments
(``url``, ``request_timeout``, ``raise_on_error``) are no longer supported.
Class-I compatibility wrappers retain their 8–11 residue default windows.
"""

from .netmhc_cons import NetMHCcons
from .netmhc_pan41 import NetMHCpan41_BA
from .netmhcii_pan import NetMHCIIpan43_BA as IedbNetMHCIIpan
from .pred import Kind
from .smm import SMM


class IedbNetMHCcons(NetMHCcons):
    """Local NetMHCcons with the historical IEDB window defaults."""

    def __init__(self, alleles=None, default_peptide_lengths=(8, 9, 10, 11),
                 program_name="netMHCcons", process_limit=0):
        super().__init__(alleles=alleles, default_peptide_lengths=default_peptide_lengths,
                         program_name=program_name, process_limit=process_limit)


class IedbNetMHCpan(NetMHCpan41_BA):
    """Local NetMHCpan 4.1 with the historical affinity-only IEDB contract."""

    def __init__(self, alleles=None, default_peptide_lengths=(8, 9, 10, 11),
                 program_name="netMHCpan-4.1", process_limit=-1, extra_flags=None):
        super().__init__(alleles=alleles, default_peptide_lengths=default_peptide_lengths,
                         program_name=program_name, process_limit=process_limit,
                         extra_flags=extra_flags or [])
        # The native binding parser preserves affinity-only output across
        # predict(), predict_pairs(), and the legacy BindingPrediction API.
        self.parse_to_preds_fn = None

    def kind_support(self):
        return {Kind.pMHC_affinity: super().kind_support()[Kind.pMHC_affinity]}


class IedbSMM(SMM):
    """Local SMM with the historical IEDB window defaults."""

    def __init__(self, alleles=None, default_peptide_lengths=(8, 9, 10, 11),
                 program_name=None, timeout=120):
        super().__init__(alleles=alleles, default_peptide_lengths=default_peptide_lengths,
                         program_name=program_name, timeout=timeout)


class IedbSMM_PMBEC(IedbSMM):
    """Local SMM-PMBEC with the historical IEDB window defaults."""

    prediction_method = "smmpmbec"


__all__ = ["IedbNetMHCcons", "IedbNetMHCpan", "IedbNetMHCIIpan", "IedbSMM", "IedbSMM_PMBEC"]
