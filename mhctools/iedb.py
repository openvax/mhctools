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

"""Historical IEDB names now resolve to local predictors; no HTTP fallback.

NetMHCpan uses 4.1 BA; NetMHCIIpan uses 4.3 BA. HTTP-only arguments
(``url``, ``request_timeout``, ``raise_on_error``) are no longer supported.
"""

from .netmhc_cons import NetMHCcons as IedbNetMHCcons
from .netmhc_pan41 import NetMHCpan41_BA as IedbNetMHCpan
from .netmhcii_pan import NetMHCIIpan43_BA as IedbNetMHCIIpan
from .smm import SMM as IedbSMM, SMMPMBEC as IedbSMM_PMBEC

__all__ = ["IedbNetMHCcons", "IedbNetMHCpan", "IedbNetMHCIIpan", "IedbSMM", "IedbSMM_PMBEC"]
