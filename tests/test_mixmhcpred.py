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

from tempfile import NamedTemporaryFile
from mhctools.mixmhcpred import parse_mixmhcpred_results
from .common import eq_

example_output =  """Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_A0201\t%Rank_A0201
MLDDFSAGA\t0.182093\tA0201\t0.3\t0.182093\t0.3
SPEGEETII\t-0.655341\tA0201\t51.0\t-0.655341\t51.0
ILDRIITNA\t0.203906\tA0201\t0.3\t0.203906\t0.3"""

def test_parse_mixmhcpred_results():
    with NamedTemporaryFile(mode="r+") as f:
        f.write(example_output)
        f.flush()
        binding_results = parse_mixmhcpred_results(f.name)

        eq_(len(binding_results), 3)
        eq_(binding_results[0].peptide, "MLDDFSAGA")
        eq_(binding_results[1].peptide, "SPEGEETII")
        eq_(binding_results[2].peptide, "ILDRIITNA")
