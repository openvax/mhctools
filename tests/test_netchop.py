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

import pytest 
from numpy import testing
from mhctools import NetChop
from .arch import apple_silicon

# Peptides from http://tools.iedb.org/netchop/example/
peptides = """
MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGIVHLILWIIDRLFSKSIYRIFKHGLKRGPSTEGVPESMREEYREEQQNAVDADDGHFVSIELE
MDSHTVSSFQVDCFLWHVRKQVADQDLGDAPFLDRLRRDQKSLKGRGSTLGLNIETATCVGKQIVERILKEESDEAFKMTMASALASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCVRMDQAIMDKNIILKANFSVIFDRLENLTLLRAFTEEGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSETLQRFTWRSSNETGGPPFTPTQKRKMAGTIRSEV
MDSHTVSSFQDILMRMSKMQLGSSSGDLNGMITQFESLKLYRDSLGEAVMRLGDLHSLQHRNGKWREQLGQKFEEIRWLIEEVRHKLKTTENSFEQITFMQALQLLFEVEQEIRTFSFQLI
""".strip().split()


@pytest.mark.skipif(apple_silicon, reason="Can't run netChop on arm64 architecture")
def test_simple():
    obj = NetChop()
    result = obj.predict(peptides)
    assert len(result) == 3
    assert len(result[0]) == len(peptides[0])
    assert len(result[1]) == len(peptides[1])
    assert len(result[2]) == len(peptides[2])

    # These numbers are from running http://tools.iedb.org/netchop
    # via the web interface on 12/19/2016.
    testing.assert_almost_equal(result[0][95], 0.976629)
    testing.assert_almost_equal(result[0][22], 0.022000)
    testing.assert_almost_equal(result[1][146], 0.977417)
    testing.assert_almost_equal(result[1][84], 0.285210)
    testing.assert_almost_equal(result[2][0], 0.547588)
    testing.assert_almost_equal(result[2][84], 0.104684)
