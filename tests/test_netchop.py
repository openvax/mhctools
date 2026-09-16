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
import subprocess
from numpy import testing
from mhctools import NetChop
from mhctools.netchop import NETCHOP_CONTAINER_IMAGE
from mhctools.proteasome_predictor import ProteasomePredictor
from mhctools.pred import Kind, PeptideResult

# Peptides from http://tools.iedb.org/netchop/example/
peptides = """
MSLLTEVETPIRNEWGCRCNDSSDPLVVAASIIGIVHLILWIIDRLFSKSIYRIFKHGLKRGPSTEGVPESMREEYREEQQNAVDADDGHFVSIELE
MDSHTVSSFQVDCFLWHVRKQVADQDLGDAPFLDRLRRDQKSLKGRGSTLGLNIETATCVGKQIVERILKEESDEAFKMTMASALASRYLTDMTIEEMSRDWFMLMPKQKVAGPLCVRMDQAIMDKNIILKANFSVIFDRLENLTLLRAFTEEGAIVGEISPLPSLPGHTNEDVKNAIGVLIGGLEWNDNTVRVSETLQRFTWRSSNETGGPPFTPTQKRKMAGTIRSEV
MDSHTVSSFQDILMRMSKMQLGSSSGDLNGMITQFESLKLYRDSLGEAVMRLGDLHSLQHRNGKWREQLGQKFEEIRWLIEEVRHKLKTTENSFEQITFMQALQLLFEVEQEIRTFSFQLI
""".strip().split()


def test_is_proteasome_predictor():
    assert issubclass(NetChop, ProteasomePredictor)


def _netchop_or_skip():
    try:
        predictor = NetChop(model_variant=0)
    except FileNotFoundError as error:
        pytest.skip(str(error))
    if predictor.execution == "container":
        try:
            inspect = subprocess.run(
                ["docker", "image", "inspect", NETCHOP_CONTAINER_IMAGE],
                capture_output=True,
                timeout=30,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired) as error:
            pytest.skip("NetChop container runtime unavailable: %s" % error)
        if inspect.returncode:
            diagnostic = inspect.stderr.decode(errors="replace").strip()
            pytest.skip(
                "Pinned NetChop compatibility image is unavailable: %s"
                % diagnostic[-500:])
    return predictor


def test_cleavage_probs():
    obj = _netchop_or_skip()
    for pep in peptides:
        probs = obj.cleavage_probs(pep)
        assert len(probs) == len(pep)

    # Spot-check values from http://tools.iedb.org/netchop (12/19/2016)
    probs0 = obj.cleavage_probs(peptides[0])
    probs1 = obj.cleavage_probs(peptides[1])
    probs2 = obj.cleavage_probs(peptides[2])
    testing.assert_almost_equal(probs0[95], 0.976629)
    testing.assert_almost_equal(probs0[22], 0.022000)
    testing.assert_almost_equal(probs1[146], 0.977417)
    testing.assert_almost_equal(probs1[84], 0.285210)
    testing.assert_almost_equal(probs2[0], 0.547588)
    testing.assert_almost_equal(probs2[84], 0.104684)


def test_predict_proteins():
    obj = _netchop_or_skip()
    result = obj.predict_proteins({"pep0": peptides[0]})
    assert "pep0" in result
    pp_list = result["pep0"]
    assert isinstance(pp_list, list)
    assert all(isinstance(pp, PeptideResult) for pp in pp_list)
    for pp in pp_list:
        pred = pp.preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert pred.source_sequence_name == "pep0"
        assert pred.predictor_name == "netchop"
        assert 0.0 <= pred.score <= 1.0
