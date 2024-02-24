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

from mhctools.binding_prediction import BindingPrediction
from .common import eq_

def test_binding_prediction_fields():
    bp = BindingPrediction(
            source_sequence_name="seq",
            offset=0,
            peptide="SIINFEKL",
            allele="H-2-K-d",
            affinity=200.0,
            percentile_rank=0.3)
    eq_(bp.source_sequence_name, "seq")
    eq_(bp.offset, 0)
    eq_(bp.peptide, "SIINFEKL")
    eq_(bp.allele, "H-2-K-d")
    eq_(bp.affinity, 200.0)
    eq_(bp.percentile_rank, 0.3)

def test_binding_prediction_str_repr():
    bp = BindingPrediction(
            source_sequence_name="seq",
            offset=0,
            peptide="SIINFEKL",
            allele="H-2-K-d",
            affinity=200.0,
            percentile_rank=0.3)
    eq_(str(bp), repr(bp))
    assert "SIINFEKL" in str(bp)
    assert "200.0" in str(bp)

def test_binding_predicton_eq():
    bp1 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKL",
        allele="H-2-K-d",
        affinity=200.0,
        percentile_rank=0.3)
    bp2 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKLY",
        allele="H-2-K-d",
        affinity=5000.0,
        percentile_rank=9.3)
    eq_(bp1, bp1)
    eq_(bp2, bp2)
    assert bp1 != bp2

def test_binding_predicton_hash():
    bp1 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKL",
        allele="H-2-K-d",
        affinity=200.0,
        percentile_rank=0.3)
    bp2 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKLY",
        allele="H-2-K-d",
        affinity=5000.0,
        percentile_rank=9.3)
    eq_(hash(bp1), hash(bp1))
    assert hash(bp1) != hash(bp2)
