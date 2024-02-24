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

from mhctools import BindingPrediction, BindingPredictionCollection
from .common import eq_

def test_collection_to_dataframe():
    bp = BindingPrediction(
        peptide="SIINFEKL",
        allele="A0201",
        affinity=1.5,
        percentile_rank=0.1)
    collection = BindingPredictionCollection([bp])
    df = collection.to_dataframe()
    eq_(df.peptide.iloc[0], "SIINFEKL")
    eq_(df.affinity.iloc[0], 1.5)
    eq_(df.allele.iloc[0], "A0201")
    eq_(df.percentile_rank.iloc[0], 0.1)
