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

import sys

from .common  import eq_
from numpy import testing

from mhcflurry import Class1AffinityPredictor
from mhctools import MHCflurry

DEFAULT_ALLELE = "HLA-A*02:01"

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}


def test_mhcflurry():
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    eq_(4, len(binding_predictions),
        "Expected 4 binding predictions from %s" % (binding_predictions,))

    prediction_scores = {
        (x.peptide, x.allele): x.affinity for x in binding_predictions
    }

    predictor = Class1AffinityPredictor.load()
    # test one prediction at a time to make sure there's no peptide/allele mixup
    for (peptide, allele), affinity in prediction_scores.items():
        prediction = predictor.predict([peptide], allele=allele)
        assert len(prediction) == 1
        # we've seen results differ a bit so doing an approximate check, not an error condition
        testing.assert_almost_equal(prediction[0], affinity, decimal=0)
