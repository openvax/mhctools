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

from numpy.testing import assert_allclose
from mhctools import NetMHCstabpan


DEFAULT_ALLELE = 'HLA-A*02:01'

# Protein sequences to test with against known values of netMHCstabpan web server.
protein_sequences = [
    "ALDKNLHQL",
    "ALHEEVVCV",
    "ALPPTVYEV",
    "AVLGSFSYV",
    "EMASVLFKA",
]

web_server_predictions = [
    4.6393,
    10.6011,
    11.8201,
    4.6722,
    1.8404,
]

# REMINDER: a program called "netMHCstabpan" MUST be installed and working for
# this test suite to succeeed. Also all peptides must be the same length.


def test_netmhc_stabpan_accuracy():    
    # Check that the netMHCstabpan program is working and returning th eexpected outputs.
    predictor = NetMHCstabpan(
        alleles=[DEFAULT_ALLELE], program_name='netMHCstabpan')

    binding_predictions = predictor.predict_peptides(protein_sequences)
    stability_predictions = [p.score for p in binding_predictions]
    rank_predictions = [p.percentile_rank for p in binding_predictions]

    assert len(web_server_predictions) == len(binding_predictions)
    assert len(stability_predictions) == len(binding_predictions)

    for prank in rank_predictions:
        # Make sure that correct mapping is done by checking percentiles aren't above 100.
        assert prank < 100
    
    for i, (expected, actual) in enumerate(zip(web_server_predictions, stability_predictions)):
        # Check to make sure that the stability predictions are within 0.01 of the webserver values.
        # This could be the result of different versions of dependencies or the nature of the ANN itself.
        assert_allclose(expected, actual, atol=0.01, err_msg="Peptide %d: expected %f but got %f" % (i, expected, actual))

    