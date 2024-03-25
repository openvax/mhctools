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
from mhctools import IedbNetMHCpan
from .common import assert_raises


DEFAULT_ALLELE = 'HLA-A*02:01'
UNSUPPORTED_ALLELE = 'HLA-A*24:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

@pytest.mark.xfail(reason="IEDB server giving 403 errors from GitHub actions runners")
def test_netmhcpan_iedb():
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[8,9])
    assert len(binding_predictions) == 10, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

@pytest.mark.xfail(reason="IEDB server giving 403 errors from GitHub actions runners")
def test_netmhcpan_iedb_unsupported_allele():
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE, UNSUPPORTED_ALLELE], raise_on_error=False)
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

    # check that the error is raised when raise_on_error is left at default (True)
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE, UNSUPPORTED_ALLELE])
    with assert_raises(ValueError):
        binding_predictions = predictor.predict_subsequences(
            protein_sequence_dict,
            peptide_lengths=[9])
