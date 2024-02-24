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


"""
Class I predictors typically allow 8mer or longer peptides, whereas
NetMHCIIpan allows 9mer or longer. So far only MHCflurry has a max
length (of 15mer).
"""

from mhctools import NetMHCIIpan, NetMHC
from .common import eq_, assert_raises


def test_class2_9mer_success():
    ii_pan_predictor = NetMHCIIpan(alleles=["HLA-DRB1*01:01"])
    predictions = ii_pan_predictor.predict_peptides(["A" * 9])
    eq_(len(predictions), 1)

def test_class2_8mer_fails():
    ii_pan_predictor = NetMHCIIpan(alleles=["HLA-DRB1*01:01"])
    with assert_raises(ValueError):
        ii_pan_predictor.predict_peptides(["A" * 8])

def test_class1_8mer_success():
    netmhc = NetMHC(alleles=["HLA-A0201"])
    predictions = netmhc.predict_peptides(["A" * 8])
    eq_(len(predictions), 1)

def test_class1_7mer_failure():
    netmhc = NetMHC(alleles=["HLA-A0201"])
    with assert_raises(ValueError):
        netmhc.predict_peptides(["A" * 7])
