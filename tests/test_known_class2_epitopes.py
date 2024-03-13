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
Make sure all class II binding predictors give a high IC50 and percentile rank.
"""
import mhctools

mhc_classes = [
    mhctools.NetMHCIIpan,
    mhctools.IedbNetMHCIIpan,
]

def expect_binder(mhc_model, peptide):
    prediction = mhc_model.predict(peptide)[0]
    if prediction.value:
        assert prediction.value < 500, "Expected %s to have IC50 < 500nM, got %s" % (
            peptide, prediction)
    else:
        assert prediction.percentile_rank < 5, "Expected %s to have percent rank < 3, got %s" % (
            peptide, prediction)

def test_Gag_epitope():
    # Test the DRB1*04:01 HIV epitope LERFAVNPGLLETSE from
    #   Dendritic Cell Mediated Delivery of Plasmid DNA Encoding LAMP/HIV-1 Gag
    #   Fusion Immunogen Enhances T Cell Epitope Responses in HLA DR4
    #   Transgenic Mice
    for mhc_class in mhc_classes:
        mhc_model = mhc_class("HLA-DRB1*04:01", default_peptide_lengths=15)
        expect_binder(mhc_model, "LERFAVNPGLLETSE")
