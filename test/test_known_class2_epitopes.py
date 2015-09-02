# Copyright (c) 2014. Mount Sinai School of Medicine
#
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
    assert prediction.value < 500, "Expected %s to have IC50 < 500nM, got %s" % (
        peptide, prediction)

def test_Gag171_epitope():
    # Test the DRB*07:01 HIV (Gag 171-185) epitope QGQMVHQAISPRTLN from
    #   Identification and antigenicity of broadly cross-reactive and conserved
    #   human immunodeficiency virus type 1-derived helper T-lymphocyte epitopes
    for mhc_class in mhc_classes:
        mhc_model = mhc_class("HLA-DRB*07:01", epitope_lengths=9)
        yield (expect_binder, mhc_model, "QGQMVHQAISPRTLN")
