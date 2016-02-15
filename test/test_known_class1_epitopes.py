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
Make sure all class I binding predictors give a high IC50 and percentile rank.
"""
import mhctools

mhc_classes = [
    mhctools.NetMHCcons,
    mhctools.NetMHCpan,
    mhctools.NetMHCcons,
    mhctools.NetMHC,
    mhctools.IedbNetMHCcons,
    mhctools.IedbNetMHCpan,
    mhctools.IedbSMM,
    mhctools.IedbSMM_PMBEC,
    mhctools.MHCFlurry,
]

def expect_binder(mhc_model, peptide):
    prediction = mhc_model.predict(peptide)[0]
    assert prediction.value < 500, "Expected %s to have IC50 < 500nM, got %s" % (
        peptide, prediction)

def test_MAGE_epitope():
    # Test the A1 MAGE epitope ESDPIVAQY from
    #   Identification of a Titin-Derived HLA-A1-Presented Peptide
    #   as a Cross-Reactive Target for Engineered MAGE A3-Directed
    #   T Cells
    for mhc_class in mhc_classes:
        mhc_model = mhc_class("HLA-A*01:01", epitope_lengths=9)
        yield (expect_binder, mhc_model, "ESDPIVAQY")

def test_HIV_epitope():
    # Test the A2 HIV epitope SLYNTVATL from
    #    The HIV-1 HLA-A2-SLYNTVATL Is a Help-Independent CTL Epitope
    for mhc_class in mhc_classes:
        mhc_model = mhc_class("HLA-A*02:01", epitope_lengths=9)
        yield (expect_binder, mhc_model, "SLYNTVATL")
