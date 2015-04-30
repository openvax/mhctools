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
Make sure all binding predictors give a high IC50 and percentile rank.
"""
import mhctools

mhc_classes = [
    mhctools.NetMHCpan,
    mhctools.IedbNetMHCcons,
    mhctools.IedbNetMHCpan,
    mhctools.IedbSMM,
    mhctools.IedbSMM_PMBEC,
]

def test_MAGE_epitope():
    """
    Test the A1 MAGE epitope from
        Identification of a Titin-Derived HLA-A1-Presented Peptide
        as a Cross-Reactive Target for Engineered MAGE A3-Directed
        T Cells
    """
    for mhc_class in mhc_classes:
        mhc_model = mhc_class("HLA-A*01:01", epitope_lengths=9)
        epitope = mhc_model.predict("ESDPIVAQY")[0]
        assert epitope.value < 500, "Expected %s to have IC50 < 500nM" % epitope
        assert epitope.percentile_rank < 1.0, \
            "Expected %s to have percentile rank <= 1.0" % epitope