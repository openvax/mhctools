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
Regression tests for Class II MHC allele name formats found in the wild.

Covers the formats listed in issue #31 and other common variants.
"""

from mhctools.allele_normalization import normalize_allele_name


# -- Formats from issue #31 --

def test_dra1_drb1_slash_no_colons():
    """HLA-DRA1*0101/HLA-DRB1*0401 (PMID 25342727)"""
    assert normalize_allele_name("HLA-DRA1*0101/HLA-DRB1*0401") == \
        "HLA-DRA1*01:01-DRB1*04:01"


def test_dra_drb1_slash_with_colons():
    """HLA-DRA*01:01/DRB1*03:01 (IEDB MHCalleleId/631)"""
    assert normalize_allele_name("HLA-DRA*01:01/DRB1*03:01") == \
        "HLA-DRA1*01:01-DRB1*03:01"


def test_dqa1_drb1_semicolon():
    """HLA-DQA1*01:02;DRB1*15:01 (JBC 2015)"""
    assert normalize_allele_name("HLA-DQA1*01:02;DRB1*15:01") == \
        "HLA-DQA1*01:02-DRB1*15:01"


# -- Other common Class II formats --

def test_drb1_beta_only_infers_alpha():
    """Beta-only DRB1 should infer DRA1*01:01 alpha."""
    assert normalize_allele_name("HLA-DRB1*07:01") == \
        "HLA-DRA1*01:01-DRB1*07:01"


def test_dpa1_dpb1_slash():
    assert normalize_allele_name("HLA-DPA1*01:03/DPB1*04:01") == \
        "HLA-DPA1*01:03-DPB1*04:01"


def test_dqa1_dqb1_slash():
    assert normalize_allele_name("HLA-DQA1*05:01/DQB1*02:01") == \
        "HLA-DQA1*05:01-DQB1*02:01"


def test_drb1_no_colons():
    """Four-digit format without colons."""
    assert normalize_allele_name("HLA-DRB1*0401") == \
        "HLA-DRA1*01:01-DRB1*04:01"


# -- Mouse Class II --

def test_mouse_ia():
    assert normalize_allele_name("H-2-IAb") == "H-2-IAb"


def test_mouse_ia_no_dash():
    assert normalize_allele_name("H2-IAd") == "H-2-IAd"


def test_mouse_ie():
    assert normalize_allele_name("H-2-IEk") == "H-2-IEk"


# -- Class I sanity (no regressions) --

def test_class1_standard():
    assert normalize_allele_name("HLA-A*02:01") == "HLA-A*02:01"


def test_class1_compact():
    assert normalize_allele_name("A0201") == "HLA-A*02:01"


def test_class1_mouse():
    assert normalize_allele_name("H-2-Kb") == "H-2-Kb"
