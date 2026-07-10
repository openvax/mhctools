# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the Calis (IEDB) class-I immunogenicity predictor.

Fully binary-free — the model is self-contained. Reference scores below were
produced by the canonical IEDB ``predict_immunogenicity.py`` algorithm (Calis
et al. 2013) with the default P1/P2/C-terminus mask.
"""

import pytest

from mhctools import Calis, Kind
from mhctools.calis import (
    IMMUNOSCALE,
    IMMUNOWEIGHT,
    immunogenicity_score,
    position_weights,
)


# peptide -> canonical IEDB score (default anchor mask, rounded to 5 dp).
_REFERENCE = {
    "SIINFEKL": 0.06294,     # 8-mer (H-2Kb)
    "GILGFVFTL": 0.30484,    # 9-mer (influenza M1)
    "ELAGIGILTV": 0.348,     # 10-mer (MART-1)
    "NLVPMVATV": -0.0742,    # 9-mer (CMV pp65)
    "LLFGYPVYV": 0.09074,    # 9-mer (HTLV-1 Tax)
    "FLRGRAYGL": 0.15481,    # 9-mer (EBV)
}


# --- parameters ------------------------------------------------------------

def test_scale_and_weights_shape():
    assert len(IMMUNOSCALE) == 20
    # W most immunogenic, K least — the paper's headline residues.
    assert max(IMMUNOSCALE, key=IMMUNOSCALE.get) == "W"
    assert min(IMMUNOSCALE, key=IMMUNOSCALE.get) == "K"
    assert IMMUNOWEIGHT == [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]


def test_position_weights_stretch_for_long_peptides():
    assert position_weights(9) == IMMUNOWEIGHT
    # 11-mer: two extra 0.30 interior weights inserted after P5.
    w11 = position_weights(11)
    assert len(w11) == 11
    assert w11 == [0.00, 0.00, 0.10, 0.31, 0.30, 0.30, 0.30, 0.29, 0.26, 0.18, 0.00]


# --- scoring against the canonical IEDB output -----------------------------

@pytest.mark.parametrize("peptide,expected", sorted(_REFERENCE.items()))
def test_immunogenicity_score_matches_iedb(peptide, expected):
    assert immunogenicity_score(peptide) == pytest.approx(expected, abs=1e-5)


def test_case_insensitive():
    assert immunogenicity_score("gilgfvftl") == immunogenicity_score("GILGFVFTL")


# --- predictor API ---------------------------------------------------------

def test_predict_shape_and_kind():
    predictor = Calis()
    peptides = ["GILGFVFTL", "NLVPMVATV"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.immunogenicity
        assert pred.allele == ""
        assert pred.predictor_name == "calis"
        assert pred.score == pytest.approx(_REFERENCE[peptide], abs=1e-5)
        # convenience accessor resolves to the same prediction
        assert result.immunogenicity is pred


def test_predict_accepts_single_string():
    (result,) = Calis().predict("GILGFVFTL")
    assert result.preds[0].score == pytest.approx(0.30484, abs=1e-5)


def test_default_kind_and_support():
    predictor = Calis()
    assert predictor._default_pred_kind() == Kind.immunogenicity
    support = predictor.kind_support()[Kind.immunogenicity]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "I"
    assert predictor.supported_kinds == (Kind.immunogenicity,)


def test_predict_dataframe():
    df = Calis().predict_dataframe(["GILGFVFTL"], sample_name="s1")
    assert len(df) == 1
    row = df.iloc[0]
    assert row["peptide"] == "GILGFVFTL"
    assert row["kind"] == Kind.immunogenicity
    assert row["sample_name"] == "s1"
    assert row["score"] == pytest.approx(0.30484, abs=1e-5)


# --- validation ------------------------------------------------------------

def test_rejects_out_of_range_length():
    predictor = Calis()  # default 8-11
    with pytest.raises(ValueError, match="length 8-11"):
        predictor.predict(["SIINFEK"])       # 7-mer, too short
    with pytest.raises(ValueError, match="length 8-11"):
        predictor.predict(["A" * 12])         # 12-mer, too long


def test_rejects_non_standard_amino_acids():
    with pytest.raises(ValueError, match="non-standard amino acid"):
        Calis().predict(["GILGFVFTX"])


def test_custom_length_bounds():
    # A 10/11-mer window opened up explicitly still scores like IEDB.
    predictor = Calis(min_peptide_length=10, max_peptide_length=10)
    (result,) = predictor.predict(["ELAGIGILTV"])
    assert result.preds[0].score == pytest.approx(0.348, abs=1e-5)
    with pytest.raises(ValueError):
        predictor.predict(["GILGFVFTL"])      # 9-mer now out of range
