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

import pytest

from mhctools.pepsickle import Pepsickle
from mhctools.pred import Kind, PeptidePreds, COLUMNS

pepsickle = pytest.importorskip("pepsickle")


@pytest.fixture(scope="module")
def predictor():
    return Pepsickle()


@pytest.fixture(scope="module")
def predictor_invitro():
    return Pepsickle(model_type="in-vitro", proteasome_type="C")


PROTEIN = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVY"
PEPTIDES = ["SIINFEKL", "GILGFVFTL", "NLVPMVATV"]


def test_init_defaults():
    p = Pepsickle()
    assert p.model_type == "epitope"
    assert p.proteasome_type == "C"
    assert p.threshold == 0.5
    assert p.human_only is False
    assert p.alleles == []
    assert p.default_peptide_lengths == [9]


def test_init_invalid_model_type():
    with pytest.raises(ValueError, match="model_type"):
        Pepsickle(model_type="bad")


def test_init_invalid_proteasome_type():
    with pytest.raises(ValueError, match="proteasome_type"):
        Pepsickle(proteasome_type="X")


def test_str():
    p = Pepsickle()
    s = str(p)
    assert "Pepsickle" in s
    assert "epitope" in s


def test_predict_returns_peptide_preds(predictor):
    results = predictor.predict(PEPTIDES)
    assert isinstance(results, list)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert isinstance(pp, PeptidePreds)
        assert len(pp.preds) == 1
        pred = pp.preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert 0.0 <= pred.score <= 1.0
        assert pred.predictor_name == "pepsickle"


def test_predict_proteins_scores_vary(predictor):
    # Isolated peptides lack flanking context so C-terminal scores are edge-
    # affected.  predict_proteins on a full protein should produce varying
    # scores across sub-peptides.
    result = predictor.predict_proteins({"spike": PROTEIN})
    scores = [pp.preds[0].score for pp in result["spike"]]
    assert len(set(round(s, 4) for s in scores)) > 1


def test_predict_dataframe(predictor):
    df = predictor.predict_dataframe(PEPTIDES)
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == len(PEPTIDES)
    assert (df["kind"] == "proteasome_cleavage").all()
    assert (df["predictor_name"] == "pepsickle").all()


def test_predict_proteins(predictor):
    result = predictor.predict_proteins({"spike": PROTEIN})
    assert "spike" in result
    preds_list = result["spike"]
    assert isinstance(preds_list, list)
    assert all(isinstance(pp, PeptidePreds) for pp in preds_list)
    expected_count = len(PROTEIN) - 9 + 1
    assert len(preds_list) == expected_count
    for pp in preds_list:
        pred = pp.preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert pred.source_sequence_name == "spike"
        assert 0.0 <= pred.score <= 1.0


def test_predict_proteins_dataframe(predictor):
    df = predictor.predict_proteins_dataframe(
        {"spike": PROTEIN}, sample_name="sample1")
    assert list(df.columns) == list(COLUMNS)
    assert (df["source_sequence_name"] == "spike").all()
    assert (df["sample_name"] == "sample1").all()


def test_predict_proteins_string_input(predictor):
    result = predictor.predict_proteins(PROTEIN)
    assert "seq" in result


def test_predict_proteins_multiple_lengths(predictor):
    p = Pepsickle(default_peptide_lengths=[8, 9, 10])
    result = p.predict_proteins({"s": PROTEIN})
    preds_list = result["s"]
    expected = sum(len(PROTEIN) - k + 1 for k in [8, 9, 10])
    assert len(preds_list) == expected


def test_predict_cleavage_sites(predictor):
    result = predictor.predict_cleavage_sites({"spike": PROTEIN})
    assert "spike" in result
    sites = result["spike"]
    assert len(sites) == len(PROTEIN)
    for pos, residue, prob in sites:
        assert isinstance(pos, int)
        assert pos >= 1
        assert isinstance(residue, str)
        assert 0.0 <= prob <= 1.0


def test_predict_cleavage_sites_string(predictor):
    result = predictor.predict_cleavage_sites(PROTEIN)
    assert "seq" in result


@pytest.mark.xfail(
    reason="pepsickle in-vitro GB model requires older sklearn pickle format",
    raises=ModuleNotFoundError,
)
def test_invitro_model(predictor_invitro):
    results = predictor_invitro.predict(PEPTIDES)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert pp.preds[0].kind == Kind.proteasome_cleavage
        assert 0.0 <= pp.preds[0].score <= 1.0


def test_default_pred_kind():
    p = Pepsickle()
    assert p._default_pred_kind() == Kind.proteasome_cleavage
