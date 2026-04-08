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


# -- init --

def test_init_defaults():
    p = Pepsickle()
    assert p.model_type == "epitope"
    assert p.proteasome_type == "C"
    assert p.threshold == 0.5
    assert p.human_only is False
    assert p.default_peptide_lengths == [9]
    assert p.scoring == "nterm_cterm_max_internal"
    assert not hasattr(p, "alleles")


def test_init_invalid_model_type():
    with pytest.raises(ValueError, match="model_type"):
        Pepsickle(model_type="bad")


def test_init_invalid_proteasome_type():
    with pytest.raises(ValueError, match="proteasome_type"):
        Pepsickle(proteasome_type="X")


def test_init_invalid_scoring():
    with pytest.raises(ValueError, match="scoring"):
        Pepsickle(scoring="bad")


def test_str():
    p = Pepsickle()
    s = str(p)
    assert "Pepsickle" in s
    assert "epitope" in s


# -- cleavage_probs --

def test_cleavage_probs(predictor):
    probs = predictor.cleavage_probs(PROTEIN)
    assert len(probs) == len(PROTEIN)
    assert all(0.0 <= p <= 1.0 for p in probs)


# -- predict --

def test_predict_returns_peptide_preds(predictor):
    results = predictor.predict(PEPTIDES)
    assert isinstance(results, list)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert isinstance(pp, PeptidePreds)
        assert len(pp.preds) == 1
        pred = pp.preds[0]
        assert pred.kind == Kind.antigen_processing
        assert 0.0 <= pred.score <= 1.0
        assert pred.predictor_name == "pepsickle"


def test_predict_dataframe(predictor):
    df = predictor.predict_dataframe(PEPTIDES)
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == len(PEPTIDES)
    assert (df["kind"] == "antigen_processing").all()
    assert (df["predictor_name"] == "pepsickle").all()


# -- predict_proteins --

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
        assert pred.kind == Kind.antigen_processing
        assert pred.source_sequence_name == "spike"
        assert 0.0 <= pred.score <= 1.0


def test_predict_proteins_scores_vary(predictor):
    result = predictor.predict_proteins({"spike": PROTEIN})
    scores = [pp.preds[0].score for pp in result["spike"]]
    assert len(set(round(s, 4) for s in scores)) > 1


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


# -- predict_cleavage_sites --

def test_predict_cleavage_sites(predictor):
    result = predictor.predict_cleavage_sites({"spike": PROTEIN})
    assert "spike" in result
    probs = result["spike"]
    assert len(probs) == len(PROTEIN)
    assert all(0.0 <= p <= 1.0 for p in probs)


def test_predict_cleavage_sites_string(predictor):
    result = predictor.predict_cleavage_sites(PROTEIN)
    assert "seq" in result


# -- scoring methods --

def test_scoring_cterm_only():
    p = Pepsickle(scoring="cterm")
    result = p.predict_proteins({"s": PROTEIN})
    probs = p.cleavage_probs(PROTEIN)
    pp = result["s"][0]
    # First 9-mer: offset=0, c_term = probs[8]
    assert pp.preds[0].score == pytest.approx(probs[8])


def test_scoring_methods_produce_different_scores():
    methods = ["cterm", "nterm_cterm", "cterm_max_internal",
               "nterm_cterm_max_internal"]
    scores_by_method = {}
    for method in methods:
        p = Pepsickle(scoring=method)
        result = p.predict_proteins({"s": PROTEIN})
        # Use a peptide with offset > 0 so n_term is available
        scores = [pp.preds[0].score for pp in result["s"]]
        scores_by_method[method] = scores
    # At least some methods should give different results
    unique_score_lists = set(
        tuple(round(s, 6) for s in v) for v in scores_by_method.values()
    )
    assert len(unique_score_lists) > 1


# -- in-vitro model --

@pytest.mark.xfail(
    reason="pepsickle in-vitro GB model requires older sklearn pickle format",
    raises=ModuleNotFoundError,
)
def test_invitro_model(predictor_invitro):
    results = predictor_invitro.predict(PEPTIDES)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert pp.preds[0].kind == Kind.antigen_processing
        assert 0.0 <= pp.preds[0].score <= 1.0
