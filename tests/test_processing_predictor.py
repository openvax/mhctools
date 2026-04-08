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

from mhctools.processing_predictor import (
    ProcessingPredictor,
    score_cterm,
    score_nterm_cterm,
    score_cterm_anti_max_internal,
    score_cterm_anti_mean_internal,
    score_nterm_cterm_anti_max_internal,
    score_nterm_cterm_anti_mean_internal,
    _geomean,
)
from mhctools.proteasome_predictor import ProteasomePredictor
from mhctools.pred import Kind, PeptidePreds, COLUMNS


class StubPredictor(ProcessingPredictor):
    """Predictor with fixed per-position cleavage probs for testing."""

    def __init__(self, probs_map, **kwargs):
        super().__init__(**kwargs)
        self.probs_map = probs_map

    def cleavage_probs(self, sequence):
        if sequence in self.probs_map:
            return self.probs_map[sequence]
        return [0.5] * len(sequence)


class StubProteasome(ProteasomePredictor):
    """ProteasomePredictor with fixed cleavage probs."""

    def __init__(self, probs_map, **kwargs):
        super().__init__(**kwargs)
        self.probs_map = probs_map

    def cleavage_probs(self, sequence):
        if sequence in self.probs_map:
            return self.probs_map[sequence]
        return [0.5] * len(sequence)


# A 12-residue "protein" with known cleavage probs
PROTEIN = "ABCDEFGHIJKL"
#          0    1    2    3    4    5    6    7    8    9   10   11
PROBS = [0.1, 0.2, 0.8, 0.3, 0.05, 0.4, 0.1, 0.9, 0.2, 0.6, 0.3, 0.7]


# -- init --

def test_init_defaults():
    p = StubPredictor({})
    assert p.default_peptide_lengths == [9]
    assert p.scoring is score_nterm_cterm_anti_max_internal
    assert not hasattr(p, "alleles")


def test_init_custom_scoring():
    p = StubPredictor({}, scoring=score_cterm)
    assert p.scoring is score_cterm


def test_init_lambda_scoring():
    fn = lambda c, n, i: c
    p = StubPredictor({}, scoring=fn)
    assert p.scoring is fn


def test_init_non_callable_scoring_raises():
    with pytest.raises(TypeError, match="callable"):
        StubPredictor({}, scoring="not_a_function")


# -- component helpers --

def test_c_term_prob():
    assert ProcessingPredictor.c_term_prob(PROBS, 2, 4) == 0.4


def test_n_term_prob():
    assert ProcessingPredictor.n_term_prob(PROBS, 2, 4) == 0.2


def test_n_term_prob_offset_zero():
    assert ProcessingPredictor.n_term_prob(PROBS, 0, 4) is None


def test_internal_probs():
    assert ProcessingPredictor.internal_probs(PROBS, 2, 4) == [0.8, 0.3, 0.05]


def test_max_internal_prob():
    assert ProcessingPredictor.max_internal_prob(PROBS, 2, 4) == 0.8


def test_mean_internal_prob():
    expected = (0.8 + 0.3 + 0.05) / 3
    assert ProcessingPredictor.mean_internal_prob(PROBS, 2, 4) == pytest.approx(expected)


def test_max_internal_empty():
    assert ProcessingPredictor.max_internal_prob(PROBS, 0, 1) == 0.0


def test_mean_internal_empty():
    assert ProcessingPredictor.mean_internal_prob(PROBS, 0, 1) == 0.0


# -- scoring functions --

class TestScoringFunctions:

    def test_score_cterm(self):
        assert score_cterm(0.4, 0.2, [0.8, 0.3, 0.05]) == 0.4

    def test_score_nterm_cterm(self):
        expected = _geomean(0.4, 0.2)
        assert score_nterm_cterm(0.4, 0.2, [0.8]) == pytest.approx(expected)

    def test_score_nterm_cterm_no_nterm(self):
        assert score_nterm_cterm(0.4, None, [0.8]) == 0.4

    def test_score_cterm_anti_max_internal(self):
        assert score_cterm_anti_max_internal(0.4, 0.2, [0.8, 0.3, 0.05]) == \
            pytest.approx(0.4 * (1.0 - 0.8))

    def test_score_cterm_anti_mean_internal(self):
        mean_i = (0.8 + 0.3 + 0.05) / 3
        assert score_cterm_anti_mean_internal(0.4, 0.2, [0.8, 0.3, 0.05]) == \
            pytest.approx(0.4 * (1.0 - mean_i))

    def test_score_nterm_cterm_anti_max_internal(self):
        expected = _geomean(0.4, 0.2, 1.0 - 0.8)
        assert score_nterm_cterm_anti_max_internal(0.4, 0.2, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_score_nterm_cterm_anti_max_internal_no_nterm(self):
        expected = _geomean(0.4, 1.0 - 0.8)
        assert score_nterm_cterm_anti_max_internal(0.4, None, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_score_nterm_cterm_anti_mean_internal(self):
        mean_i = (0.8 + 0.3) / 2
        expected = _geomean(0.4, 0.2, 1.0 - mean_i)
        assert score_nterm_cterm_anti_mean_internal(0.4, 0.2, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_empty_internal(self):
        # length-1 peptide: no internal positions
        assert score_cterm_anti_max_internal(0.4, 0.2, []) == \
            pytest.approx(0.4 * 1.0)

    def test_zero_cterm_propagates(self):
        assert score_cterm_anti_max_internal(0.0, 0.5, [0.3]) == 0.0


# -- _peptide_score integration --

class TestPeptideScore:

    def _make(self, scoring):
        return StubPredictor({PROTEIN: PROBS}, scoring=scoring)

    def test_cterm(self):
        p = self._make(score_cterm)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(0.4)

    def test_nterm_cterm(self):
        p = self._make(score_nterm_cterm)
        expected = _geomean(0.4, 0.2)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_nterm_cterm_offset_zero(self):
        p = self._make(score_nterm_cterm)
        # offset=0 → n_term is None → falls back to cterm
        assert p._peptide_score(PROBS, 0, 4) == pytest.approx(PROBS[3])

    def test_custom_scoring(self):
        def my_scoring(c, n, internal):
            return c * 2
        p = self._make(my_scoring)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(0.8)


# -- predict with flanking --

def test_predict_no_flanks():
    stub = StubPredictor({})
    results = stub.predict(["ABCD", "EFGH"])
    assert len(results) == 2
    for pp in results:
        assert isinstance(pp, PeptidePreds)
        assert pp.preds[0].kind == Kind.antigen_processing
        assert pp.preds[0].n_flank == ""
        assert pp.preds[0].c_flank == ""


def test_predict_with_flanks():
    # With flanking, the model sees the concatenated sequence
    full_seq = "XXABCDYY"
    probs_full = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    stub = StubPredictor({full_seq: probs_full}, scoring=score_cterm)
    results = stub.predict(["ABCD"], n_flanks=["XX"], c_flanks=["YY"])
    assert len(results) == 1
    pred = results[0].preds[0]
    assert pred.n_flank == "XX"
    assert pred.c_flank == "YY"
    # c_term for peptide at offset=2, length=4 → probs_full[5] = 0.6
    assert pred.score == pytest.approx(0.6)


def test_predict_dataframe_with_flanks():
    stub = StubPredictor({})
    df = stub.predict_dataframe(
        ["ABCD"], n_flanks=["XX"], c_flanks=["YY"], sample_name="s1")
    assert list(df.columns) == list(COLUMNS)
    assert df.iloc[0]["n_flank"] == "XX"
    assert df.iloc[0]["c_flank"] == "YY"


# -- predict_proteins with flank_length --

def test_predict_proteins():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins({"p": PROTEIN})
    assert "p" in result
    pp_list = result["p"]
    expected_count = len(PROTEIN) - 4 + 1
    assert len(pp_list) == expected_count


def test_predict_proteins_string():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins(PROTEIN)
    assert "seq" in result


def test_predict_proteins_flank_length():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins({"p": PROTEIN}, flank_length=2)
    pp_list = result["p"]
    # First peptide at offset=0: n_flank="" (nothing before), c_flank="EF"
    first = pp_list[0].preds[0]
    assert first.n_flank == ""
    assert first.c_flank == "EF"
    # Second peptide at offset=1: n_flank="A", c_flank="FG"
    second = pp_list[1].preds[0]
    assert second.n_flank == "A"
    assert second.c_flank == "FG"
    # Last peptide: c_flank truncated at protein end
    last = pp_list[-1].preds[0]
    assert last.c_flank == ""
    assert last.n_flank == "GH"


def test_predict_proteins_no_flank_length():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins({"p": PROTEIN})
    for pp in result["p"]:
        assert pp.preds[0].n_flank == ""
        assert pp.preds[0].c_flank == ""


def test_predict_proteins_dataframe():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    df = stub.predict_proteins_dataframe({"p": PROTEIN}, sample_name="s1")
    assert list(df.columns) == list(COLUMNS)
    assert (df["source_sequence_name"] == "p").all()


# -- predict_cleavage_sites --

def test_predict_cleavage_sites():
    stub = StubPredictor({PROTEIN: PROBS})
    result = stub.predict_cleavage_sites({"p": PROTEIN})
    assert result["p"] == PROBS


# -- ProteasomePredictor --

def test_proteasome_default_scoring():
    p = StubProteasome({PROTEIN: PROBS})
    assert p.scoring is score_cterm_anti_max_internal


def test_proteasome_pred_kind():
    p = StubProteasome({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = p.predict_proteins({"p": PROTEIN})
    for pp in result["p"]:
        assert pp.preds[0].kind == Kind.proteasome_cleavage


def test_proteasome_custom_scoring():
    p = StubProteasome({PROTEIN: PROBS}, scoring=score_cterm)
    assert p.scoring is score_cterm
    result = p.predict_proteins({"p": PROTEIN}, peptide_lengths=[4])
    first = result["p"][0].preds[0]
    # offset=0, length=4, cterm = PROBS[3] = 0.3
    assert first.score == pytest.approx(0.3)
