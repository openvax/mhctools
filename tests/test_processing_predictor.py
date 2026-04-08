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

"""
Comprehensive unit tests for the antigen-processing prediction framework:
ProcessingPredictor, ProteasomePredictor, scoring functions, component
helpers, flanking support, and edge cases.
"""

import math

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


# ======================================================================
# Test fixtures / stubs
# ======================================================================

class StubPredictor(ProcessingPredictor):
    """ProcessingPredictor with fixed per-position cleavage probs."""

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

# A second protein for multi-protein tests
PROTEIN2 = "MNOPQR"
PROBS2 = [0.3, 0.5, 0.7, 0.1, 0.9, 0.4]

ALL_SCORING_FNS = [
    score_cterm,
    score_nterm_cterm,
    score_cterm_anti_max_internal,
    score_cterm_anti_mean_internal,
    score_nterm_cterm_anti_max_internal,
    score_nterm_cterm_anti_mean_internal,
]


# ======================================================================
# _geomean
# ======================================================================

class TestGeomean:

    def test_single_value(self):
        assert _geomean(0.5) == pytest.approx(0.5)

    def test_two_values(self):
        assert _geomean(0.25, 1.0) == pytest.approx(math.sqrt(0.25))

    def test_three_values(self):
        assert _geomean(0.8, 0.2, 0.5) == pytest.approx(
            (0.8 * 0.2 * 0.5) ** (1.0 / 3))

    def test_zero_propagates(self):
        assert _geomean(0.0, 0.5, 0.9) == 0.0

    def test_all_ones(self):
        assert _geomean(1.0, 1.0, 1.0) == pytest.approx(1.0)


# ======================================================================
# ProcessingPredictor.__init__
# ======================================================================

class TestInit:

    def test_defaults(self):
        p = StubPredictor({})
        assert p.default_peptide_lengths == [9]
        assert p.scoring is score_nterm_cterm_anti_max_internal
        assert not hasattr(p, "alleles")

    def test_custom_peptide_lengths_list(self):
        p = StubPredictor({}, default_peptide_lengths=[8, 10, 11])
        assert p.default_peptide_lengths == [8, 10, 11]

    def test_peptide_lengths_int_coerced_to_list(self):
        p = StubPredictor({}, default_peptide_lengths=10)
        assert p.default_peptide_lengths == [10]

    def test_custom_scoring(self):
        p = StubPredictor({}, scoring=score_cterm)
        assert p.scoring is score_cterm

    def test_lambda_scoring(self):
        fn = lambda c, n, i: c
        p = StubPredictor({}, scoring=fn)
        assert p.scoring is fn

    def test_non_callable_scoring_raises(self):
        with pytest.raises(TypeError, match="callable"):
            StubPredictor({}, scoring="not_a_function")

    def test_non_callable_scoring_int_raises(self):
        with pytest.raises(TypeError, match="callable"):
            StubPredictor({}, scoring=42)


# ======================================================================
# __str__ / __repr__
# ======================================================================

class TestRepr:

    def test_str_contains_class_name(self):
        p = StubPredictor({})
        assert "StubPredictor" in str(p)

    def test_str_contains_scoring_name(self):
        p = StubPredictor({}, scoring=score_cterm)
        assert "score_cterm" in str(p)

    def test_str_lambda_scoring(self):
        p = StubPredictor({}, scoring=lambda c, n, i: c)
        s = str(p)
        assert "StubPredictor" in s
        # lambda has no __name__ in the usual sense, should still work
        assert "lambda" in s or "function" in s or "<" in s

    def test_repr_equals_str(self):
        p = StubPredictor({})
        assert repr(p) == str(p)


# ======================================================================
# Abstract method
# ======================================================================

def test_cleavage_probs_abstract():
    with pytest.raises(NotImplementedError, match="ProcessingPredictor"):
        ProcessingPredictor().cleavage_probs("ABC")


# ======================================================================
# _pred_kind / _predictor_name
# ======================================================================

def test_processing_predictor_pred_kind():
    p = StubPredictor({})
    assert p._pred_kind() == Kind.antigen_processing


def test_predictor_name_default():
    p = StubPredictor({})
    assert p._predictor_name() == "stubpredictor"


# ======================================================================
# Component helpers
# ======================================================================

class TestComponentHelpers:
    """Test c_term_prob, n_term_prob, internal_probs, max/mean variants."""

    # -- c_term_prob --

    def test_c_term_first_peptide(self):
        # offset=0, length=4 → probs[3]
        assert ProcessingPredictor.c_term_prob(PROBS, 0, 4) == 0.3

    def test_c_term_middle_peptide(self):
        # offset=2, length=4 → probs[5]
        assert ProcessingPredictor.c_term_prob(PROBS, 2, 4) == 0.4

    def test_c_term_last_peptide(self):
        # offset=8, length=4 → probs[11]
        assert ProcessingPredictor.c_term_prob(PROBS, 8, 4) == 0.7

    def test_c_term_length_1(self):
        assert ProcessingPredictor.c_term_prob(PROBS, 5, 1) == 0.4

    # -- n_term_prob --

    def test_n_term_offset_zero(self):
        assert ProcessingPredictor.n_term_prob(PROBS, 0, 4) is None

    def test_n_term_offset_one(self):
        # offset=1 → probs[0]
        assert ProcessingPredictor.n_term_prob(PROBS, 1, 4) == 0.1

    def test_n_term_middle(self):
        # offset=5 → probs[4]
        assert ProcessingPredictor.n_term_prob(PROBS, 5, 3) == 0.05

    # -- internal_probs --

    def test_internal_probs_normal(self):
        # offset=2, length=4 → probs[2:5] = [0.8, 0.3, 0.05]
        assert ProcessingPredictor.internal_probs(PROBS, 2, 4) == [0.8, 0.3, 0.05]

    def test_internal_probs_length_1(self):
        # length=1 → empty (no positions before c_term)
        assert ProcessingPredictor.internal_probs(PROBS, 3, 1) == []

    def test_internal_probs_length_2(self):
        # offset=3, length=2 → probs[3:4] = [0.3]
        assert ProcessingPredictor.internal_probs(PROBS, 3, 2) == [0.3]

    def test_internal_probs_long_peptide(self):
        # offset=0, length=6 → probs[0:5]
        assert ProcessingPredictor.internal_probs(PROBS, 0, 6) == \
            [0.1, 0.2, 0.8, 0.3, 0.05]

    # -- max_internal_prob --

    def test_max_internal_normal(self):
        assert ProcessingPredictor.max_internal_prob(PROBS, 2, 4) == 0.8

    def test_max_internal_empty(self):
        assert ProcessingPredictor.max_internal_prob(PROBS, 0, 1) == 0.0

    def test_max_internal_single(self):
        assert ProcessingPredictor.max_internal_prob(PROBS, 3, 2) == 0.3

    # -- mean_internal_prob --

    def test_mean_internal_normal(self):
        expected = (0.8 + 0.3 + 0.05) / 3
        assert ProcessingPredictor.mean_internal_prob(PROBS, 2, 4) == \
            pytest.approx(expected)

    def test_mean_internal_empty(self):
        assert ProcessingPredictor.mean_internal_prob(PROBS, 0, 1) == 0.0

    def test_mean_internal_single(self):
        assert ProcessingPredictor.mean_internal_prob(PROBS, 3, 2) == \
            pytest.approx(0.3)


# ======================================================================
# Built-in scoring functions
# ======================================================================

class TestScoringFunctions:
    """Test each scoring function's formula directly."""

    # -- score_cterm --

    def test_cterm_ignores_nterm_and_internal(self):
        assert score_cterm(0.7, 0.3, [0.5, 0.6]) == 0.7

    def test_cterm_zero(self):
        assert score_cterm(0.0, 0.9, [0.1]) == 0.0

    def test_cterm_one(self):
        assert score_cterm(1.0, 0.0, []) == 1.0

    # -- score_nterm_cterm --

    def test_nterm_cterm_both_present(self):
        expected = _geomean(0.4, 0.2)
        assert score_nterm_cterm(0.4, 0.2, [0.8]) == pytest.approx(expected)

    def test_nterm_cterm_none_nterm(self):
        assert score_nterm_cterm(0.4, None, [0.8]) == 0.4

    def test_nterm_cterm_equal(self):
        # geomean(x, x) == x
        assert score_nterm_cterm(0.6, 0.6, []) == pytest.approx(0.6)

    # -- score_cterm_anti_max_internal --

    def test_cterm_anti_max_internal(self):
        assert score_cterm_anti_max_internal(0.4, 0.2, [0.8, 0.3, 0.05]) == \
            pytest.approx(0.4 * 0.2)

    def test_cterm_anti_max_internal_empty(self):
        assert score_cterm_anti_max_internal(0.4, 0.2, []) == \
            pytest.approx(0.4)

    def test_cterm_anti_max_internal_single(self):
        assert score_cterm_anti_max_internal(0.5, None, [0.3]) == \
            pytest.approx(0.5 * 0.7)

    def test_cterm_anti_max_internal_all_one_internal(self):
        # max internal = 1.0 → anti = 0 → score = 0
        assert score_cterm_anti_max_internal(0.9, None, [1.0, 0.5]) == 0.0

    # -- score_cterm_anti_mean_internal --

    def test_cterm_anti_mean_internal(self):
        mean_i = (0.8 + 0.3 + 0.05) / 3
        assert score_cterm_anti_mean_internal(0.4, 0.2, [0.8, 0.3, 0.05]) == \
            pytest.approx(0.4 * (1.0 - mean_i))

    def test_cterm_anti_mean_internal_empty(self):
        assert score_cterm_anti_mean_internal(0.5, 0.3, []) == \
            pytest.approx(0.5)

    def test_cterm_anti_mean_internal_uniform(self):
        # All internal probs equal → mean == that value
        assert score_cterm_anti_mean_internal(0.6, None, [0.2, 0.2, 0.2]) == \
            pytest.approx(0.6 * 0.8)

    # -- score_nterm_cterm_anti_max_internal --

    def test_nterm_cterm_anti_max_full(self):
        expected = _geomean(0.4, 0.2, 1.0 - 0.8)
        assert score_nterm_cterm_anti_max_internal(0.4, 0.2, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_nterm_cterm_anti_max_no_nterm(self):
        expected = _geomean(0.4, 1.0 - 0.8)
        assert score_nterm_cterm_anti_max_internal(0.4, None, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_nterm_cterm_anti_max_empty_internal(self):
        expected = _geomean(0.4, 0.2, 1.0)  # anti = 1.0
        assert score_nterm_cterm_anti_max_internal(0.4, 0.2, []) == \
            pytest.approx(expected)

    # -- score_nterm_cterm_anti_mean_internal --

    def test_nterm_cterm_anti_mean_full(self):
        mean_i = (0.8 + 0.3) / 2
        expected = _geomean(0.4, 0.2, 1.0 - mean_i)
        assert score_nterm_cterm_anti_mean_internal(0.4, 0.2, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_nterm_cterm_anti_mean_no_nterm(self):
        mean_i = (0.8 + 0.3) / 2
        expected = _geomean(0.4, 1.0 - mean_i)
        assert score_nterm_cterm_anti_mean_internal(0.4, None, [0.8, 0.3]) == \
            pytest.approx(expected)

    def test_nterm_cterm_anti_mean_empty_internal(self):
        expected = _geomean(0.4, 0.2, 1.0)
        assert score_nterm_cterm_anti_mean_internal(0.4, 0.2, []) == \
            pytest.approx(expected)

    # -- cross-function edge cases --

    def test_zero_cterm_propagates_through_all(self):
        for fn in ALL_SCORING_FNS:
            assert fn(0.0, 0.5, [0.3]) == 0.0, fn.__name__

    def test_all_zeros(self):
        for fn in ALL_SCORING_FNS:
            assert fn(0.0, 0.0, [0.0]) == 0.0, fn.__name__

    def test_perfect_scores(self):
        # c=1, n=1, no internal → should be 1.0 for all
        for fn in ALL_SCORING_FNS:
            assert fn(1.0, 1.0, []) == pytest.approx(1.0), fn.__name__


# ======================================================================
# Scoring monotonicity
# ======================================================================

class TestScoringMonotonicity:
    """Higher c_term should raise scores; higher internal should lower scores."""

    @pytest.mark.parametrize("fn", ALL_SCORING_FNS, ids=lambda f: f.__name__)
    def test_higher_cterm_gives_higher_score(self, fn):
        low = fn(0.2, 0.5, [0.3])
        high = fn(0.8, 0.5, [0.3])
        assert high >= low

    @pytest.mark.parametrize("fn", [
        score_cterm_anti_max_internal,
        score_cterm_anti_mean_internal,
        score_nterm_cterm_anti_max_internal,
        score_nterm_cterm_anti_mean_internal,
    ], ids=lambda f: f.__name__)
    def test_higher_internal_gives_lower_score(self, fn):
        low_internal = fn(0.5, 0.5, [0.1])
        high_internal = fn(0.5, 0.5, [0.9])
        assert low_internal >= high_internal


# ======================================================================
# _peptide_score integration
# ======================================================================

class TestPeptideScore:
    """Verify _peptide_score extracts components and delegates correctly."""

    def _make(self, scoring):
        return StubPredictor({PROTEIN: PROBS}, scoring=scoring)

    @pytest.mark.parametrize("fn", ALL_SCORING_FNS, ids=lambda f: f.__name__)
    def test_all_scoring_functions_at_offset_2_length_4(self, fn):
        p = self._make(fn)
        # offset=2, length=4
        # c=PROBS[5]=0.4, n=PROBS[1]=0.2, internal=PROBS[2:5]=[0.8,0.3,0.05]
        actual = p._peptide_score(PROBS, 2, 4)
        expected = fn(0.4, 0.2, [0.8, 0.3, 0.05])
        assert actual == pytest.approx(expected)

    @pytest.mark.parametrize("fn", ALL_SCORING_FNS, ids=lambda f: f.__name__)
    def test_all_scoring_functions_at_offset_0(self, fn):
        p = self._make(fn)
        # offset=0, length=4 → c=PROBS[3]=0.3, n=None, internal=[0.1,0.2,0.8]
        actual = p._peptide_score(PROBS, 0, 4)
        expected = fn(0.3, None, [0.1, 0.2, 0.8])
        assert actual == pytest.approx(expected)

    def test_custom_scoring(self):
        def my_scoring(c, n, internal):
            return c + (n or 0) + sum(internal)
        p = self._make(my_scoring)
        # offset=2, length=4: c=0.4, n=0.2, internal=[0.8,0.3,0.05]
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(
            0.4 + 0.2 + 0.8 + 0.3 + 0.05)

    def test_length_1_peptide(self):
        p = self._make(score_cterm_anti_max_internal)
        # offset=5, length=1 → c=0.4, internal=[]
        # score = 0.4 * (1 - 0) = 0.4
        assert p._peptide_score(PROBS, 5, 1) == pytest.approx(0.4)


# ======================================================================
# predict()
# ======================================================================

class TestPredict:

    def test_basic(self):
        stub = StubPredictor({})
        results = stub.predict(["ABCD", "EFGH"])
        assert len(results) == 2
        for pp in results:
            assert isinstance(pp, PeptidePreds)
            assert len(pp.preds) == 1
            assert pp.preds[0].kind == Kind.antigen_processing
            assert pp.preds[0].predictor_name == "stubpredictor"

    def test_empty_list(self):
        stub = StubPredictor({})
        results = stub.predict([])
        assert results == []

    def test_no_flanks_records_empty_strings(self):
        stub = StubPredictor({})
        results = stub.predict(["ABCD"])
        pred = results[0].preds[0]
        assert pred.n_flank == ""
        assert pred.c_flank == ""

    def test_with_both_flanks(self):
        full_seq = "XXABCDYY"
        probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        stub = StubPredictor({full_seq: probs}, scoring=score_cterm)
        results = stub.predict(["ABCD"], n_flanks=["XX"], c_flanks=["YY"])
        pred = results[0].preds[0]
        assert pred.n_flank == "XX"
        assert pred.c_flank == "YY"
        assert pred.peptide == "ABCD"
        # c_term for peptide at offset=2, length=4 → probs[5] = 0.6
        assert pred.score == pytest.approx(0.6)

    def test_with_only_n_flanks(self):
        full_seq = "XXABCD"
        probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        stub = StubPredictor({full_seq: probs}, scoring=score_cterm)
        results = stub.predict(["ABCD"], n_flanks=["XX"])
        pred = results[0].preds[0]
        assert pred.n_flank == "XX"
        assert pred.c_flank == ""
        # c_term at offset=2, length=4 → probs[5] = 0.6
        assert pred.score == pytest.approx(0.6)

    def test_with_only_c_flanks(self):
        full_seq = "ABCDYY"
        probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        stub = StubPredictor({full_seq: probs}, scoring=score_cterm)
        results = stub.predict(["ABCD"], c_flanks=["YY"])
        pred = results[0].preds[0]
        assert pred.n_flank == ""
        assert pred.c_flank == "YY"
        # c_term at offset=0, length=4 → probs[3] = 0.4
        assert pred.score == pytest.approx(0.4)

    def test_flanks_change_nterm_availability(self):
        """With n_flank, the scoring function receives a non-None n_term."""
        calls = []

        def spy_scoring(c, n, internal):
            calls.append((c, n, internal))
            return c

        # Without flank
        stub = StubPredictor({}, scoring=spy_scoring)
        stub.predict(["ABCD"])
        assert calls[-1][1] is None  # n_term is None

        # With flank
        stub.predict(["ABCD"], n_flanks=["X"])
        assert calls[-1][1] is not None  # n_term is set

    def test_multiple_peptides_with_flanks(self):
        stub = StubPredictor({}, scoring=score_cterm)
        results = stub.predict(
            ["ABC", "DEF"],
            n_flanks=["X", "YY"],
            c_flanks=["Z", "WW"],
        )
        assert len(results) == 2
        assert results[0].preds[0].n_flank == "X"
        assert results[0].preds[0].c_flank == "Z"
        assert results[1].preds[0].n_flank == "YY"
        assert results[1].preds[0].c_flank == "WW"

    def test_score_uses_scoring_function(self):
        """Verify the full chain: probs → components → scoring → Pred."""
        probs = [0.1, 0.2, 0.3, 0.4, 0.5]
        stub = StubPredictor({"ABCDE": probs}, scoring=score_cterm_anti_max_internal)
        results = stub.predict(["ABCDE"])
        pred = results[0].preds[0]
        # offset=0, length=5 → c=probs[4]=0.5, internal=probs[0:4]=[0.1,0.2,0.3,0.4]
        # max_internal = 0.4, score = 0.5 * (1 - 0.4) = 0.3
        assert pred.score == pytest.approx(0.3)


# ======================================================================
# predict_dataframe()
# ======================================================================

class TestPredictDataframe:

    def test_basic(self):
        stub = StubPredictor({})
        df = stub.predict_dataframe(["ABCD"])
        assert list(df.columns) == list(COLUMNS)
        assert len(df) == 1
        assert df.iloc[0]["kind"] == "antigen_processing"

    def test_sample_name(self):
        stub = StubPredictor({})
        df = stub.predict_dataframe(["ABCD"], sample_name="s1")
        assert df.iloc[0]["sample_name"] == "s1"

    def test_with_flanks_recorded(self):
        stub = StubPredictor({})
        df = stub.predict_dataframe(
            ["ABCD"], n_flanks=["XX"], c_flanks=["YY"], sample_name="s1")
        assert df.iloc[0]["n_flank"] == "XX"
        assert df.iloc[0]["c_flank"] == "YY"

    def test_empty_input(self):
        stub = StubPredictor({})
        df = stub.predict_dataframe([])
        assert list(df.columns) == list(COLUMNS)
        assert len(df) == 0

    def test_multiple_peptides(self):
        stub = StubPredictor({})
        df = stub.predict_dataframe(["ABC", "DEF", "GHI"])
        assert len(df) == 3
        assert list(df["peptide"]) == ["ABC", "DEF", "GHI"]


# ======================================================================
# predict_proteins()
# ======================================================================

class TestPredictProteins:

    def test_basic(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"p": PROTEIN})
        assert "p" in result
        pp_list = result["p"]
        expected_count = len(PROTEIN) - 4 + 1  # 9
        assert len(pp_list) == expected_count

    def test_peptide_strings_correct(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"p": PROTEIN})
        peptides = [pp.preds[0].peptide for pp in result["p"]]
        expected = [PROTEIN[i:i+4] for i in range(len(PROTEIN) - 3)]
        assert peptides == expected

    def test_offsets_correct(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"p": PROTEIN})
        offsets = [pp.preds[0].offset for pp in result["p"]]
        assert offsets == list(range(len(PROTEIN) - 3))

    def test_source_sequence_name(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"myprotein": PROTEIN})
        for pp in result["myprotein"]:
            assert pp.preds[0].source_sequence_name == "myprotein"

    def test_scores_match_manual_computation(self):
        stub = StubPredictor({PROTEIN: PROBS},
                             default_peptide_lengths=[4],
                             scoring=score_cterm)
        result = stub.predict_proteins({"p": PROTEIN})
        for pp in result["p"]:
            pred = pp.preds[0]
            # cterm score = PROBS[offset + 4 - 1]
            assert pred.score == pytest.approx(PROBS[pred.offset + 3])

    def test_string_input(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins(PROTEIN)
        assert "seq" in result
        assert len(result["seq"]) == len(PROTEIN) - 3

    def test_list_input(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins([PROTEIN])
        assert PROTEIN in result

    def test_multiple_proteins(self):
        stub = StubPredictor(
            {PROTEIN: PROBS, PROTEIN2: PROBS2},
            default_peptide_lengths=[3])
        result = stub.predict_proteins({"a": PROTEIN, "b": PROTEIN2})
        assert len(result["a"]) == len(PROTEIN) - 2
        assert len(result["b"]) == len(PROTEIN2) - 2

    def test_multiple_peptide_lengths(self):
        stub = StubPredictor({PROTEIN: PROBS},
                             default_peptide_lengths=[3, 4, 5])
        result = stub.predict_proteins({"p": PROTEIN})
        expected = sum(len(PROTEIN) - k + 1 for k in [3, 4, 5])
        assert len(result["p"]) == expected

    def test_peptide_lengths_override(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[9])
        result = stub.predict_proteins({"p": PROTEIN}, peptide_lengths=[3])
        assert len(result["p"]) == len(PROTEIN) - 2

    # -- flank_length --

    def test_flank_length_zero_no_flanks(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"p": PROTEIN}, flank_length=0)
        for pp in result["p"]:
            assert pp.preds[0].n_flank == ""
            assert pp.preds[0].c_flank == ""

    def test_flank_length_records_flanks(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = stub.predict_proteins({"p": PROTEIN}, flank_length=2)
        pp_list = result["p"]

        # First peptide (offset=0): no n_flank, c_flank="EF"
        first = pp_list[0].preds[0]
        assert first.n_flank == ""
        assert first.c_flank == "EF"

        # Second peptide (offset=1): n_flank="A", c_flank="FG"
        second = pp_list[1].preds[0]
        assert second.n_flank == "A"
        assert second.c_flank == "FG"

        # Third peptide (offset=2): n_flank="AB", c_flank="GH"
        third = pp_list[2].preds[0]
        assert third.n_flank == "AB"
        assert third.c_flank == "GH"

        # Last peptide (offset=8): n_flank="GH", c_flank=""
        last = pp_list[-1].preds[0]
        assert last.n_flank == "GH"
        assert last.c_flank == ""

    def test_flank_length_does_not_change_scores(self):
        """flank_length is cosmetic — scores come from the full protein."""
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        without = stub.predict_proteins({"p": PROTEIN}, flank_length=0)
        with_flanks = stub.predict_proteins({"p": PROTEIN}, flank_length=5)
        scores_a = [pp.preds[0].score for pp in without["p"]]
        scores_b = [pp.preds[0].score for pp in with_flanks["p"]]
        for a, b in zip(scores_a, scores_b):
            assert a == pytest.approx(b)


# ======================================================================
# predict_proteins_dataframe()
# ======================================================================

class TestPredictProteinsDataframe:

    def test_basic(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        df = stub.predict_proteins_dataframe({"p": PROTEIN}, sample_name="s1")
        assert list(df.columns) == list(COLUMNS)
        assert (df["source_sequence_name"] == "p").all()
        assert (df["sample_name"] == "s1").all()

    def test_with_flank_length(self):
        stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
        df = stub.predict_proteins_dataframe(
            {"p": PROTEIN}, flank_length=3, sample_name="s1")
        # Middle peptides should have non-empty flanks
        mid = df.iloc[len(df) // 2]
        assert len(mid["n_flank"]) > 0
        assert len(mid["c_flank"]) > 0

    def test_empty_protein(self):
        stub = StubPredictor({}, default_peptide_lengths=[4])
        df = stub.predict_proteins_dataframe({"p": "AB"})
        assert list(df.columns) == list(COLUMNS)
        assert len(df) == 0


# ======================================================================
# predict_cleavage_sites()
# ======================================================================

class TestPredictCleavageSites:

    def test_dict_input(self):
        stub = StubPredictor({PROTEIN: PROBS})
        result = stub.predict_cleavage_sites({"p": PROTEIN})
        assert result["p"] == PROBS

    def test_string_input(self):
        stub = StubPredictor({PROTEIN: PROBS})
        result = stub.predict_cleavage_sites(PROTEIN)
        assert "seq" in result
        assert result["seq"] == PROBS

    def test_multiple_sequences(self):
        stub = StubPredictor({PROTEIN: PROBS, PROTEIN2: PROBS2})
        result = stub.predict_cleavage_sites({"a": PROTEIN, "b": PROTEIN2})
        assert result["a"] == PROBS
        assert result["b"] == PROBS2


# ======================================================================
# _resolve_peptide_lengths()
# ======================================================================

class TestResolvePeptideLengths:

    def test_none_uses_default(self):
        stub = StubPredictor({}, default_peptide_lengths=[8, 10])
        assert stub._resolve_peptide_lengths(None) == [8, 10]

    def test_int_coerced_to_list(self):
        stub = StubPredictor({})
        assert stub._resolve_peptide_lengths(11) == [11]

    def test_list_passed_through(self):
        stub = StubPredictor({})
        assert stub._resolve_peptide_lengths([7, 8]) == [7, 8]

    def test_empty_raises(self):
        stub = StubPredictor({})
        with pytest.raises(ValueError, match="non-empty"):
            stub._resolve_peptide_lengths([])


# ======================================================================
# ProteasomePredictor
# ======================================================================

class TestProteasomePredictor:

    def test_default_scoring(self):
        p = StubProteasome({})
        assert p.scoring is score_cterm_anti_max_internal

    def test_custom_scoring(self):
        p = StubProteasome({}, scoring=score_cterm)
        assert p.scoring is score_cterm

    def test_pred_kind(self):
        p = StubProteasome({})
        assert p._pred_kind() == Kind.proteasome_cleavage

    def test_predict_uses_proteasome_kind(self):
        p = StubProteasome({})
        results = p.predict(["ABCD"])
        assert results[0].preds[0].kind == Kind.proteasome_cleavage

    def test_predict_proteins_uses_proteasome_kind(self):
        p = StubProteasome({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = p.predict_proteins({"p": PROTEIN})
        for pp in result["p"]:
            assert pp.preds[0].kind == Kind.proteasome_cleavage

    def test_inherits_from_processing_predictor(self):
        assert issubclass(ProteasomePredictor, ProcessingPredictor)

    def test_scoring_formula_matches(self):
        """Verify the default scoring matches c_term * (1 - max_internal)."""
        p = StubProteasome({PROTEIN: PROBS}, default_peptide_lengths=[4])
        result = p.predict_proteins({"p": PROTEIN})
        for pp in result["p"]:
            pred = pp.preds[0]
            c = PROBS[pred.offset + 3]
            internal = PROBS[pred.offset:pred.offset + 3]
            max_i = max(internal) if internal else 0.0
            expected = c * (1.0 - max_i)
            assert pred.score == pytest.approx(expected)

    def test_predict_with_flanks(self):
        """Flanking works through the ProteasomePredictor layer."""
        full_seq = "XXABCDYY"
        probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        p = StubProteasome({full_seq: probs}, scoring=score_cterm)
        results = p.predict(["ABCD"], n_flanks=["XX"], c_flanks=["YY"])
        pred = results[0].preds[0]
        assert pred.kind == Kind.proteasome_cleavage
        assert pred.n_flank == "XX"
        assert pred.c_flank == "YY"


# ======================================================================
# Integration: scoring consistency across predict() and predict_proteins()
# ======================================================================

class TestScoringConsistency:
    """
    When predict() is given flanks that reconstruct the full protein,
    the score at that position should match predict_proteins().
    """

    def test_flanked_predict_matches_predict_proteins(self):
        stub = StubPredictor(
            {PROTEIN: PROBS},
            default_peptide_lengths=[4],
            scoring=score_cterm_anti_max_internal,
        )

        # Get scores from predict_proteins
        protein_result = stub.predict_proteins({"p": PROTEIN})
        protein_scores = {
            pp.preds[0].offset: pp.preds[0].score
            for pp in protein_result["p"]
        }

        # Reconstruct the same predictions via predict() with flanks
        for offset in range(len(PROTEIN) - 3):
            peptide = PROTEIN[offset:offset + 4]
            n_flank = PROTEIN[:offset]
            c_flank = PROTEIN[offset + 4:]
            # The full_seq is the original protein
            results = stub.predict(
                [peptide], n_flanks=[n_flank], c_flanks=[c_flank])
            pred_score = results[0].preds[0].score
            assert pred_score == pytest.approx(protein_scores[offset]), \
                f"Mismatch at offset {offset}"

    @pytest.mark.parametrize("fn", ALL_SCORING_FNS, ids=lambda f: f.__name__)
    def test_predict_proteins_all_scoring_functions(self, fn):
        """Every scoring function should produce valid results."""
        stub = StubPredictor(
            {PROTEIN: PROBS}, default_peptide_lengths=[4], scoring=fn)
        result = stub.predict_proteins({"p": PROTEIN})
        for pp in result["p"]:
            # Should be a finite number
            assert math.isfinite(pp.preds[0].score), fn.__name__
