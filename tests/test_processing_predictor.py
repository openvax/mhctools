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

import math
import pytest

from mhctools.processing_predictor import ProcessingPredictor
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


# A 12-residue "protein" with known cleavage probs
PROTEIN = "ABCDEFGHIJKL"
#          0    1    2    3    4    5    6    7    8    9   10   11
PROBS = [0.1, 0.2, 0.8, 0.3, 0.05, 0.4, 0.1, 0.9, 0.2, 0.6, 0.3, 0.7]


def _geomean(*values):
    product = 1.0
    for v in values:
        product *= v
    return product ** (1.0 / len(values))


# -- init --

def test_init_defaults():
    p = StubPredictor({})
    assert p.default_peptide_lengths == [9]
    assert p.scoring == "nterm_cterm_max_internal"
    assert not hasattr(p, "alleles")


def test_init_invalid_scoring():
    with pytest.raises(ValueError, match="scoring"):
        StubPredictor({}, scoring="garbage")


# -- _peptide_score math --

class TestPeptideScore:
    """Test aggregation formulas on a known cleavage profile."""

    def _make(self, scoring):
        return StubPredictor({PROTEIN: PROBS}, scoring=scoring)

    def test_cterm(self):
        p = self._make("cterm")
        # peptide at offset=2, length=4 → CDEF → c_term = probs[5] = 0.4
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(0.4)

    def test_nterm_cterm(self):
        p = self._make("nterm_cterm")
        # offset=2, length=4 → n_term=probs[1]=0.2, c_term=probs[5]=0.4
        expected = _geomean(0.4, 0.2)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_nterm_cterm_offset_zero(self):
        p = self._make("nterm_cterm")
        # offset=0 → no n_term, falls back to cterm only
        # c_term = probs[3] = 0.3
        assert p._peptide_score(PROBS, 0, 4) == pytest.approx(0.3)

    def test_cterm_max_internal(self):
        p = self._make("cterm_max_internal")
        # offset=2, length=4 → c_term=0.4, internal=[0.8, 0.3, 0.05]
        # max_internal=0.8, anti=0.2
        expected = _geomean(0.4, 1.0 - 0.8)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_cterm_mean_internal(self):
        p = self._make("cterm_mean_internal")
        # internal=[0.8, 0.3, 0.05], mean=0.383333...
        expected = _geomean(0.4, 1.0 - (0.8 + 0.3 + 0.05) / 3)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_nterm_cterm_max_internal(self):
        p = self._make("nterm_cterm_max_internal")
        # n=0.2, c=0.4, max_internal=0.8, anti=0.2
        expected = _geomean(0.4, 0.2, 1.0 - 0.8)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_nterm_cterm_mean_internal(self):
        p = self._make("nterm_cterm_mean_internal")
        mean_i = (0.8 + 0.3 + 0.05) / 3
        expected = _geomean(0.4, 0.2, 1.0 - mean_i)
        assert p._peptide_score(PROBS, 2, 4) == pytest.approx(expected)

    def test_length_2_no_internal(self):
        # length=2 → internal has only 1 position (offset to offset+length-2)
        p = self._make("cterm_max_internal")
        # offset=3, length=2 → DEFG? No, peptide DE → internal=[probs[3]]=0.3
        # c_term=probs[4]=0.05
        expected = _geomean(0.05, 1.0 - 0.3)
        assert p._peptide_score(PROBS, 3, 2) == pytest.approx(expected)

    def test_zero_cterm(self):
        p = self._make("nterm_cterm_max_internal")
        # If c_term is 0, score should be 0
        probs = [0.5, 0.5, 0.0]  # c_term at pos 2 = 0
        score = p._peptide_score(probs, 0, 3)
        assert score == pytest.approx(0.0)


# -- predict --

def test_predict():
    stub = StubPredictor({})
    results = stub.predict(["ABCD", "EFGH"])
    assert len(results) == 2
    for pp in results:
        assert isinstance(pp, PeptidePreds)
        assert len(pp.preds) == 1
        assert pp.preds[0].kind == Kind.antigen_processing


def test_predict_dataframe():
    stub = StubPredictor({})
    df = stub.predict_dataframe(["ABCD"], sample_name="s1")
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 1
    assert df.iloc[0]["kind"] == "antigen_processing"
    assert df.iloc[0]["sample_name"] == "s1"


# -- predict_proteins --

def test_predict_proteins():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins({"p": PROTEIN})
    assert "p" in result
    pp_list = result["p"]
    expected_count = len(PROTEIN) - 4 + 1
    assert len(pp_list) == expected_count
    for pp in pp_list:
        pred = pp.preds[0]
        assert pred.kind == Kind.antigen_processing
        assert pred.source_sequence_name == "p"


def test_predict_proteins_string():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    result = stub.predict_proteins(PROTEIN)
    assert "seq" in result


def test_predict_proteins_dataframe():
    stub = StubPredictor({PROTEIN: PROBS}, default_peptide_lengths=[4])
    df = stub.predict_proteins_dataframe({"p": PROTEIN}, sample_name="s1")
    assert list(df.columns) == list(COLUMNS)
    assert (df["source_sequence_name"] == "p").all()
    assert (df["sample_name"] == "s1").all()


# -- predict_cleavage_sites --

def test_predict_cleavage_sites():
    stub = StubPredictor({PROTEIN: PROBS})
    result = stub.predict_cleavage_sites({"p": PROTEIN})
    assert result["p"] == PROBS


def test_predict_cleavage_sites_string():
    stub = StubPredictor({PROTEIN: PROBS})
    result = stub.predict_cleavage_sites(PROTEIN)
    assert "seq" in result
