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

import os

import numpy as np
import pytest

from mhctools import TCR
from mhctools.nettcr import _BLOSUM50, _encode_feature
from mhctools.pred import COLUMNS, Kind, PeptideResult


# ---------------------------------------------------------------------------
# Encoding unit tests — no NetTCR install or TFLite runtime required.
# These pin the exact encoding reproduced from NetTCR-2.2 src/predict.py.
# ---------------------------------------------------------------------------

def test_encode_shape():
    arr = _encode_feature(["SIINFEKL"], "pep")
    assert arr.shape == (1, 12, 20)   # pep max length is 12
    assert arr.dtype == np.float32


def test_encode_residues_are_blosum_over_five():
    # First peptide residue 'S' -> BLOSUM50['S'] / 5
    arr = _encode_feature(["SIINFEKL"], "pep")
    np.testing.assert_allclose(arr[0, 0], _BLOSUM50["S"] / 5.0, rtol=1e-6)


def test_encode_padding_is_minus_one():
    # An 8-mer in a length-12 feature leaves positions 8..11 padded.
    arr = _encode_feature(["SIINFEKL"], "pep")
    assert np.all(arr[0, 8:] == -1.0)


def test_encode_is_case_insensitive():
    lower = _encode_feature(["siinfekl"], "pep")
    upper = _encode_feature(["SIINFEKL"], "pep")
    np.testing.assert_array_equal(lower, upper)


def test_encode_rejects_too_long():
    with pytest.raises(ValueError, match="max length"):
        _encode_feature(["A" * 13], "pep")   # pep max is 12


def test_encode_rejects_unknown_amino_acid():
    with pytest.raises(ValueError, match="Unknown amino acid"):
        _encode_feature(["SIINFEKB"], "pep")


def test_encode_batch():
    arr = _encode_feature(["SIINFEKL", "GILGFVFTL"], "pep")
    assert arr.shape == (2, 12, 20)


# ---------------------------------------------------------------------------
# Model tests — require a cloned NetTCR-2.2 (weights) and a TFLite runtime.
# ---------------------------------------------------------------------------

NETTCR_DIR = None
for candidate in [
    os.environ.get("NETTCR_DIR", ""),
    os.path.join(os.path.expanduser("~"), "NetTCR-2.2"),
    os.path.join(os.path.expanduser("~"), "code", "NetTCR-2.2"),
]:
    if candidate and os.path.isdir(
            os.path.join(candidate, "models", "nettcr_2_2_pan", "checkpoint")):
        NETTCR_DIR = candidate
        break

pytestmark = pytest.mark.skipif(
    NETTCR_DIR is None,
    reason="NetTCR-2.2 not installed (set NETTCR_DIR or clone to ~/NetTCR-2.2)")


# Three peptides paired with their cognate-ish TCRs. LLWNGPMAV is a strong
# recognised pair; RAKFKQLL here is a mismatch (should score low).
PAIRS = [
    ("RAKFKQLL", TCR(cdr1a="NSAFQY", cdr2a="TYSSGN", cdr3a="AMSGDGGSQGNLI",
                     cdr1b="LNHDA", cdr2b="SQIVND", cdr3b="ASSIRAAYEQY")),
    ("GILGFVFTL", TCR(cdr1a="TSGFYG", cdr2a="NALDGL", cdr3a="AVRPTSGGSYIPT",
                      cdr1b="SGHRS", cdr2b="YFSETQ", cdr3b="ASSIRSSYEQY")),
    ("LLWNGPMAV", TCR(cdr1a="NSASQS", cdr2a="VYSSG", cdr3a="VVEGDKVI",
                      cdr1b="MGHRA", cdr2b="YSYEKL", cdr3b="ASSHSGYEQF")),
]


@pytest.fixture(scope="module")
def predictor():
    from mhctools import NetTCR
    return NetTCR(nettcr_path=NETTCR_DIR)


def test_init_lazy(predictor):
    assert predictor._interpreters is None
    assert len(predictor._model_paths) > 0
    assert "not loaded" in str(predictor)


def test_predict_pairs_shape(predictor):
    results = predictor.predict_pairs(PAIRS)
    assert len(results) == len(PAIRS)
    for pp in results:
        assert isinstance(pp, PeptideResult)
        assert len(pp.preds) == 1


def test_predict_pairs_kind_and_fields(predictor):
    results = predictor.predict_pairs(PAIRS)
    for (pep, tcr), pp in zip(PAIRS, results):
        pred = pp.preds[0]
        assert pred.kind == Kind.pMHC_TCR_binding
        assert pred.peptide == pep
        assert pred.tcr == tcr.identifier
        assert pred.allele == ""
        assert pred.predictor_name == "nettcr"


def test_predict_pairs_scores_are_probabilities(predictor):
    results = predictor.predict_pairs(PAIRS)
    for pp in results:
        assert 0.0 <= pp.preds[0].score <= 1.0


def test_predict_recovers_binding_signal(predictor):
    """The cognate LLWNGPMAV pair should score well above the RAKFKQLL
    mismatch — a coarse check that encoding + ensemble are wired correctly."""
    scores = {pep: pp.preds[0].score
              for (pep, _), pp in zip(PAIRS, predictor.predict_pairs(PAIRS))}
    assert scores["LLWNGPMAV"] > 0.5
    assert scores["LLWNGPMAV"] > scores["RAKFKQLL"]


def test_predict_cross_product(predictor):
    peptides = ["LLWNGPMAV", "GILGFVFTL"]
    tcrs = [PAIRS[0][1], PAIRS[2][1]]
    results = predictor.predict(peptides, tcrs)
    assert len(results) == len(peptides)
    for pp in results:
        assert len(pp.preds) == len(tcrs)
        assert pp.tcrs == {t.identifier for t in tcrs}


def test_predict_single_peptide_single_tcr(predictor):
    results = predictor.predict("LLWNGPMAV", PAIRS[2][1])
    assert len(results) == 1
    assert len(results[0].preds) == 1


def test_predict_models_stay_loaded(predictor):
    predictor.predict_pairs(PAIRS[:1])
    assert predictor._interpreters is not None
    assert "loaded" in str(predictor)


def test_predict_repeated_calls_consistent(predictor):
    r1 = predictor.predict_pairs(PAIRS)
    r2 = predictor.predict_pairs(PAIRS)
    for a, b in zip(r1, r2):
        assert a.preds[0].score == b.preds[0].score


def test_predict_dataframe_schema(predictor):
    df = predictor.predict_dataframe(["LLWNGPMAV"], [PAIRS[2][1]])
    assert list(df.columns) == list(COLUMNS)
    assert df["tcr"].iloc[0] == PAIRS[2][1].identifier
    assert df["kind"].iloc[0] == Kind.pMHC_TCR_binding


def test_bad_tcr_type_raises(predictor):
    with pytest.raises(TypeError, match="TCR"):
        predictor.predict_pairs([("SIINFEKL", "not-a-tcr")])
