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

import os

import pytest

from mhctools.bigmhc import BigMHC
from mhctools.pred import Kind, PeptidePreds, COLUMNS


# Skip entire module if BigMHC is not installed
BIGMHC_DIR = None
for candidate in [
    os.environ.get("BIGMHC_DIR", ""),
    os.path.join(os.path.expanduser("~"), "bigmhc"),
    os.path.join(os.path.expanduser("~"), "code", "bigmhc"),
]:
    if candidate and os.path.isdir(candidate):
        src_dir = os.path.join(candidate, "src", "bigmhc.py")
        if os.path.isfile(src_dir):
            BIGMHC_DIR = candidate
            break

pytestmark = pytest.mark.skipif(
    BIGMHC_DIR is None,
    reason="BigMHC not installed (set BIGMHC_DIR or clone to ~/bigmhc)")

ALLELES = ["HLA-A*02:01"]
PEPTIDES = ["SIINFEKL", "GILGFVFTL", "NLVPMVATV"]


# -- init --

def test_init_el():
    p = BigMHC(ALLELES, mode="el", bigmhc_path=BIGMHC_DIR)
    assert p.mode == "el"
    assert p.alleles == ALLELES
    assert p._models is None  # lazy, not loaded yet


def test_init_im():
    p = BigMHC(ALLELES, mode="im", bigmhc_path=BIGMHC_DIR)
    assert p.mode == "im"


def test_init_string_allele():
    p = BigMHC("HLA-A*02:01", mode="el", bigmhc_path=BIGMHC_DIR)
    assert p.alleles == ["HLA-A*02:01"]


def test_init_invalid_mode():
    with pytest.raises(ValueError, match="mode"):
        BigMHC(ALLELES, mode="bad", bigmhc_path=BIGMHC_DIR)


def test_init_missing_path():
    with pytest.raises(FileNotFoundError):
        BigMHC(ALLELES, bigmhc_path="/nonexistent/path")


def test_str_before_load():
    p = BigMHC(ALLELES, mode="el", bigmhc_path=BIGMHC_DIR)
    s = str(p)
    assert "BigMHC" in s
    assert "not loaded" in s


# -- predict (EL) --

@pytest.fixture(scope="module")
def el_predictor():
    return BigMHC(ALLELES, mode="el", bigmhc_path=BIGMHC_DIR)


def test_predict_el_returns_peptide_preds(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    assert isinstance(results, list)
    assert len(results) == len(PEPTIDES)
    for pp in results:
        assert isinstance(pp, PeptidePreds)
        assert len(pp.preds) == len(ALLELES)


def test_predict_el_scores_are_probabilities(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    for pp in results:
        for pred in pp.preds:
            assert 0.0 <= pred.score <= 1.0


def test_predict_el_kind(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    for pp in results:
        for pred in pp.preds:
            assert pred.kind == Kind.pMHC_presentation


def test_predict_el_predictor_name(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    assert results[0].preds[0].predictor_name == "bigmhc_el"


def test_predict_el_peptides_correct(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    for pep, pp in zip(PEPTIDES, results):
        assert pp.preds[0].peptide == pep


def test_predict_el_alleles_correct(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    for pp in results:
        assert pp.preds[0].allele == ALLELES[0]


def test_predict_el_scores_vary(el_predictor):
    results = el_predictor.predict(PEPTIDES)
    scores = [pp.preds[0].score for pp in results]
    assert len(set(round(s, 4) for s in scores)) > 1


def test_predict_models_stay_loaded(el_predictor):
    """After first predict(), models should stay in memory."""
    el_predictor.predict(["SIINFEKL"])
    assert el_predictor._models is not None
    assert "loaded" in str(el_predictor)


def test_predict_repeated_calls_consistent(el_predictor):
    """Multiple predict() calls return identical scores."""
    r1 = el_predictor.predict(PEPTIDES)
    r2 = el_predictor.predict(PEPTIDES)
    for pp1, pp2 in zip(r1, r2):
        assert abs(pp1.preds[0].score - pp2.preds[0].score) < 1e-6


# -- predict (IM) --

@pytest.fixture(scope="module")
def im_predictor():
    return BigMHC(ALLELES, mode="im", bigmhc_path=BIGMHC_DIR)


def test_predict_im_kind(im_predictor):
    results = im_predictor.predict(PEPTIDES)
    for pp in results:
        for pred in pp.preds:
            assert pred.kind == Kind.immunogenicity


def test_predict_im_predictor_name(im_predictor):
    results = im_predictor.predict(PEPTIDES)
    assert results[0].preds[0].predictor_name == "bigmhc_im"


def test_predict_im_scores_are_probabilities(im_predictor):
    results = im_predictor.predict(PEPTIDES)
    for pp in results:
        for pred in pp.preds:
            assert 0.0 <= pred.score <= 1.0


def test_predict_el_im_scores_differ(el_predictor, im_predictor):
    """EL and IM models should produce different scores."""
    el = el_predictor.predict(PEPTIDES)
    im = im_predictor.predict(PEPTIDES)
    el_scores = [pp.preds[0].score for pp in el]
    im_scores = [pp.preds[0].score for pp in im]
    assert el_scores != im_scores


# -- ordering --

def test_predict_mixed_lengths_ordering(el_predictor):
    """Peptides of different lengths are batched separately;
    verify output order matches input order."""
    peptides = ["YLQPRTFL", "GILGFVFTL", "YLQPRTFLLK"]  # 8, 9, 10
    results = el_predictor.predict(peptides)
    for pep, pp in zip(peptides, results):
        assert pp.preds[0].peptide == pep


def test_predict_duplicate_peptides(el_predictor):
    """Duplicate peptides should produce duplicate results in order."""
    peptides = ["SIINFEKL", "GILGFVFTL", "SIINFEKL"]
    results = el_predictor.predict(peptides)
    assert len(results) == 3
    assert results[0].preds[0].peptide == "SIINFEKL"
    assert results[1].preds[0].peptide == "GILGFVFTL"
    assert results[2].preds[0].peptide == "SIINFEKL"
    assert abs(results[0].preds[0].score - results[2].preds[0].score) < 1e-6


# -- multiple alleles --

def test_predict_multiple_alleles():
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    p = BigMHC(alleles, mode="el", bigmhc_path=BIGMHC_DIR)
    results = p.predict(["SIINFEKL"])
    assert len(results) == 1
    pp = results[0]
    assert len(pp.preds) == 2
    assert pp.preds[0].allele == "HLA-A*02:01"
    assert pp.preds[1].allele == "HLA-B*07:02"


def test_predict_multiple_alleles_mixed_lengths():
    """Multi-allele + mixed-length peptides: ordering stress test."""
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    peptides = ["YLQPRTFL", "GILGFVFTL", "YLQPRTFLLK"]
    p = BigMHC(alleles, mode="el", bigmhc_path=BIGMHC_DIR)
    results = p.predict(peptides)
    assert len(results) == 3
    for pep, pp in zip(peptides, results):
        assert len(pp.preds) == 2
        assert pp.preds[0].peptide == pep
        assert pp.preds[1].peptide == pep
        assert pp.preds[0].allele == "HLA-A*02:01"
        assert pp.preds[1].allele == "HLA-B*07:02"


# -- batch performance --

def test_predict_long_list(el_predictor):
    """Run on a realistic protein-length input to verify batching works."""
    protein = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAI"
    peptides = [protein[i:i + 9] for i in range(len(protein) - 9 + 1)]
    assert len(peptides) > 50
    results = el_predictor.predict(peptides)
    assert len(results) == len(peptides)
    scores = [pp.preds[0].score for pp in results]
    assert all(0.0 <= s <= 1.0 for s in scores)
    # Scores should not all be identical
    assert len(set(round(s, 4) for s in scores)) > 1


# -- dataframe --

def test_predict_dataframe(el_predictor):
    df = el_predictor.predict_dataframe(PEPTIDES)
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == len(PEPTIDES) * len(ALLELES)
    assert (df["kind"] == "pMHC_presentation").all()
    assert (df["predictor_name"] == "bigmhc_el").all()


def test_predict_dataframe_multiple_alleles():
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    p = BigMHC(alleles, mode="el", bigmhc_path=BIGMHC_DIR)
    df = p.predict_dataframe(PEPTIDES)
    assert len(df) == len(PEPTIDES) * len(alleles)


# -- single string input --

def test_predict_single_string(el_predictor):
    results = el_predictor.predict("SIINFEKL")
    assert len(results) == 1
    assert results[0].preds[0].peptide == "SIINFEKL"
