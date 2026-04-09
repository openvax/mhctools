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
        predict_py = os.path.join(candidate, "src", "predict.py")
        if os.path.isfile(predict_py):
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


def test_init_im():
    p = BigMHC(ALLELES, mode="im", bigmhc_path=BIGMHC_DIR)
    assert p.mode == "im"


def test_init_invalid_mode():
    with pytest.raises(ValueError, match="mode"):
        BigMHC(ALLELES, mode="bad", bigmhc_path=BIGMHC_DIR)


def test_str():
    p = BigMHC(ALLELES, mode="el", bigmhc_path=BIGMHC_DIR)
    s = str(p)
    assert "BigMHC" in s
    assert "el" in s


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


# -- dataframe --

def test_predict_dataframe(el_predictor):
    df = el_predictor.predict_dataframe(PEPTIDES)
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == len(PEPTIDES) * len(ALLELES)
    assert (df["kind"] == "pMHC_presentation").all()
    assert (df["predictor_name"] == "bigmhc_el").all()


# -- single string input --

def test_predict_single_string(el_predictor):
    results = el_predictor.predict("SIINFEKL")
    assert len(results) == 1
    assert results[0].preds[0].peptide == "SIINFEKL"
