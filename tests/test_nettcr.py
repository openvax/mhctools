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
import pandas as pd
import pytest

from mhctools import TCR
from mhctools.nettcr import (
    NetTCR,
    _BLOSUM50,
    _encode_feature,
    _suppress_native_stderr,
)
from mhctools.pred import COLUMNS, Kind


# ---------------------------------------------------------------------------
# Encoding unit tests — no NetTCR install or TFLite runtime required.
# These pin the exact encoding reproduced from NetTCR-2.2 src/predict.py.
# ---------------------------------------------------------------------------

def test_suppress_native_stderr_hides_fd2_writes_and_restores(capfd):
    # Emulate a native library writing to fd 2 (like TFLite's XNNPACK line):
    # writes inside the block are dropped, and stderr works again afterward.
    with _suppress_native_stderr():
        os.write(2, b"NATIVE-XNNPACK-NOISE\n")
    os.write(2, b"real-stderr-after\n")
    err = capfd.readouterr().err
    assert "NATIVE-XNNPACK-NOISE" not in err
    assert "real-stderr-after" in err


def test_suppress_native_stderr_restores_on_exception(capfd):
    # fd 2 must be restored even if the body raises.
    with pytest.raises(ValueError):
        with _suppress_native_stderr():
            raise ValueError("boom")
    os.write(2, b"still-works\n")  # would raise if fd 2 were left closed
    assert "still-works" in capfd.readouterr().err


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
# Constructor error paths — no NetTCR install required.
# ---------------------------------------------------------------------------

def test_init_missing_path_raises():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        NetTCR(nettcr_path="/nonexistent/nettcr")


def test_init_no_models_raises(tmp_path):
    # Directory exists but contains no *.tflite ensemble.
    with pytest.raises(FileNotFoundError, match="No NetTCR"):
        NetTCR(nettcr_path=str(tmp_path))


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

requires_nettcr = pytest.mark.skipif(
    NETTCR_DIR is None,
    reason="NetTCR-2.2 not installed (set NETTCR_DIR or clone to ~/NetTCR-2.2)")


class _FakeNetTCR(NetTCR):
    def __init__(self):
        self.calls = []

    def _predict_raw(self, peptides, tcrs):
        self.calls.append((
            list(peptides),
            [tcr.cdr_dict() for tcr in tcrs],
        ))
        return np.array(
            [0.25 + 0.1 * i for i in range(len(peptides))],
            dtype=np.float32)


def _nettcr_dataframe(column_names=None):
    values = {
        "peptide": ["SPRWYFYYL", "AVFDRKSDAK"],
        "cdr1a": ["KALYS", "VGISA"],
        "cdr2a": ["LLKGGEQ", "LSSGK"],
        "cdr3a": ["GTEIGGGTSYGKLT", "AVFNTGNQFY"],
        "cdr1b": ["MNHEY", "SGDLS"],
        "cdr2b": ["SMNVEV", "YYNGEE"],
        "cdr3b": ["ASGTETQY", "ASTPWGRGTDTQY"],
    }
    if column_names:
        values = {column_names.get(key, key): value
                  for key, value in values.items()}
    return pd.DataFrame(values)


def test_predict_dataframe_builds_tcrs_from_canonical_columns():
    predictor = _FakeNetTCR()
    df = _nettcr_dataframe()

    out = predictor.predict_dataframe(df, sample_name="sample1")

    assert predictor.calls == [(
        ["SPRWYFYYL", "AVFDRKSDAK"],
        [
            {
                "a1": "KALYS",
                "a2": "LLKGGEQ",
                "a3": "GTEIGGGTSYGKLT",
                "b1": "MNHEY",
                "b2": "SMNVEV",
                "b3": "ASGTETQY",
            },
            {
                "a1": "VGISA",
                "a2": "LSSGK",
                "a3": "AVFNTGNQFY",
                "b1": "SGDLS",
                "b2": "YYNGEE",
                "b3": "ASTPWGRGTDTQY",
            },
        ],
    )]
    assert list(out.columns) == list(COLUMNS)
    assert out["sample_name"].tolist() == ["sample1", "sample1"]
    assert out["peptide"].tolist() == ["SPRWYFYYL", "AVFDRKSDAK"]
    assert out["tcr"].tolist() == [
        "GTEIGGGTSYGKLT/ASGTETQY",
        "AVFNTGNQFY/ASTPWGRGTDTQY",
    ]
    assert out["score"].tolist() == pytest.approx([0.25, 0.35])


def test_predict_dataframe_accepts_custom_column_mapping():
    predictor = _FakeNetTCR()
    df = _nettcr_dataframe({
        "peptide": "epitope",
        "cdr1a": "A1",
        "cdr2a": "A2",
        "cdr3a": "A3",
        "cdr1b": "B1",
        "cdr2b": "B2",
        "cdr3b": "B3",
    })

    out = predictor.predict_dataframe(
        df,
        peptide_col="epitope",
        cdr_cols={
            "cdr1a": "A1",
            "cdr2a": "A2",
            "cdr3a": "A3",
            "cdr1b": "B1",
            "cdr2b": "B2",
            "cdr3b": "B3",
        })

    assert predictor.calls[0][0] == ["SPRWYFYYL", "AVFDRKSDAK"]
    assert predictor.calls[0][1][0]["a3"] == "GTEIGGGTSYGKLT"
    assert out["peptide"].tolist() == ["SPRWYFYYL", "AVFDRKSDAK"]


def test_predict_dataframe_preserves_existing_cross_product_api():
    predictor = _FakeNetTCR()
    tcr = TCR(cdr3a="CAVR", cdr3b="CASS")

    out = predictor.predict_dataframe(["SPRWYFYYL"], [tcr])

    assert predictor.calls == [(
        ["SPRWYFYYL"],
        [{
            "a1": "",
            "a2": "",
            "a3": "CAVR",
            "b1": "",
            "b2": "",
            "b3": "CASS",
        }],
    )]
    assert out["tcr"].tolist() == ["CAVR/CASS"]


def test_predict_dataframe_missing_column_raises():
    predictor = _FakeNetTCR()
    df = _nettcr_dataframe().drop(columns=["cdr2b"])

    with pytest.raises(ValueError, match="cdr2b"):
        predictor.predict_dataframe(df)


def test_predict_dataframe_unknown_cdr_mapping_key_raises():
    predictor = _FakeNetTCR()
    df = _nettcr_dataframe()

    with pytest.raises(ValueError, match="cdr4a"):
        predictor.predict_dataframe(df, cdr_cols={"cdr4a": "cdr4a"})


# Reference ensemble predictions produced by running NetTCR-2.2's OWN
# `src/predict.py` over all 20 pan-model checkpoints and averaging the
# outputs -- the canonical ensemble defined in `src/make_webserver_prediction.py`
# (`avg_prediction / 20`). These come from upstream code, not this wrapper, so
# the test is not circular. Inputs are real TCRs from NetTCR's own
# `data/nettcr_2_2_limited_dataset.csv`. Fields: (peptide, A1, A2, A3, B1, B2, B3).
PUBLISHED_ENSEMBLE = [
    (("SPRWYFYYL", "KALYS", "LLKGGEQ", "GTEIGGGTSYGKLT", "MNHEY", "SMNVEV", "ASGTETQY"), 0.060078),
    (("KSKRTPMGF", "DSAIYN", "IQSSQRE", "AVRNYGGATNKLI", "PRHDT", "FYEKMQ", "ASSLTTGGRNEQF"), 0.062686),
    (("AVFDRKSDAK", "VGISA", "LSSGK", "AVFNTGNQFY", "SGDLS", "YYNGEE", "ASTPWGRGTDTQY"), 0.566735),
    (("RPPIFIRRL", "TTLSN", "LVKSGEV", "AGADAGNNRKLI", "SGHRS", "YFSETQ", "ASSLDQGAYEQY"), 0.011056),
    (("KLGGALQAK", "DSAIYN", "IQSSQRE", "AVRPHSGGGADGLT", "SGHDY", "FNNNVP", "ASSPGDYGYT"), 0.410461),
    (("GILGFVFTL", "VSGLRG", "LYSAGEE", "AVPTILTGGGNKLT", "LNHNV", "YYDKDF", "ATSRVQETQY"), 0.024195),
]
# AVFDRKSDAK is the one true binder (binder=1) in this sample.
PUBLISHED_BINDER = "AVFDRKSDAK"


def _row_to_pair(row):
    peptide, a1, a2, a3, b1, b2, b3 = row
    return peptide, TCR(cdr1a=a1, cdr2a=a2, cdr3a=a3,
                        cdr1b=b1, cdr2b=b2, cdr3b=b3)


@pytest.fixture(scope="module")
def predictor():
    return NetTCR(nettcr_path=NETTCR_DIR)


@requires_nettcr
def test_init_lazy(predictor):
    assert predictor._interpreters is None
    assert len(predictor._model_paths) > 0
    assert "not loaded" in str(predictor)


@requires_nettcr
def test_reproduces_published_ensemble(predictor):
    """The wrapper must reproduce NetTCR's own 20-model ensemble output."""
    pairs = [_row_to_pair(row) for row, _ in PUBLISHED_ENSEMBLE]
    results = predictor.predict_pairs(pairs)
    for (row, expected), pp in zip(PUBLISHED_ENSEMBLE, results):
        got = pp.preds[0].score
        assert got == pytest.approx(expected, abs=1e-3), (
            "%s: expected %.6f (upstream ensemble), got %.6f"
            % (row[0], expected, got))


@requires_nettcr
def test_published_binder_ranks_top(predictor):
    """The labeled binder should outscore every non-binder in the sample."""
    pairs = [_row_to_pair(row) for row, _ in PUBLISHED_ENSEMBLE]
    scores = {row[0]: pp.preds[0].score
              for (row, _), pp in zip(PUBLISHED_ENSEMBLE,
                                      predictor.predict_pairs(pairs))}
    binder = scores[PUBLISHED_BINDER]
    others = [s for pep, s in scores.items() if pep != PUBLISHED_BINDER]
    assert binder > max(others)


@requires_nettcr
def test_predict_pairs_kind_and_fields(predictor):
    pairs = [_row_to_pair(row) for row, _ in PUBLISHED_ENSEMBLE]
    for (row, _), pp in zip(PUBLISHED_ENSEMBLE, predictor.predict_pairs(pairs)):
        pep, tcr = _row_to_pair(row)
        pred = pp.preds[0]
        assert pred.kind == Kind.pMHC_TCR_binding
        assert pred.peptide == pep
        assert pred.tcr == tcr.identifier
        assert pred.allele == ""
        assert pred.predictor_name == "nettcr"
        assert pred.value is None          # no native units for TCR binding
        assert 0.0 <= pred.score <= 1.0


@requires_nettcr
def test_predict_cross_product(predictor):
    peptides = ["AVFDRKSDAK", "GILGFVFTL"]
    tcrs = [_row_to_pair(PUBLISHED_ENSEMBLE[2][0])[1],
            _row_to_pair(PUBLISHED_ENSEMBLE[5][0])[1]]
    results = predictor.predict(peptides, tcrs)
    assert len(results) == len(peptides)
    for pp in results:
        assert len(pp.preds) == len(tcrs)
        assert pp.tcrs == {t.identifier for t in tcrs}


@requires_nettcr
def test_predict_single_peptide_single_tcr(predictor):
    _, tcr = _row_to_pair(PUBLISHED_ENSEMBLE[2][0])
    results = predictor.predict("AVFDRKSDAK", tcr)
    assert len(results) == 1
    assert len(results[0].preds) == 1


@requires_nettcr
def test_batch_matches_single(predictor):
    """A batched call and per-pair calls must give identical scores
    (guards the allocation-caching path against batch-size bugs)."""
    pairs = [_row_to_pair(row) for row, _ in PUBLISHED_ENSEMBLE]
    batched = [pp.preds[0].score for pp in predictor.predict_pairs(pairs)]
    singly = [predictor.predict_pairs([p])[0].preds[0].score for p in pairs]
    np.testing.assert_allclose(batched, singly, atol=1e-6)


@requires_nettcr
def test_predict_models_stay_loaded(predictor):
    predictor.predict_pairs([_row_to_pair(PUBLISHED_ENSEMBLE[0][0])])
    assert predictor._interpreters is not None
    assert "loaded" in str(predictor)


@requires_nettcr
def test_predict_repeated_calls_consistent(predictor):
    pairs = [_row_to_pair(row) for row, _ in PUBLISHED_ENSEMBLE]
    r1 = predictor.predict_pairs(pairs)
    r2 = predictor.predict_pairs(pairs)
    for a, b in zip(r1, r2):
        assert a.preds[0].score == b.preds[0].score


@requires_nettcr
def test_predict_empty_tcrs(predictor):
    results = predictor.predict(["GILGFVFTL"], [])
    assert len(results) == 1
    assert results[0].preds == ()


@requires_nettcr
def test_predict_dataframe_schema(predictor):
    _, tcr = _row_to_pair(PUBLISHED_ENSEMBLE[2][0])
    df = predictor.predict_dataframe(["AVFDRKSDAK"], [tcr])
    assert list(df.columns) == list(COLUMNS)
    assert df["tcr"].iloc[0] == tcr.identifier
    assert df["kind"].iloc[0] == Kind.pMHC_TCR_binding


@requires_nettcr
def test_bad_tcr_type_raises(predictor):
    with pytest.raises(TypeError, match="TCR"):
        predictor.predict_pairs([("SIINFEKL", "not-a-tcr")])
