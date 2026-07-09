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

"""Tests for the TULIP-TCR wrapper (mhctools/tulip.py).

The wrapper shells out to a user-provided TULIP checkout via a separate
interpreter, so most logic (input-CSV construction, per-peptide output
parsing, order/score mapping, result construction) is tested with the
subprocess mocked -- no torch, transformers, or TULIP install required.

The end-to-end tests near the bottom run only when both TULIP_HOME and
TULIP_PYTHON point at a working install + isolated env; otherwise they skip.
"""

import os
import subprocess

import pandas as pd
import pytest

from mhctools import TCR
from mhctools.pred import COLUMNS, Kind
from mhctools.tulip import (
    Tulip,
    MISSING_TOKEN,
    _find_tulip_home,
    _find_tulip_python,
)


# ---------------------------------------------------------------------------
# Resolution / constructor error paths -- no TULIP install required.
# ---------------------------------------------------------------------------

def test_find_tulip_home_missing_raises():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        _find_tulip_home("/nonexistent/tulip")


def test_find_tulip_home_without_predict_py_raises(tmp_path):
    # Directory exists but isn't a TULIP checkout (no predict.py).
    with pytest.raises(FileNotFoundError, match="no predict.py"):
        _find_tulip_home(str(tmp_path))


def test_find_tulip_python_not_executable_raises(tmp_path):
    not_exe = tmp_path / "python"
    not_exe.write_text("#!/bin/sh\n")   # exists but not chmod +x
    with pytest.raises(FileNotFoundError, match="not an executable"):
        _find_tulip_python(str(not_exe))


def test_clean_maps_empty_to_missing_token():
    assert Tulip._clean("") == MISSING_TOKEN
    assert Tulip._clean(None) == MISSING_TOKEN
    assert Tulip._clean("  CASSIR  ") == "CASSIR"


# ---------------------------------------------------------------------------
# Prediction plumbing with the TULIP subprocess mocked.
# ---------------------------------------------------------------------------

def _stub_tulip(tmp_path):
    """A Tulip whose __init__ filesystem checks are bypassed; only the
    attributes _run_predict reads are set. The subprocess itself is mocked."""
    t = Tulip.__new__(Tulip)
    t.tulip_home = str(tmp_path)
    t.tulip_python = "/usr/bin/python3"
    t.model_config = str(tmp_path / "cfg.json")
    t.weights = str(tmp_path / "w.bin")
    t.batch_size = 512
    return t


def _install_fake_predict(monkeypatch, fail=False):
    """Replace subprocess.run with a stand-in for TULIP's predict.py.

    It reads the input CSV the wrapper wrote, and for each unique peptide
    writes ``<output><peptide>.csv`` in the SAME row order, assigning each row
    a score equal to its 0-based position in the input CSV. That lets tests
    assert the wrapper maps scores back to the right (peptide, mhc, cdr3) key
    purely by position -- exactly how it aligns with real predict.py output.
    """
    def fake_run(cmd, **kwargs):
        test_dir = cmd[cmd.index("--test_dir") + 1]
        output = cmd[cmd.index("--output") + 1]
        if fail:
            return subprocess.CompletedProcess(cmd, 1, "", "boom: TULIP failed")
        df = pd.read_csv(test_dir).reset_index(drop=True)
        df["__score"] = df.index.astype(float)
        for peptide, grp in df.groupby("peptide", sort=False):
            out = pd.DataFrame({
                "CDR3a": grp["CDR3a"].values,
                "CDR3b": grp["CDR3b"].values,
                "peptide": peptide,
                "score": grp["__score"].values,
                "rank": range(len(grp)),
            })
            out.to_csv(output + str(peptide) + ".csv")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")

    monkeypatch.setattr(subprocess, "run", fake_run)


TCR1 = TCR(cdr3a="CAAA", cdr3b="CBBB", name="clone1")
TCR2 = TCR(cdr3a="CAAAAA", cdr3b="CBBBBB", name="clone2")


def test_predict_pairs_maps_scores_by_position(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    pairs = [
        ("GILGFVFTL", TCR1, "HLA-A*02:01"),
        ("GILGFVFTL", TCR2, "HLA-A*02:01"),
        ("NLVPMVATV", TCR1, "HLA-A*02:01"),
    ]
    results = t.predict_pairs(pairs)
    # Score equals the row's position in the (deduped) input CSV: 0, 1, 2.
    assert [r.preds[0].score for r in results] == [0.0, 1.0, 2.0]
    for r, (pep, tcr, mhc) in zip(results, pairs):
        p = r.preds[0]
        assert p.kind == Kind.pMHC_TCR_binding
        assert p.peptide == pep
        assert p.tcr == tcr.identifier
        assert p.allele == mhc
        assert p.predictor_name == "tulip"


def test_predict_pairs_same_peptide_different_mhc_stay_distinct(tmp_path, monkeypatch):
    # Same peptide + same TCR, two MHCs => two distinct keys, mapped by the
    # position they occupy in the per-peptide output (MHC isn't echoed back).
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    pairs = [
        ("GILGFVFTL", TCR1, "HLA-A*02:01"),
        ("GILGFVFTL", TCR1, "HLA-B*07:02"),
    ]
    results = t.predict_pairs(pairs)
    assert [r.preds[0].score for r in results] == [0.0, 1.0]
    assert [r.preds[0].allele for r in results] == ["HLA-A*02:01", "HLA-B*07:02"]


def test_predict_pairs_missing_mhc_gives_empty_allele(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    results = t.predict_pairs([("GILGFVFTL", TCR1)])   # no MHC
    p = results[0].preds[0]
    assert p.allele == ""          # missing token surfaced as empty allele
    # And the input the wrapper wrote used the missing token, not a blank.
    # (verified indirectly: it produced a score, i.e. a well-formed CSV row)
    assert p.score == 0.0


def test_predict_pairs_dedupes_identical_items(tmp_path, monkeypatch):
    # Two identical (peptide, TCR, MHC) items => one subprocess row, and both
    # results get that row's score.
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    pair = ("GILGFVFTL", TCR1, "HLA-A*02:01")
    results = t.predict_pairs([pair, pair])
    assert [r.preds[0].score for r in results] == [0.0, 0.0]


def test_predict_grid(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    results = t.predict(["GILGFVFTL", "NLVPMVATV"], [TCR1, TCR2], mhc="HLA-A*02:01")
    assert len(results) == 2               # one PeptideResult per peptide
    for r in results:
        assert len(r.preds) == 2           # one Prediction per TCR
        assert {p.tcr for p in r.preds} == {"clone1", "clone2"}
        assert all(p.allele == "HLA-A*02:01" for p in r.preds)


def test_predict_grid_parallel_mhc_list(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    results = t.predict(
        ["GILGFVFTL", "NLVPMVATV"], [TCR1], mhc=["HLA-A*02:01", "HLA-B*07:02"])
    assert results[0].preds[0].allele == "HLA-A*02:01"
    assert results[1].preds[0].allele == "HLA-B*07:02"


def test_predict_mhc_list_length_mismatch_raises(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    with pytest.raises(ValueError, match="mhc list length"):
        t.predict(["GILGFVFTL", "NLVPMVATV"], [TCR1], mhc=["HLA-A*02:01"])


def test_predict_failure_raises_with_stderr(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch, fail=True)
    t = _stub_tulip(tmp_path)
    with pytest.raises(RuntimeError, match="boom: TULIP failed"):
        t.predict_pairs([("GILGFVFTL", TCR1, "HLA-A*02:01")])


def test_bad_tcr_type_raises(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    with pytest.raises(TypeError, match="mhctools.TCR"):
        t.predict_pairs([("GILGFVFTL", "not-a-tcr", "HLA-A*02:01")])


def test_predict_dataframe_schema(tmp_path, monkeypatch):
    _install_fake_predict(monkeypatch)
    t = _stub_tulip(tmp_path)
    df = t.predict_dataframe(["GILGFVFTL"], [TCR1, TCR2], mhc="HLA-A*02:01")
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 2
    assert set(df["kind"]) == {Kind.pMHC_TCR_binding}


# ---------------------------------------------------------------------------
# End-to-end -- requires a real TULIP checkout + isolated interpreter.
# ---------------------------------------------------------------------------

TULIP_HOME = os.environ.get("TULIP_HOME")
TULIP_PYTHON = os.environ.get("TULIP_PYTHON")

requires_tulip = pytest.mark.skipif(
    not (TULIP_HOME and os.path.isfile(os.path.join(TULIP_HOME or "", "predict.py"))
         and TULIP_PYTHON and os.access(TULIP_PYTHON or "", os.X_OK)),
    reason="TULIP not installed (set TULIP_HOME + TULIP_PYTHON; "
           "see scripts/setup_tulip_env.sh)")


# Two real TCRs from TULIP's own data/VDJ_test_2.csv; clone1 is a known
# GILGFVFTL (influenza M1) binder.
_E2E_TCR1 = TCR(cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF", name="clone1")
_E2E_TCR2 = TCR(cdr3a="CALSGETSGSRLTF", cdr3b="CASGLVPGGLVYEQYF", name="clone2")


@pytest.fixture(scope="module")
def predictor():
    return Tulip()


@requires_tulip
def test_e2e_predict_grid(predictor):
    results = predictor.predict(
        ["GILGFVFTL", "NLVPMVATV"], [_E2E_TCR1, _E2E_TCR2], mhc="HLA-A*02:01")
    assert len(results) == 2
    for r in results:
        assert len(r.preds) == 2
        for p in r.preds:
            assert p.kind == Kind.pMHC_TCR_binding
            assert p.allele == "HLA-A*02:01"
            assert isinstance(p.score, float)


@requires_tulip
def test_e2e_known_binder_scores_higher(predictor):
    # For the influenza epitope GILGFVFTL, the known binder (clone1) should
    # score higher than an unrelated receptor (clone2).
    results = predictor.predict("GILGFVFTL", [_E2E_TCR1, _E2E_TCR2], mhc="HLA-A*02:01")
    by_tcr = {p.tcr: p.score for p in results[0].preds}
    assert by_tcr["clone1"] > by_tcr["clone2"]


@requires_tulip
def test_e2e_mhc_changes_score(predictor):
    # Supplying the MHC allele should change the score vs. MHC-agnostic.
    with_mhc = predictor.predict_pairs(
        [("GILGFVFTL", _E2E_TCR1, "HLA-A*02:01")])[0].preds[0].score
    without_mhc = predictor.predict_pairs(
        [("GILGFVFTL", _E2E_TCR1)])[0].preds[0].score
    assert with_mhc != without_mhc
