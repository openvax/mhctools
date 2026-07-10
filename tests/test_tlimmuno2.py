# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

"""Tests for the TLimmuno2 class-II immunogenicity wrapper.

Parser, allele-mapping, and validation tests are binary-free. The end-to-end
tests run only when a TLimmuno2 checkout is available (``TLIMMUNO2_HOME``
pointing at a clone, with ``TLIMMUNO2_PYTHON`` optionally naming an interpreter
that has TensorFlow with Keras 2, or newer TensorFlow plus the ``tf-keras``
shim). They are slow: TLimmuno2 scores ~90k background peptides per allele.
"""

import os
from os.path import isfile, join
from tempfile import NamedTemporaryFile

import pytest

from mhctools import TLimmuno2, Kind
from mhctools.tlimmuno2 import _tlimmuno2_allele, parse_tlimmuno2_results


# Captured TLimmuno2 result.csv (DRB1_0803, two 15-mers).
_OUTPUT = (
    ",pep,HLA,sequence,prediction,Rank\n"
    "0,FHTMWHVTRGAVLMY,DRB1_0803,QEFFIASGAA,0.98740566,0.001933311852090558\n"
    "1,GLLFRRLTSREVLLL,DRB1_0803,QEFFIASGAA,0.9916319,0.001033321851979463\n")

_VALID = {"DRB1_0803", "DRB3_0101", "H-2-IAb",
          "HLA-DPA10103-DPB10101", "HLA-DQA10501-DQB10201"}


def _write(text):
    f = NamedTemporaryFile("w", suffix="_tlimmuno2.csv", delete=False)
    f.write(text)
    f.close()
    return f.name


# --- parser (binary-free) ---------------------------------------------------

def test_parse_tlimmuno2_results():
    path = _write(_OUTPUT)
    try:
        by_pair = parse_tlimmuno2_results(path)
    finally:
        os.remove(path)

    assert set(by_pair) == {
        ("FHTMWHVTRGAVLMY", "DRB1_0803"),
        ("GLLFRRLTSREVLLL", "DRB1_0803"),
    }
    score, pct_rank = by_pair[("FHTMWHVTRGAVLMY", "DRB1_0803")]
    assert score == pytest.approx(0.98740566)
    # Rank is a 0-1 fraction rescaled to a 0-100 percentile (lower = stronger).
    assert pct_rank == pytest.approx(0.1933311852)


def test_parse_tlimmuno2_rejects_missing_column():
    path = _write(",pep,HLA\n0,FHTMWHVTRGAVLMY,DRB1_0803\n")
    try:
        with pytest.raises(ValueError, match="prediction"):
            parse_tlimmuno2_results(path)
    finally:
        os.remove(path)


# --- allele mapping (binary-free) -------------------------------------------

def test_allele_mapping_native_and_dr_forms():
    # Native keys pass through.
    assert _tlimmuno2_allele("DRB1_0803", _VALID) == "DRB1_0803"
    assert _tlimmuno2_allele("H-2-IAb", _VALID) == "H-2-IAb"
    # Common DR short forms convert.
    assert _tlimmuno2_allele("HLA-DRB1*08:03", _VALID) == "DRB1_0803"
    assert _tlimmuno2_allele("DRB1*08:03", _VALID) == "DRB1_0803"
    # mhctools canonicalizes DR to an alpha/beta pair; keep the beta chain.
    assert _tlimmuno2_allele(
        "HLA-DRA1*01:01-DRB1*08:03", _VALID) == "DRB1_0803"
    # DP/DQ heterodimers: strip punctuation.
    assert _tlimmuno2_allele(
        "HLA-DPA1*01:03-DPB1*01:01", _VALID) == "HLA-DPA10103-DPB10101"


def test_allele_mapping_rejects_unknown():
    with pytest.raises(ValueError, match="does not recognize"):
        _tlimmuno2_allele("DRB1_9999", _VALID)


# --- construction / validation (no binary needed) ---------------------------

def _stub_home(tmp_path):
    """A minimal fake TLimmuno2 checkout (just a Python/TLimmuno2.py marker)."""
    (tmp_path / "Python").mkdir()
    (tmp_path / "Python" / "TLimmuno2.py").write_text("# stub\n")
    return str(tmp_path)


def test_default_kind_and_support(tmp_path):
    predictor = TLimmuno2(
        alleles=["DRB1_0803"], tlimmuno2_home=_stub_home(tmp_path))
    assert predictor._default_pred_kind() == Kind.immunogenicity
    support = predictor.kind_support()[Kind.immunogenicity]
    assert support["mhc_dependence"] == "single_allele"
    assert support["mhc_class"] == "II"
    assert predictor.supported_kinds == (Kind.immunogenicity,)


def test_rejects_missing_home(tmp_path):
    with pytest.raises(FileNotFoundError, match="TLimmuno2.py"):
        TLimmuno2(alleles=["DRB1_0803"], tlimmuno2_home=str(tmp_path))


def test_rejects_no_alleles(tmp_path):
    with pytest.raises(ValueError, match="at least one allele"):
        TLimmuno2(alleles=[], tlimmuno2_home=_stub_home(tmp_path))


def test_peptide_validation(tmp_path):
    predictor = TLimmuno2(
        alleles=["DRB1_0803"], tlimmuno2_home=_stub_home(tmp_path))
    # These raise before any allele lookup or subprocess.
    with pytest.raises(ValueError, match="up to 21"):
        predictor.predict(["A" * 22])
    with pytest.raises(ValueError, match="Empty peptide"):
        predictor.predict([""])
    assert predictor.predict([]) == []


# --- end-to-end (requires a TLimmuno2 checkout) -----------------------------

TLIMMUNO2_HOME = os.environ.get("TLIMMUNO2_HOME")
_has_tlimmuno2 = bool(TLIMMUNO2_HOME) and isfile(
    join(TLIMMUNO2_HOME, "Python", "TLimmuno2.py"))

requires_tlimmuno2 = pytest.mark.skipif(
    not _has_tlimmuno2,
    reason="TLimmuno2 not installed (set TLIMMUNO2_HOME to a clone; optionally "
           "TLIMMUNO2_PYTHON to an interpreter with TensorFlow and Keras 2 / "
           "tf-keras)")


@requires_tlimmuno2
def test_tlimmuno2_end_to_end():
    predictor = TLimmuno2(
        alleles=["DRB1_0803"], tlimmuno2_home=TLIMMUNO2_HOME)
    peptides = ["FHTMWHVTRGAVLMY", "GLLFRRLTSREVLLL"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.immunogenicity
        assert pred.allele == "DRB1_0803"
        assert pred.predictor_name == "tlimmuno2"
        assert 0.0 <= pred.score <= 1.0
        assert pred.percentile_rank is not None
        assert result.immunogenicity is pred

    # TLimmuno2 is deterministic; the score matches a local reference run.
    by_peptide = {r.peptide: r.preds[0].score for r in results}
    assert by_peptide["FHTMWHVTRGAVLMY"] == pytest.approx(0.9874, abs=1e-3)
