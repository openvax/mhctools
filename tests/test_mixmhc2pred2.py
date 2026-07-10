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

"""Tests for the MixMHC2pred class-II presentation wrapper.

Parser and allele-mapping tests are binary-free. The end-to-end tests run only
when a MixMHC2pred executable is available (``MIXMHC2PRED_EXECUTABLE`` env var
or ``MixMHC2pred`` / ``MixMHC2pred_unix`` on ``PATH``).
"""

import os
import shutil
from tempfile import NamedTemporaryFile

import pytest

from mhctools import MixMHC2pred
from mhctools.mixmhc2pred import parse_mixmhc2pred_results, to_mixmhc2pred_allele
from mhctools.pred import Kind


# Captured MixMHC2pred (v2.0.2, --extra_out) output for two alleles. Note the
# columns are NOT grouped per allele: every allele's %Rank/CoreP1/SubSpec come
# first, then every allele's Score/%RankPWM/ScorePWM — which is why the parser
# must look columns up by name, not by position.
_MM2_OUTPUT = (
    "####################\n"
    "# Output from MixMHC2pred (v2.0.2)\n"
    "# Alleles: DRB1_15_01, DQA1_01_02__DQB1_06_02\n"
    "####################\n"
    "Peptide\tContext\tBestAllele\t%Rank_best\tCore_best\tCoreP1_best\t"
    "SubSpec_best\t%Rank_DRB1_15_01\tCoreP1_DRB1_15_01\tSubSpec_DRB1_15_01\t"
    "%Rank_DQA1_01_02__DQB1_06_02\tCoreP1_DQA1_01_02__DQB1_06_02\t"
    "SubSpec_DQA1_01_02__DQB1_06_02\tScore_DRB1_15_01\t%RankPWM_DRB1_15_01\t"
    "ScorePWM_DRB1_15_01\tScore_DQA1_01_02__DQB1_06_02\t"
    "%RankPWM_DQA1_01_02__DQB1_06_02\tScorePWM_DQA1_01_02__DQB1_06_02\n"
    "PKYVKQNTLKLAT\t\tDQA1_01_02__DQB1_06_02\t11.3\tVKQNTLKLA\t4\t1\t74.9\t4\t1\t"
    "11.3\t4\t1\t0.0791\t43.8\t0.00912\t0.218\t6.45\t2.26\n"
    "GELIGTLNAAKVPAD\t\tDRB1_15_01\t2.52\tIGTLNAAKV\t4\t1\t2.52\t4\t1\t"
    "18.2\t7\t1\t0.408\t3.07\t3.25\t0.166\t10.6\t0.895\n"
)


def _write_output():
    f = NamedTemporaryFile("w", suffix="_mm2.txt", delete=False)
    f.write(_MM2_OUTPUT)
    f.close()
    return f.name


# --- allele nomenclature mapping --------------------------------------------

@pytest.mark.parametrize("allele,expected", [
    ("HLA-DRB1*15:01", "DRB1_15_01"),
    ("HLA-DRA1*01:01-DRB1*15:01", "DRB1_15_01"),   # mhcgnomes-canonical
    ("DRB1_15_01", "DRB1_15_01"),
    ("HLA-DQA1*01:02-DQB1*06:02", "DQA1_01_02__DQB1_06_02"),
    ("DQA1_01_02__DQB1_06_02", "DQA1_01_02__DQB1_06_02"),
    ("DPA1_02_01__DPB1_01_01", "DPA1_02_01__DPB1_01_01"),
])
def test_to_mixmhc2pred_allele(allele, expected):
    assert to_mixmhc2pred_allele(allele) == expected


# --- parser (binary-free) ---------------------------------------------------

def test_parse_mixmhc2pred_shape():
    alleles = ["HLA-DRA1*01:01-DRB1*15:01", "HLA-DQA1*01:02-DQB1*06:02"]
    cli_names = ["DRB1_15_01", "DQA1_01_02__DQB1_06_02"]
    path = _write_output()
    try:
        by_peptide = parse_mixmhc2pred_results(path, alleles, cli_names)
    finally:
        os.remove(path)

    assert set(by_peptide) == {"PKYVKQNTLKLAT", "GELIGTLNAAKVPAD"}
    for peptide, preds in by_peptide.items():
        assert len(preds) == 2
        assert {p.allele for p in preds} == set(alleles)
        for p in preds:
            assert p.kind == Kind.pMHC_presentation
            assert p.predictor_name == "mixmhc2pred"


def test_parse_mixmhc2pred_values_by_name_not_position():
    # DRB1_15_01 columns come before DQA... for %Rank but the Score_* block is
    # ordered separately; name-based lookup must still pair them correctly.
    alleles = ["HLA-DRA1*01:01-DRB1*15:01", "HLA-DQA1*01:02-DQB1*06:02"]
    cli_names = ["DRB1_15_01", "DQA1_01_02__DQB1_06_02"]
    path = _write_output()
    try:
        by_peptide = parse_mixmhc2pred_results(path, alleles, cli_names)
    finally:
        os.remove(path)

    by_allele = {p.allele: p for p in by_peptide["PKYVKQNTLKLAT"]}
    drb = by_allele["HLA-DRA1*01:01-DRB1*15:01"]
    dq = by_allele["HLA-DQA1*01:02-DQB1*06:02"]
    assert drb.score == 0.0791 and drb.percentile_rank == 74.9
    assert dq.score == 0.218 and dq.percentile_rank == 11.3


def test_parse_mixmhc2pred_missing_column_raises():
    path = _write_output()
    try:
        with pytest.raises(ValueError, match="missing column"):
            parse_mixmhc2pred_results(path, ["A_1"], ["NOT_A_REAL_ALLELE"])
    finally:
        os.remove(path)


def test_mixmhc2pred_default_kind_and_support():
    predictor = MixMHC2pred(alleles=["DRB1_15_01"], program_name="MixMHC2pred")
    assert predictor._default_pred_kind() == Kind.pMHC_presentation
    assert predictor.mhc_class == "II"
    support = predictor.kind_support()[Kind.pMHC_presentation]
    assert support["mhc_dependence"] == "single_allele"
    assert support["mhc_class"] == "II"


# --- end-to-end (requires a MixMHC2pred install) ----------------------------

MIXMHC2PRED = (
    os.environ.get("MIXMHC2PRED_EXECUTABLE")
    or shutil.which("MixMHC2pred_unix")
    or shutil.which("MixMHC2pred"))

requires_mixmhc2pred = pytest.mark.skipif(
    not MIXMHC2PRED,
    reason="MixMHC2pred not installed (set MIXMHC2PRED_EXECUTABLE or put "
           "MixMHC2pred / MixMHC2pred_unix on PATH)")


@requires_mixmhc2pred
def test_mixmhc2pred_end_to_end():
    predictor = MixMHC2pred(
        alleles=["DRB1_15_01", "HLA-DQA1*01:02-DQB1*06:02"],
        program_name=MIXMHC2PRED)
    peptides = ["PKYVKQNTLKLAT", "GELIGTLNAAKVPAD"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 2
        for pred in result.preds:
            assert pred.kind == Kind.pMHC_presentation
            assert pred.score is not None
            assert pred.percentile_rank is not None
        assert result.presentation is not None


def test_predict_dataframe_inherits_flank_signature():
    import inspect
    params = inspect.signature(MixMHC2pred.predict_dataframe).parameters
    assert "n_flanks" in params and "c_flanks" in params
