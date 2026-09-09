# Copyright (c) 2026 Mount Sinai School of Medicine
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

"""Tests for the PlifePred2 whole-blood half-life wrapper.

Parser, unit-conversion, kind and validation tests need no model. The
end-to-end tests run only when both ``PLIFEPRED2_HOME`` (an installed
``plifepred2`` package) and ``PFEATURE_HOME`` (Pfeature's ``Standalone``
directory) are set.
"""

import math
import os
from pathlib import Path
import tempfile

import pytest

from mhctools import Kind, PlifePred2
from mhctools.plifepred2 import (
    PLIFEPRED2_MAX_PEPTIDE_LENGTH,
    PLIFEPRED2_MIN_PEPTIDE_LENGTH,
    half_life_hours,
    parse_plifepred2_results,
)
from mhctools.pred import VALUE_BEST_DIRECTIONS, best_direction


_OUTPUT = (
    "__mhctools_id,peptide,log10_seconds\n"
    "0,SIINFEKLGGALQAKKY,3.157266\n"
    "1,GILGFVFTLAAAKKWWWQ,3.793472\n")


def _write(text):
    handle = tempfile.NamedTemporaryFile(
        "w", suffix="_plifepred2.csv", delete=False)
    handle.write(text)
    handle.close()
    return handle.name


# --- units ------------------------------------------------------------------

def test_half_life_hours_inverts_log10_seconds():
    # One hour is 3600 s, so log10(3600) must come back as exactly 1 hour.
    assert half_life_hours(math.log10(3600.0)) == pytest.approx(1.0)
    assert half_life_hours(math.log10(86400.0)) == pytest.approx(24.0)


def test_units_reproduce_the_documented_training_floor():
    # The lineage paper (PLoS ONE 2018) discarded peptides with a half-life
    # below 20 seconds, and the shipped forests' lowest leaf value inverts to
    # that floor under log10-seconds. This is the evidence the transform rests
    # on, so it is pinned here: if someone "fixes" the conversion to log2 or
    # ln, the training range stops making sense.
    lowest_leaf_value = 1.30535
    assert half_life_hours(lowest_leaf_value) * 3600.0 == pytest.approx(
        20.2, abs=0.1)
    # The natural model's highest leaf value is exactly seven days.
    assert half_life_hours(5.78161) == pytest.approx(24.0 * 7, rel=1e-4)


def test_log2_and_ln_readings_contradict_the_documented_dataset():
    # Guard the reasoning, not just the result. The paper's dataset spans
    # 20 seconds to at least 24 hours, and the shipped forests' extreme leaf
    # values are 1.30535 and 6.91424. Under log2 or ln seconds the whole
    # training range collapses to a couple of seconds up to ~2-17 minutes:
    # the floor falls under the documented 20 s and the ceiling nowhere near
    # the documented 24 h. Only log10 spans the documented range.
    lowest, highest = 1.30535, 6.91424
    assert 2.0 ** lowest < 20.0                  # below the documented floor
    assert 2.0 ** highest < 24 * 3600.0          # far below the documented cap
    assert math.exp(lowest) < 20.0
    assert math.exp(highest) < 24 * 3600.0
    # log10 covers it: floor at 20 s, ceiling well past 24 h.
    assert 10.0 ** lowest == pytest.approx(20.2, abs=0.1)
    assert 10.0 ** highest > 24 * 3600.0


# --- kind semantics ---------------------------------------------------------

def test_blood_half_life_is_distinct_from_serum_and_pmhc_stability():
    assert Kind.blood_half_life == "blood_half_life"
    assert Kind.blood_half_life != Kind.serum_half_life
    assert Kind.blood_half_life != Kind.pMHC_stability


def test_blood_half_life_value_direction_is_max():
    assert VALUE_BEST_DIRECTIONS[Kind.blood_half_life] == "max"
    assert best_direction(Kind.blood_half_life, "value") == "max"


def test_annotate_exposes_blood_half_life_separately_from_serum():
    from mhctools.annotate import _OUTPUT_FIELDS, output_field_tokens
    assert "blood_half_life" in output_field_tokens()
    assert _OUTPUT_FIELDS["blood_half_life"] == (Kind.blood_half_life, "value")
    assert _OUTPUT_FIELDS["serum_half_life"] != _OUTPUT_FIELDS["blood_half_life"]


def test_accessors_do_not_mix_matrices():
    from mhctools.pred import PeptideResult, Prediction
    result = PeptideResult(preds=(
        Prediction(kind=Kind.blood_half_life, score=3.0, value=3.0,
                   peptide="SIINFEKLGGALQAKKY"),
        Prediction(kind=Kind.serum_half_life, score=8.0, value=8.0,
                   peptide="SIINFEKLGGALQAKKY"),
    ))
    assert result.blood_half_life.value == 3.0
    assert result.serum_half_life.value == 8.0


# --- parser -----------------------------------------------------------------

def test_parse_results_converts_to_hours():
    path = _write(_OUTPUT)
    try:
        frame = parse_plifepred2_results(
            path, ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    finally:
        os.remove(path)
    assert frame["log10_seconds"].tolist() == [3.157266, 3.793472]
    assert frame["hours"].tolist() == pytest.approx([0.39899, 1.72651], rel=1e-4)


def test_parse_results_restores_input_order():
    shuffled = (
        "__mhctools_id,peptide,log10_seconds\n"
        "1,GILGFVFTLAAAKKWWWQ,3.793472\n"
        "0,SIINFEKLGGALQAKKY,3.157266\n")
    path = _write(shuffled)
    try:
        frame = parse_plifepred2_results(
            path, ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    finally:
        os.remove(path)
    assert frame["peptide"].tolist() == [
        "SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"]
    assert frame["log10_seconds"].tolist() == [3.157266, 3.793472]


def test_parse_results_rejects_shifted_rows():
    # Pfeature links features to inputs by position only, so a shifted row set
    # would silently attach each score to the wrong peptide.
    shifted = (
        "__mhctools_id,peptide,log10_seconds\n"
        "0,GILGFVFTLAAAKKWWWQ,3.793472\n"
        "1,SIINFEKLGGALQAKKY,3.157266\n")
    path = _write(shifted)
    try:
        with pytest.raises(RuntimeError, match="different peptide"):
            parse_plifepred2_results(
                path, ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    finally:
        os.remove(path)


def test_parse_results_rejects_dropped_peptide():
    path = _write("__mhctools_id,peptide,log10_seconds\n0,SIINFEKLGGALQAKKY,3.1\n")
    try:
        with pytest.raises(RuntimeError, match="did not preserve"):
            parse_plifepred2_results(
                path, ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    finally:
        os.remove(path)


def test_parse_results_rejects_non_finite_prediction():
    path = _write("__mhctools_id,peptide,log10_seconds\n0,SIINFEKLGGALQAKKY,nan\n")
    try:
        with pytest.raises(RuntimeError, match="non-finite"):
            parse_plifepred2_results(path, ["SIINFEKLGGALQAKKY"])
    finally:
        os.remove(path)


def test_parse_results_rejects_missing_column():
    path = _write("__mhctools_id,peptide\n0,SIINFEKLGGALQAKKY\n")
    try:
        with pytest.raises(ValueError, match="missing column 'log10_seconds'"):
            parse_plifepred2_results(path, ["SIINFEKLGGALQAKKY"])
    finally:
        os.remove(path)


# --- construction and validation --------------------------------------------

def _fake_homes(tmp_path):
    plifepred2 = tmp_path / "plifepred2"
    (plifepred2 / "models").mkdir(parents=True)
    (plifepred2 / "models" / "plifepred2_natural_model.sav").write_text("")
    pfeature = tmp_path / "Standalone"
    (pfeature / "Data").mkdir(parents=True)
    (pfeature / "pfeature_comp.py").write_text("")
    for name in ("Schneider-Wrede.csv", "Grantham.csv"):
        (pfeature / "Data" / name).write_text("")
    return str(plifepred2), str(pfeature)


def test_missing_plifepred2_home_is_reported(tmp_path, monkeypatch):
    monkeypatch.delenv("PLIFEPRED2_HOME", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    with pytest.raises(FileNotFoundError, match="PlifePred2 not found"):
        PlifePred2()


def test_missing_model_file_is_reported(tmp_path, monkeypatch):
    monkeypatch.delenv("PLIFEPRED2_HOME", raising=False)
    with pytest.raises(FileNotFoundError, match="natural_model.sav not found"):
        PlifePred2(plifepred2_home=str(tmp_path))


def test_missing_pfeature_home_is_reported(tmp_path, monkeypatch):
    plifepred2_home, _ = _fake_homes(tmp_path)
    monkeypatch.delenv("PFEATURE_HOME", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    with pytest.raises(FileNotFoundError, match="Pfeature not found"):
        PlifePred2(plifepred2_home=plifepred2_home)


def test_pfeature_home_without_distance_matrices_is_reported(tmp_path):
    plifepred2_home, pfeature_home = _fake_homes(tmp_path)
    os.remove(os.path.join(pfeature_home, "Data", "Grantham.csv"))
    with pytest.raises(FileNotFoundError, match="Data/Grantham.csv not found"):
        PlifePred2(
            plifepred2_home=plifepred2_home, pfeature_home=pfeature_home)


def _predictor(tmp_path):
    plifepred2_home, pfeature_home = _fake_homes(tmp_path)
    return PlifePred2(
        plifepred2_home=plifepred2_home, pfeature_home=pfeature_home)


def test_supported_kinds_and_mhc_context(tmp_path):
    predictor = _predictor(tmp_path)
    assert predictor.supported_kinds == (Kind.blood_half_life,)
    support = predictor.kind_support()[Kind.blood_half_life]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "none"


def test_empty_peptide_list_returns_nothing(tmp_path):
    assert _predictor(tmp_path).predict([]) == []


def test_short_peptide_is_rejected_rather_than_silently_dropped(tmp_path):
    # Upstream's CLI writes these to eliminated_sequences.csv and carries on
    # with a shorter result set; that must not reach a caller unannounced.
    with pytest.raises(ValueError, match="12-100 residues"):
        _predictor(tmp_path).predict(["SIINFEKL"])


def test_long_peptide_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="12-100 residues"):
        _predictor(tmp_path).predict(["A" * (PLIFEPRED2_MAX_PEPTIDE_LENGTH + 1)])


def test_modified_residues_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="non-standard residues"):
        _predictor(tmp_path).predict(["SIINFEKLGGALQAKKB"])


def test_empty_peptide_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="Empty peptide"):
        _predictor(tmp_path).predict([""])


def test_missing_python_is_reported(tmp_path):
    plifepred2_home, pfeature_home = _fake_homes(tmp_path)
    with pytest.raises(FileNotFoundError, match="Python does not exist"):
        PlifePred2(
            plifepred2_home=plifepred2_home,
            pfeature_home=pfeature_home,
            plifepred2_python="/nonexistent/python")


def test_predictor_version_names_the_natural_model(tmp_path):
    assert _predictor(tmp_path).predictor_version == "1.0:natural"


# --- end-to-end (opt-in) ----------------------------------------------------

requires_plifepred2 = pytest.mark.skipif(
    not (os.environ.get("PLIFEPRED2_HOME") and os.environ.get("PFEATURE_HOME")),
    reason="set PLIFEPRED2_HOME and PFEATURE_HOME to run this")


@requires_plifepred2
def test_end_to_end_matches_the_reference_values():
    # Computed by running Pfeature's QSO extractor and the shipped forest
    # directly, outside this wrapper.
    predictor = PlifePred2()
    results = predictor.predict(["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    assert len(results) == 2
    raw = predictor.last_qc["log10_seconds"].tolist()
    assert raw == pytest.approx([3.1572664, 3.79347232], rel=1e-6)
    for peptide, result in zip(
            ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"], results):
        pred = result.preds[0]
        assert pred.kind == Kind.blood_half_life
        assert pred.peptide == peptide
        assert pred.allele == ""
        assert pred.value > 0
        assert pred.value == pred.score


@requires_plifepred2
def test_end_to_end_scores_follow_their_peptides_when_reordered():
    # Pfeature associates feature rows with inputs by position only, so this
    # is the test that the ordering contract actually holds through it.
    predictor = PlifePred2()
    forward = predictor.predict(["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    reverse = predictor.predict(["GILGFVFTLAAAKKWWWQ", "SIINFEKLGGALQAKKY"])
    assert forward[0].preds[0].value == pytest.approx(reverse[1].preds[0].value)
    assert forward[1].preds[0].value == pytest.approx(reverse[0].preds[0].value)


@requires_plifepred2
def test_end_to_end_dataframe_has_the_standard_columns():
    from mhctools.pred import COLUMNS
    frame = PlifePred2().predict_dataframe(["SIINFEKLGGALQAKKY"])
    assert list(frame.columns) == list(COLUMNS)
    assert frame["kind"].tolist() == [Kind.blood_half_life]


@requires_plifepred2
def test_end_to_end_minimum_peptide_length_is_accepted():
    peptide = "SIINFEKLGGAL"
    assert len(peptide) == PLIFEPRED2_MIN_PEPTIDE_LENGTH
    results = PlifePred2().predict([peptide])
    assert results[0].preds[0].value > 0
