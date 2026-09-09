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


def test_inferred_transform_is_self_consistent_with_the_forest_extrema():
    # NOT a validation of the transform -- upstream documents no target, and
    # agreement with our own conversion is not external evidence. This pins
    # the arithmetic behind the inference recorded in the module docstring, so
    # that changing `half_life_hours` without revisiting that reasoning fails.
    lowest_leaf_value = 1.30535
    assert half_life_hours(lowest_leaf_value) * 3600.0 == pytest.approx(
        20.2, abs=0.1)
    # The natural model's highest leaf value comes out at exactly seven days.
    assert half_life_hours(5.78161) == pytest.approx(24.0 * 7, rel=1e-4)


def test_log2_and_ln_readings_give_implausible_durations():
    # Guard the reasoning, not just the result. Under log2 or ln seconds the
    # forests' extreme leaf values (1.30535, 6.91424) invert to a range of a
    # couple of seconds up to about two minutes, which no peptide half-life
    # dataset would span. This is what makes log10 the strongest reading; it
    # is not proof, since the training targets themselves are undocumented.
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

def test_parse_results_names_the_derived_column_for_its_assumption():
    path = _write(_OUTPUT)
    try:
        frame = parse_plifepred2_results(
            path, ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    finally:
        os.remove(path)
    assert frame["log10_seconds"].tolist() == [3.157266, 3.793472]
    # Not "hours": the conversion rests on an inferred transform, and the
    # column name has to carry that so a QC frame cannot be read as measured.
    assert "hours" not in frame.columns
    assert frame["hours_if_log10_seconds"].tolist() == pytest.approx(
        [0.39899, 1.72651], rel=1e-4)


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
    for peptide, result, native in zip(
            ["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"], results, raw):
        pred = result.preds[0]
        assert pred.kind == Kind.blood_half_life
        assert pred.peptide == peptide
        assert pred.allele == ""
        # score is the native model output, and by default no duration is
        # claimed at all.
        assert pred.score == pytest.approx(native)
        assert pred.value is None


@requires_plifepred2
def test_end_to_end_scores_follow_their_peptides_when_reordered():
    # Pfeature associates feature rows with inputs by position only, so this
    # is the test that the ordering contract actually holds through it.
    predictor = PlifePred2()
    forward = predictor.predict(["SIINFEKLGGALQAKKY", "GILGFVFTLAAAKKWWWQ"])
    reverse = predictor.predict(["GILGFVFTLAAAKKWWWQ", "SIINFEKLGGALQAKKY"])
    assert forward[0].preds[0].score == pytest.approx(reverse[1].preds[0].score)
    assert forward[1].preds[0].score == pytest.approx(reverse[0].preds[0].score)


@requires_plifepred2
def test_end_to_end_hours_require_explicit_opt_in():
    # The duration is gated because the transform is inferred, not documented.
    default = PlifePred2().predict(["SIINFEKLGGALQAKKY"])[0].preds[0]
    assert default.value is None

    opted_in = PlifePred2(assume_log10_seconds=True)
    pred = opted_in.predict(["SIINFEKLGGALQAKKY"])[0].preds[0]
    assert pred.value == pytest.approx(half_life_hours(default.score))
    assert pred.value == pytest.approx(0.39899, rel=1e-4)
    # score stays the native output either way, so ranking is unaffected.
    assert pred.score == pytest.approx(default.score)


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
    assert results[0].preds[0].score is not None


# --- Pfeature workspace isolation (offline, uses a stub extractor) ----------

# Upstream's real pfeature_comp.py reads Data/ by relative path, scatters
# intermediates into the working directory, and ends with an unscoped
# glob.glob("sam_allcomp*") + os.remove. This stub reproduces exactly those
# three behaviours so the isolation can be tested without a real install.
_DESTRUCTIVE_STUB = '''\
import argparse, glob, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-i"); parser.add_argument("-o"); parser.add_argument("-j")
args = parser.parse_args()

# 1. Data must be reachable by RELATIVE path, as upstream reads it.
open(os.path.join("Data", "Schneider-Wrede.csv")).close()

# 2. Upstream scatters intermediates into the working directory.
for name in ("sam_allcomp.qso_st", "sam_input.csv"):
    with open(name, "w") as handle:
        handle.write("scratch\\n")

# 3. Upstream's unscoped wildcard cleanup.
for path in glob.glob("sam_allcomp*"):
    os.remove(path)

with open(args.o, "w") as handle:
    handle.write("QSO1_SC_A\\n1.0\\n")
'''


def _stub_installation(tmp_path):
    """A Pfeature-shaped installation whose extractor is the stub above."""
    home = tmp_path / "Standalone"
    (home / "Data").mkdir(parents=True)
    for name in ("Schneider-Wrede.csv", "Grantham.csv"):
        (home / "Data" / name).write_text("Name,A\nA,0\n")
    (home / "pfeature_comp.py").write_text(_DESTRUCTIVE_STUB)
    return home


def test_pfeature_cleanup_cannot_delete_files_in_the_installation(tmp_path):
    # Regression for the wildcard cleanup deleting unrelated user files.
    from mhctools.plifepred2_sidecar import _run_pfeature

    home = _stub_installation(tmp_path)
    sentinel = home / "sam_allcomp.review_sentinel"
    sentinel.write_text("do not delete me")

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">0\nSIINFEKLGGAL\n")
    out = tmp_path / "qso.csv"
    _run_pfeature(str(home), fasta, out)

    assert out.exists(), "extraction still has to produce its output"
    assert sentinel.exists(), "upstream cleanup reached into the installation"
    assert sentinel.read_text() == "do not delete me"


def test_pfeature_leaves_no_intermediates_in_the_installation(tmp_path):
    from mhctools.plifepred2_sidecar import _run_pfeature

    home = _stub_installation(tmp_path)
    before = {path.name for path in home.iterdir()}

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">0\nSIINFEKLGGAL\n")
    _run_pfeature(str(home), fasta, tmp_path / "qso.csv")

    assert {path.name for path in home.iterdir()} == before


def test_pfeature_works_with_a_read_only_installation(tmp_path):
    import stat
    from mhctools.plifepred2_sidecar import _run_pfeature

    home = _stub_installation(tmp_path)
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">0\nSIINFEKLGGAL\n")
    out = tmp_path / "qso.csv"

    read_only = stat.S_IRUSR | stat.S_IXUSR
    paths = [home / "Data", home]
    original = [path.stat().st_mode for path in paths]
    for path in paths:
        path.chmod(read_only)
    try:
        _run_pfeature(str(home), fasta, out)
    finally:
        for path, mode in zip(paths, original):
            path.chmod(mode)
    assert out.exists()


def test_concurrent_extractions_do_not_collide_on_one_installation(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from mhctools.plifepred2_sidecar import _run_pfeature

    home = _stub_installation(tmp_path)

    def run(index):
        fasta = tmp_path / ("in_%d.fasta" % index)
        fasta.write_text(">0\nSIINFEKLGGAL\n")
        out = tmp_path / ("qso_%d.csv" % index)
        _run_pfeature(str(home), fasta, out)
        return out.exists()

    with ThreadPoolExecutor(max_workers=4) as pool:
        assert all(pool.map(run, range(4)))


def test_missing_data_directory_is_reported(tmp_path):
    from mhctools.plifepred2_sidecar import _run_pfeature

    home = tmp_path / "Standalone"
    home.mkdir()
    (home / "pfeature_comp.py").write_text(_DESTRUCTIVE_STUB)
    with pytest.raises(RuntimeError, match="no Data directory"):
        _run_pfeature(str(home), tmp_path / "in.fasta", tmp_path / "out.csv")


def test_value_is_withheld_by_default(tmp_path):
    # The units gate is a constructor-level contract, checkable without a model.
    assert _predictor(tmp_path).assume_log10_seconds is False


@pytest.mark.parametrize("via_environment", [False, True])
def test_relative_asset_homes_survive_directory_changes(
        tmp_path, monkeypatch, via_environment):
    from mhctools.plifepred2_sidecar import _run_pfeature

    assets = tmp_path / "model assets"
    assets.mkdir()
    model_home, feature_home = _fake_homes(assets)
    Path(feature_home, "pfeature_comp.py").write_text(_DESTRUCTIVE_STUB)
    (tmp_path / "nested").mkdir()
    monkeypatch.chdir(tmp_path)
    relative_model = "nested/../model assets/plifepred2"
    relative_features = "nested/../model assets/Standalone"
    if via_environment:
        monkeypatch.setenv("PLIFEPRED2_HOME", relative_model)
        monkeypatch.setenv("PFEATURE_HOME", relative_features)
        predictor = PlifePred2()
    else:
        predictor = PlifePred2(
            plifepred2_home=relative_model, pfeature_home=relative_features)
    monkeypatch.chdir(tmp_path / "nested")
    assert predictor.plifepred2_home == str(Path(model_home).resolve())
    assert predictor.pfeature_home == str(Path(feature_home).resolve())
    fasta = tmp_path / "input.fa"
    fasta.write_text(">0\nSIINFEKLGGAL\n")
    output = tmp_path / "output.csv"
    _run_pfeature(predictor.pfeature_home, fasta, output)
    assert "QSO1_SC_A" in output.read_text()


def test_opt_in_is_visible_in_the_repr(tmp_path):
    plifepred2_home, pfeature_home = _fake_homes(tmp_path)
    predictor = PlifePred2(
        plifepred2_home=plifepred2_home,
        pfeature_home=pfeature_home,
        assume_log10_seconds=True)
    assert "assume_log10_seconds=True" in str(predictor)
