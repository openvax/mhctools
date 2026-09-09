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

"""Tests for the PeptiVerse serum half-life wrapper.

Parser, manifest, kind and validation tests need no model and no network. The
end-to-end test runs only when ``PEPTIVERSE_HOME`` points at a snapshot with its
half-life weights pulled, optionally with ``PEPTIVERSE_PYTHON`` naming an
interpreter that has torch + transformers.
"""

import os
from io import StringIO
from pathlib import Path
import sys
import tempfile

import pytest

from mhctools import Kind, PeptiVerse
from mhctools.peptiverse import (
    PEPTIVERSE_MAX_PEPTIDE_LENGTH,
    _HALF_LIFE_MANIFEST,
    parse_peptiverse_results,
)
from mhctools.pred import VALUE_BEST_DIRECTIONS, best_direction


# Sidecar output for three peptides, in the shape peptiverse_sidecar.py writes.
_OUTPUT = (
    "__mhctools_id,peptide,hours,emb_tag,uncertainty,uncertainty_type\n"
    "0,SIINFEKL,1.4832,wt,,\n"
    "1,GILGFVFTL,0.6271,wt,,\n"
    "2,KLGGALQAK,12.9044,wt,,\n")


def _write(text):
    handle = tempfile.NamedTemporaryFile(
        "w", suffix="_peptiverse.csv", delete=False)
    handle.write(text)
    handle.close()
    return handle.name


# --- kind semantics (no model, no network) ----------------------------------

def test_serum_half_life_is_not_pmhc_stability():
    # The whole point of the new kind: a serum half-life must never land in the
    # field NetMHCstabpan writes, which is pMHC complex dissociation.
    assert Kind.serum_half_life != Kind.pMHC_stability
    assert Kind.serum_half_life == "serum_half_life"


def test_serum_half_life_value_direction_is_max():
    # Hours in serum: longer-lived is "better", same as pMHC stability but for
    # an unrelated reason, so it needs its own registered direction.
    assert VALUE_BEST_DIRECTIONS[Kind.serum_half_life] == "max"
    assert best_direction(Kind.serum_half_life, "value") == "max"
    assert best_direction(Kind.serum_half_life, "score") == "max"


def test_annotate_exposes_serum_half_life_as_a_units_bearing_field():
    from mhctools.annotate import _OUTPUT_FIELDS, output_field_tokens
    assert "serum_half_life" in output_field_tokens()
    # `value` and not `score`, because hours are the meaningful quantity.
    assert _OUTPUT_FIELDS["serum_half_life"] == (Kind.serum_half_life, "value")


def test_peptide_result_accessor():
    from mhctools.pred import PeptideResult, Prediction
    result = PeptideResult(preds=(
        Prediction(kind=Kind.serum_half_life, score=2.0, value=2.0,
                   peptide="SIINFEKL"),
        Prediction(kind=Kind.serum_half_life, score=9.0, value=9.0,
                   peptide="SIINFEKL"),
        Prediction(kind=Kind.pMHC_stability, score=99.0, peptide="SIINFEKL"),
    ))
    assert result.serum_half_life.value == 9.0
    # The pMHC stability prediction must not leak into the serum accessor.
    assert result.stability.score == 99.0


# --- manifest ---------------------------------------------------------------

def test_manifest_selects_only_the_half_life_sequence_model():
    lines = [
        line for line in _HALF_LIFE_MANIFEST.strip().split("\n") if line.strip()
    ]
    assert len(lines) == 2, "header plus exactly one property row"
    header, row = lines
    assert "Best_Model_WT" in header and "Best_Model_SMILES" in header
    fields = [field.strip() for field in row.split(",")]
    assert fields[0] == "Halflife"
    # "Transformer" resolves upstream to transformer_wt_log, the only half-life
    # variant whose output is expm1'd into hours.
    assert fields[1] == "Transformer"
    # SMILES model blanked: its output is not on the hours scale.
    assert fields[2] == "-"
    assert fields[3] == "Regression"


# --- parser (no model) ------------------------------------------------------

def test_parse_results():
    path = _write(_OUTPUT)
    try:
        frame = parse_peptiverse_results(
            path, ["SIINFEKL", "GILGFVFTL", "KLGGALQAK"])
    finally:
        os.remove(path)
    assert frame["hours"].tolist() == [1.4832, 0.6271, 12.9044]
    assert frame["peptide"].tolist() == ["SIINFEKL", "GILGFVFTL", "KLGGALQAK"]


def test_parse_results_restores_input_order():
    shuffled = (
        "__mhctools_id,peptide,hours,emb_tag,uncertainty,uncertainty_type\n"
        "2,KLGGALQAK,12.9044,wt,,\n"
        "0,SIINFEKL,1.4832,wt,,\n"
        "1,GILGFVFTL,0.6271,wt,,\n")
    path = _write(shuffled)
    try:
        frame = parse_peptiverse_results(
            path, ["SIINFEKL", "GILGFVFTL", "KLGGALQAK"])
    finally:
        os.remove(path)
    assert frame["peptide"].tolist() == ["SIINFEKL", "GILGFVFTL", "KLGGALQAK"]
    assert frame["hours"].tolist() == [1.4832, 0.6271, 12.9044]


def test_parse_results_keeps_duplicate_peptides_separate():
    duplicated = (
        "__mhctools_id,peptide,hours,emb_tag,uncertainty,uncertainty_type\n"
        "0,SIINFEKL,1.4832,wt,,\n"
        "1,SIINFEKL,1.4832,wt,,\n")
    path = _write(duplicated)
    try:
        frame = parse_peptiverse_results(path, ["SIINFEKL", "SIINFEKL"])
    finally:
        os.remove(path)
    assert len(frame) == 2


def test_parse_results_rejects_dropped_peptide():
    truncated = (
        "__mhctools_id,peptide,hours,emb_tag,uncertainty,uncertainty_type\n"
        "0,SIINFEKL,1.4832,wt,,\n")
    path = _write(truncated)
    try:
        with pytest.raises(RuntimeError, match="did not preserve"):
            parse_peptiverse_results(path, ["SIINFEKL", "GILGFVFTL"])
    finally:
        os.remove(path)


def test_parse_results_rejects_altered_peptide():
    # Upstream's embedders truncate over-long input rather than failing, so a
    # peptide coming back changed means the score is for a different molecule.
    altered = (
        "__mhctools_id,peptide,hours,emb_tag,uncertainty,uncertainty_type\n"
        "0,SIINFEK,1.4832,wt,,\n")
    path = _write(altered)
    try:
        with pytest.raises(RuntimeError, match="different peptide"):
            parse_peptiverse_results(path, ["SIINFEKL"])
    finally:
        os.remove(path)


def test_parse_results_rejects_missing_column():
    path = _write("__mhctools_id,peptide\n0,SIINFEKL\n")
    try:
        with pytest.raises(ValueError, match="missing column 'hours'"):
            parse_peptiverse_results(path, ["SIINFEKL"])
    finally:
        os.remove(path)


@pytest.mark.parametrize("hours", ["nan", "inf", "-inf", "-1.0"])
def test_parse_rejects_invalid_duration_with_row_context(hours):
    data = StringIO("__mhctools_id,peptide,hours\n0,SIINFEKL,%s\n" % hours)
    with pytest.raises(RuntimeError, match="invalid half-lives.*SIINFEKL"):
        parse_peptiverse_results(data, ["SIINFEKL"])


def test_zero_duration_is_preserved():
    data = StringIO("__mhctools_id,peptide,hours\n0,SIINFEKL,0\n")
    assert parse_peptiverse_results(data, ["SIINFEKL"])["hours"].tolist() == [0.0]


def _stub_inference(tmp_path, hours=1.25, retained_limit=1020):
    """Run the actual sidecar against a dependency-free upstream-shaped stub."""
    home = tmp_path / "Pepti Verse"
    home.mkdir()
    _fake_home(home)
    (home / "inference.py").write_text(f'''
from pathlib import Path

class Mask:
    def __init__(self, n): self.n = n
    def sum(self): return self
    def item(self): return self.n

class Embedder:
    def _tokenize(self, peptides):
        return {{"input_ids": peptides[0], "attention_mask": None}}
    def _valid_mask(self, ids, mask):
        return Mask(min(len(ids), {retained_limit}))

class PeptiVersePredictor:
    def __init__(self, manifest_path, classifier_weight_root, device):
        assert Path(classifier_weight_root) == Path(__file__).resolve().parent
        self.wt_embedder = Embedder()
    def predict_property(self, prop, col, peptide, uncertainty):
        return {{"score": {hours!r}, "emb_tag": "wt"}}
''')
    return home


@pytest.mark.parametrize("via_environment", [False, True])
def test_relative_home_survives_caller_and_sidecar_directory_changes(
        tmp_path, monkeypatch, via_environment):
    home = _stub_inference(tmp_path)
    (tmp_path / "nested").mkdir()
    relative = "nested/../Pepti Verse"
    monkeypatch.chdir(tmp_path)
    if via_environment:
        monkeypatch.setenv("PEPTIVERSE_HOME", relative)
        predictor = PeptiVerse(peptiverse_python=sys.executable)
    else:
        predictor = PeptiVerse(
            peptiverse_home=relative, peptiverse_python=sys.executable)
    monkeypatch.chdir(tmp_path / "nested")
    assert predictor.peptiverse_home == str(home.resolve())
    peptides = ["A" * 1019 + "K", "A" * 1019 + "R"]
    results = predictor.predict(peptides)
    assert [r.peptide for r in results] == peptides
    assert [r.serum_half_life.value for r in results] == [1.25, 1.25]


@pytest.mark.parametrize("hours", ["nan", "inf", "-inf", -1.0])
def test_public_predict_rejects_invalid_sidecar_duration(tmp_path, hours):
    home = _stub_inference(tmp_path, hours=hours)
    predictor = PeptiVerse(peptiverse_home=home, peptiverse_python=sys.executable)
    with pytest.raises(RuntimeError, match="invalid half-lives"):
        predictor.predict(["SIINFEKL"])


def test_sidecar_checks_actual_tokenizer_residue_count(tmp_path):
    home = _stub_inference(tmp_path, retained_limit=7)
    predictor = PeptiVerse(peptiverse_home=home, peptiverse_python=sys.executable)
    with pytest.raises(RuntimeError, match="retained 7 of 8 residues"):
        predictor.predict(["SIINFEKL"])


@pytest.mark.parametrize("length", [1021, 1022])
def test_rejects_lengths_that_fit_only_without_special_tokens(tmp_path, length):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    with pytest.raises(ValueError, match="up to 1020 residues"):
        predictor.predict(["A" * (length - 1) + "K"])


# --- construction and validation (no model) ---------------------------------

def _fake_home(tmp_path):
    """A directory shaped enough like a snapshot to pass resolution."""
    (tmp_path / "inference.py").write_text("")
    (tmp_path / "training_classifiers" / "half_life").mkdir(parents=True)
    return str(tmp_path)


def test_missing_home_is_reported(tmp_path, monkeypatch):
    monkeypatch.delenv("PEPTIVERSE_HOME", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    with pytest.raises(FileNotFoundError, match="PeptiVerse not found"):
        PeptiVerse()


def test_home_without_inference_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="inference.py not found"):
        PeptiVerse(peptiverse_home=str(tmp_path))


def test_home_without_weights_is_reported(tmp_path):
    (tmp_path / "inference.py").write_text("")
    with pytest.raises(FileNotFoundError, match="half_life not found"):
        PeptiVerse(peptiverse_home=str(tmp_path))


def test_supported_kinds_and_mhc_context(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    assert predictor.supported_kinds == (Kind.serum_half_life,)
    support = predictor.kind_support()[Kind.serum_half_life]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "none"


def test_empty_peptide_list_returns_nothing(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    assert predictor.predict([]) == []


def test_modified_residues_are_rejected_not_silently_stripped(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    with pytest.raises(ValueError, match="unmodified sequence"):
        predictor.predict(["SIINFEKLB"])


def test_empty_peptide_is_rejected(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    with pytest.raises(ValueError, match="Empty peptide"):
        predictor.predict([""])


def test_over_long_peptide_is_rejected_rather_than_truncated(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    with pytest.raises(ValueError, match="up to 1020 residues"):
        predictor.predict(["A" * (PEPTIVERSE_MAX_PEPTIDE_LENGTH + 1)])


def test_max_peptide_length_cannot_exceed_model_limit(tmp_path):
    with pytest.raises(ValueError, match="at most 1020 residues"):
        PeptiVerse(
            peptiverse_home=_fake_home(tmp_path),
            max_peptide_length=PEPTIVERSE_MAX_PEPTIDE_LENGTH + 1)


def test_missing_python_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="Python does not exist"):
        PeptiVerse(
            peptiverse_home=_fake_home(tmp_path),
            peptiverse_python="/nonexistent/python")


def test_predictor_version_pins_the_upstream_revision(tmp_path):
    predictor = PeptiVerse(peptiverse_home=_fake_home(tmp_path))
    assert predictor.predictor_version.startswith("8cf0b21")
    assert predictor.predictor_version.endswith("transformer_wt_log")


# --- end-to-end (opt-in; needs a real snapshot) -----------------------------

PEPTIVERSE_HOME = os.environ.get("PEPTIVERSE_HOME")

requires_peptiverse = pytest.mark.skipif(
    not PEPTIVERSE_HOME,
    reason="set PEPTIVERSE_HOME to a PeptiVerse snapshot to run this")


@requires_peptiverse
def test_end_to_end_returns_hours():
    predictor = PeptiVerse(device="cpu")
    peptides = ["SIINFEKL", "KLGGALQAK"]
    results = predictor.predict(peptides)
    assert len(results) == len(peptides)
    for peptide, result in zip(peptides, results):
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.serum_half_life
        assert pred.peptide == peptide
        assert pred.allele == ""
        assert pred.predictor_name == "peptiverse"
        # expm1 of a log1p-trained target: a duration in hours, so positive and
        # not a 0-1 probability.
        assert pred.value > 0
        assert pred.value == pred.score
    assert len(predictor.last_qc) == len(peptides)


@requires_peptiverse
def test_end_to_end_is_deterministic():
    predictor = PeptiVerse(device="cpu")
    first = predictor.predict(["SIINFEKL"])[0].preds[0].value
    second = predictor.predict(["SIINFEKL"])[0].preds[0].value
    assert first == second


@requires_peptiverse
def test_end_to_end_dataframe_has_the_standard_columns():
    from mhctools.pred import COLUMNS
    predictor = PeptiVerse(device="cpu")
    frame = predictor.predict_dataframe(["SIINFEKL"])
    assert list(frame.columns) == list(COLUMNS)
    assert frame["kind"].tolist() == [Kind.serum_half_life]


@requires_peptiverse
def test_end_to_end_uncertainty_stays_out_of_the_prediction():
    # Upstream's half-life conformal bundle references a class defined in its
    # training script's __main__, so it cannot be unpickled at inference time.
    # Whatever upstream reports must stay in last_qc, never in a field that
    # would read as a calibrated bound on the hours estimate.
    predictor = PeptiVerse(device="cpu", uncertainty=True)
    pred = predictor.predict(["SIINFEKL"])[0].preds[0]
    assert pred.percentile_rank is None
    assert pred.value == pred.score
    assert "uncertainty_type" in predictor.last_qc.columns
