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

from mhctools import Kind, PeptideContext, PeptideInput, PeptiVerse
from mhctools.peptiverse import (
    PEPTIVERSE_MAX_PEPTIDE_LENGTH,
    _HALF_LIFE_MANIFEST,
    _MODEL_DIRECTORY,
    parse_peptiverse_results,
)
from mhctools.pred import VALUE_BEST_DIRECTIONS, best_direction


@pytest.fixture(autouse=True)
def peptiverse_environment(monkeypatch):
    """Keep synthetic snapshots independent of the user's real ESM install."""
    configured = {}
    for argument, variable in (
            ("peptiverse_home", "PEPTIVERSE_HOME"),
            ("peptiverse_python", "PEPTIVERSE_PYTHON"),
            ("peptiverse_esm_home", "PEPTIVERSE_ESM_HOME")):
        configured[argument] = os.environ.get(variable)
        monkeypatch.delenv(variable, raising=False)
    return configured


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

def test_peptide_half_life_is_not_pmhc_stability():
    assert Kind.peptide_half_life != Kind.pMHC_stability
    assert Kind.peptide_half_life == "peptide_half_life"


def test_peptide_half_life_has_no_context_free_best_direction():
    assert VALUE_BEST_DIRECTIONS[Kind.peptide_half_life] == "max"
    with pytest.raises(ValueError, match="context-dependent"):
        best_direction(Kind.peptide_half_life, "value")


def test_annotate_exposes_serum_half_life_as_a_units_bearing_field():
    from mhctools.annotate import _OUTPUT_FIELDS, output_field_tokens
    assert "serum_half_life" in output_field_tokens()
    # `value` and not `score`, because hours are the meaningful quantity.
    assert _OUTPUT_FIELDS["serum_half_life"] == (
        Kind.peptide_half_life, "value")


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
    # The exact name prevents upstream's generic Transformer alias from falling
    # back to transformer_wt, whose output has different provenance and units.
    assert fields[1] == "Transformer_WT_Log"
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


def _stub_inference(
        tmp_path,
        hours=1.25,
        retained_limit=1020,
        model_name="transformer_wt_log"):
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
    def __init__(self, manifest_path, classifier_weight_root, esm_name, device):
        assert Path(classifier_weight_root) == Path(__file__).resolve().parent
        assert Path(esm_name) == Path(classifier_weight_root) / "esm2_t33_650M_UR50D"
        self.wt_embedder = Embedder()
        artifact = Path(classifier_weight_root) / {str(_MODEL_DIRECTORY)!r} / "best_model.pt"
        self.meta = {{("halflife", "wt"): {{
            "artifact": str(artifact),
            "model_name": {model_name!r},
            "emb_tag": "wt",
            "kind": "torch_ckpt",
        }}}}
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
        predictor = PeptiVerse(
            peptiverse_python=sys.executable,
            allow_unverified_assets=True)
    else:
        predictor = PeptiVerse(
            peptiverse_home=relative,
            peptiverse_python=sys.executable,
            allow_unverified_assets=True)
    monkeypatch.chdir(tmp_path / "nested")
    assert predictor.peptiverse_home == str(home.resolve())
    peptides = ["A" * 1019 + "K", "A" * 1019 + "R"]
    results = predictor.predict(peptides)
    assert [r.peptide for r in results] == peptides
    assert [r.serum_half_life.value for r in results] == [1.25, 1.25]
    assert predictor.artifact_inventory.inference_status == "reproduced"
    assert all(
        result.serum_half_life.predictor_version == predictor.predictor_version
        for result in results)
    assert (results[0].peptide_half_life.measurement_context is
            results[1].peptide_half_life.measurement_context)


def test_contextual_inputs_preserve_order_duplicates_and_partial_failure(
        tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_stub_inference(tmp_path),
        peptiverse_python=sys.executable,
        allow_unverified_assets=True)
    context = PeptideContext(matrix="serum", assay_species="Homo sapiens")
    inputs = [
        PeptideInput("SIINFEKL", occurrence_id="first", context=context),
        PeptideInput(
            "SIINFEKL", c_term="amidated", occurrence_id="unsupported",
            context=context),
        PeptideInput("SIINFEKL", occurrence_id="second", context=context),
    ]

    results = predictor.predict(inputs, on_unsupported="record")
    predictions = [result.preds[0] for result in results]

    assert [pred.peptide_input.occurrence_id for pred in predictions] == [
        "first", "unsupported", "second"]
    assert [pred.measurement_context.status for pred in predictions] == [
        "available", "unsupported", "available"]
    assert predictions[1].score is None
    assert predictions[1].value is None
    assert "C-terminal chemistry" in predictions[1].measurement_context.detail
    assert predictions[0].cache_key == predictions[2].cache_key
    assert predictions[1].cache_key != predictions[0].cache_key
    assert predictions[0].peptide_input.record_identity_sha256 != \
        predictions[2].peptide_input.record_identity_sha256
    assert predictions[0].measurement_context.matrix == "human serum"


def test_contextual_modified_input_is_rejected_by_default(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    with pytest.raises(ValueError, match="unsupported C-terminal chemistry"):
        predictor.predict([PeptideInput("SIINFEKL", c_term="amidated")])


@pytest.mark.parametrize("hours", ["nan", "inf", "-inf", -1.0])
def test_public_predict_rejects_invalid_sidecar_duration(tmp_path, hours):
    home = _stub_inference(tmp_path, hours=hours)
    predictor = PeptiVerse(
        peptiverse_home=home,
        peptiverse_python=sys.executable,
        allow_unverified_assets=True)
    with pytest.raises(RuntimeError, match="invalid half-lives"):
        predictor.predict(["SIINFEKL"])


def test_sidecar_checks_actual_tokenizer_residue_count(tmp_path):
    home = _stub_inference(tmp_path, retained_limit=7)
    predictor = PeptiVerse(
        peptiverse_home=home,
        peptiverse_python=sys.executable,
        allow_unverified_assets=True)
    with pytest.raises(RuntimeError, match="retained 7 of 8 residues"):
        predictor.predict(["SIINFEKL"])


def test_sidecar_rejects_unexpected_loaded_model_metadata(tmp_path):
    home = _stub_inference(tmp_path, model_name="transformer_wt")
    predictor = PeptiVerse(
        peptiverse_home=home,
        peptiverse_python=sys.executable,
        allow_unverified_assets=True)
    with pytest.raises(RuntimeError, match="Unexpected PeptiVerse.*metadata"):
        predictor.predict(["SIINFEKL"])


@pytest.mark.parametrize("length", [1021, 1022])
def test_rejects_lengths_that_fit_only_without_special_tokens(tmp_path, length):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    with pytest.raises(ValueError, match="up to 1020 residues"):
        predictor.predict(["A" * (length - 1) + "K"])


# --- construction and validation (no model) ---------------------------------

def _fake_home(tmp_path):
    """A directory shaped enough like a snapshot to pass resolution."""
    (tmp_path / "inference.py").write_text("")
    model_dir = tmp_path / _MODEL_DIRECTORY
    model_dir.mkdir(parents=True)
    for name in ("best_model.pt", "best_params.json", "mapie_calibration.joblib"):
        (model_dir / name).write_text("fake " + name)
    esm_dir = tmp_path / "esm2_t33_650M_UR50D"
    esm_dir.mkdir()
    for name in (
            "config.json", "tokenizer_config.json", "special_tokens_map.json",
            "vocab.txt", "model.safetensors"):
        (esm_dir / name).write_text("fake " + name)
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
    with pytest.raises(FileNotFoundError, match="transformer_wt_log not found"):
        PeptiVerse(peptiverse_home=str(tmp_path))


def test_generic_transformer_fallback_is_not_accepted(tmp_path):
    (tmp_path / "inference.py").write_text("")
    fallback = tmp_path / "training_classifiers" / "half_life" / "transformer_wt"
    fallback.mkdir(parents=True)
    (fallback / "best_model.pt").write_text("different model")
    with pytest.raises(FileNotFoundError, match="will not fall back"):
        PeptiVerse(peptiverse_home=str(tmp_path))


def test_supported_kinds_and_mhc_context(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    assert predictor.supported_kinds == (Kind.peptide_half_life,)
    support = predictor.kind_support()[Kind.peptide_half_life]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "none"


def test_empty_peptide_list_returns_nothing(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    assert predictor.predict([]) == []


def test_modified_residues_are_rejected_not_silently_stripped(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    with pytest.raises(ValueError, match="unmodified sequence"):
        predictor.predict(["SIINFEKLB"])


def test_empty_peptide_is_rejected(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    with pytest.raises(ValueError, match="Empty peptide"):
        predictor.predict([""])


def test_over_long_peptide_is_rejected_rather_than_truncated(tmp_path):
    predictor = PeptiVerse(
        peptiverse_home=_fake_home(tmp_path),
        allow_unverified_assets=True)
    with pytest.raises(ValueError, match="up to 1020 residues"):
        predictor.predict(["A" * (PEPTIVERSE_MAX_PEPTIDE_LENGTH + 1)])


def test_max_peptide_length_cannot_exceed_model_limit(tmp_path):
    with pytest.raises(ValueError, match="at most 1020 residues"):
        PeptiVerse(
            peptiverse_home=_fake_home(tmp_path),
            max_peptide_length=PEPTIVERSE_MAX_PEPTIDE_LENGTH + 1,
            allow_unverified_assets=True)


def test_missing_python_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="Python does not exist"):
        PeptiVerse(
            peptiverse_home=_fake_home(tmp_path),
            peptiverse_python="/nonexistent/python",
            allow_unverified_assets=True)


def test_unverified_assets_are_rejected_by_default(tmp_path):
    with pytest.raises(RuntimeError, match="mismatch artifacts"):
        PeptiVerse(peptiverse_home=_fake_home(tmp_path))


def test_predictor_version_carries_actual_asset_identity(tmp_path):
    home = Path(_fake_home(tmp_path))
    predictor = PeptiVerse(
        peptiverse_home=home,
        allow_unverified_assets=True)
    first_version = predictor.predictor_version
    assert "developed-against=peptiverse@8cf0b21" in first_version
    assert ";assets-sha256=" in first_version
    assert first_version.endswith(";status=mismatch")

    (home / _MODEL_DIRECTORY / "best_model.pt").write_text("replacement")
    replacement = PeptiVerse(
        peptiverse_home=home,
        allow_unverified_assets=True)
    assert replacement.predictor_version != first_version


def test_output_affecting_device_setting_changes_backend_identity(tmp_path):
    home = _fake_home(tmp_path)
    cpu = PeptiVerse(
        peptiverse_home=home,
        device="cpu",
        allow_unverified_assets=True)
    cuda = PeptiVerse(
        peptiverse_home=home,
        device="cuda",
        allow_unverified_assets=True)
    assert cpu.artifact_inventory.identity_sha256 != \
        cuda.artifact_inventory.identity_sha256


def test_missing_esm_snapshot_is_reported(tmp_path, monkeypatch):
    home = Path(_fake_home(tmp_path))
    esm = home / "esm2_t33_650M_UR50D"
    esm.rename(home / "removed-esm")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "empty"))
    with pytest.raises(FileNotFoundError, match="will not download"):
        PeptiVerse(peptiverse_home=home)


# --- end-to-end (opt-in; needs a real snapshot) -----------------------------

PEPTIVERSE_HOME = os.environ.get("PEPTIVERSE_HOME")

requires_peptiverse = pytest.mark.skipif(
    not PEPTIVERSE_HOME,
    reason="set PEPTIVERSE_HOME to a PeptiVerse snapshot to run this")


@requires_peptiverse
def test_end_to_end_returns_hours(peptiverse_environment):
    predictor = PeptiVerse(device="cpu", **peptiverse_environment)
    peptides = ["SIINFEKL", "KLGGALQAK"]
    results = predictor.predict(peptides)
    assert len(results) == len(peptides)
    for peptide, result in zip(peptides, results):
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.peptide_half_life
        assert pred.peptide == peptide
        assert pred.allele == ""
        assert pred.predictor_name == "peptiverse"
        # expm1 of a log1p-trained target: a duration in hours, so positive and
        # not a 0-1 probability.
        assert pred.value > 0
        assert pred.value == pred.score
    assert len(predictor.last_qc) == len(peptides)


@requires_peptiverse
def test_end_to_end_is_deterministic(peptiverse_environment):
    predictor = PeptiVerse(device="cpu", **peptiverse_environment)
    first = predictor.predict(["SIINFEKL"])[0].preds[0].value
    second = predictor.predict(["SIINFEKL"])[0].preds[0].value
    assert first == second


@requires_peptiverse
def test_end_to_end_dataframe_has_the_standard_columns(peptiverse_environment):
    from mhctools.pred import COLUMNS
    predictor = PeptiVerse(device="cpu", **peptiverse_environment)
    frame = predictor.predict_dataframe(["SIINFEKL"])
    assert list(frame.columns) == list(COLUMNS)
    assert frame["kind"].tolist() == [Kind.peptide_half_life]


@requires_peptiverse
def test_end_to_end_uncertainty_stays_out_of_the_prediction(peptiverse_environment):
    # Upstream's half-life conformal bundle references a class defined in its
    # training script's __main__, so it cannot be unpickled at inference time.
    # Whatever upstream reports must stay in last_qc, never in a field that
    # would read as a calibrated bound on the hours estimate.
    predictor = PeptiVerse(device="cpu", uncertainty=True, **peptiverse_environment)
    pred = predictor.predict(["SIINFEKL"])[0].preds[0]
    assert pred.percentile_rank is None
    assert pred.value == pred.score
    assert "uncertainty_type" in predictor.last_qc.columns
