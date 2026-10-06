"""Published-equation and independent author-app descriptor conformance."""

import hashlib
from importlib import resources
import json
import math
from pathlib import Path
import socket

import pandas as pd
import pytest

from mhctools import CavacoHalfLife, Kind, PeptideContext, PeptideInput, Prediction
from mhctools.cavaco import MODEL_SHA256, _model
from mhctools.cli.annotate_table import main as table_main
from mhctools.cli.args import mhc_predictors


REFERENCE = json.loads((Path(__file__).parent / "data" /
                        "cavaco_source_reference.json").read_text())


@pytest.mark.parametrize("row", REFERENCE["records"], ids=lambda row: row["peptide"])
def test_unmodified_js_descriptors_and_published_equation(row):
    predictor = CavacoHalfLife()
    pred = predictor.predict(row["peptide"])[0].peptide_half_life
    qc = predictor.last_qc.iloc[0]
    assert qc.pI == row["pi"]  # raw author-app result, before display rounding
    assert qc.nonpolar_percent == row["np_percent"]
    assert pred.score == pytest.approx(row["ln_published"], abs=1e-14)
    assert pred.value == pytest.approx(math.exp(row["ln_published"]) / 60)
    assert pred.kind == Kind.peptide_half_life
    assert pred.allele == ""
    assert pred.measurement_context.matrix is None
    assert pred.measurement_context.unit == "hours"
    assert pred.measurement_context.transform == "linear"
    assert "ln(half-life in minutes)" in pred.measurement_context.score_semantics
    assert MODEL_SHA256 in pred.predictor_version


def _log_r_squared(observed, predicted):
    x = [math.log(value) for value in observed]
    y = [math.log(value) for value in predicted]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return (sum((a - mx) * (b - my) for a, b in zip(x, y)) ** 2 /
            (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)))


def test_source_panel_reproduction_is_distinct_from_bugged_app_and_external_validation():
    panel = [row for row in REFERENCE["records"] if "source_panel_id" in row]
    assert len(panel) == 16
    assert all(row["observed_chemistry"] == "C-terminal carboxamide" for row in panel)
    assert all(row["observed_matrix"] == "50% human serum" for row in panel)
    assert all(row["observed_replicates"] == 3 for row in panel)
    assert panel[-1]["observed_minutes"] == 15
    observed = [row["observed_minutes"] for row in panel]
    predictions = CavacoHalfLife().predict([row["peptide"] for row in panel])
    minutes = [result.peptide_half_life.value * 60 for result in predictions]
    assert _log_r_squared(observed, minutes) == pytest.approx(0.7601643279437353)
    app_minutes = [math.exp(row["ln_app"]) for row in panel]
    assert _log_r_squared(observed, app_minutes) == pytest.approx(0.5568722659936396)
    rmse = math.sqrt(sum(math.log(p / o) ** 2 for p, o in zip(minutes, observed)) / 16)
    assert rmse == pytest.approx(0.9206305109002861)
    assert minutes[3] == pytest.approx(568.31, abs=0.01)
    assert observed[3] == 89


@pytest.mark.parametrize("residue", list("ACDEFGHIKLMNPQRSTVWY"))
def test_nonpolar_descriptor_matches_deposited_classification(residue):
    predictor = CavacoHalfLife()
    predictor.predict(residue, isoelectric_points=[6])
    assert predictor.last_qc.iloc[0].nonpolar_percent == (
        100 if residue in "ACFILMPVWY" else 0)


def test_independent_coefficients_and_binary_boundaries():
    predictor = CavacoHalfLife()
    seqs = ["AAA", "WAA", "WWA", "YAA", "YYA", "YYY", "GGG", "AAA", "AAA"]
    preds = [r.peptide_half_life for r in predictor.predict(
        seqs, isoelectric_points=[6] * 7 + [9.999999, 10])]
    assert preds[0].score == pytest.approx(2.226 + 5.3)
    assert preds[6].score == pytest.approx(2.226)
    assert preds[1].score - preds[0].score == pytest.approx(-1.515)
    assert preds[2].score == preds[1].score  # W presence, not count
    assert preds[3].score == preds[0].score  # one Y has no binary bonus
    assert preds[4].score - preds[0].score == pytest.approx(1.290)
    assert preds[5].score == preds[4].score
    assert preds[8].score - preds[7].score == pytest.approx(-1.052)
    assert preds[8].value / preds[7].value == pytest.approx(math.exp(-1.052))


def test_provided_pi_is_per_row_and_affects_identity_and_serialization():
    predictor = CavacoHalfLife()
    result = predictor.predict(["AAA"] * 3, isoelectric_points=[None, 6, 10])
    preds = [r.peptide_half_life for r in result]
    assert len({p.cache_key for p in preds}) == 3
    assert preds[0].predictor_version != preds[1].predictor_version
    assert "provided pI=10.0" in preds[2].measurement_context.detail
    assert predictor.last_qc.pI_source.tolist() == [
        "author-app-free-terminal", "provided", "provided"]
    for pred in preds:
        assert Prediction.from_dict(json.loads(json.dumps(pred.to_dict()))) == pred
    assert predictor.predict("AAA", isoelectric_points=[10])[0].peptide_half_life == preds[2]


def test_occurrences_duplicates_context_and_unavailable_chemistry_are_retained():
    inputs = [
        PeptideInput("SIINFEKL", occurrence_id="one", source_gene="first"),
        PeptideInput("SIINFEKL", c_term="amidated", occurrence_id="amide"),
        PeptideInput("SIINFEKL", occurrence_id="two", source_gene="second"),
        PeptideInput("SIINFEKL", n_term="acetylated", occurrence_id="acetyl"),
        PeptideInput("SIINFEKL", attachments=(("1", "phosphorylation"),)),
        PeptideInput("SIINFEKL", context=PeptideContext(matrix="human serum")),
    ]
    predictor = CavacoHalfLife()
    with pytest.raises(ValueError, match="Input 1 is unsupported"):
        predictor.predict(inputs)
    # Convenience accessors select available estimates; unavailable records
    # remain in preds and must not disappear from the original input order.
    results = predictor.predict(inputs, on_unsupported="record")
    preds = [r.preds[0] for r in results]
    assert results[1].peptide_half_life is None
    assert [p.peptide_input for p in preds] == inputs
    assert [p.measurement_context.status for p in preds] == [
        "available", "unsupported", "available", "unsupported", "unsupported", "available"]
    assert preds[0].cache_key == preds[2].cache_key  # occurrence-independent inference
    assert preds[0].cache_key != preds[5].cache_key  # context is part of input identity
    assert preds[1].cache_key != preds[0].cache_key  # chemistry is part of input identity
    for pred in preds:
        assert Prediction.from_dict(json.loads(json.dumps(pred.to_dict()))) == pred
    for index in (1, 3, 4):
        assert preds[index].score is None and preds[index].value is None
        assert preds[index].measurement_context.unit is None
        assert "Unsupported:" in preds[index].measurement_context.detail
    assert all(p.measurement_context.matrix is None for p in preds)
    assert predictor.last_qc.input_index.tolist() == list(range(6))


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -0.1, 14.1, 10 ** 1000, True, "6"])
def test_invalid_pi_rejected(value):
    with pytest.raises(ValueError, match="finite number"):
        CavacoHalfLife().predict("AAA", isoelectric_points=[value])


def test_pi_length_validation_empty_and_input_validation():
    predictor = CavacoHalfLife()
    with pytest.raises(ValueError, match="one entry per input"):
        predictor.predict(["AAA", "WWW"], isoelectric_points=[6])
    with pytest.raises(ValueError, match="on_unsupported"):
        predictor.predict("AAA", on_unsupported="skip")
    for invalid in ("", "XAA", "A-D"):
        with pytest.raises(ValueError):
            predictor.predict(invalid, on_unsupported="record")
    assert predictor.predict([]) == []
    assert predictor.last_qc.empty
    assert predictor.predict(" siinFekl ")[0].peptide == "SIINFEKL"


def test_long_formula_inputs_do_not_acquire_validation_or_serum_semantics():
    results = CavacoHalfLife().predict(["G" * 80, "G"])
    assert results[0].peptide_half_life.value == results[1].peptide_half_life.value
    pred = results[0].peptide_half_life
    assert "long vaccine peptide accuracy unestablished" in pred.measurement_context.detail
    assert results[0].serum_half_life is None
    assert results[0].blood_half_life is None


def test_offline_cli_dataframe_and_registry(tmp_path, monkeypatch):
    def no_network(*args, **kwargs):
        raise AssertionError("Inference must be offline")
    monkeypatch.setattr(socket, "create_connection", no_network)
    assert mhc_predictors["cavaco"] is CavacoHalfLife
    source, output = tmp_path / "input.csv", tmp_path / "output.csv"
    source.write_text("peptide,label\nSIINFEKL,a\nAAA,b\nSIINFEKL,c\n")
    table_main(["--input", str(source), "--out", str(output),
                "--predictor", "cavaco:half_life_hours:peptide_half_life"])
    frame = pd.read_csv(output)
    assert frame.label.tolist() == ["a", "b", "c"]
    predictor = CavacoHalfLife()
    expected = [r.peptide_half_life.value for r in predictor.predict(frame.peptide)]
    assert frame.half_life_hours.tolist() == pytest.approx(expected, rel=1e-5)
    assert len(predictor.predict_dataframe(frame.peptide)) == 3
    assert predictor.kind_support()[Kind.peptide_half_life]["mhc_dependence"] == "none"


def test_bundled_resource_identity_integrity_and_immutability(monkeypatch):
    raw = resources.files("mhctools.data").joinpath("cavaco_published_model.json").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == MODEL_SHA256
    with pytest.raises(TypeError):
        _model()["coefficients"]["intercept"] = 0
    with pytest.raises(TypeError):
        _model()["pka"]["A"]["alpha_amino"] = 0
    class BadResource:
        def joinpath(self, name):
            return self
        def read_bytes(self):
            return b"{}"
    _model.cache_clear()
    try:
        with monkeypatch.context() as patched:
            patched.setattr("mhctools.cavaco.resources.files", lambda package: BadResource())
            with pytest.raises(ValueError, match="checksum mismatch"):
                CavacoHalfLife()
    finally:
        _model.cache_clear()
