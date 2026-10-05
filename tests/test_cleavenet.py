"""Whole-substrate contracts and pinned-source CleaveNet conformance."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from mhctools import CleaveNet, PeptideInput
from mhctools import cleavenet

REFERENCE = json.loads((Path(__file__).parent / "data/cleavenet_source_reference.json").read_text())


@pytest.fixture
def predictor(tmp_path, monkeypatch):
    # Keep an actual checksum verification path without creating model stand-ins.
    from mhctools.optional_backend import sha256_file
    source = tmp_path / "source.txt"
    source.write_text("reviewed runtime")
    monkeypatch.setattr(cleavenet, "_ARTIFACTS", {
        "source.txt": ("fixture", sha256_file(source), "text")})
    return CleaveNet(cleavenet_home=tmp_path, cleavenet_python=sys.executable)


def mock_sidecar(monkeypatch, transform=lambda output: output):
    def run(name, python, script, args, **kwargs):
        assert name == "CleaveNet"
        assert Path(script).name == "cleavenet_sidecar.py"
        assert "TF_USE_LEGACY_KERAS" not in kwargs["environment"]
        home, input_path, output_path = args
        sequences = json.loads(Path(input_path).read_text())
        means, deviations = [], []
        for sequence in sequences:
            index = REFERENCE["sequences"].index(sequence)
            means.append(list(REFERENCE["means"][index]))
            deviations.append(list(REFERENCE["ensemble_sd"][index]))
        output = dict(sequences=sequences, enzymes=list(cleavenet.ENZYMES),
                      means=means, ensemble_sd=deviations,
                      runtime={"tensorflow": "2.18.0", "python": "3.12.6"})
        Path(output_path).write_text(json.dumps(transform(output)))
    monkeypatch.setattr(cleavenet, "run_python_sidecar", run)


def test_centered_padding_native_scores_duplicates_and_provenance(predictor, monkeypatch):
    mock_sidecar(monkeypatch)
    monkeypatch.setenv("TF_USE_LEGACY_KERAS", "1")
    inputs = [PeptideInput("LRVFL", source_sequence_name="vaccine", source_start=17),
              PeptideInput("LRVFL", source_sequence_name="vaccine", source_start=32),
              PeptideInput("PRVFQLRVFL")]
    results = predictor.predict(inputs)
    assert [r.peptide_input for r in results] == inputs
    assert [r.padded_sequence for r in results] == ["--LRVFL---", "--LRVFL---", "PRVFQLRVFL"]
    assert results[0].cache_key == results[1].cache_key
    assert results[0].scores == results[1].scores
    assert [s.enzyme for s in results[0].scores] == list(cleavenet.ENZYMES)
    assert "MMP24" in [s.enzyme for s in results[0].scores]
    assert [s.z_score for s in results[0].scores] == REFERENCE["means"][1]
    assert [s.ensemble_sd for s in results[0].scores] == REFERENCE["ensemble_sd"][1]
    record = results[0].to_dict()
    assert record["peptide_input"]["source_start"] == 17
    assert record["endpoint"] == cleavenet.ENDPOINT
    assert record["inventory"]["inference_status"] == "reproduced"
    assert "bond" not in record
    assert "limited validation" in record["applicability"]
    json.dumps(record, allow_nan=False)


@pytest.mark.parametrize("input_", ["aHA", "AHA ", "A-X", "", "ACDEFGHIKLM",
                                    PeptideInput("LRVFL", c_term="amidated"),
                                    PeptideInput("LRVFL", attachments=(("N", "cargo"),))])
def test_unsupported_inputs_rejected_before_inference(predictor, monkeypatch, input_):
    monkeypatch.setattr(predictor, "_run_sidecar", lambda *a: pytest.fail("inference was attempted"))
    with pytest.raises(ValueError):
        predictor.predict(["LRVFL", input_])


def test_empty_batch_does_not_invoke_runtime(predictor, monkeypatch):
    monkeypatch.setattr(predictor, "_run_sidecar", lambda *a: pytest.fail("inference was attempted"))
    assert predictor.predict([]) == []
    assert predictor.predict_dataframe([]).empty


def test_windows_keep_coordinates_and_whole_substrate_endpoint(predictor, monkeypatch):
    requested = []
    def predict(inputs):
        requested.extend(inputs)
        return inputs
    monkeypatch.setattr(predictor, "predict", predict)
    original = PeptideInput("ACDEFGHIKLMN", source_sequence_name="vaccine", source_start=40)
    predictor.predict_windows(original)
    assert [i.sequence for i in requested] == ["ACDEFGHIKL", "CDEFGHIKLM", "DEFGHIKLMN"]
    assert [i.source_start for i in requested] == [40, 41, 42]
    assert all(i.source_sequence_name == "vaccine" for i in requested)
    with pytest.raises(ValueError, match="at least 10"):
        predictor.predict_windows("LRVFL")
    with pytest.raises(ValueError, match="chemistry"):
        predictor.predict_windows(PeptideInput("ACDEFGHIKLMN", c_term="amidated"))


def test_dataframe_keeps_enzymes_and_source_intervals(predictor, monkeypatch):
    mock_sidecar(monkeypatch)
    frame = predictor.predict_dataframe([PeptideInput("LRVFL", source_start=8)])
    assert len(frame) == 18
    assert frame.enzyme.tolist() == list(cleavenet.ENZYMES)
    assert set(frame.source_start) == {8}
    assert set(frame.source_end) == {13}
    assert set(frame.endpoint) == {cleavenet.ENDPOINT}
    assert "bond" not in frame.columns
    assert frame.z_score.tolist() == REFERENCE["means"][1]


@pytest.mark.parametrize("fault", ["identity", "enzymes", "shape", "nonfinite", "negative_sd"])
def test_bad_sidecar_outputs_fail_without_partial_results(predictor, monkeypatch, fault):
    def alter(output):
        if fault == "identity":
            output["sequences"] = ["AAAAAAAAAA"]
        elif fault == "enzymes":
            output["enzymes"] = output["enzymes"][::-1]
        elif fault == "shape":
            output["means"] = []
        elif fault == "nonfinite":
            output["means"][0][0] = float("nan")
        else:
            output["ensemble_sd"][0][0] = -1
        return output
    mock_sidecar(monkeypatch, alter)
    with pytest.raises((ValueError, RuntimeError)):
        predictor.predict("LRVFL")
    assert predictor.artifact_inventory.inference_status == "not_run"


def test_asset_tampering_prevents_construction(predictor):
    (predictor.home / "source.txt").write_text("changed runtime")
    with pytest.raises(RuntimeError, match="mismatch"):
        CleaveNet(cleavenet_home=predictor.home)


def test_host_import_does_not_load_tensorflow():
    subprocess.run([sys.executable, "-c", "import sys; from mhctools import CleaveNet; "
                    "assert 'tensorflow' not in sys.modules"], check=True)


@pytest.mark.requires_external_tool
def test_real_models_match_official_source():
    if not os.environ.get("CLEAVENET_HOME") or not os.environ.get("CLEAVENET_PYTHON"):
        pytest.skip("Set CLEAVENET_HOME and CLEAVENET_PYTHON to the provisioned runtime")
    predictor = CleaveNet()
    results = predictor.predict(["PRVFQLRVFL", "LRVFL", "PRVFQLRVFL"])
    for result, index in zip(results, (0, 1, 0)):
        assert [s.z_score for s in result.scores] == pytest.approx(REFERENCE["means"][index], abs=2e-5)
        assert [s.ensemble_sd for s in result.scores] == pytest.approx(REFERENCE["ensemble_sd"][index], abs=2e-5)
        assert result.inventory.status == "verified"
        assert result.inventory.inference_status == "reproduced"
        assert dict(result.runtime)["tensorflow"] == "2.18.0"
