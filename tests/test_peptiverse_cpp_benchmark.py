"""Source-label denominators, provenance gaps and native CPP integration."""

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from mhctools import Kind, MeasurementContext, PeptideResult, Prediction
from mhctools.cli.benchmark import main
from mhctools import peptiverse_cpp_benchmark as cpp


SYNTHETIC = ("sequence,label,id,split\n"
             "ACD,1,seq_0,train\nACDE,0,seq_1,train\n"
             "ACD,0,seq_2,val\nKKKK,1,seq_3,val\nRRRRR,1,seq_4,val\n")
METADATA = (Path(os.environ["PEPTIVERSE_HOME"]) / cpp.CPP_METADATA_PATH
            if os.environ.get("PEPTIVERSE_HOME") else None)


@pytest.fixture
def synthetic(monkeypatch):
    rows = cpp._parse_metadata(SYNTHETIC)
    monkeypatch.setattr(cpp, "load_cpp_metadata", lambda path: rows)
    return rows


class StubCPP:
    def __init__(self, scores=(0.7, 0.5493, 0.1), failure=None):
        self.scores = scores
        self.failure = failure
        self.artifact_inventory = SimpleNamespace(to_dict=lambda: {"synthetic": True})

    def predict(self, inputs, on_unsupported):
        assert on_unsupported == "record"
        self.inputs = inputs
        if self.failure:
            raise RuntimeError(self.failure)
        return [PeptideResult(preds=(Prediction(
            kind=Kind.cpp_classification, peptide=item.sequence,
            peptide_input=item, score=score,
            measurement_context=MeasurementContext(
                estimate_type="ml_predicted", status="available", score_semantics="native SVC positive-class P(CPP)",
                class_label="CPP" if score >= cpp.CPP_THRESHOLD else "non-CPP")),))
            for item, score in zip(inputs, self.scores)]


@pytest.mark.parametrize("replacement", [
    ("sequence,label,id,split", "id,sequence,label,split"),
    ("seq_4", "seq_3"), ("seq_4", "source-4"), ("RRRRR", "RRXRR"),
    ("RRRRR", "rrrrr"), ("RRRRR,1", "RRRRR,2"), ("seq_4,val", "seq_4,test"),
    ("RRRRR,1,seq_4,val", "RRRRR,1,seq_4,val,extra"),
    ("RRRRR,1,seq_4,val", "RRRRR,1,seq_4"),
])
def test_parser_rejects_corrupt_source_rows(replacement):
    with pytest.raises(ValueError):
        cpp._parse_metadata(SYNTHETIC.replace(*replacement))


def test_no_validation_split_and_wrong_hash_fail_before_model_load(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="nonempty"):
        cpp._parse_metadata(SYNTHETIC.replace(",val", ",train"))
    path = tmp_path / "metadata.csv"
    path.write_text(SYNTHETIC)
    monkeypatch.setattr(cpp, "_score_rows", lambda *args: pytest.fail("Model must not load"))
    with pytest.raises(ValueError, match="SHA-256"):
        cpp.evaluate_cpp_metadata(path, predict=True)


def test_audit_only_keeps_unknown_chemistry_and_all_denominators(synthetic, monkeypatch):
    monkeypatch.setattr(cpp, "_score_rows", lambda selected, predictor:
                        ([], None) if not selected else pytest.fail("No native runtime in audit"))
    report = cpp.evaluate_cpp_metadata("synthetic")
    audit = report["source_audit"]
    assert audit["records"] == 5
    assert audit["duplicate_sequences"] == 1
    assert audit["conflicting_label_sequences"] == 1
    assert audit["exact_sequence_overlap_count"] == 1
    assert audit["study_overlap"] == audit["family_overlap"] == "unknown"
    assert audit["splits"]["val"]["records_with_unknown_metadata"]["chemical_form"] == 3
    assert report["cohort"]["status_counts"] == {"not_assessed": 3}
    assert report["cohort"]["native_threshold_confusion"] is None
    assert not report["cohort"]["full_validation_cohort"]
    assert report["artifact_inventory"] is None
    records = report["benchmark"]["records"]
    assert [r["measurement"]["source_measurement_id"] for r in records] == ["seq_2", "seq_3", "seq_4"]
    assert all(r["measurement"]["chemistry"] == "unknown" for r in records)
    assert all(not r["independence_established"] for r in records)
    assert all(d["evidence"] == "no_evidence_in_supplied_data"
               for d in report["benchmark"]["requested_domains"])


def test_native_probability_and_threshold_reproduction_not_external_validation(synthetic):
    predictor = StubCPP()
    report = cpp.evaluate_cpp_metadata("synthetic", predict=True, predictor=predictor)
    assert [item.occurrence_id for item in predictor.inputs] == ["seq_2", "seq_3", "seq_4"]
    cohort = report["cohort"]
    assert cohort["full_validation_cohort"] and cohort["scored_records"] == 3
    assert cohort["native_threshold_confusion"] == dict(tp=1, tn=0, fp=1, fn=1)
    group = report["benchmark"]["groups"][0]
    assert group["measurement_count"] == group["comparable_count"] == 3
    assert group["descriptive_metrics"]["brier"] == pytest.approx((0.7**2 + (1-.5493)**2 + .9**2)/3)
    assert group["claim"] == "reproduction" and group["independent_measurement_count"] == 0
    assert group["prediction_interval_coverage"] == []
    assert "sequence" in report["benchmark"]["records"][0]["training_overlap"]


def test_partial_cohort_and_runtime_failure_preserve_every_source_id(synthetic):
    report = cpp.evaluate_cpp_metadata("synthetic", predict=True, source_ids=["seq_4", "seq_2"],
                                       predictor=StubCPP(scores=(.7, .1)))
    assert report["cohort"]["selected_source_ids"] == ["seq_2", "seq_4"]
    assert report["cohort"]["validation_records"] == 3
    assert not report["cohort"]["full_validation_cohort"]
    assert report["cohort"]["status_counts"] == dict(scored=2, not_assessed=1)
    failed = cpp.evaluate_cpp_metadata("synthetic", predict=True, source_ids=["seq_3"],
                                       predictor=StubCPP(failure="missing CPU weights"))
    assert failed["cohort"]["status_counts"] == dict(failed=1, not_assessed=2)
    assert len(failed["benchmark"]["records"]) == 3
    row = next(r for r in failed["benchmark"]["records"] if r["measurement"]["measurement_id"] == "seq_3")
    assert row["prediction"]["reason"] == "missing CPU weights"
    assert failed["cohort"]["native_threshold_confusion"] is None


@pytest.mark.parametrize("ids", [[], ["seq_0"], ["seq_9"], ["seq_2", "seq_2"]])
def test_invalid_cohort_never_loads_predictor(synthetic, ids, monkeypatch):
    monkeypatch.setattr(cpp, "_score_rows", lambda *args: pytest.fail("Invalid selection"))
    with pytest.raises(ValueError, match="validation source IDs"):
        cpp.evaluate_cpp_metadata("synthetic", predict=True, source_ids=ids)


def test_prediction_is_opt_in_and_cli_mode_flags_are_checked(synthetic, tmp_path, capsys):
    with pytest.raises(ValueError, match="explicit prediction"):
        cpp.evaluate_cpp_metadata("synthetic", source_ids=["seq_2"])
    with pytest.raises(ValueError, match="explicit prediction"):
        cpp.evaluate_cpp_metadata("synthetic", predictor=StubCPP())
    out = tmp_path / "report.json"
    main(["--peptiverse-cpp-metadata", "synthetic", "--out", str(out)])
    assert json.loads(out.read_text())["cohort"]["validation_records"] == 3
    for args in (["--lineage-inventory", "--predict-cpp"],
                 ["--lineage-inventory", "--source-id", "seq_2"],
                 ["--peptiverse-cpp-metadata", "synthetic", "--model", "dpp4-qpisa"],
                 ["--peptiverse-cpp-metadata", "synthetic", "--source-id", "seq_2"]):
        with pytest.raises(SystemExit) as error:
            main(args)
        assert error.value.code == 2
    assert "error:" in capsys.readouterr().err


@pytest.mark.parametrize("mode", ["lost_record", "wrong_endpoint", "wrong_identity", "wrong_scale"])
def test_adapter_corruption_invalidates_batch_without_dropping_records(synthetic, mode):
    predictor = StubCPP()
    original = predictor.predict
    def corrupt(inputs, on_unsupported):
        results = original(inputs, on_unsupported)
        if mode == "lost_record":
            return results[:-1]
        original_pred = results[0].preds[0]
        pred = SimpleNamespace(**{key: getattr(original_pred, key) for key in
                                  ("kind", "peptide_input", "peptide", "value", "score", "measurement_context")})
        if mode == "wrong_endpoint":
            pred.kind = Kind.serum_half_life
        elif mode == "wrong_identity":
            pred.peptide_input = inputs[1]
        else:
            pred.measurement_context = MeasurementContext(estimate_type="ml_predicted", status="available", class_label="CPP",
                                                          score_semantics="fraction entering cells")
        results[0] = SimpleNamespace(preds=(pred,))
        return results
    predictor.predict = corrupt
    report = cpp.evaluate_cpp_metadata("synthetic", predict=True, predictor=predictor)
    assert report["cohort"]["status_counts"] == {"failed": 3}
    assert report["benchmark"]["groups"][0]["descriptive_metrics"] is None
    reason = report["benchmark"]["records"][0]["prediction"]["reason"]
    assert reason.startswith("CPP benchmark")


def test_unsupported_native_input_retains_source_label(synthetic):
    predictor = StubCPP()
    def unsupported(inputs, on_unsupported):
        return [SimpleNamespace(preds=(SimpleNamespace(kind=Kind.cpp_classification,
            peptide_input=item, peptide=item.sequence, value=None, score=None,
            measurement_context=MeasurementContext(estimate_type="ml_predicted",
                status="unsupported", detail="Synthetic computational limit")),)) for item in inputs]
    predictor.predict = unsupported
    report = cpp.evaluate_cpp_metadata("synthetic", predict=True, predictor=predictor)
    assert report["cohort"]["status_counts"] == {"unsupported": 3}
    assert report["benchmark"]["records"][0]["measurement"]["value"] == 0
    assert report["benchmark"]["records"][0]["prediction"]["reason"] == "Synthetic computational limit"


@pytest.mark.requires_external_tool
@pytest.mark.skipif(METADATA is None, reason="provision PeptiVerse CPP metadata and CPU runtime")
def test_actual_pinned_metadata_and_partial_native_source_reproduction():
    from mhctools.peptiverse_cpp import PeptiVerseCPP
    rows = cpp.load_cpp_metadata(METADATA)
    audit = cpp._audit(rows)
    summary = json.loads((Path(__file__).parent / "data" /
                          "peptiverse_cpp_validation_summary.json").read_text())
    assert json.loads(json.dumps(audit)) == summary["source_audit"]
    assert summary["cohort"]["status_counts"] == {"scored": 465}
    confusion = summary["cohort"]["native_threshold_confusion"]
    assert confusion["tp"] + confusion["fn"] == audit["splits"]["val"]["labels"]["1"]
    assert confusion["tn"] + confusion["fp"] == audit["splits"]["val"]["labels"]["0"]
    assert audit["records"] == 2324
    assert audit["duplicate_source_ids"] == audit["duplicate_sequences"] == 0
    assert audit["conflicting_label_sequences"] == audit["exact_sequence_overlap_count"] == 0
    train, val = audit["splits"]["train"], audit["splits"]["val"]
    assert (train["records"], train["min_length"], train["max_length"]) == (1859, 3, 61)
    assert (val["records"], val["min_length"], val["max_length"]) == (465, 3, 52)
    assert train["labels"] == {"0": 918, "1": 941}
    assert val["labels"] == {"0": 244, "1": 221}
    assert sum(n for length, n in val["length_counts"].items() if length >= 30) == 56
    ids = [r["id"] for r in rows if r["split"] == "val"][:3]
    report = cpp.evaluate_cpp_metadata(METADATA, predict=True, source_ids=ids,
                                       predictor=PeptiVerseCPP(device="cpu"))
    assert report["cohort"]["status_counts"] == {"scored": 3, "not_assessed": 462}
    assert not report["cohort"]["full_validation_cohort"]
    assert report["artifact_inventory"]["artifact_status"] == "verified"
    assert report["artifact_inventory"]["inference_status"] == "reproduced"
    assert len(report["benchmark"]["records"]) == 465
    assert report["benchmark"]["groups"][0]["claim"] == "reproduction"
    assert all(not r["independence_established"] for r in report["benchmark"]["records"])
