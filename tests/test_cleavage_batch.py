"""Composed ingestion/inference/overlay/save/reload and CLI workflows."""

from copy import deepcopy
import json
import subprocess
import sys

import pytest

from mhctools import (
    CleavageInput, CleavageResult, get_cleavage_model, load_cleavage_batch,
    normalize_cleavage_input, predict_cleavage_batch, write_cleavage_batch,
)
from mhctools.cli.script import main


def test_catalog_and_motif_batch_do_not_import_optional_runtimes():
    completed = subprocess.run([sys.executable, "-c", """
import sys
from mhctools import cleavage_models, predict_cleavage_batch
assert any(m.name.startswith('pepsickle-') for m in cleavage_models(include_optional=True))
report = predict_cleavage_batch(
    [dict(id='vaccine', scope='construct', sequence='RPPGFSPFR', n_term='free', c_term='free')],
    [dict(id='serum', context='extracellular', models=['cpn-basic'])])
assert report['assessments'][0]['status'] == 'assessed'
assert not {'torch', 'tensorflow', 'pepsickle.model_functions'} & sys.modules.keys()
"""], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr


def request():
    return dict(schema_version=1, inputs=[
        dict(id="native-a", peptide="RPPGFSPFR", n_flank="AA", c_flank="GG",
             source_id="protein-A", source_start=12, evidence=[
                 dict(source="imported-table", kind="binding_affinity", value=12.3, units="nM")]),
        dict(id="vaccine", scope="construct", sequence="AARPPGFSPFRGG",
             n_term="free", c_term="free", epitopes=[
                 dict(id="target", start=2, end=11, sequence="RPPGFSPFR")],
             fragments=[dict(id="trimmed", start=2, end=11, n_term="free", c_term="free",
                             assumption="Only if boundary cuts release the target")]),
    ], scenarios=[
        dict(id="tumor-processing", context="tumor", models=["app1-xp"]),
        dict(id="apc-processing", context="apc", compartments=["endosome"], models=["lnpep-observed"]),
        dict(id="serum-incubation", context="extracellular", compartments=["serum"],
             models=["dpp4-qpisa", "cpn-basic", "cpb2-basic"]),
    ])


def run_request(value):
    return predict_cleavage_batch(value["inputs"], value["scenarios"])


def test_mixed_context_evidence_round_trip_and_cli_agree(tmp_path, monkeypatch):
    value = request()
    report = run_request(value)
    path = tmp_path / "evidence.json"
    html_path = tmp_path / "evidence.html"
    write_cleavage_batch(report, path, html_path=html_path)
    assert load_cleavage_batch(path) == report
    assert report["inputs"][0]["evidence"] == value["inputs"][0]["evidence"]
    input_path = tmp_path / "request.json"
    input_path.write_text(json.dumps(value))
    cli_path = tmp_path / "cli.json"
    main(["cleavage", "--input", str(input_path), "--out", str(cli_path)])
    assert load_cleavage_batch(cli_path) == report
    # Reload/report must not execute a predictor.
    monkeypatch.setattr("mhctools.cleavage_batch.get_cleavage_model",
                        lambda *a, **k: pytest.fail("Unexpected inference on reload"))
    assert load_cleavage_batch(path) == report
    text = html_path.read_text()
    assert "Conditional fragment" in text and "matched" in text and "unassessed" in text
    assert "imported-table" in text and "lnpep-observed" in text


def test_native_window_does_not_become_a_free_peptide():
    report = run_request(request())
    native = [r for r in report["assessments"] if r["input_id"] == "native-a"]
    assert all(r["status"] == "unsupported" for r in native)
    assert all(not r["result"]["sites"] for r in native)
    parent = next(r for r in report["assessments"] if r["input_id"] == "vaccine"
                  and r["model"] == "app1-xp" and r["fragment"] is None)
    fragment = next(r for r in report["assessments"] if r["input_id"] == "vaccine"
                    and r["model"] == "app1-xp" and r["fragment"] is not None)
    assert parent["result"]["sites"][0]["status"] == "not_matched"
    assert fragment["result"]["sites"][0]["status"] == "matched"
    overlay = fragment["overlays"][0]
    assert overlay["conditional_on"]["assumption"]
    assert overlay["internal"][0]["bond"] == 3
    assert overlay["internal"][0]["status"] == "matched"
    assert overlay["n_boundary"]["status"] == "unassessed"
    assert overlay["c_boundary"]["status"] == "unassessed"


def test_duplicate_inference_retains_occurrences_and_offsets():
    class Counting:
        def __init__(self):
            self.predictor = get_cleavage_model("cpn-basic")
            self.model = self.predictor.model
            self.calls = []

        def predict(self, peptide):
            self.calls.append(peptide)
            return self.predictor.predict(peptide)

    predictor = Counting()
    inputs = [dict(id=name, scope="construct", sequence="RPPGFSPFR",
                   n_term="free", c_term="free", source_start=offset,
                   epitopes=[dict(id="target", start=0, end=9)])
              for name, offset in (("a", 0), ("b", 30))]
    report = predict_cleavage_batch(inputs, [dict(
        id="serum", context="extracellular", models=["cpn-basic"])],
        predictors={"cpn-basic": predictor})
    assert len(predictor.calls) == 1
    rows = report["assessments"]
    assert [r["result"]["sites"][0]["source_bond"] for r in rows] == [8, 38]
    assert [r["input_id"] for r in rows] == ["a", "b"]
    assert rows[0]["overlays"][0]["n_boundary"]["status"] == "sequence_endpoint"
    assert rows[0]["overlays"][0]["internal"][-1]["score"] is None


def test_failed_backend_is_visible_and_can_raise():
    class Broken:
        model = get_cleavage_model("cpn-basic").model

        def predict(self, peptide):
            raise RuntimeError("weights missing")

    kwargs = dict(inputs=[dict(id="a", sequence="RPPGFSPFR")], scenarios=[dict(
        id="serum", context="extracellular", models=["cpn-basic"])],
        predictors={"cpn-basic": Broken()})
    row = predict_cleavage_batch(**kwargs)["assessments"][0]
    assert row["status"] == "failed" and "weights missing" in row["error"]
    assert row["result"] is None and row["overlays"] == []
    with pytest.raises(RuntimeError, match="weights missing"):
        predict_cleavage_batch(**kwargs, raise_on_error=True)


def test_saved_evidence_corruption_is_rejected(tmp_path):
    report = run_request(request())
    path = tmp_path / "evidence.json"
    write_cleavage_batch(report, path)
    report["assessments"].pop()
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="missing requested"):
        load_cleavage_batch(path)
    result = get_cleavage_model("cpn-basic").predict(CleavageInput("RPPGFSPFR", source_start=12))
    assert CleavageResult.from_dict(result.to_dict()) == result
    corrupted = result.to_dict()
    corrupted["sites"][0]["source_bond"] = 0
    with pytest.raises(ValueError, match="source_bond"):
        CleavageResult.from_dict(corrupted)
    with pytest.raises(ValueError, match="distinct"):
        write_cleavage_batch(report, path, html_path=path)


@pytest.mark.parametrize("change", [
    dict(n_term="free"),
    dict(sequence="RPPGFSPFR", epitopes=[dict(id="x", start=0, end=9, sequence="AAAAAAAAA")]),
    dict(sequence="RPPGFSPFR", fragments=[dict(id="x", start=1, end=8, n_term="free", c_term="free")]),
])
def test_invalid_input_refused_before_prediction(change):
    value = dict(id="a", sequence="RPPGFSPFR")
    value.update(change)
    with pytest.raises(ValueError):
        normalize_cleavage_input(value)


def test_changed_flanks_remain_distinct_inputs():
    a = request()["inputs"][0]
    b = deepcopy(a)
    b.update(id="native-b", c_flank="AA")
    result = predict_cleavage_batch([a, b], [dict(
        id="tumor", context="tumor", models=["app1-xp"])])
    assert result["inputs"][0]["sequence"] != result["inputs"][1]["sequence"]
    assert len(result["assessments"]) == 2


def test_imported_assay_panel_retains_sources_and_abstains_on_unseen_sequences(tmp_path):
    from mhctools._resources import load_json_resource

    panel = deepcopy(load_json_resource("intracellular_substrate_evidence.json")["models"][0])
    panel["model"]["name"] = "imported-assay-panel"
    case = panel["cases"][0]
    inputs = [dict(id="observed", scope="construct", sequence=case["sequence"],
                   n_term=case["n_term"], c_term=case["c_term"]),
              dict(id="unseen", scope="construct", sequence="AAAAAAAAAAAAAA",
                   n_term="free", c_term="free")]
    report = predict_cleavage_batch(inputs, [dict(
        id="assay", context="tumor", models=["imported-assay-panel"])], reference_panels=[panel])
    observed, unseen = report["assessments"]
    assert observed["result"]["substrate_observation"] == case["substrate_observation"]
    assert dict(observed["result"]["conditions"])["source_measurement_id"] == case["source_measurement_id"]
    assert unseen["status"] == "unsupported"
    path = tmp_path / "imported.json"
    write_cleavage_batch(report, path)
    assert load_cleavage_batch(path) == report
    assert report["reference_panels"] == [panel]


def test_batch_predictor_deduplicates_and_rejects_truncated_results():
    class Batched:
        model = get_cleavage_model("cpn-basic").model

        def __init__(self):
            self.calls = []

        def predict_many(self, peptides):
            self.calls.append(peptides)
            return [get_cleavage_model("cpn-basic").predict(p) for p in peptides]

    predictor = Batched()
    inputs = [dict(id=name, sequence="RPPGFSPFR", scope="construct", n_term="free", c_term="free")
              for name in ("a", "b")]
    scenarios = [dict(id="serum", context="extracellular", models=["cpn-basic"])]
    result = predict_cleavage_batch(inputs, scenarios, predictors={"cpn-basic": predictor})
    assert len(predictor.calls) == 1 and len(predictor.calls[0]) == 1
    assert len(result["assessments"]) == 2
    predictor.predict_many = lambda peptides: []
    result = predict_cleavage_batch(inputs, scenarios, predictors={"cpn-basic": predictor})
    assert all(r["status"] == "failed" and "number of results" in r["error"]
               for r in result["assessments"])


def test_reference_panel_name_does_not_change_chemical_form():
    from mhctools._resources import load_json_resource

    panel = deepcopy(load_json_resource("intracellular_substrate_evidence.json")["models"][0])
    panel["model"]["name"] = "pepsickle-source-control"
    case = panel["cases"][0]
    case["n_term"] = case["c_term"] = "unknown"
    panel["cases"] = [case]
    report = predict_cleavage_batch(
        [dict(id="source", sequence=case["sequence"])],
        [dict(id="assay", context="tumor", models=[panel["model"]["name"]])],
        reference_panels=[panel], raise_on_error=True)
    row = report["assessments"][0]
    assert row["status"] == "assessed"
    assert row["result"]["substrate_observation"] == case["substrate_observation"]
    assert row["result"]["peptide"]["n_term"] == "unknown"
    assert "input_scope" not in dict(row["result"]["conditions"])
