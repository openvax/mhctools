"""Assay separation, leakage, missingness and report regressions."""

from dataclasses import replace
import json

import pytest

from mhctools.benchmark import (
    AssayMeasurement, BenchmarkPrediction, ModelLineage,
    evaluate_benchmark, model_lineage_inventory, predict_cleavage_measurements,
)
from mhctools.cli.script import main


def measurement(mid="one", **kwargs):
    values = dict(measurement_id=mid, source_measurement_id=mid,
        source="https://example.org/study", dataset="fixture", study="study",
        assay="serum incubation", sequence="HAEGT", chemistry="linear_L_free",
        endpoint="serum_half_life", units="hours", species="Homo sapiens",
        matrix="serum", value=2.0, family="glp1")
    values.update(kwargs)
    return AssayMeasurement(**values)


def prediction(mid="one", **kwargs):
    values = dict(measurement_id=mid, model="model", endpoint="serum_half_life",
                  units="hours", value=3.0)
    values.update(kwargs)
    return BenchmarkPrediction(**values)


def test_native_metrics_and_repeats_are_not_independent_peptides():
    rows = [measurement(), measurement("two", value=4)]
    pred = [prediction(), prediction("two", value=2)]
    report = evaluate_benchmark(rows, pred)
    group = report["groups"][0]
    assert group["measurement_count"] == 2
    assert group["unique_sequences"] == 1
    assert group["unique_chemical_forms"] == 1
    assert group["descriptive_metrics"] == pytest.approx({"mae":1.5,"rmse":2.5**0.5,"bias":-0.5})
    assert group["claim"] == "unverified_external"
    assert group["independent_measurement_count"] == 0


@pytest.mark.parametrize("changed", [
    {"endpoint":"systemic_half_life"}, {"units":"minutes"},
    {"matrix":"plasma"}, {"species":"Mus musculus"},
    {"assay":"another incubation"}, {"study":"another study"},
    {"conditions":{"pH":"6.0"}}, {"cell_type":"DC"},
])
def test_assays_and_domains_are_not_pooled(changed):
    rows = [measurement(), measurement("two", **changed)]
    result = evaluate_benchmark(rows, [prediction(), prediction("two")])
    assert len(result["groups"]) == 2
    assert all(g["measurement_count"] == 1 for g in result["groups"])


def test_no_unit_or_endpoint_conversion():
    for changes in ({"units":"minutes"}, {"endpoint":"fluorescence_uptake"}):
        result = evaluate_benchmark([measurement()], [prediction(**changes)])
        assert result["records"][0]["exclusion"] == "incompatible_endpoint_or_units"
        assert result["groups"][0]["descriptive_metrics"] is None


def test_censored_missing_failed_and_absent_predictions_remain_visible():
    rows = [measurement(), measurement("two", censoring="greater_than"),
            measurement("three"), measurement("four")]
    preds = [prediction(status="unsupported", value=None, reason="chemistry"),
             prediction("two"), prediction("three", status="failed", value=None, reason="runtime")]
    report = evaluate_benchmark(rows, preds)
    assert [r["exclusion"] for r in report["records"]] == [
        "unsupported", "unknown_or_censored_observation", "failed", "missing_prediction"]
    assert sum(g["measurement_count"] for g in report["groups"]) == 4
    assert all(g["descriptive_metrics"] is None for g in report["groups"])


def test_complete_provenance_is_required_and_overlap_is_reported():
    lineage = ModelLineage("model", "https://example.org/training", "complete",
        studies=("training-study",), sequences=("AAAA",), families=("training-family",),
        chemical_forms=(("AAAA","linear_L_free"),), cell_types=("HeLa",))
    clean = evaluate_benchmark([measurement()], [prediction()], [lineage])
    assert clean["groups"][0]["claim"] == "audited_external"
    row = measurement(sequence="AAAA", chemistry="linear_L_C_amidated", family="training-family", cell_type="HeLa")
    dirty = evaluate_benchmark([row], [prediction()], [lineage])
    assert set(dirty["records"][0]["training_overlap"]) == {"sequence","family","cell_type"}
    assert "chemical_form" not in dirty["records"][0]["training_overlap"]
    assert dirty["groups"][0]["claim"] == "unverified_external"
    no_family = evaluate_benchmark([measurement(family=None)], [prediction()], [lineage])
    assert no_family["groups"][0]["claim"] == "unverified_external"
    with pytest.raises(ValueError, match="inventories"):
        ModelLineage("model", "source", "complete")


def test_split_overlap_and_reference_reproduction():
    train = measurement("train", split="train")
    test = measurement()
    report = evaluate_benchmark([train,test], [prediction()])
    assert {"sequence","study","assay","family","chemical_form"} <= set(report["records"][0]["partition_overlap"])
    ref = measurement(split="reference")
    external = evaluate_benchmark([ref], [prediction()])
    assert external["records"][0]["exclusion"] == "reference_record"
    reproduced = evaluate_benchmark([ref], [prediction()], evaluation="reproduction")
    assert reproduced["groups"][0]["claim"] == "reproduction"
    assert reproduced["groups"][0]["comparable_count"] == 1


def test_native_site_score_never_becomes_a_probability():
    row = measurement(endpoint="site_cleavage", units="binary", enzyme="DPP4", bond=2, value=1)
    predictions = predict_cleavage_measurements([row], ["dpp4-qpisa"])
    assert predictions[0].value == pytest.approx(2.1694)
    report = evaluate_benchmark([row], predictions)
    assert report["groups"][0]["descriptive_metrics"] is None
    assert report["records"][0]["exclusion"] == "incompatible_endpoint_or_units"


def test_binary_decisions_and_probabilities_keep_distinct_metrics():
    rows = [measurement(endpoint="site_cleavage", units="binary", enzyme="CPN1", bond=4, value=1),
            measurement("two", endpoint="site_cleavage", units="binary", enzyme="CPN1", bond=4, value=0)]
    decisions = [prediction(endpoint="site_cleavage", units="binary", scale="decision", value=1),
                 prediction("two", endpoint="site_cleavage", units="binary", scale="decision", value=1)]
    assert evaluate_benchmark(rows,decisions)["groups"][0]["descriptive_metrics"] == {"tp":1,"fp":1,"tn":0,"fn":0}
    probabilities = [replace(p,scale="probability",value=0.8) for p in decisions]
    assert evaluate_benchmark(rows,probabilities)["groups"][0]["descriptive_metrics"]["brier"] == pytest.approx(0.34)


def test_only_explicit_measurement_intervals_are_evaluated():
    p = prediction(interval=(1,4), interval_level=0.9, interval_target="individual_measurement")
    group = evaluate_benchmark([measurement()],[p])["groups"][0]
    assert group["prediction_interval_coverage"] == [{"nominal_level":0.9,"n":1,"coverage":1.0}]
    with pytest.raises(ValueError, match="target"):
        prediction(interval=(1,4), interval_level=0.9, interval_target="model_disagreement")


def test_unknown_chemistry_and_unassessed_bond_do_not_become_negatives():
    row = measurement(sequence="RPPGFSPFR", chemistry="cyclic",endpoint="site_cleavage",units="binary",enzyme="CPN1",bond=8,value=1)
    assert predict_cleavage_measurements([row],["cpn-basic"])[0].status == "unsupported"
    row = replace(row,chemistry="linear_L_free",bond=4)
    p = predict_cleavage_measurements([row],["cpn-basic"])[0]
    assert p.status == "not_assessed" and p.value is None


def test_identity_errors_are_rejected():
    with pytest.raises(ValueError, match="Duplicate measurement"):
        evaluate_benchmark([measurement(),measurement()],[prediction()])
    with pytest.raises(ValueError, match="Duplicate measurement/model"):
        evaluate_benchmark([measurement()],[prediction(),prediction()])
    with pytest.raises(ValueError, match="unknown measurement"):
        evaluate_benchmark([measurement()],[prediction("missing")])
    for value in (float("nan"),float("inf"),True):
        with pytest.raises(ValueError, match="finite"):
            prediction(value=value)


def test_reference_cli_and_lineage_inventory(capsys):
    main(["benchmark","--reference-cleavage"])
    data = json.loads(capsys.readouterr().out)
    assert data["evaluation"] == "reproduction"
    assert len(data["records"]) == 12
    assert sum(g["comparable_count"] for g in data["groups"]) == 4
    assert all(g["claim"] == "reproduction" for g in data["groups"])
    main(["benchmark","--lineage-inventory"])
    inventory = json.loads(capsys.readouterr().out)
    assert inventory == model_lineage_inventory()
    assert any(m["model"] == "pepADMET" and m["provenance"] == "unknown" for m in inventory["models"])
    assert all(m["source"] and m["unresolved"] for m in inventory["models"])


def test_requested_domains_report_absent_evidence_without_threshold_invention():
    domains = [{"name":"human serum","matrix":"serum","species":"Homo sapiens"},
               {"name":"requested long peptides","min_length":25},
               {"name":"primary DC","cell_type":"primary dendritic cell"}]
    report = evaluate_benchmark([measurement()],[prediction()],requested_domains=domains)
    assert [r["measurement_count"] for r in report["requested_domains"]] == [1,0,0]
    assert report["requested_domains"][1]["evidence"] == "no_evidence_in_supplied_data"
    assert report["requested_domains"][1]["requested"]["min_length"] == 25
