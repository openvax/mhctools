"""Source-observation, abstention and provenance regressions for intracellular candidates."""

from dataclasses import replace
from importlib.resources import files
import json

import pytest

from mhctools import CleavageInput, cleavage_models, get_cleavage_model, predict_cleavage
from mhctools.benchmark import AssayMeasurement, evaluate_benchmark, predict_cleavage_measurements
from mhctools.cli.script import main


def _panel():
    return json.loads(files("mhctools").joinpath(
        "data/intracellular_cleavage_reference.json").read_text())


def test_source_reference_returns_observations_without_extrapolating():
    thop1 = get_cleavage_model("thop1-observed")
    result = thop1.predict(CleavageInput("GPLGPLGPL", source_start=4))
    assert [(s.bond, s.status) for s in result.sites] == [(3, "reported"), (6, "reported")]
    assert result.substrate_observation == "cleavage_reported"
    assert result.model.evidence == "substrate_reference"
    assert [s["source_bond"] for s in result.to_dict()["sites"]] == [7, 10]
    assert dict(result.conditions)["source_measurement_id"] == "knight1995-thop1-gpl3"
    # Bradykinin is cleaved at the Phe-Ser bond, giving RPPGF and SPFR.
    assert [s.bond for s in thop1.predict("RPPGFSPFR").sites] == [5]
    # One residue longer is a different chemical form, so the catalog abstains.
    unseen = thop1.predict("GPLGPLGPLG")
    assert unseen.unsupported_reason and not unseen.sites
    assert unseen.substrate_observation is None
    # Modified termini are a different chemical form too.
    assert thop1.predict(CleavageInput("RPPGFSPFR", n_term="acetylated")).unsupported_reason


def test_reported_non_cleavage_never_becomes_a_site_label():
    nln = get_cleavage_model("nln-observed")
    resistant = nln.predict("YGGFLRRIRPKLK")
    assert resistant.substrate_observation == "no_cleavage_detected"
    assert not resistant.sites
    assert resistant.unsupported_reason is None
    assert [s.bond for s in nln.predict("YGGFLRRI").sites] == [5]
    # Degradation was reported without a bond being pinned down, so no bond is invented.
    degraded = nln.predict("YGGFLRRIR")
    assert degraded.substrate_observation == "cleavage_reported" and not degraded.sites
    lnpep = get_cleavage_model("lnpep-observed")
    assert [s.bond for s in lnpep.predict("KSLYNTVATL").sites] == [1]
    assert lnpep.predict("SLYNTVATL").substrate_observation == "no_cleavage_detected"


def test_disputed_source_sequence_is_excluded_from_curation():
    # Georgiadou 2010 prints two different resistant precursors; see mhctools issue 332.
    catalog = json.loads(files("mhctools").joinpath(
        "data/intracellular_substrate_evidence.json").read_text())
    panel = _panel()
    curated = {c["sequence"] for entry in catalog["models"] for c in entry["cases"]}
    curated |= {m["sequence"] for m in panel["measurements"]}
    assert not curated & {"DIRSSVQNKL", "DIRSSQVNKL"}
    # Both conflicting source strings stay on the record as an explicit exclusion.
    lnpep = next(e for e in catalog["models"] if e["model"]["enzyme"] == "LNPEP")
    for disputed in ("DIRSSVQNKL", "DIRSSQVNKL"):
        assert disputed in lnpep["model"]["limitations"] or disputed in panel["notice"]
    lnpep = get_cleavage_model("lnpep-observed")
    assert lnpep.predict("DIRSSVQNKL").unsupported_reason


def test_compartments_separate_cytosolic_and_endosomal_enzymes():
    cytosolic = {r.model.enzyme for r in predict_cleavage("VPYGSFKHV", compartment="cytosol")}
    assert {"THOP1", "NLN"} <= cytosolic and "LNPEP" not in cytosolic
    endosomal = predict_cleavage("KSLYNTVATL", compartment="endosome")
    assert [r.model.enzyme for r in endosomal] == ["LNPEP"]
    assert [s.bond for s in endosomal[0].sites] == [1]
    with pytest.raises(ValueError, match="not annotated"):
        predict_cleavage("KSLYNTVATL", models="thop1-observed", compartment="endosome")


def test_source_observations_refuse_external_validation_claims():
    measurement = next(m for m in _panel()["measurements"]
                       if m["measurement_id"] == "knight1995-thop1-gpl3-bond3")
    held_out = AssayMeasurement(**dict(measurement, split="test",
                                       measurement_id="held-out"))
    prediction = predict_cleavage_measurements([held_out], ["thop1-observed"])[0]
    assert prediction.value is None
    assert "reproduction" in prediction.reason


def test_intracellular_reference_reproduces_every_source_observation():
    data = _panel()
    measurements = [AssayMeasurement(**m) for m in data["measurements"]]
    predictions = predict_cleavage_measurements(measurements, data["models"])
    report = evaluate_benchmark(measurements, predictions, evaluation="reproduction")
    comparable = [r for r in report["records"] if r["exclusion"] is None]
    assert len(measurements) == 42 and len(comparable) == 41
    assert all(r["measurement"]["value"] == r["prediction"]["value"] for r in comparable)
    assert sum(r["measurement"]["value"] == 0 for r in comparable) == 6
    assert all(g["claim"] == "reproduction" for g in report["groups"])
    # Amidated substance P is outside the aminopeptidase P rule's documented input domain.
    amidated = [r for r in report["records"]
                if r["measurement"]["measurement_id"] == "cottrell2000-xpnpep1-substance-p-bond1"
                and r["model"] == "app1-xp"]
    assert len(amidated) == 1 and amidated[0]["exclusion"] == "not_assessed"
    # Whole-substrate labels must stay binary detection calls, not fractional depletion.
    depletion = next(r for r in comparable
                     if r["measurement"]["endpoint"] == "substrate_depletion")
    assert depletion["prediction"]["scale"] == "decision"
    assert "not fractional depletion" in depletion["prediction"]["reason"]


def test_motif_rules_grade_strictness_and_cite_sources():
    graded = {m.name: m.motif_strictness for m in cleavage_models(include_optional=True)
              if m.evidence == "motif_rule"}
    assert graded["tpp2-tripeptidyl"] == "permissive"
    assert graded["npepps-n-terminal"] == "permissive"
    assert graded["app1-xp"] == "required"
    assert graded["mme-hydrophobic"] == "preferred"
    for model in cleavage_models(include_optional=True):
        assert model.references, model.name
        if model.evidence == "motif_rule":
            assert model.motif_strictness in ("required", "preferred", "permissive"), model.name
            assert len(model.strictness_basis) > 40, model.name
        else:
            assert model.motif_strictness is None and model.strictness_basis is None, model.name


def test_intracellular_cli_reports_observations_and_missing_domains(capsys):
    main(["cleavage", "--sequence", "YGGFLRRI", "--model", "nln-observed"])
    result = json.loads(capsys.readouterr().out)["results"][0]
    assert result["substrate_observation"] == "cleavage_reported"
    assert result["sites"][0]["bond"] == 5
    assert result["model"]["motif_strictness"] is None
    main(["cleavage", "--list-models"])
    listed = json.loads(capsys.readouterr().out)["models"]
    assert {m["name"] for m in listed} >= {"thop1-observed", "nln-observed", "lnpep-observed"}
    assert all(m["strictness_basis"] for m in listed if m["evidence"] == "motif_rule")
    main(["benchmark", "--reference-cleavage", "intracellular"])
    report = json.loads(capsys.readouterr().out)
    assert report["evaluation"] == "reproduction"
    assert "issue 332" in report["dataset_notice"]
    absent = [d for d in report["requested_domains"] if d["measurement_count"] == 0]
    assert {d["requested"]["name"] for d in absent} == {
        "antigen presentation in primary dendritic cells",
        "cleavage measured in cytosol rather than purified enzyme"}


def test_reference_models_abstain_rather_than_reporting_absence():
    for model in cleavage_models():
        if model.evidence != "substrate_reference":
            continue
        predictor = get_cleavage_model(model.name)
        result = predictor.predict("WWWWWWWW")
        assert result.unsupported_reason and not result.sites
        assert result.substrate_observation is None
        with pytest.raises(ValueError, match="Explicit enzyme state"):
            get_cleavage_model(model.name, enzyme_state="active")
    with pytest.raises(ValueError, match="Substrate observations require"):
        replace(get_cleavage_model("thop1-observed").predict("RPPGFSPFR"),
                model=get_cleavage_model("prep-pro").model, sites=())
