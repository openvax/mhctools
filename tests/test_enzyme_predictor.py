"""One enzyme entry point without mixing motifs, scores, and assay fractions."""

import json

import pytest

from mhctools import DPP4, FAP, EnzymePredictor, PeptideInput, enzyme_predictors, enzymes


PEPTIDE = "HAE" + "G" * 27


def test_dpp4_catalog_and_distinct_predictions():
    predictor = DPP4()
    assert predictor.enzyme == "DPP4"
    methods = predictor.describe()["methods"]
    assert [method["name"] for method in methods] == [
        "dpp4-qpisa", "dpp4-p1-pro-ala", "dpp4-p1prime-nonpro",
        "dpp4-qpisa-keane2011-parent-loss"]
    assert [method["evidence_type"] for method in methods] == [
        "quantitative_model", "motif_rule", "motif_rule", "assay_estimate"]
    assert methods[1]["activity"] == "Removes 2 N-terminal residue(s)"
    assert methods[1]["strictness"] == "preferred"
    assert "10.1515/bc.2011.028" in methods[1]["references"][0]
    assert "10.1038/s44320-024-00071-4" in methods[2]["references"][0]

    result = predictor.predict(PEPTIDE, minutes=10)
    assert result.by_model["dpp4-qpisa"].sites[0].score == pytest.approx(2.1694)
    assert result.by_model["dpp4-p1-pro-ala"].pattern_mask[1] == 1
    assert result.by_model["dpp4-p1prime-nonpro"].pattern_mask[1] == 1
    assert result.by_model["dpp4-qpisa-keane2011-parent-loss"]["parent_lost_fraction"] == pytest.approx(
        0.3325, abs=0.001)
    assert [row["endpoint"] for row in result.summary_rows()] == [
        "substrate_depletion", "recognition_pattern", "recognition_pattern",
        "parent_peptide_loss_fraction"]
    assert "33.2% parent lost" in result.summary()
    assert "X-(Pro/Ala)|X preference" in result.summary()
    serialized = result.to_dict()
    assert serialized["schema_version"] == 1
    assert len(serialized["cleavage_results"]) == 3
    assert len(serialized["assay_fractions"]) == 1
    json.dumps(serialized)


def test_dpp4_preference_nonmatch_does_not_gate_scored_model():
    result = DPP4().predict("HSE" + "G" * 27, minutes=10)
    assert result.by_model["dpp4-p1-pro-ala"].sites[0].match == 0
    assert result.by_model["dpp4-p1prime-nonpro"].sites[0].match == 1
    assert result.by_model["dpp4-qpisa"].sites[0].score == pytest.approx(0.6744)
    assert result.assay_fractions[0]["status"] == "unassessed"
    assert "parent_lost_fraction" not in result.assay_fractions[0]
    pro_prime = DPP4().predict("HAP" + "G" * 27)
    assert pro_prime.by_model["dpp4-p1prime-nonpro"].sites[0].match == 0
    assert pro_prime.by_model["dpp4-qpisa"].sites[0].score == pytest.approx(-0.2018)


def test_fap_multiple_activities_and_assay_estimate():
    predictor = FAP()
    assert [method["activity"] for method in predictor.describe()["methods"][:2]] == [
        "Removes 2 N-terminal residue(s)", "Cuts an internal peptide bond"]
    result = predictor.predict(PEPTIDE, minutes=10)
    assert result.by_model["fap-dipeptidyl"].pattern_mask[1] == 1
    assert sum(result.by_model["fap-endo-gp"].pattern_mask) == 0
    assert result.by_model["fap-keane2011-parent-loss"]["parent_lost_fraction"] == pytest.approx(
        0.005508, rel=1e-3)


def test_unsupported_chemistry_and_enzymes_without_assay():
    modified = DPP4().predict(
        PeptideInput(PEPTIDE, n_term="acetylated"), minutes=10)
    assert all(result.unsupported_reason for result in modified.cleavage_results)
    assert modified.assay_fractions[0]["status"] == "unassessed"
    assert all(row["status"] == "unassessed" for row in modified.summary_rows())

    ace = EnzymePredictor.from_name("ACE").predict(PEPTIDE, minutes=10)
    assert not ace.assay_fractions
    assert "No assay-calibrated parent fraction" in ace.summary()
    assert "ACE" in EnzymePredictor.available()
    with pytest.raises(ValueError, match="Unknown enzyme"):
        EnzymePredictor.from_name("not-an-enzyme")
    with pytest.raises(ValueError, match="minutes"):
        EnzymePredictor.from_name("ACE").predict(PEPTIDE, minutes=-1)


def test_activation_state_is_forwarded_to_existing_model():
    predictor = EnzymePredictor.from_name("CPB2")
    unselected = predictor.predict("YFPGQFAFSK")
    assert unselected.cleavage_results[0].unsupported_reason
    active = predictor.predict("YFPGQFAFSK", enzyme_state="active")
    assert not active.cleavage_results[0].unsupported_reason
    assert active.cleavage_results[0].conditions == (("enzyme_state", "active"),)


def test_concrete_classes_and_configuration_registry():
    assert enzyme_predictors["dpp4"] is DPP4
    assert enzyme_predictors["fap"] is FAP
    assert isinstance(EnzymePredictor.from_name("DPP4"), DPP4)
    assert isinstance(EnzymePredictor.from_name("fap"), FAP)
    assert EnzymePredictor(DPP4).enzyme == "DPP4"
    assert EnzymePredictor(DPP4()).enzyme == "DPP4"
    assert EnzymePredictor("DPP4").enzyme == "DPP4"  # 4.0.5 compatibility
    assert EnzymePredictor.from_name("ace").enzyme == "ACE"
    assert isinstance(enzymes.DPP4, DPP4)
    assert enzymes.ACE.enzyme == "ACE"
    assert enzymes.CathepsinS.enzyme == "Cathepsin S"
    assert "CathepsinS" in dir(enzymes)
    with pytest.raises(AttributeError, match="No cataloged enzyme"):
        _ = enzymes.NotAnEnzyme
    with pytest.raises(TypeError, match="already selects"):
        DPP4("FAP")
