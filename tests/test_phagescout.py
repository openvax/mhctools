"""Deposited native-score conformance, sparse coverage and bond coordinates."""

import json
from pathlib import Path

import pytest

from mhctools import (
    CleavageInput, CleavageResult, PhageScout, cleavage_models,
    get_cleavage_model, predict_cleavage, predict_cleavage_batch,
)
from mhctools.cli.cleavage import main
from mhctools.phagescout import PHAGESCOUT_MODELS


REFERENCE = json.loads((Path(__file__).parent / "data/phagescout_source_reference.json").read_text())
FEATURES = {
    "pwm-deseq2": "pwm_deseq2",
    "pwm-relaxed-unaligned": "pwm_relaxed_unaligned",
    "pwm-relaxed-aligned": "pwm_relaxed_aligned",
    "peptide-relaxed-unaligned": "pep_relaxed_unaligned",
    "peptide-relaxed-aligned": "pep_relaxed_aligned",
}


@pytest.mark.parametrize("name,settings", PHAGESCOUT_MODELS.items())
def test_reproduces_deposited_scores_and_missing_values(name, settings):
    predictor = get_cleavage_model(name)
    source = "elastase" if settings["enzyme"] == "ELANE" else "cathepsin G"
    feature = FEATURES[settings["profile"]]
    for record in REFERENCE["proteins"]:
        peptide = CleavageInput(record["sequence"], source_id=record["accession"],
                               source_start=record["source_start"])
        result = predictor.predict(peptide)
        rows = record["expected"][source]
        expected = {row["pos"]: row[feature] for row in rows if row[feature] is not None}
        actual = {site["source_bond"]: site["score"] for site in result.to_dict()["sites"]}
        assert actual == pytest.approx(expected, abs=1e-10)
        gaps = json.loads(dict(result.conditions)["unassessed_bonds"])
        assert {int(bond) + record["source_start"] for bond in gaps} == {
            row["pos"] for row in rows if row[feature] is None}
        assert CleavageResult.from_dict(result.to_dict()) == result
        assert result.model.scored_endpoint == "site_cleavage"
        assert dict(result.conditions)["normalization"] == "none"


def test_five_mer_feature_spans_four_candidate_bonds_without_assigning_a_cut():
    result = PhageScout("ELANE", "peptide-relaxed-unaligned").predict("AAFIF")
    assert [s.bond for s in result.sites] == [1, 2, 3, 4]
    assert [s.score for s in result.sites] == pytest.approx([2.95272845823985] * 4)
    assert "does not experimentally identify" in result.model.limitations


def test_aligned_wildcard_ties_preserve_source_order_across_masks():
    # Released AASFV---- row 2 wins over -ASFVA--- row 159.
    result = PhageScout("ELANE", "peptide-relaxed-aligned").predict("AASFVAAAA")
    assert [(s.bond, s.score) for s in result.sites] == [(5, 1.9929172627577)]
    gaps = json.loads(dict(result.conditions)["unassessed_bonds"])
    assert set(gaps) == {"1", "2", "3", "4", "6", "7", "8"}
    assert set(gaps.values()) == {"Requires complete nine-mer context"}


def test_cached_source_weights_cannot_be_mutated_across_predictors():
    pwm = PhageScout()
    with pytest.raises(TypeError):
        pwm.weights["A"] = (0,) * 9
    profile = PhageScout("ELANE", "peptide-relaxed-unaligned")
    with pytest.raises(TypeError):
        profile.weights[(0, 1, 2, 3, 4)]["AAFIF"] = (0, 0)


@pytest.mark.parametrize("profile", ["peptide-relaxed-unaligned", "peptide-relaxed-aligned"])
def test_missing_profile_matches_do_not_become_zero_or_resistance(profile):
    result = PhageScout("ELANE", profile).predict("PPPPPPPPP")
    assert result.sites == ()
    assert result.unsupported_reason
    assert "No matching released peptide profile" in json.loads(
        dict(result.conditions)["unassessed_bonds"]).values()


@pytest.mark.parametrize("name", PHAGESCOUT_MODELS)
@pytest.mark.parametrize("peptide", [
    CleavageInput("AAFLDAAFF", n_term="unknown"),
    CleavageInput("AAFLDAAFF", n_term="acetylated"),
    CleavageInput("AAFLDAAFF", c_term="amidated"), "A",
])
def test_unestablished_chemistry_and_missing_bonds_abstain(name, peptide):
    result = get_cleavage_model(name).predict(peptide)
    assert result.unsupported_reason and not result.sites


def test_short_contexts_and_truncated_pwm_flanks_have_distinct_reasons():
    pwm = PhageScout("CTSG").predict("AA")
    assert [s.bond for s in pwm.sites] == [1]
    assert "missing flanks omitted" in pwm.sites[0].reason
    five_mer = PhageScout("CTSG", "pwm-deseq2").predict("AAAA")
    assert five_mer.unsupported_reason and not five_mer.sites
    assert set(json.loads(dict(five_mer.conditions)["unassessed_bonds"]).values()) == {
        "Requires a complete five-mer window"}
    sparse = PhageScout("CTSG", "peptide-relaxed-unaligned").predict("AAAA")
    assert set(json.loads(dict(sparse.conditions)["unassessed_bonds"]).values()) == {
        "Requires a complete five-mer window"}


def test_discovery_selection_and_cli_preserve_native_endpoint(capsys):
    models = {m.name: m for m in cleavage_models() if m.name in PHAGESCOUT_MODELS}
    assert set(models) == set(PHAGESCOUT_MODELS)
    assert len({m.version for m in models.values()}) == 10
    name = "phagescout-elane-peptide-relaxed-unaligned"
    result, = predict_cleavage("AAFIF", models=[name], compartment="extracellular")
    main(["--sequence", "AAFIF", "--model", name])
    report = json.loads(capsys.readouterr().out)
    assert report["results"][0]["sites"][0]["score"] == pytest.approx(result.sites[0].score, rel=1e-5)
    assert report["models"][name]["score_units"] == result.model.score_units
    with pytest.raises(ValueError, match="enzyme state"):
        get_cleavage_model(name, enzyme_state="active")


def test_batch_epitope_overlay_uses_same_bonds_and_keeps_sparse_gaps():
    name = "phagescout-elane-peptide-relaxed-aligned"
    report = predict_cleavage_batch(
        [dict(id="synthetic", scope="construct", sequence="AASFVAAAA", n_term="free", c_term="free",
              epitopes=[dict(id="core", start=1, end=8, sequence="ASFVAAA")])],
        [dict(id="inflammation", context="extracellular", models=[name])])
    assessment, = report["assessments"]
    scored = [site for site in assessment["overlays"][0]["internal"] if site["status"] == "scored"]
    assert [(s["bond"], s["score"]) for s in scored] == [(5, 1.9929172627577)]
    assert assessment["overlays"][0]["n_boundary"]["status"] == "unassessed"


@pytest.mark.parametrize("enzyme,profile", [("ELANE", "all"), ("elastase", "pwm-relaxed-aligned"),
                                           ("PRTN3", "pwm-relaxed-aligned")])
def test_no_silent_enzyme_or_classifier_substitution(enzyme, profile):
    with pytest.raises(ValueError):
        PhageScout(enzyme, profile)
