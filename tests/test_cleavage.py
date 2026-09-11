"""Coordinate, chemistry and evidence contract regressions."""

from dataclasses import FrozenInstanceError, replace
import json

import pytest

from mhctools import CleavageInput, CleavageResult, CleavageSite, DPP4qPISA


@pytest.mark.parametrize("sequence", ["", "hap", " HAP", "HAP\n", "HAX", "HAP-amide", None])
def test_exact_canonical_identity(sequence):
    with pytest.raises(ValueError, match="canonical"):
        CleavageInput(sequence)


def test_immutable_identity_and_explicit_fragment_chemistry():
    parent = CleavageInput("AAHAEG", n_term="acetylated", source_id="construct", source_start=40)
    with pytest.raises(FrozenInstanceError):
        parent.sequence = "HAEG"
    with pytest.raises(TypeError):
        parent.fragment(2, 6)
    with pytest.raises(ValueError, match="preserve chemistry"):
        parent.fragment(0, 6, n_term="free", c_term="free")
    fragment = parent.fragment(2, 6, n_term="free", c_term="free")
    result = DPP4qPISA().predict(fragment)
    assert result.peptide.sequence == "HAEG"
    assert result.peptide.source_start == 42
    assert result.to_dict()["sites"][0]["source_bond"] == 44
    assert result.sites[0].bond == 2
    assert json.loads(json.dumps(result.to_dict()))["peptide"]["source_id"] == "construct"


@pytest.mark.parametrize("bond", [0, -1, True, 1.5, 3, 4])
def test_only_real_peptide_bonds(bond):
    with pytest.raises(ValueError):
        CleavageResult(CleavageInput("HAP"), DPP4qPISA.model,
                       (CleavageSite(bond, "scored", "fixture", 2.0),))


def test_evidence_types_do_not_convert_to_probabilities():
    with pytest.raises(ValueError, match="numerical"):
        CleavageSite(1, "matched", "rule", 1.0)
    for value in (float("nan"), float("inf"), True, None):
        with pytest.raises(ValueError, match="finite"):
            CleavageSite(1, "scored", "model", value)
    with pytest.raises(ValueError, match="semantics"):
        CleavageResult(CleavageInput("HAP"), DPP4qPISA.model,
                       (CleavageSite(1, "matched", "rule"),))
    motif = replace(DPP4qPISA.model, evidence="motif_rule", score_name=None, score_units=None,
                    scored_endpoint=None, motif_strictness="required", strictness_basis="fixture")
    with pytest.raises(ValueError, match="semantics"):
        CleavageResult(CleavageInput("HAP"), motif,
                       (CleavageSite(1, "scored", "model", 2.0),))


def test_duplicate_and_unsupported_sites_rejected():
    site = CleavageSite(2, "scored", "fixture", 2.0)
    with pytest.raises(ValueError, match="distinct"):
        CleavageResult(CleavageInput("HAP"), DPP4qPISA.model, [site, site])
    with pytest.raises(ValueError, match="Unsupported"):
        CleavageResult(CleavageInput("HAP"), DPP4qPISA.model, [site], "unknown chemistry")


def test_motif_strictness_is_graded_and_attributable():
    motif = replace(DPP4qPISA.model, evidence="motif_rule", score_name=None, score_units=None,
                    scored_endpoint=None, motif_strictness="preferred", strictness_basis="source preference")
    assert motif.motif_strictness == "preferred"
    with pytest.raises(ValueError, match="required, preferred or permissive"):
        replace(motif, motif_strictness="strict")
    with pytest.raises(ValueError, match="required, preferred or permissive"):
        replace(motif, motif_strictness=None)
    with pytest.raises(ValueError, match="observation behind it"):
        replace(motif, strictness_basis="")
    # A grade would be meaningless on evidence that is not a recognition pattern.
    for evidence in ("quantitative_model", "substrate_reference"):
        with pytest.raises(ValueError, match="Only motif rules"):
            replace(DPP4qPISA.model, evidence=evidence,
                    score_name=DPP4qPISA.model.score_name if evidence == "quantitative_model" else None,
                    score_units=DPP4qPISA.model.score_units if evidence == "quantitative_model" else None,
                    scored_endpoint=DPP4qPISA.model.scored_endpoint if evidence == "quantitative_model" else None,
                    motif_strictness="required", strictness_basis="fixture")


def test_every_model_cites_a_source():
    with pytest.raises(ValueError, match="cite at least one source"):
        replace(DPP4qPISA.model, references=())


def test_scored_endpoint_names_the_benchmark_endpoint_for_quantitative_models():
    assert DPP4qPISA.model.scored_endpoint == "substrate_depletion"
    with pytest.raises(ValueError, match="benchmark endpoint their score answers"):
        replace(DPP4qPISA.model, scored_endpoint=None)
    with pytest.raises(ValueError, match="benchmark endpoint their score answers"):
        replace(DPP4qPISA.model, scored_endpoint="")
    # Only a quantitative model may carry one; motif rules and source references cannot.
    with pytest.raises(ValueError, match="Only quantitative models carry a benchmark endpoint"):
        replace(DPP4qPISA.model, evidence="substrate_reference", score_name=None, score_units=None,
                motif_strictness=None, strictness_basis=None)


def test_coerce_peptide_is_shared_by_every_predictor():
    # Every predictor's predict() routes string/CleavageInput coercion through
    # one shared helper, so the accepted-input contract is identical everywhere.
    from mhctools.cleavage import coerce_peptide
    peptide = coerce_peptide("HAP")
    assert isinstance(peptide, CleavageInput) and peptide.sequence == "HAP"
    same = coerce_peptide(peptide)
    assert same is peptide
    with pytest.raises(TypeError, match="CleavageInput or canonical peptide string"):
        coerce_peptide(("H", "A", "P"))
    with pytest.raises(TypeError, match="CleavageInput or canonical peptide string"):
        coerce_peptide(None)
    for predictor in (DPP4qPISA(),):
        assert predictor.predict("HAE").sites[0].bond == 2


def test_non_quantitative_evidence_rejects_numerical_scores():
    # Pre-existing CleavageModel.__post_init__ branch, previously untested:
    # a motif-rule model must not carry score fields either.
    from mhctools import get_cleavage_model
    motif_model = get_cleavage_model("app1-xp").model
    with pytest.raises(ValueError, match="does not have numerical scores"):
        replace(motif_model, score_name="x")
    with pytest.raises(ValueError, match="does not have numerical scores"):
        replace(motif_model, score_units="x")


def test_substrate_observation_invariants_are_enforced():
    from mhctools import get_cleavage_model
    reference_model = get_cleavage_model("thop1-observed").model
    peptide = CleavageInput("RPPGFSPFR")
    # A whole-substrate non-cleavage cannot also carry a per-bond site label.
    with pytest.raises(ValueError, match="Whole-substrate non-cleavage must not create site labels"):
        CleavageResult(peptide, reference_model, (CleavageSite(5, "reported", "fixture"),),
                       substrate_observation="no_cleavage_detected")
    # An unsupported result (no exact chemical-form match) cannot carry an observation either.
    with pytest.raises(ValueError, match="Unsupported inputs cannot carry substrate observations"):
        CleavageResult(peptide, reference_model, unsupported_reason="no match",
                       substrate_observation="cleavage_reported")
