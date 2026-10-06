"""Author-source conformance and terminal topology of bundled real profiles."""

import json
from pathlib import Path

import pytest

from mhctools import CleavageInput, CleavageResult, ITCellCleavage, cleavage_models, get_cleavage_model
from mhctools.itcell_cleavage import ITCELL_MODELS


REFERENCE = json.loads((Path(__file__).parent / "data/itcell_source_reference.json").read_text())


@pytest.mark.parametrize("name", ITCELL_MODELS)
def test_all_real_profiles_reproduce_released_perl_scores(name):
    predictor = get_cleavage_model(name)
    for record in REFERENCE["profiles"][name]:
        result = predictor.predict(record["sequence"])
        assert [s.score for s in result.sites] == pytest.approx(record["scores"], abs=1e-12)
        assert CleavageResult.from_dict(result.to_dict()) == result
        assert result.model.scored_endpoint == "site_cleavage"
        assert result.model.version.startswith("matrix-sha256:")
        assert dict(result.conditions)["profile_minutes"] == str(predictor.minutes)


def test_cathepsin_h_only_scores_initial_n_terminal_bond():
    result = ITCellCleavage("H").predict(CleavageInput("MALWMRLLPLL", source_start=30))
    assert [s.bond for s in result.sites] == [1]
    assert result.to_dict()["sites"][0]["source_bond"] == 31
    assert dict(result.conditions)["topology"] == "n_terminal"
    assert "norleucine surrogate" in result.sites[0].reason
    assert "cascade" in result.model.limitations


@pytest.mark.parametrize("enzyme", ["B", "S"])
def test_internal_profiles_preserve_edge_and_unrepresented_residue_limits(enzyme):
    result = ITCellCleavage(enzyme).predict("ACMMMM")
    assert [s.bond for s in result.sites] == list(range(1, 6))
    assert "absent flanks contribute zero" in result.sites[0].reason
    assert "C unrepresented" in result.sites[0].reason
    assert "norleucine surrogate" in result.sites[0].reason
    assert dict(result.conditions)["pH"] == "6.5"
    assert result.model.score_units.startswith("sum of log2")


@pytest.mark.parametrize("input_", [
    CleavageInput("AACG", n_term="unknown"), CleavageInput("AACG", n_term="acetylated"),
    CleavageInput("AACG", c_term="amidated"), "A",
])
def test_unestablished_or_modified_chemistry_and_missing_bonds_abstain(input_):
    for settings in ITCELL_MODELS.values():
        result = ITCellCleavage(**settings).predict(input_)
        assert result.unsupported_reason
        assert result.sites == ()


@pytest.mark.parametrize("enzyme,minutes", [("L", 240), ("S", True), ("S", 240.0), ("B", 30)])
def test_no_silent_enzyme_or_time_profile_substitution(enzyme, minutes):
    with pytest.raises(ValueError):
        ITCellCleavage(enzyme, minutes)


def test_all_profiles_discoverable_with_separate_artifact_identities():
    models = {m.name: m for m in cleavage_models() if m.name in ITCELL_MODELS}
    assert set(models) == set(ITCELL_MODELS)
    assert len({m.version for m in models.values()}) == 9
    with pytest.raises(ValueError, match="enzyme state"):
        get_cleavage_model("itcell-cats-240", enzyme_state="active")
