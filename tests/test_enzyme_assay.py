"""Keane 2011 purified-enzyme parent loss, not vaccine exposure."""

from importlib.resources import files
import json
import math

import pytest

from mhctools import AssayConditions, PeptideInput, estimate_assay_parent_fraction


def test_fap_source_values_and_a_vs_p_strata():
    data = json.loads(files("mhctools").joinpath("data/fap_keane2011.json").read_text())
    assert data["enzyme_concentration_um"] == 0.2
    assert len(data["rows"]) == 9
    assert len({r["triplet"] for r in data["rows"]}) == 8
    assert next(r for r in data["rows"] if r["name"] == "PYY") == {
        "name": "PYY", "triplet": "YPI", "length": 36,
        "half_life_minutes": 60.2, "sd_minutes": 16.9, "n": 4,
    }
    ala = estimate_assay_parent_fraction("HAE" + "G" * 27, enzyme="FAP", minutes=10)
    pro = estimate_assay_parent_fraction("YPS" + "G" * 33, enzyme="FAP", minutes=10)
    assert ala["status"] == pro["status"] == "assay_estimate"
    assert ala["estimated_half_life_minutes"] == pytest.approx(1254.7913, rel=1e-5)
    assert pro["estimated_half_life_minutes"] == pytest.approx(11.56536, rel=1e-5)
    assert ala["parent_lost_fraction"] == pytest.approx(0.005508, rel=1e-3)
    assert pro["parent_lost_fraction"] == pytest.approx(0.4508, rel=1e-3)
    assert ala["validation"]["leave_one_substrate_out_typical_factor_error"] == pytest.approx(1.396, rel=1e-3)
    assert pro["validation"]["leave_one_substrate_out_typical_factor_error"] == pytest.approx(3.003, rel=1e-3)
    assert ala["parent_lost_fraction"] + ala["parent_remaining_fraction"] == pytest.approx(1)
    assert pro["parent_lost_fraction"] + pro["parent_remaining_fraction"] == pytest.approx(1)


def test_fap_abstains_outside_source_domain():
    candidates = (
        "HAE" + "G" * 7, "HAE" + "G" * 40, "AAA" + "G" * 27,
        "HSD" + "G" * 27, PeptideInput("HAE" + "G" * 27, n_term="acetylated"),
        PeptideInput("HAE" + "G" * 27, c_term="unknown"),
    )
    for peptide in candidates:
        result = estimate_assay_parent_fraction(peptide, enzyme="FAP", minutes=10)
        assert result["status"] == "unassessed"
        assert result["reason"]
        assert "parent_lost_fraction" not in result
    for minutes in (-1, math.nan, math.inf, True):
        with pytest.raises(ValueError, match="minutes"):
            estimate_assay_parent_fraction("HAE" + "G" * 27, enzyme="FAP", minutes=minutes)
    with pytest.raises(ValueError, match="Supported enzymes"):
        estimate_assay_parent_fraction("HAE" + "G" * 27, enzyme="ACE", minutes=10)


def test_shared_entry_preserves_dpp4_assay_estimate():
    source = json.loads(files("mhctools").joinpath("data/dpp4_keane2011.json").read_text())
    assert source["enzyme_concentration_um"] == 0.112
    assert "0.112 micromolar" in source["assay"]
    peptide = "HAE" + "G" * 27
    result = estimate_assay_parent_fraction(peptide, enzyme="DPP4", minutes=10)
    assert result["status"] == "assay_estimate"
    assert result["model"] == "dpp4-qpisa-keane2011-parent-loss"
    assert result["parent_lost_fraction"] == pytest.approx(0.3325, abs=0.001)
    assert result["parent_remaining_fraction"] + result["parent_lost_fraction"] == pytest.approx(1)


def test_default_exposure_and_concentration_sensitivity_share_fraction_units():
    peptide = "HAE" + "G" * 27
    dpp4 = estimate_assay_parent_fraction(peptide, enzyme="DPP4")
    fap = estimate_assay_parent_fraction(peptide, enzyme="FAP")
    assert dpp4["exposure_minutes"] == fap["exposure_minutes"] == 60
    assert dpp4["reference_enzyme_concentration_um"] == 0.112
    assert fap["reference_enzyme_concentration_um"] == 0.2
    assert dpp4["parent_remaining_fraction"] == pytest.approx(
        2 ** (-60 / dpp4["estimated_half_life_minutes"]))
    assert fap["parent_remaining_fraction"] == pytest.approx(
        2 ** (-60 / fap["estimated_half_life_minutes"]))
    assert 0 <= dpp4["parent_remaining_fraction"] <= 1
    assert 0 <= fap["parent_remaining_fraction"] <= 1

    doubled = estimate_assay_parent_fraction(
        peptide, enzyme="DPP4",
        conditions=AssayConditions(exposure_minutes=60, enzyme_concentration_um=0.224))
    assert doubled["reference_half_life_minutes"] == pytest.approx(
        dpp4["estimated_half_life_minutes"])
    assert doubled["estimated_half_life_minutes"] == pytest.approx(
        dpp4["estimated_half_life_minutes"] / 2)
    assert doubled["parent_remaining_fraction"] == pytest.approx(
        dpp4["parent_remaining_fraction"] ** 2)


def test_assay_conditions_are_explicit_and_validated():
    for value in (-1, math.nan, math.inf, True):
        with pytest.raises(ValueError, match="exposure_minutes"):
            AssayConditions(exposure_minutes=value)
    for value in (0, -1, math.nan, math.inf, True):
        with pytest.raises(ValueError, match="enzyme_concentration_um"):
            AssayConditions(enzyme_concentration_um=value)
    with pytest.raises(ValueError, match="either minutes or conditions"):
        estimate_assay_parent_fraction(
            "HAE" + "G" * 27, enzyme="DPP4", minutes=10,
            conditions=AssayConditions())
    with pytest.raises(TypeError, match="AssayConditions"):
        estimate_assay_parent_fraction("HAE" + "G" * 27, enzyme="DPP4", conditions={})
