"""Reference values from the published DPP4 qPISA Dataset EV2."""

import hashlib
from importlib.resources import files
from itertools import product
import json
import math

import pytest

from mhctools import PeptideInput, DPP4qPISA


@pytest.mark.parametrize("triplet,expected", [
    # Independent sums of the three published EV2 cells, not fitted fixtures.
    ("HAE", 2.1694),  # 2.6329 - 0.2231 - 0.2404
    ("DAE", 0.9836),  # 2.6329 - 1.4089 - 0.2404
    ("HAP", -0.2018),  # 2.6329 - 0.2231 - 2.6116; no clipping
    ("AAA", 3.1206),
    ("HSE", 0.6744),  # Ser is supported too; do not impose an XP/XA gate.
    ("HPE", 2.7760),
])
def test_published_native_scores(triplet, expected):
    result = DPP4qPISA().predict(triplet + "GTFTSDYSK")
    assert result.unsupported_reason is None
    assert len(result.sites) == 1
    assert result.sites[0].bond == 2
    assert result.sites[0].status == "scored"
    assert result.sites[0].score == pytest.approx(expected, abs=1e-12)


def test_only_exposed_terminus_is_assessed():
    result = DPP4qPISA().predict("HAEHAEHAE")
    assert [site.bond for site in result.sites] == [2]


@pytest.mark.parametrize("peptide", [
    PeptideInput("HAE", n_term="acetylated"),
    PeptideInput("HAE", n_term="unknown"),
    PeptideInput("HAE", c_term="amidated"),
    PeptideInput("HAE", c_term="unknown"),
    PeptideInput("H"), PeptideInput("HA"),
])
def test_unsupported_inputs_abstain(peptide):
    result = DPP4qPISA().predict(peptide)
    assert result.unsupported_reason
    assert not result.sites


def test_missing_coefficients_abstain_without_imputation():
    result = DPP4qPISA().predict("KCA")
    assert "P2:K" in result.unsupported_reason
    assert not result.sites


def test_published_table_identity_and_complete_triplet_coverage():
    raw = files("mhctools").joinpath("data/dpp4_qpisa.json").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "aee694d2efbd4de14dddb72059ebe468a11bf73c97be1d3cc447ae59ca484f78"
    data = json.loads(raw)
    assert data["source_sha256"] == "ee449da13b5ec66fd6fb08c16203da8c44f2e9c10572a355abdcb985694c6ddc"
    assert data["license"] == "CC0-1.0"
    counts = {"scored": 0, "unsupported": 0}
    for triplet in product("ACDEFGHIKLMNPQRSTVWY", repeat=3):
        result = DPP4qPISA().predict("".join(triplet))
        counts["unsupported" if result.unsupported_reason else "scored"] += 1
    assert counts == {"scored": 6420, "unsupported": 1580}


def test_keane_assay_uses_published_parent_half_lives_and_held_out_triplets():
    data = json.loads(files("mhctools").joinpath("data/dpp4_keane2011.json").read_text())
    rows = data["rows"]
    assert len(rows) == 11
    assert len({row["triplet"] for row in rows}) == 8
    assert next(row for row in rows if row["name"] == "GLP-1-amide") == {
        "name": "GLP-1-amide", "triplet": "HAE", "length": 30,
        "half_life_minutes": 8.63, "sd_minutes": 0.92, "n": 3,
    }
    assert next(row for row in rows if row["name"] == "VIP")["half_life_minutes"] == 173.3
    estimate = DPP4qPISA().estimate_keane_fraction("HAE" + "G" * 27, minutes=10)
    assert estimate["status"] == "assay_estimate"
    assert estimate["bond"] == 2
    assert estimate["estimated_half_life_minutes"] == pytest.approx(17.14987, rel=1e-5)
    assert estimate["parent_lost_fraction"] == pytest.approx(0.3325, abs=0.001)
    assert estimate["parent_remaining_fraction"] + estimate["parent_lost_fraction"] == pytest.approx(1)
    assert estimate["validation"]["leave_one_triplet_out_typical_factor_error"] == pytest.approx(2.81, abs=0.01)


def test_keane_assay_abstains_outside_measured_domain():
    predictor = DPP4qPISA()
    for peptide in (
            "HAE", "HAE" + "G" * 23, "HAE" + "G" * 40,
            "RPK" + "G" * 27, "KPV" + "G" * 27,
            PeptideInput("HAE" + "G" * 27, n_term="acetylated"),
            PeptideInput("HAE" + "G" * 27, c_term="amidated")):
        result = predictor.estimate_keane_fraction(peptide, minutes=10)
        assert result["status"] == "unassessed"
        assert result["reason"]
        assert "parent_lost_fraction" not in result
    for minutes in (-1, math.nan, math.inf, True):
        with pytest.raises(ValueError, match="minutes"):
            predictor.estimate_keane_fraction("HAE" + "G" * 27, minutes=minutes)


def test_keane_assay_fraction_is_bounded_and_changes_with_time():
    predictor = DPP4qPISA()
    peptide = "HSD" + "G" * 27
    at_zero = predictor.estimate_keane_fraction(peptide, minutes=0)
    at_ten = predictor.estimate_keane_fraction(peptide, minutes=10)
    at_sixty = predictor.estimate_keane_fraction(peptide, minutes=60)
    assert at_zero["parent_lost_fraction"] == 0
    assert 0 < at_ten["parent_lost_fraction"] < at_sixty["parent_lost_fraction"] < 1
