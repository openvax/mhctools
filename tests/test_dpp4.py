"""Reference values from the published DPP4 qPISA Dataset EV2."""

import hashlib
from importlib.resources import files
from itertools import product
import json

import pytest

from mhctools import CleavageInput, DPP4qPISA


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
    CleavageInput("HAE", n_term="acetylated"),
    CleavageInput("HAE", n_term="unknown"),
    CleavageInput("HAE", c_term="amidated"),
    CleavageInput("HAE", c_term="unknown"),
    CleavageInput("H"), CleavageInput("HA"),
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
