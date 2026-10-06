"""Real optional source tables versus deposited native pep_deseq2 scores."""

import json
from pathlib import Path

import pytest

from mhctools import CleavageInput, CleavageResult, get_cleavage_model
from mhctools import phagescout_artifacts as assets
from mhctools.phagescout import PHAGESCOUT_OPTIONAL_MODELS


pytestmark = pytest.mark.requires_external_tool
REFERENCE = json.loads((Path(__file__).parent / "data/phagescout_source_reference.json").read_text())


@pytest.mark.parametrize("name,settings", PHAGESCOUT_OPTIONAL_MODELS.items())
def test_full_tables_reproduce_native_scores_and_missing_masks(name, settings):
    path = assets.profile_directory()
    if not (path / assets.asset(settings["enzyme"])[0]).exists():
        pytest.skip("Run mhctools fetch phagescout")
    predictor = get_cleavage_model(name)
    source = "elastase" if settings["enzyme"] == "ELANE" else "cathepsin G"
    for record in REFERENCE["proteins"]:
        result = predictor.predict(CleavageInput(
            record["sequence"], source_start=record["source_start"], source_id=record["accession"]))
        rows = record["expected"][source]
        expected = {row["pos"]: row["pep_deseq2"] for row in rows if row["pep_deseq2"] is not None}
        actual = {site["source_bond"]: site["score"] for site in result.to_dict()["sites"]}
        assert actual == pytest.approx(expected, abs=1e-10)
        gaps = json.loads(dict(result.conditions)["unassessed_bonds"])
        assert {int(bond) + record["source_start"] for bond in gaps} == {
            row["pos"] for row in rows if row["pep_deseq2"] is None}
        assert CleavageResult.from_dict(result.to_dict()) == result
    # Negative native enrichment remains a score, never a probability.
    expected_first = -1.15831462085427 if settings["enzyme"] == "ELANE" else -1.52650628148886
    assert [site.score for site in predictor.predict("AAACD").sites] == pytest.approx([expected_first] * 4)
