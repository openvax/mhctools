"""Source-data integrity and product/epitope coordinates for the Wada case study."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

from mhctools import load_cleavage_batch, normalize_cleavage_input


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("evaluate_wada_cleavage", ROOT / "scripts/evaluate_wada_cleavage.py")
evaluation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluation)


def test_complete_source_panels_preserve_products_times_and_parent_ends():
    request = json.loads(evaluation.FIXTURE.read_text())
    a, c = request["inputs"]
    assert [len(v["observed_products"]) for v in (a, c)] == [28, 19]
    for value in (a, c):
        normalized = normalize_cleavage_input(value)
        assert normalized["n_term"] == normalized["c_term"] == "unknown"
        assert normalized["observed_products"] == value["observed_products"]
        assert len({p["id"] for p in value["observed_products"]}) == len(value["observed_products"])
        assert all(p["first_detected_hours"] in (1, 2, 4) for p in value["observed_products"])
        bonds = evaluation.product_bonds(value)
        assert 0 not in bonds and len(value["sequence"]) not in bonds
        assert all(0 < b < len(value["sequence"]) for b in bonds)
    # Same SART3 epitope: detected only at 4h in A, already at 1h in C.
    assert a["observed_products"][6]["sequence"] == c["observed_products"][15]["sequence"] == "LLQAEAPRL"
    assert a["observed_products"][6]["first_detected_hours"] == 4
    assert c["observed_products"][15]["first_detected_hours"] == 1
    # C's N-terminal SART2 product requires only its internal C-terminal cut.
    assert c["observed_products"][0]["start"] == 0
    assert any(p["product_id"].endswith("fragment1") for p in evaluation.product_bonds(c)[9])


def test_product_alignment_corruption_and_unpinned_training_maps_are_rejected(tmp_path):
    value = deepcopy(json.loads(evaluation.FIXTURE.read_text())["inputs"][0])
    value["observed_products"][0]["start"] = 1
    with pytest.raises(ValueError, match="Invalid source product"):
        evaluation.product_bonds(value)
    (tmp_path / "data/raw/digestion_map_files").mkdir(parents=True)
    with pytest.raises(ValueError, match="pinned author snapshot"):
        evaluation.audit_training([value], tmp_path)


def test_context_alignment_uses_p1_and_explicit_padding():
    assert evaluation.context_window("ACDEFGHIK", 4) == "ACDEFGH"
    assert evaluation.context_window("ACDEFGHIK", 1) == "***ACDE"
    assert evaluation.context_window("ACDEFGHIK", 9) == "GHIK***"


def test_recorded_experimental_report_retains_partial_audit_and_all_observations():
    report = load_cleavage_batch(evaluation.FIXTURE.with_name("pepsickle-gb-evaluation.json"))
    assert sum(len(row["observed_bonds"]) for row in report["experimental_product_overlays"]) == 24
    assert sum(len(row["observed_products"]) for row in report["inputs"]) == 47
    audit = report["training_audit"]
    assert audit["inventory_sha256"] == evaluation.TRAINING_MAPS_SHA256
    assert audit["files"] == 79 and audit["source_sequences"] == 58
    assert not audit["wada_study_present"]
    assert audit["unresolved_study_files"]  # Do not silently declare complete lineage.
    assert all(not row["full_sequence_overlap"] and not row["observed_bond_context_overlaps"]
               for row in audit["inputs"])
    for row in report["experimental_product_overlays"]:
        value = next(v for v in report["inputs"] if v["id"] == row["input_id"])
        assert {entry["bond"]: entry["products"] for entry in row["observed_bonds"]} == evaluation.product_bonds(value)
