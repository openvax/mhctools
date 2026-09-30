"""Live legacy-artifact tests; setup_test_backends.py pepsickle provisions these."""

import json
import os
import subprocess

import pytest

from mhctools import load_cleavage_batch, predict_cleavage_batch, write_cleavage_batch
from mhctools.pepsickle import Pepsickle


pytestmark = pytest.mark.skipif(not os.environ.get("PEPSICKLE_GB_PYTHON"),
                                reason="Configure PEPSICKLE_GB_PYTHON for legacy-artifact inference")
SEQUENCE = "LSRKVAELVHFLLLKYRAR"


def test_legacy_ci_profiles_match_direct_upstream_and_round_trip(tmp_path):
    # Conformance sequence only: this published MAGE-A3 substrate occurs in
    # upstream's training inputs and is deliberately not labeled held out.
    script = """
import json
from pepsickle.model_functions import initialize_digestion_gb_model, predict_protein_cleavage_locations
model = initialize_digestion_gb_model()
profiles = {kind: [float(row[2]) for row in predict_protein_cleavage_locations(
    'LSRKVAELVHFLLLKYRAR', model, mod_type='in-vitro', proteasome_type=kind)] for kind in ('C', 'I')}
print(json.dumps(profiles))
"""
    completed = subprocess.run([os.environ["PEPSICKLE_GB_PYTHON"], "-c", script],
                               text=True, capture_output=True, check=True, timeout=120)
    expected = json.loads(completed.stdout)
    models = ["pepsickle-in-vitro-all-mammal-" + kind
              for kind in ("constitutive", "immunoproteasome")]
    report = predict_cleavage_batch(
        [dict(id="native", peptide=SEQUENCE[3:16], n_flank=SEQUENCE[:3], c_flank=SEQUENCE[16:],
              source_start=40, evidence=[dict(kind="conformance", held_out=False)])],
        [dict(id="tumor", context="tumor", models=models)], raise_on_error=True)
    for row, kind in zip(report["assessments"], ("C", "I")):
        assert row["status"] == "assessed"
        result = row["result"]
        assert [s["score"] for s in result["sites"]] == pytest.approx(expected[kind][:-1])
        assert result["sites"][0]["source_bond"] == 41
        runtime = json.loads(dict(result["conditions"])["runtime"])
        assert runtime["python"].startswith("3.8.")
        assert runtime["packages"]["scikit-learn"] == "0.23.2"
        assert "runtime-sha256:" in result["model"]["version"]
    assert expected["C"] != expected["I"]
    path = tmp_path / "legacy.json"
    write_cleavage_batch(report, path)
    assert load_cleavage_batch(path) == report


def test_legacy_runtime_needs_no_host_pepsickle_install(monkeypatch):
    monkeypatch.setattr("mhctools.pepsickle.runtime_identity",
                        lambda *a: pytest.fail("Host Pepsickle must not be inspected"))
    result = Pepsickle(model_type="in-vitro", proteasome_type="C").predict_cleavage(SEQUENCE)
    assert len(result.sites) == len(SEQUENCE) - 1
    assert all(0 <= site.score <= 1 for site in result.sites)
