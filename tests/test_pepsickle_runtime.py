"""Runtime selection/provenance without requiring host Pepsickle assets."""

import json
import subprocess
import sys

import pytest

from mhctools import CleavageInput
from mhctools.pepsickle import Pepsickle


def foreign_identity():
    return dict(package_version="foreign-package", model_key="gradient_boosting",
                weights_sha256="a" * 64, inference_sha256="b" * 64,
                features_sha256="c" * 64, python="3.8.20",
                executable="/foreign/python", packages={"scikit-learn": "0.23.2"})


def test_external_runtime_provenance_and_changed_asset_rejection(monkeypatch):
    monkeypatch.setattr("mhctools.pepsickle._identity_cache", {})
    monkeypatch.setattr("mhctools.pepsickle.runtime_identity",
                        lambda *a: pytest.fail("Host assets must not identify a foreign runtime"))
    calls = []
    changed = False

    def run(args, input, **kwargs):
        request = json.loads(input)
        calls.append(request)
        identity = foreign_identity()
        if changed:
            identity["weights_sha256"] = "d" * 64
        response = {"identity": identity}
        if request.get("operation") != "identity":
            response["results"] = {seq: [0.2] * (len(seq) - 1) + [0.0]
                                   for seq in request["sequences"]}
        return subprocess.CompletedProcess(args, 0, json.dumps(response), "")

    monkeypatch.setattr("mhctools.pepsickle.subprocess.run", run)
    predictor = Pepsickle(model_type="in-vitro", proteasome_type="I", python_executable=sys.executable)
    assert predictor.isolate_subprocess
    result = predictor.predict_cleavage(CleavageInput("SIINFEKL"))
    assert "package:foreign-package" in result.model.version
    assert json.loads(dict(result.conditions)["runtime"]) == foreign_identity()
    assert len(result.sites) == 7
    assert len(calls) == 2  # inspect once, infer once
    changed = True
    with pytest.raises(RuntimeError, match="changed between inspection"):
        predictor.predict_cleavage("SIINFEKL")


def test_legacy_environment_selects_only_gradient_boosting(monkeypatch):
    monkeypatch.delenv("PEPSICKLE_PYTHON", raising=False)
    monkeypatch.setenv("PEPSICKLE_GB_PYTHON", sys.executable)
    assert Pepsickle(model_type="in-vitro", proteasome_type="C").python_executable == sys.executable
    assert Pepsickle().python_executable is None
    assert Pepsickle(model_type="in-vitro-2", proteasome_type="I").python_executable is None
    monkeypatch.setenv("PEPSICKLE_GB_PYTHON", "/nonexistent/pepsickle-python")
    with pytest.raises(RuntimeError, match="not runnable"):
        Pepsickle(model_type="in-vitro", proteasome_type="C")


def test_catalog_does_not_start_configured_external_runtime(monkeypatch):
    monkeypatch.setenv("PEPSICKLE_GB_PYTHON", sys.executable)
    monkeypatch.setattr("mhctools.pepsickle.subprocess.run",
                        lambda *a, **k: pytest.fail("Catalog must not start the legacy runtime"))
    model = Pepsickle.catalog_cleavage_model(model_type="in-vitro", proteasome_type="C")
    assert model.version.startswith("unresolved:")


def test_runtime_identity_inspection_error_is_explicit(monkeypatch):
    monkeypatch.setattr("mhctools.pepsickle._identity_cache", {})
    monkeypatch.setattr("mhctools.pepsickle.subprocess.run", lambda args, **kwargs:
                        subprocess.CompletedProcess(args, 1, "", "missing runtime assets"))
    predictor = Pepsickle(python_executable=sys.executable, human_only=True)
    with pytest.raises(RuntimeError, match="missing runtime assets"):
        predictor.cleavage_model()


def test_sequence_context_preserves_partial_terminal_knowledge_and_rejects_modifications(monkeypatch):
    from mhctools import predict_cleavage_batch
    from mhctools.pepsickle import PepsickleCleavage

    facade = PepsickleCleavage(model_type="in-vitro", proteasome_type="I")
    monkeypatch.setattr(facade.predictor, "_identity", foreign_identity)
    calls = []

    def score(sequences):
        calls.append(sequences)
        return {s: [0.25] * (len(s) - 1) + [0.0] for s in sequences}

    monkeypatch.setattr(facade.predictor, "cleavage_probs_many", score)
    inputs = [CleavageInput("SIINFEKL", n_term=n, c_term=c)
              for n, c in (("free", "unknown"), ("unknown", "free"),
                           ("unknown", "unknown"), ("acetylated", "unknown"),
                           ("unknown", "amidated"))]
    results = facade.predict_many(inputs)
    assert calls == [["SIINFEKL"]]
    assert [r.peptide for r in results] == inputs
    for result in results[:3]:
        assert len(result.sites) == 7
        assert "sequence_context_only" in dict(result.conditions)["input_scope"]
    assert all(r.unsupported_reason and not r.sites for r in results[3:])
    report = predict_cleavage_batch(
        [dict(id="parent", sequence="AASIINFEKL", scope="protein", fragments=[
            dict(id="conditional", start=2, end=10, n_term="free", c_term="unknown",
                 assumption="N-terminal cut occurs; retained C-terminal chemistry is unknown")])],
        [dict(id="apc", context="apc", models=[facade.model.name])],
        predictors={facade.model.name: facade}, raise_on_error=True)
    assert all(row["status"] == "assessed" for row in report["assessments"])
    assert report["assessments"][1]["result"]["peptide"]["c_term"] == "unknown"
