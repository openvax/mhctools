"""Native CPP scores, isolated estimator guards and exact-input semantics."""

from io import StringIO
import json
import math
import os
from pathlib import Path
import sys

import pandas as pd
import pytest

from mhctools import Kind, PeptideContext, PeptideInput, PeptiVerseCPP, Prediction
from mhctools.cli.annotate_table import main as table_main
from mhctools.cli.args import mhc_predictors
from mhctools.optional_backend import run_python_sidecar
from mhctools.peptiverse_cpp import (
    CPP_THRESHOLD, _CPP_DIRECTORY, _CPP_MANIFEST, parse_cpp_results,
)


REFERENCE = json.loads((Path(__file__).parent / "data" /
                        "peptiverse_cpp_source_reference.json").read_text())
CONFIGURED = {key: os.environ.get(variable) for key, variable in (
    ("peptiverse_home", "PEPTIVERSE_HOME"),
    ("peptiverse_python", "PEPTIVERSE_PYTHON"),
    ("peptiverse_esm_home", "PEPTIVERSE_ESM_HOME"))}


@pytest.fixture(autouse=True)
def isolated_environment(monkeypatch):
    for variable in ("PEPTIVERSE_HOME", "PEPTIVERSE_PYTHON", "PEPTIVERSE_ESM_HOME"):
        monkeypatch.delenv(variable, raising=False)


def _row(sequence="SIINFEKL", index=0, score=0.9):
    return dict(__mhctools_id=index, peptide=sequence, score=score,
                label=int(score >= CPP_THRESHOLD), threshold=CPP_THRESHOLD,
                emb_tag="wt", uncertainty="", uncertainty_type="")


def _parse(rows, sequences):
    return parse_cpp_results(StringIO(pd.DataFrame(rows).to_csv(index=False)), sequences)


def test_parser_restores_occurrences_and_native_threshold_boundary():
    rows = [_row(index=2, score=1), _row(index=0, score=0),
            _row(index=1, score=CPP_THRESHOLD)]
    frame = _parse(rows, ["SIINFEKL"] * 3)
    assert frame.score.tolist() == [0, CPP_THRESHOLD, 1]
    assert frame.label.tolist() == [0, 1, 1]
    assert frame.predictive_entropy_nats.tolist() == [None] * 3


@pytest.mark.parametrize("field,value,error", [
    ("__mhctools_id", 1, "every input ID"),
    ("__mhctools_id", 0.5, "every input ID"),
    ("peptide", "SIINFEK", "different peptide"),
    ("emb_tag", "chemberta", "different embedding"),
    ("score", -0.1, "scores, labels"),
    ("score", 1.1, "scores, labels"),
    ("score", float("nan"), "scores, labels"),
    ("score", float("inf"), "scores, labels"),
    ("label", 0, "scores, labels"),
    ("label", 1.5, "scores, labels"),
    ("threshold", 0.5, "native threshold"),
    ("uncertainty_type", "ensemble_std", "type without a value"),
    ("uncertainty", 0.5, "predictive entropy"),
])
def test_parser_rejects_identity_and_endpoint_corruption(field, value, error):
    row = _row()
    row[field] = value
    with pytest.raises(RuntimeError, match=error):
        _parse([row], [row["peptide"] if field != "peptide" else "SIINFEKL"])


def test_parser_missing_rows_duplicates_and_columns():
    with pytest.raises(RuntimeError, match="every input ID"):
        _parse([_row()], ["SIINFEKL"] * 2)
    with pytest.raises(RuntimeError, match="every input ID"):
        _parse([_row(), _row()], ["SIINFEKL"] * 2)
    row = _row()
    del row["score"]
    with pytest.raises(ValueError, match="missing columns"):
        _parse([row], ["SIINFEKL"])


@pytest.mark.parametrize("score", [0, 1e-10, 0.5, 1 - 1e-10, 1])
def test_native_entropy_uses_author_clipping_and_natural_logs(score):
    row = _row(score=score)
    p = min(max(score, 1e-9), 1 - 1e-9)
    row.update(uncertainty=-(p * math.log(p) + (1-p) * math.log(1-p)),
               uncertainty_type="binary_predictive_entropy_single_model")
    assert _parse([row], ["SIINFEKL"]).iloc[0].predictive_entropy_nats == pytest.approx(
        row["uncertainty"], rel=1e-14)
    row["uncertainty"] += 1e-6
    with pytest.raises(RuntimeError, match="predictive entropy"):
        _parse([row], ["SIINFEKL"])


def _stub_home(tmp_path, sklearn_version="1.7.2", retained_limit=1020,
               threshold=CPP_THRESHOLD, classes=(0, 1), variant="svm_gpu"):
    """Actual standalone bridge, with an upstream-shaped dependency-free source."""
    home = tmp_path / "CPP source"
    home.mkdir()
    model = home / _CPP_DIRECTORY / "best_model.joblib"
    model.parent.mkdir(parents=True)
    model.write_text("synthetic estimator")
    (home / "best_models.txt").write_text(_CPP_MANIFEST)
    esm = home / "esm2_t33_650M_UR50D"
    esm.mkdir()
    for name in ("config.json", "tokenizer_config.json", "special_tokens_map.json",
                 "vocab.txt", "model.safetensors"):
        (esm / name).write_text("synthetic " + name)
    sklearn = home / "sklearn"
    sklearn.mkdir()
    (sklearn / "__init__.py").write_text("__version__ = %r\n" % sklearn_version)
    (sklearn / "svm.py").write_text(f'''
class Classes:
    def tolist(self): return {list(classes)!r}
class SVC:
    classes_ = Classes()
    probability = True
''')
    (home / "inference.py").write_text(f'''
from pathlib import Path
from sklearn.svm import SVC
class Mask:
    def __init__(self, n): self.n = n
    def sum(self): return self
    def item(self): return self.n
class Embedder:
    def _tokenize(self, peptides):
        return {{"input_ids": peptides[0], "attention_mask": None}}
    def _valid_mask(self, ids, mask): return Mask(min(len(ids), {retained_limit}))
class PeptiVersePredictor:
    def __init__(self, manifest_path, classifier_weight_root, esm_name, device):
        assert "Halflife" not in Path(manifest_path).read_text()
        assert Path(esm_name).is_absolute()
        self.wt_embedder = Embedder()
        self.meta = {{("permeability_penetrance", "wt"): {{
            "artifact": str(Path(classifier_weight_root) / {str(_CPP_DIRECTORY)!r} / "best_model.joblib"),
            "model_name": {variant!r}, "emb_tag": "wt", "kind": "joblib",
            "threshold": {threshold}, "task_type": "Classifier"}}}}
        self.models = {{("permeability_penetrance", "wt"): SVC()}}
    def predict_property(self, prop, col, peptide, uncertainty):
        assert prop == "permeability_penetrance" and col == "wt"
        return {{"score": 0.9, "label": 1, "threshold": {threshold}, "emb_tag": "wt"}}
''')
    return home


def _predictor(home, **kwargs):
    return PeptiVerseCPP(peptiverse_home=home, peptiverse_python=sys.executable,
                        allow_unverified_assets=True, **kwargs)


def test_exact_input_context_identity_serialization_and_unsupported_records(tmp_path):
    predictor = _predictor(_stub_home(tmp_path))
    inputs = [PeptideInput("SIINFEKL", occurrence_id="first"),
              PeptideInput("SIINFEKL", c_term="amidated", occurrence_id="amide"),
              PeptideInput("SIINFEKL", occurrence_id="second"),
              PeptideInput("SIINFEKL", context=PeptideContext(matrix="serum")),
              PeptideInput("A" * 80), PeptideInput("A" * 1021)]
    with pytest.raises(ValueError, match="Input 1 is unsupported"):
        predictor.predict(inputs)
    preds = [r.preds[0] for r in predictor.predict(inputs, on_unsupported="record")]
    assert [p.peptide_input for p in preds] == inputs
    assert [p.score for p in preds] == [0.9, None, 0.9, 0.9, 0.9, None]
    assert all(p.value is None and p.kind == Kind.cpp_classification for p in preds)
    assert all(p.measurement_context.matrix is None for p in preds)
    assert preds[0].measurement_context.class_label == "CPP"
    assert preds[1].measurement_context.class_label is None
    assert preds[0].cache_key == preds[2].cache_key
    assert len({preds[i].cache_key for i in (0, 1, 3)}) == 3
    for pred in preds:
        assert Prediction.from_dict(json.loads(json.dumps(pred.to_dict()))) == pred
    assert predictor.last_qc.outside_source_training_lengths.tolist() == [False] * 4 + [True, True]
    assert predictor.last_qc.input_index.tolist() == list(range(6))
    assert predictor.artifact_inventory.inference_status == "reproduced"
    first_key = preds[0].cache_key
    assert predictor.predict(inputs[:1])[0].preds[0].cache_key == first_key
    assert predictor.predict([]) == [] and predictor.last_qc.empty
    assert predictor.supported_kinds == (Kind.cpp_classification,)
    assert predictor.kind_support()[Kind.cpp_classification]["mhc_dependence"] == "none"


@pytest.mark.parametrize("maximum", [0, 1021, 3.5, True])
def test_invalid_capacity_is_rejected(maximum):
    with pytest.raises(ValueError, match="integer in"):
        PeptiVerseCPP(max_peptide_length=maximum)


@pytest.mark.parametrize("length", [1, 2, 3, 61, 62, 1020])
def test_source_domain_length_is_a_flag_not_a_biological_cutoff(tmp_path, length):
    predictor = _predictor(_stub_home(tmp_path))
    pred = predictor.predict("A" * length)[0].preds[0]
    assert pred.score == 0.9
    assert bool(predictor.last_qc.iloc[0].outside_source_training_lengths) == (not 3 <= length <= 61)


@pytest.mark.parametrize("kwargs,error", [
    ({"sklearn_version": "1.9.1"}, "requires scikit-learn==1.7.2"),
    ({"retained_limit": 7}, "retained 7 of 8 residues"),
    ({"threshold": 0.5}, "estimator/classes/threshold"),
    ({"classes": (1, 0)}, "estimator/classes/threshold"),
    ({"variant": "svm_cpu"}, "Unexpected PeptiVerse model metadata"),
])
def test_bridge_rejects_runtime_or_model_drift(tmp_path, kwargs, error):
    home = _stub_home(tmp_path, **kwargs)
    if "sklearn_version" in kwargs:
        # The version guard must run before source import or deserialization.
        (home / "inference.py").write_text("raise AssertionError('loaded too early')\n")
    with pytest.raises(RuntimeError, match=error):
        _predictor(home).predict("SIINFEKL")


def test_assets_are_verified_and_cpp_does_not_require_half_life_weights(tmp_path):
    home = _stub_home(tmp_path)
    with pytest.raises(RuntimeError, match="mismatch artifacts"):
        PeptiVerseCPP(peptiverse_home=home)
    predictor = _predictor(home)
    before = predictor.predictor_version
    (home / _CPP_DIRECTORY / "best_model.joblib").write_text("replacement")
    assert _predictor(home).predictor_version != before
    (home / _CPP_DIRECTORY / "best_model.joblib").unlink()
    with pytest.raises(FileNotFoundError, match="missing required artifacts"):
        _predictor(home)


def test_cli_preserves_rows_and_exposes_class_score(tmp_path, monkeypatch):
    predictor = _predictor(_stub_home(tmp_path))
    assert mhc_predictors["peptiverse-cpp"] is PeptiVerseCPP
    monkeypatch.setitem(mhc_predictors, "peptiverse-cpp", lambda **kwargs: predictor)
    source, output = tmp_path / "input.csv", tmp_path / "output.csv"
    source.write_text("peptide,label\nSIINFEKL,a\nRRRRRRRR,b\nSIINFEKL,c\n")
    table_main(["--input", str(source), "--out", str(output),
                "--predictor", "peptiverse-cpp:cpp_probability:score"])
    frame = pd.read_csv(output)
    assert frame.label.tolist() == ["a", "b", "c"]
    assert frame.cpp_probability.tolist() == [0.9] * 3


@pytest.mark.requires_external_tool
@pytest.mark.skipif(not CONFIGURED["peptiverse_home"], reason="set PEPTIVERSE_HOME for native CPP tests")
def test_native_source_and_fixed_reference_conformance(tmp_path):
    host_ml = {name: sys.modules.get(name) for name in ("torch", "sklearn")}
    host_offline = {name: os.environ.get(name) for name in (
        "HF_HUB_OFFLINE", "TRANSFORMERS_OFFLINE", "PYTHONNOUSERSITE", "WANDB_MODE")}
    # Independent original API call, without the mhctools endpoint bridge.
    audit = tmp_path / "native.py"
    output = tmp_path / "native.json"
    sequences = [r["peptide"] for r in REFERENCE["records"]] + ["SIINFEKL"]
    audit.write_text('''
import json,sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import inference
class Unused:
    def __init__(self,*args,**kwargs): pass
inference.SMILESEmbedder = Unused
inference.ChemBERTaEmbedder = Unused
predictor = inference.PeptiVersePredictor(
    manifest_path=sys.argv[2], classifier_weight_root=sys.argv[1],
    esm_name=sys.argv[3], device="cpu")
rows = [predictor.predict_property("permeability_penetrance", "wt", sequence,
    uncertainty=True) for sequence in json.loads(sys.argv[4])]
Path(sys.argv[5]).write_text(json.dumps(rows))
''')
    manifest = tmp_path / "native-manifest.txt"
    manifest.write_text(_CPP_MANIFEST)
    run_python_sidecar(backend_name="native CPP oracle", python=CONFIGURED["peptiverse_python"],
                       sidecar=audit, cwd=CONFIGURED["peptiverse_home"], timeout=3600,
                       args=[CONFIGURED["peptiverse_home"], str(manifest),
                             CONFIGURED["peptiverse_esm_home"], json.dumps(sequences), str(output)])
    native = json.loads(output.read_text())
    predictor = PeptiVerseCPP(device="cpu", uncertainty=True, **CONFIGURED)
    preds = [r.preds[0] for r in predictor.predict(sequences)]
    assert predictor.artifact_inventory.status == "verified"
    differences = []
    for i, (pred, direct) in enumerate(zip(preds, native)):
        fixed = REFERENCE["records"][i % len(REFERENCE["records"])]
        # Same-runtime conformance stays tight. Saved CPU controls also
        # include float32 ESM2 platform/library rounding (see #531).
        assert pred.score == pytest.approx(direct["score"], rel=0, abs=1e-10)
        differences.append((pred.peptide, abs(pred.score - fixed["score"])))
        assert direct["label"] == fixed["label"]
        assert pred.value is None and pred.percentile_rank is None
        assert pred.measurement_context.class_label == ("CPP" if direct["label"] else "non-CPP")
        qc = predictor.last_qc.iloc[i]
        assert qc.threshold == CPP_THRESHOLD
        assert qc.predictive_entropy_nats == pytest.approx(direct["uncertainty"], rel=0, abs=1e-10)
    assert all(error <= 1e-5 for _, error in differences), differences
    assert preds[0].cache_key == preds[-1].cache_key
    assert {name: sys.modules.get(name) for name in host_ml} == host_ml
    assert {name: os.environ.get(name) for name in host_offline} == host_offline
