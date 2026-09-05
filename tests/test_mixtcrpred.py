# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the MixTCRpred catalog, downloads, and wrapper."""

from io import BytesIO
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from mhctools import Kind, MixTCRpred, TCR
from mhctools import mixtcrpred
from mhctools.cli.script import main


MODEL = "A0201_GILGFVFTL"


def _make_home(tmp_path, checkpoint=True):
    home = tmp_path / "MixTCRpred"
    for relative in mixtcrpred._REQUIRED_HOME_PATHS:
        path = home / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    catalog = pd.DataFrame([{
        "MixTCRpred_model_name": MODEL,
        "Peptide": "GILGFVFTL",
        "Origin": "Influenza A virus",
        "MHC_class": "MHCI",
        "MHC": "HLA-A*02:01",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 2079,
        "AUC_5fold": 0.9432831141,
    }, {
        "MixTCRpred_model_name": "DPB10401_TFEYVSQPFLMDLE",
        "Peptide": "TFEYVSQPFLMDLE",
        "Origin": "SARS-CoV-2",
        "MHC_class": "MHCII",
        "MHC": "HLA-DPB1*04:01",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 20,
        "AUC_5fold": 0.7,
    }])
    catalog.to_csv(
        home / "pretrained_models" / "info_models.csv", index=False)
    model_path = home / "pretrained_models" / ("model_%s.ckpt" % MODEL)
    if checkpoint:
        model_path.write_bytes(b"checkpoint")
    return home


def _tcr(name="clone1"):
    return TCR(
        cdr3a="CAGASGNTGKLIF",
        cdr3b="CASSIRASYEQYF",
        name=name,
        trav="TRAV27",
        traj="TRAJ42",
        trbv="TRBV19",
        trbj="TRBJ2-6",
    )


def _sidecar_output(count=1):
    return pd.DataFrame({
        "__mhctools_id": range(count),
        "score": [3.076019] * count,
        "percentile_rank": [0.001] * count,
        "trav_corrected": ["TRAV27"] * count,
        "traj_corrected": ["TRAJ42"] * count,
        "trbv_corrected": ["TRBV19"] * count,
        "trbj_corrected": ["TRBJ2-6"] * count,
        "cdr1a_derived": ["SVFSS"] * count,
        "cdr2a_derived": ["VVTGGEV"] * count,
        "cdr1b_derived": ["LNHDA"] * count,
        "cdr2b_derived": ["SQIVND"] * count,
        "warning": ["-"] * count,
    })


def test_model_catalog_includes_target_metadata_and_weight_status(tmp_path):
    home = _make_home(tmp_path)
    models = mixtcrpred.model_catalog(home)
    assert len(models) == 2
    assert models[0].name == MODEL
    assert models[0].allele == "HLA-A*02:01"
    assert models[0].mhc_class == "I"
    assert models[0].training_tcrs == 2079
    assert models[0].high_confidence is True
    assert models[0].status == "ready"
    assert models[0].manager == "user"
    assert models[1].mhc_class == "II"
    assert models[1].high_confidence is False
    assert models[1].status == "missing"


def test_resolve_model_uses_normalized_allele(tmp_path):
    home = _make_home(tmp_path)
    model = MixTCRpred.resolve_model("GILGFVFTL", "HLA-A0201", home)
    assert model.name == MODEL


def test_constructor_requires_selected_checkpoint(tmp_path):
    home = _make_home(tmp_path, checkpoint=False)
    with pytest.raises(FileNotFoundError, match="mhctools fetch mixtcrpred"):
        MixTCRpred(MODEL, mixtcrpred_path=home)


def test_predict_maps_score_rank_target_and_version(monkeypatch, tmp_path):
    home = _make_home(tmp_path)
    predictor = MixTCRpred(MODEL, mixtcrpred_path=home)
    monkeypatch.setattr(
        predictor, "_run_sidecar", lambda tcrs: _sidecar_output(len(tcrs)))
    result = predictor.predict([_tcr()])[0].preds[0]
    assert result.kind == Kind.pMHC_TCR_binding
    assert result.peptide == "GILGFVFTL"
    assert result.allele == "HLA-A*02:01"
    assert result.tcr == "clone1"
    assert result.score == pytest.approx(3.076019)
    assert result.percentile_rank == pytest.approx(0.001)
    assert MODEL in result.predictor_version
    assert predictor.kind_support()[Kind.pMHC_TCR_binding] == {
        "mhc_dependence": "single_allele",
        "mhc_class": "I",
    }


def test_annotate_dataframe_captures_scores_metadata_and_qc(
        monkeypatch, tmp_path):
    home = _make_home(tmp_path)
    predictor = MixTCRpred(MODEL, mixtcrpred_path=home)
    monkeypatch.setattr(
        predictor, "_run_sidecar", lambda tcrs: _sidecar_output(len(tcrs)))
    input_df = pd.DataFrame([{
        "sample": "s1",
        "cdr3_TRA": "CAGASGNTGKLIF",
        "cdr3_TRB": "CASSIRASYEQYF",
        "TRAV": "TRAV27",
        "TRAJ": "TRAJ42",
        "TRBV": "TRBV19",
        "TRBJ": "TRBJ2-6",
    }])
    output = predictor.annotate_dataframe(input_df)
    assert output["sample"].tolist() == ["s1"]
    assert output["mixtcrpred_model"].tolist() == [MODEL]
    assert output["mixtcrpred_score"].tolist() == pytest.approx([3.076019])
    assert output["mixtcrpred_percentile_rank"].tolist() == [0.001]
    assert output["mixtcrpred_trav_corrected"].tolist() == ["TRAV27"]
    assert output["mixtcrpred_cdr1a_derived"].tolist() == ["SVFSS"]
    assert output["mixtcrpred_warning"].tolist() == ["-"]


@pytest.mark.parametrize("field,value,error", [
    ("cdr3a", "", "requires cdr3a"),
    ("cdr3b", "C" * 21, "limited to 20"),
    ("cdr3a", "CAV?", "unsupported residues"),
])
def test_input_validation(field, value, error):
    values = _tcr().to_dict()
    values[field] = value
    with pytest.raises(ValueError, match=error):
        mixtcrpred._validate_tcr(TCR.from_dict(values))


def test_fetch_models_downloads_atomically_and_records_hashes(
        monkeypatch, tmp_path):
    home = _make_home(tmp_path, checkpoint=False)
    content = b"downloaded checkpoint"
    md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    record = {
        "id": int(mixtcrpred.ZENODO_RECORD),
        "files": [{
            "key": "model_%s.ckpt" % MODEL,
            "size": len(content),
            "checksum": "md5:%s" % md5,
            "links": {"self": "https://zenodo.example/model"},
        }],
    }

    def fake_urlopen(url, timeout):
        if url == mixtcrpred.ZENODO_API_URL:
            return BytesIO(json.dumps(record).encode())
        assert url == "https://zenodo.example/model"
        return BytesIO(content)

    monkeypatch.setattr(mixtcrpred, "urlopen", fake_urlopen)
    selected = mixtcrpred.fetch_models(home, models=[MODEL])
    assert selected[0].status == "ready"
    assert Path(selected[0].path).read_bytes() == content
    manifest = json.loads((
        home / "pretrained_models" /
        mixtcrpred._MODEL_MANIFEST).read_text())
    assert manifest["record"] == mixtcrpred.ZENODO_RECORD
    assert manifest["models"][MODEL]["md5"] == md5
    assert len(manifest["models"][MODEL]["sha256"]) == 64


def test_model_catalog_json_is_serializable(tmp_path):
    home = _make_home(tmp_path)
    json.dumps([model.to_dict() for model in MixTCRpred.catalog(home)])


def test_ls_models_cli_reports_downloaded_weights(
        monkeypatch, tmp_path, capsys):
    home = _make_home(tmp_path)
    monkeypatch.setenv("MIXTCRPRED_HOME", str(home))
    main(["ls", "mixtcrpred", "--models", "--downloaded", "--json"])
    output = json.loads(capsys.readouterr().out)
    assert [row["name"] for row in output] == [MODEL]
    assert output[0]["manager"] == "user"
    assert output[0]["status"] == "ready"


def test_prediction_cli_preserves_input_and_appends_all_outputs(
        monkeypatch, tmp_path):
    home = _make_home(tmp_path)
    input_path = tmp_path / "tcrs.csv"
    output_path = tmp_path / "scored.csv"
    pd.DataFrame([{
        "sample": "s1",
        "cdr3_TRA": "CAGASGNTGKLIF",
        "cdr3_TRB": "CASSIRASYEQYF",
        "TRAV": "TRAV27",
        "TRAJ": "TRAJ42",
        "TRBV": "TRBV19",
        "TRBJ": "TRBJ2-6",
    }]).to_csv(input_path, index=False)
    monkeypatch.setattr(
        MixTCRpred,
        "_run_sidecar",
        lambda self, tcrs: _sidecar_output(len(tcrs)),
    )

    main([
        "mixtcrpred", "--model", MODEL,
        "--mixtcrpred-path", str(home),
        "--input", str(input_path), "--out", str(output_path),
    ])

    output = pd.read_csv(output_path, keep_default_na=False)
    assert output["sample"].tolist() == ["s1"]
    assert output["mixtcrpred_score"].tolist() == pytest.approx([3.076019])
    assert output["mixtcrpred_percentile_rank"].tolist() == [0.001]
    assert output["mixtcrpred_trav_corrected"].tolist() == ["TRAV27"]
    assert output["mixtcrpred_warning"].tolist() == ["-"]


def _integration_predictor():
    try:
        catalog = MixTCRpred.catalog()
    except FileNotFoundError:
        pytest.skip("MixTCRpred checkout is not available")
    model = next((entry for entry in catalog if entry.name == MODEL), None)
    if model is None or model.status != "ready":
        pytest.skip("MixTCRpred GILGFVFTL checkpoint is not available")
    return MixTCRpred(MODEL, batch_size=3)


def test_real_checkpoint_regression_current_pytorch():
    predictor = _integration_predictor()
    tcrs = [
        _tcr("one"),
        TCR(
            cdr3a="CALSGETSGSRLTF",
            cdr3b="CASGLVPGGLVYEQYF",
            trav="TRAV12-2", traj="TRAJ45",
            trbv="TRBV2", trbj="TRBJ2-5", name="two"),
        TCR(
            cdr3a="CAVR", cdr3b="CASS",
            trav="TRAV999", traj="TRAJ999",
            trbv="TRBV999", trbj="TRBJ999", name="bad-v"),
    ]
    output = predictor.predict_dataframe(tcrs)
    # Upstream's reference output is published to three decimal places. Keep
    # that scientifically meaningful precision while tolerating the observed
    # ~4e-5 CPU-kernel difference between macOS and Linux.
    assert output["score"].tolist() == pytest.approx([
        3.076, -1.651, -1.871], abs=5e-4)
    assert output["percentile_rank"].tolist() == pytest.approx([
        0.001, 22.142, 28.797], abs=5e-4)
    assert predictor.last_qc["warning"].tolist() == [
        "-", "-", "Error TRAV-TRBV"]
