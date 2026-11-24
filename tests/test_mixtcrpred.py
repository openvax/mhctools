# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the MixTCRpred catalog, downloads, and wrapper."""

from mhctools.artifacts import _SNAPSHOTS
from functools import lru_cache
from io import BytesIO
import hashlib
import json
from pathlib import Path
from urllib.error import HTTPError

import pandas as pd
import pytest

from mhctools import (
    AlleleIdentityError,
    Kind,
    MHC_DEPENDENCE_VALUES,
    MixTCRpred,
    MixTCRpredModel,
    MixTCRpredNonSpecificModel,
    TCR,
)
from mhctools import mixtcrpred
from mhctools.mixtcrpred import (
    classify_catalog_mhc,
    mixtcrpred_runtime_available,
    non_specific_model_catalog,
    partition_model_catalog,
    resolve_non_specific_mhc_target,
)
from mhctools.cli.script import main


MODEL = "A0201_GILGFVFTL"


def _make_home(tmp_path, checkpoint=True):
    home = tmp_path / "MixTCRpred"
    for relative in _SNAPSHOTS["mixtcrpred"].required_paths:
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
    }, {
        # The pinned upstream catalog has one such placeholder but publishes
        # no corresponding checkpoint. It is not an inference model.
        "MixTCRpred_model_name": "A0201_SLLMWITQC",
        "Peptide": "SLLMWITQC",
        "Origin": "Homo sapiens",
        "MHC_class": "MHCI",
        "MHC": "HLA-A*02:01",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 0,
        "AUC_5fold": 0,
    }])
    catalog.to_csv(
        home / "pretrained_models" / "info_models.csv", index=False)
    model_path = home / "pretrained_models" / ("model_%s.ckpt" % MODEL)
    if checkpoint:
        model_path.write_bytes(b"checkpoint")
    else:
        model_path.unlink()
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
    assert all(model.name != "A0201_SLLMWITQC" for model in models)


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
    home = _make_home(tmp_path)
    model = "DPB10401_TFEYVSQPFLMDLE"
    content = b"downloaded checkpoint"
    md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    record = {
        "id": int(mixtcrpred.ZENODO_RECORD),
        "files": [{
            "key": "model_%s.ckpt" % model,
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
    selected = mixtcrpred.fetch_models(home, models=[model])
    assert selected[0].status == "ready"
    assert Path(selected[0].path).read_bytes() == content
    manifest = json.loads((
        home / "pretrained_models" /
        mixtcrpred._MODEL_MANIFEST).read_text())
    assert manifest["record"] == mixtcrpred.ZENODO_RECORD
    assert manifest["models"][model]["md5"] == md5
    assert len(manifest["models"][model]["sha256"]) == 64

    def unavailable(url, timeout):
        raise OSError("fixture offline")

    monkeypatch.setattr(mixtcrpred, "urlopen", unavailable)
    monkeypatch.setattr(mixtcrpred, "ZENODO_RETRY_DELAYS", (0, 0))
    # A previously verified checkpoint can be used with no network access.
    assert mixtcrpred.fetch_models(home, models=[model])[0].status == "ready"
    with pytest.raises(RuntimeError, match="fixture offline"):
        mixtcrpred.fetch_models(home, models=[model], refresh_metadata=True)

    path = Path(selected[0].path)
    # Same size, different contents: size alone must never verify a cache.
    path.write_bytes(b"x" * len(content))
    with pytest.raises(RuntimeError, match="does not match its recorded"):
        mixtcrpred.fetch_models(home, models=[model])

    path.unlink()
    with pytest.raises(RuntimeError, match="fixture offline"):
        mixtcrpred.fetch_models(home, models=[model])


def test_unrecorded_cached_model_requires_remote_verification(monkeypatch, tmp_path):
    home = _make_home(tmp_path)

    def unavailable(url, timeout):
        raise OSError("fixture offline")

    monkeypatch.setattr(mixtcrpred, "urlopen", unavailable)
    monkeypatch.setattr(mixtcrpred, "ZENODO_RETRY_DELAYS", (0, 0))
    with pytest.raises(RuntimeError, match="fixture offline"):
        mixtcrpred.fetch_models(home, models=[MODEL])


def test_zenodo_metadata_retries_transient_504(monkeypatch, capsys):
    record = {"id": int(mixtcrpred.ZENODO_RECORD), "files": []}
    responses = [
        HTTPError(mixtcrpred.ZENODO_API_URL, 504, "Gateway Timeout", {}, None),
        HTTPError(mixtcrpred.ZENODO_API_URL, 503, "Unavailable", {}, None),
        BytesIO(json.dumps(record).encode()),
    ]
    calls = []
    sleeps = []

    def fake_urlopen(url, timeout):
        calls.append((url, timeout))
        response = responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response

    monkeypatch.setattr(mixtcrpred, "urlopen", fake_urlopen)
    monkeypatch.setattr(mixtcrpred, "ZENODO_RETRY_DELAYS", (2, 5))
    monkeypatch.setattr(mixtcrpred.time, "sleep", sleeps.append)

    assert mixtcrpred._zenodo_files() == {}
    assert len(calls) == 3
    assert sleeps == [2, 5]
    diagnostic = capsys.readouterr().err
    assert "attempt 1/3" in diagnostic
    assert "attempt 2/3" in diagnostic


def test_zenodo_download_retry_truncates_partial_file(monkeypatch, tmp_path):
    destination = tmp_path / "model.ckpt"
    responses = [OSError("connection reset"), BytesIO(b"complete")]
    sleeps = []

    def fake_urlopen(url, timeout):
        response = responses.pop(0)
        if isinstance(response, Exception):
            destination.write_bytes(b"partial")
            raise response
        return response

    monkeypatch.setattr(mixtcrpred, "urlopen", fake_urlopen)
    monkeypatch.setattr(mixtcrpred, "ZENODO_RETRY_DELAYS", (0,))
    monkeypatch.setattr(mixtcrpred.time, "sleep", sleeps.append)

    mixtcrpred._download_zenodo_file("https://zenodo.example/model", destination)

    assert destination.read_bytes() == b"complete"
    assert sleeps == [0]


def test_zenodo_does_not_retry_nontransient_http_error(monkeypatch):
    calls = []

    def missing():
        calls.append(True)
        raise HTTPError("https://zenodo.example/missing", 404, "Missing", {}, None)

    monkeypatch.setattr(mixtcrpred.time, "sleep", lambda delay: None)
    with pytest.raises(RuntimeError, match="after 1 attempt.*404"):
        mixtcrpred._retry_zenodo("download missing model", missing)
    assert len(calls) == 1


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


def test_prediction_cli_output_write_failure_has_no_traceback(
        monkeypatch, tmp_path, capsys):
    home = _make_home(tmp_path)
    input_path = tmp_path / "tcrs.csv"
    output_path = tmp_path / "missing" / "scored.csv"
    pd.DataFrame([{
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

    with pytest.raises(SystemExit) as raised:
        main([
            "mixtcrpred", "--model", MODEL,
            "--mixtcrpred-path", str(home),
            "--input", str(input_path), "--out", str(output_path),
        ])

    assert raised.value.code == 2
    stderr = capsys.readouterr().err
    assert str(output_path.parent) in stderr
    assert "Traceback" not in stderr
    assert not output_path.exists()


# Probed lazily and cached: it starts an interpreter and imports torch, which
# the pure-unit tests in this module must not pay for.
@lru_cache(maxsize=None)
def _runtime_ready():
    return mixtcrpred_runtime_available()


def _integration_predictor():
    try:
        catalog = MixTCRpred.catalog()
    except FileNotFoundError:
        pytest.skip("MixTCRpred checkout is not available")
    model = next((entry for entry in catalog if entry.name == MODEL), None)
    if model is None or model.status != "ready":
        pytest.skip("MixTCRpred GILGFVFTL checkpoint is not available")
    # "ready" means the checkpoint files were located, not that inference can
    # run: upstream's src/ imports torchvision and pytorch_lightning, which the
    # interpreter running the tests need not have.
    if not _runtime_ready():
        pytest.skip(
            "MixTCRpred runtime is not available (set MIXTCRPRED_PYTHON to an "
            "interpreter with torch, torchvision and pytorch-lightning; "
            "see docs/testing.md)")
    return MixTCRpred(MODEL, batch_size=3)


@pytest.mark.requires_external_tool
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



def test_catalog_classifies_ambiguous_targets_without_blocking_valid_model(tmp_path, caplog):
    home = _make_home(tmp_path)
    path = home / 'pretrained_models' / 'info_models.csv'
    catalog = pd.read_csv(path)
    extra = catalog.iloc[0].copy()
    extra['MixTCRpred_model_name'] = 'H-2q_example'
    extra['MHC'] = 'H-2q'
    pd.concat([catalog, extra.to_frame().T]).to_csv(path, index=False)
    with caplog.at_level('WARNING', logger='mhctools.mixtcrpred'):
        models = mixtcrpred.model_catalog(home)
    assert {model.name for model in models} == {MODEL, 'DPB10401_TFEYVSQPFLMDLE'}
    assert 'H-2q_example' in caplog.text
    assert 'Haplotype' in caplog.text
    assert MixTCRpred.resolve_model('GILGFVFTL', 'HLA-A*02:01', home).name == MODEL


NON_SPECIFIC_UPSTREAM_ROWS = (
    {
        "MixTCRpred_model_name": "H-2q_GPEGAQGPRGEPGTP",
        "Peptide": "GPEGAQGPRGEPGTP",
        "Origin": "Mus musculus",
        "MHC_class": "MHCII",
        "MHC": "H-2q",
        "Host_species": "MusMusculus",
        "Number_training_abTCR": 60,
        "AUC_5fold": 0.85,
        "expected_target": "H2-q",
        "expected_kind": "haplotype",
        "query_variants": ("H-2q", "H2-q"),
    },
    {
        "MixTCRpred_model_name": "B08_QIKVRVDMV",
        "Peptide": "QIKVRVDMV",
        "Origin": "CMV",
        "MHC_class": "MHCI",
        "MHC": "HLA-B*08",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 42,
        "AUC_5fold": 0.81,
        "expected_target": "HLA-B*08",
        "expected_kind": "allele_group",
        "query_variants": ("HLA-B*08",),
    },
    {
        "MixTCRpred_model_name": "DQA01_LLEFYLAMPFATP",
        "Peptide": "LLEFYLAMPFATP",
        "Origin": "Homo sapiens",
        "MHC_class": "MHCII",
        "MHC": "HLA-DQA*01",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 55,
        "AUC_5fold": 0.78,
        "expected_target": "HLA-DQA1*01",
        "expected_kind": "allele_group",
        "query_variants": ("HLA-DQA*01", "HLA-DQA1*01"),
    },
    {
        "MixTCRpred_model_name": "DQ2_PQPELPYPQPE",
        "Peptide": "PQPELPYPQPE",
        "Origin": "Triticum aestivum",
        "MHC_class": "MHCII",
        "MHC": "HLA-DQ2",
        "Host_species": "HomoSapiens",
        "Number_training_abTCR": 120,
        "AUC_5fold": 0.91,
        "expected_target": "HLA-DQ2",
        "expected_kind": "serotype",
        "query_variants": ("HLA-DQ2",),
    },
)


def make_home_with_non_specific_targets(tmp_path, checkpoints=()):
    """Create a MixTCRpred fixture containing both valid and the four #649 models."""
    home = _make_home(tmp_path)
    csv_path = home / "pretrained_models" / "info_models.csv"
    catalog = pd.read_csv(csv_path)
    extra_rows = [
        {
            key: value
            for key, value in row.items()
            if key not in ("expected_target", "expected_kind", "query_variants")
        }
        for row in NON_SPECIFIC_UPSTREAM_ROWS
    ]
    pd.concat([catalog, pd.DataFrame(extra_rows)], ignore_index=True).to_csv(
        csv_path, index=False)
    for model_name in checkpoints:
        (home / "pretrained_models" / ("model_%s.ckpt" % model_name)).write_bytes(
            b"checkpoint")
    return home


@pytest.mark.parametrize("row", NON_SPECIFIC_UPSTREAM_ROWS)
def test_classify_catalog_mhc_preserves_non_specific_targets_without_relabeling(row):
    target = classify_catalog_mhc(row["MHC"])
    assert target.is_specific_allele is False
    assert target.upstream_mhc == row["MHC"]
    assert target.canonical_target == row["expected_target"]
    assert target.target_kind == row["expected_kind"]
    assert target.reason
    assert target.mhcgnomes_version
    json.dumps(target.to_dict())


def test_classify_catalog_mhc_keeps_specific_alleles_and_pairs():
    allele = classify_catalog_mhc("HLA-A*02:01")
    assert allele.is_specific_allele is True
    assert allele.canonical_target == "HLA-A*02:01"
    assert allele.target_kind == "allele"
    assert allele.reason == ""

    pair = classify_catalog_mhc("HLA-DQA1*01:02/DQB1*06:02")
    assert pair.is_specific_allele is True
    assert pair.canonical_target == "HLA-DQA1*01:02/DQB1*06:02"
    assert pair.target_kind == "pair"


def test_partitioned_catalog_separates_all_four_non_specific_targets(tmp_path, caplog):
    home = make_home_with_non_specific_targets(tmp_path)
    partitioned = partition_model_catalog(home)
    assert [m.name for m in partitioned.models] == [MODEL, "DPB10401_TFEYVSQPFLMDLE"]
    assert all(isinstance(m, MixTCRpredModel) for m in partitioned.models)

    non_specific = partitioned.non_specific_models
    assert len(non_specific) == 4
    assert non_specific == non_specific_model_catalog(home)
    assert non_specific == MixTCRpred.non_specific_catalog(home)
    assert all(isinstance(m, MixTCRpredNonSpecificModel) for m in non_specific)
    for entry, expected in zip(non_specific, NON_SPECIFIC_UPSTREAM_ROWS):
        assert entry.name == expected["MixTCRpred_model_name"]
        assert entry.peptide == expected["Peptide"]
        assert entry.mhc_target == expected["expected_target"]
        assert entry.mhc_target_kind == expected["expected_kind"]
        assert entry.upstream_mhc == expected["MHC"]
        assert entry.upstream_revision == mixtcrpred.UPSTREAM_REVISION
        assert entry.mhcgnomes_version
        assert entry.reason
        assert not hasattr(entry, "allele")

    # Serialization works cleanly for both partitions.
    payload = json.loads(json.dumps(partitioned.to_dict()))
    assert len(payload["models"]) == 2
    assert len(payload["non_specific_models"]) == 4

    # Single-allele catalog still excludes all 4 with warnings while retaining
    # valid checkpoints.
    with caplog.at_level("WARNING", logger="mhctools.mixtcrpred"):
        single_models = MixTCRpred.catalog(home)
    assert {m.name for m in single_models} == {MODEL, "DPB10401_TFEYVSQPFLMDLE"}
    for expected in NON_SPECIFIC_UPSTREAM_ROWS:
        assert expected["MixTCRpred_model_name"] in caplog.text


@pytest.mark.parametrize("row", NON_SPECIFIC_UPSTREAM_ROWS)
def test_resolve_model_rejects_ambiguous_targets_while_non_specific_resolver_works(
        tmp_path, row):
    home = make_home_with_non_specific_targets(tmp_path)
    with pytest.raises(AlleleIdentityError):
        MixTCRpred.resolve_model(row["Peptide"], row["MHC"], home)

    for variant in row["query_variants"]:
        resolved = MixTCRpred.resolve_non_specific_model(
            row["Peptide"], variant, home)
        assert isinstance(resolved, MixTCRpredNonSpecificModel)
        assert resolved.name == row["MixTCRpred_model_name"]
        assert resolved.mhc_target == row["expected_target"]
        assert resolved.mhc_target_kind == row["expected_kind"]

    # Valid single-allele targets remain unaffected.
    assert MixTCRpred.resolve_model("GILGFVFTL", "HLA-A*02:01", home).name == MODEL


def test_resolve_non_specific_mhc_target_rejects_specific_allele():
    with pytest.raises(ValueError, match="resolve_model"):
        resolve_non_specific_mhc_target("HLA-A*02:01")


@pytest.mark.parametrize("row", NON_SPECIFIC_UPSTREAM_ROWS)
def test_constructor_rejects_non_specific_models_without_opt_in(tmp_path, row):
    model_name = row["MixTCRpred_model_name"]
    home = make_home_with_non_specific_targets(tmp_path, checkpoints=[model_name])
    with pytest.raises(AlleleIdentityError, match=row["expected_kind"]):
        MixTCRpred(model_name, mixtcrpred_path=home)


def test_non_specific_model_explicit_opt_in_contract(monkeypatch, tmp_path):
    model_name = "DQ2_PQPELPYPQPE"
    home = make_home_with_non_specific_targets(tmp_path, checkpoints=[model_name])
    predictor = MixTCRpred(
        model_name,
        mixtcrpred_path=home,
        allow_non_specific_target=True,
    )
    assert predictor.is_non_specific_target is True
    assert isinstance(predictor.model_info, MixTCRpredNonSpecificModel)
    assert "mhc_target=HLA-DQ2" in repr(predictor)
    assert "mhc_target_kind=serotype" in repr(predictor)

    support = predictor.kind_support()[Kind.pMHC_TCR_binding]
    assert support["mhc_dependence"] == "mhc_group"
    assert support["mhc_dependence"] in MHC_DEPENDENCE_VALUES
    assert support["mhc_dependence"] != "single_allele"
    assert support["mhc_class"] == "II"
    assert support["mhc_target"] == "HLA-DQ2"
    assert support["mhc_target_kind"] == "serotype"

    monkeypatch.setattr(
        predictor, "_run_sidecar", lambda tcrs: _sidecar_output(len(tcrs)))
    prediction = predictor.predict([_tcr()])[0].preds[0]
    assert prediction.peptide == "PQPELPYPQPE"
    assert prediction.allele == ""
    assert prediction.measurement_context is not None
    assert "serotype=HLA-DQ2" in prediction.measurement_context.detail
    assert "upstream=HLA-DQ2" in prediction.measurement_context.detail

    input_df = pd.DataFrame([{
        "sample": "s1",
        "cdr3_TRA": "CAGASGNTGKLIF",
        "cdr3_TRB": "CASSIRASYEQYF",
        "TRAV": "TRAV27",
        "TRAJ": "TRAJ42",
        "TRBV": "TRBV19",
        "TRBJ": "TRBJ2-6",
    }])
    annotated = predictor.annotate_dataframe(input_df)
    assert annotated["mixtcrpred_model"].tolist() == [model_name]
    assert annotated["mixtcrpred_allele"].tolist() == [""]
    assert annotated["mixtcrpred_mhc_target"].tolist() == ["HLA-DQ2"]
    assert annotated["mixtcrpred_mhc_target_kind"].tolist() == ["serotype"]
    assert annotated["mixtcrpred_upstream_mhc"].tolist() == ["HLA-DQ2"]


def test_fetch_models_enforces_opt_in_for_non_specific_targets(monkeypatch, tmp_path):
    home = make_home_with_non_specific_targets(tmp_path)
    target_model = "B08_QIKVRVDMV"
    content = b"b08 checkpoint"
    md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    record = {
        "id": int(mixtcrpred.ZENODO_RECORD),
        "files": [{
            "key": "model_%s.ckpt" % target_model,
            "size": len(content),
            "checksum": "md5:%s" % md5,
            "links": {"self": "https://zenodo.example/b08"},
        }],
    }

    def fake_urlopen(url, timeout):
        if url == mixtcrpred.ZENODO_API_URL:
            return BytesIO(json.dumps(record).encode())
        return BytesIO(content)

    monkeypatch.setattr(mixtcrpred, "urlopen", fake_urlopen)

    with pytest.raises(AlleleIdentityError, match="allele_group"):
        mixtcrpred.fetch_models(home, models=[target_model])

    # Single-allele bulk selection never pulls non-specific models.
    selected_high = mixtcrpred.select_catalog_models(
        partition_model_catalog(home), high_confidence=True)
    assert [m.name for m in selected_high] == [MODEL]

    fetched = mixtcrpred.fetch_models(
        home, models=[target_model], include_non_specific=True)
    assert len(fetched) == 1
    assert isinstance(fetched[0], MixTCRpredNonSpecificModel)
    assert fetched[0].status == "ready"
    assert Path(fetched[0].path).read_bytes() == content


def test_cli_non_specific_catalog_listing_and_opt_in_scoring(
        monkeypatch, tmp_path, capsys):
    model_name = "DQ2_PQPELPYPQPE"
    home = make_home_with_non_specific_targets(tmp_path, checkpoints=[model_name])
    monkeypatch.setenv("MIXTCRPRED_HOME", str(home))

    main(["ls", "mixtcrpred", "--models", "--non-specific", "--json"])
    listed = json.loads(capsys.readouterr().out)
    assert [entry["name"] for entry in listed] == [
        row["MixTCRpred_model_name"] for row in NON_SPECIFIC_UPSTREAM_ROWS
    ]
    assert all("allele" not in entry for entry in listed)
    assert listed[3]["mhc_target"] == "HLA-DQ2"
    assert listed[3]["mhc_target_kind"] == "serotype"

    main(["ls", "mixtcrpred", "--models", "--non-specific"])
    table = capsys.readouterr().out
    assert "MHC_TARGET" in table
    assert "TARGET_KIND" in table
    assert "H2-q" in table
    assert "HLA-DQ2" in table

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

    with pytest.raises(SystemExit) as raised:
        main([
            "mixtcrpred", "--model", model_name,
            "--mixtcrpred-path", str(home),
            "--input", str(input_path), "--out", str(output_path),
        ])
    assert raised.value.code == 2
    assert "non-specific" in capsys.readouterr().err

    main([
        "mixtcrpred", "--model", model_name,
        "--allow-non-specific",
        "--mixtcrpred-path", str(home),
        "--input", str(input_path), "--out", str(output_path),
    ])
    output = pd.read_csv(output_path, keep_default_na=False)
    assert output["mixtcrpred_allele"].tolist() == [""]
    assert output["mixtcrpred_mhc_target"].tolist() == ["HLA-DQ2"]
    assert output["mixtcrpred_mhc_target_kind"].tolist() == ["serotype"]
    assert output["mixtcrpred_upstream_mhc"].tolist() == ["HLA-DQ2"]
