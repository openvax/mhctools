"""Optional source-data acquisition, integrity, lazy discovery and inference."""

import hashlib
import io
import json
from pathlib import Path
import shutil

import pytest

from mhctools import CleavageInput, PhageScout, cleavage_models, predict_cleavage_batch
from mhctools.artifacts import artifact_status, fetch, managed_path
from mhctools.cli.artifacts import fetch_main, ls_main
from mhctools.cli.cleavage import main
from mhctools import phagescout_artifacts as assets
from mhctools.phagescout import PHAGESCOUT_OPTIONAL_MODELS, _full_weights, _peptide_weights


@pytest.fixture
def small_profiles(monkeypatch, tmp_path):
    """Small pinned test assets exercise the actual streaming/install path."""
    monkeypatch.delenv("PHAGESCOUT_HOME", raising=False)
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(tmp_path / "data"))
    payloads = {}
    files = []
    for enzyme, name, _, _ in assets.FILES:
        raw = b"peptide\tlog2FoldChange\nAAACD\t-2\nAACDE\t4\n"
        payloads[name] = raw
        files.append((enzyme, name, len(raw), hashlib.sha256(raw).hexdigest()))
    monkeypatch.setattr(assets, "FILES", tuple(files))
    requests = []

    def download(url, timeout):
        assert timeout == 30
        requests.append(url)
        name = next(name for name in payloads if assets.asset(
            next(row[0] for row in files if row[1] == name))[3] == url)
        return io.BytesIO(payloads[name])

    monkeypatch.setattr(assets, "urlopen", download)
    _full_weights.cache_clear()
    yield payloads, requests
    _full_weights.cache_clear()


def test_inventory_and_catalog_are_read_only_and_lazy(small_profiles, monkeypatch):
    def unexpected(*args, **kwargs):
        pytest.fail("Optional catalog loaded external profiles")
    monkeypatch.setattr(assets, "selected_profile", unexpected)
    base = {model.name for model in cleavage_models()}
    optional = {model.name for model in cleavage_models(include_optional=True)}
    assert not base & PHAGESCOUT_OPTIONAL_MODELS.keys()
    assert PHAGESCOUT_OPTIONAL_MODELS.keys() <= optional
    status = artifact_status("phagescout")
    assert status.status == "missing" and status.fetchable
    assert status.version == assets.RECORD
    assert not Path(status.path).parent.exists()
    assert small_profiles[1] == []


def test_fetch_installs_exact_files_license_and_provenance_and_is_noop(small_profiles):
    payloads, requests = small_profiles
    result = fetch("phagescout", version=assets.RECORD)
    root = Path(result.path)
    assert result.status == "ready" and result.manager == "mhctools"
    assert root == managed_path("phagescout")
    assert {name: (root / name).read_bytes() for name in payloads} == payloads
    manifest = json.loads((root / ".mhctools-artifact.json").read_text())
    assert manifest["license"] == "CC-BY-4.0"
    assert manifest["source"] == "https://zenodo.org/records/21387981"
    assert "PhageScout" in (root / "PHAGESCOUT_LICENSE.txt").read_text()
    assert len(requests) == 2
    assert "%20" in requests[1]
    assert fetch("phagescout") == result
    assert len(requests) == 2


@pytest.mark.parametrize("failure", ["checksum", "truncated", "oversized", "network", "stream"])
def test_second_download_failure_never_publishes_partial_install(small_profiles, failure, monkeypatch):
    payloads, _ = small_profiles
    original = assets.urlopen
    target = assets.managed_directory()

    def download(url, timeout):
        if url == assets.asset("ELANE")[3]:
            return original(url, timeout)
        if failure == "network":
            raise OSError("interrupted transfer")
        data = payloads[assets.asset("CTSG")[0]]
        if failure == "stream":
            class Interrupted(io.BytesIO):
                def read(self, size):
                    if self.tell():
                        raise OSError("interrupted after partial write")
                    return super().read(10)
            return Interrupted(data)
        data = {"checksum": data.replace(b"-2", b"-3"), "truncated": data[:-1],
                "oversized": data + b"extra"}[failure]
        return io.BytesIO(data)

    monkeypatch.setattr(assets, "urlopen", download)
    with pytest.raises((ValueError, OSError)):
        fetch("phagescout")
    assert not target.exists()
    assert list(target.parent.iterdir()) == []


def test_concurrent_complete_fetch_is_verified_and_reused(small_profiles, monkeypatch):
    target = assets.managed_directory()

    def competing_fetch(staging, destination):
        assert destination == target
        shutil.copytree(staging, target)
        raise FileExistsError(17, "concurrent fetch won")

    monkeypatch.setattr(Path, "rename", competing_fetch)
    assert fetch("phagescout").status == "ready"
    assert list(target.parent.iterdir()) == [target]


def test_existing_invalid_install_is_preserved_and_not_overwritten(small_profiles):
    target = assets.managed_directory()
    target.mkdir(parents=True)
    sentinel = target / "user-file"
    sentinel.write_text("preserve")
    with pytest.raises(RuntimeError, match="Move it aside"):
        fetch("phagescout")
    assert sentinel.read_text() == "preserve"
    assert small_profiles[1] == []


def test_invalid_requested_version_never_downloads_even_if_ready(small_profiles):
    for ready in (False, True):
        if ready:
            fetch("phagescout")
        before = len(small_profiles[1])
        with pytest.raises(ValueError, match="only Zenodo record 21387981"):
            fetch("phagescout", version="latest")
        assert len(small_profiles[1]) == before


@pytest.mark.parametrize("corruption", ["profile", "checksum", "manifest", "license"])
def test_inventory_rejects_corruption_and_inference_rechecks_cached_file(small_profiles, corruption):
    root = Path(fetch("phagescout").path)
    predictor = PhageScout(profile="peptide-deseq2")
    assert predictor.predict("AAACD").sites
    name = {"profile": assets.asset("ELANE")[0], "checksum": assets.asset("ELANE")[0],
            "manifest": ".mhctools-artifact.json",
            "license": "PHAGESCOUT_LICENSE.txt"}[corruption]
    (root / name).write_bytes((root / name).read_bytes().replace(b"-2", b"-3")
                             if corruption == "checksum" else b"invalid")
    assert artifact_status("phagescout").status == "missing"
    with pytest.raises(RuntimeError, match="Move it aside"):
        fetch("phagescout")
    if corruption in ("profile", "checksum"):
        with pytest.raises(RuntimeError, match="invalid"):
            PhageScout(profile="peptide-deseq2")


def test_explicit_data_directory_overrides_environment(small_profiles, monkeypatch, tmp_path):
    monkeypatch.setenv("PHAGESCOUT_HOME", str(tmp_path / "invalid"))
    with pytest.raises(RuntimeError, match="PHAGESCOUT_HOME"):
        fetch("phagescout")
    result = fetch("phagescout", data_dir=tmp_path / "explicit")
    assert Path(result.path) == assets.managed_directory(tmp_path / "explicit")
    assert artifact_status("phagescout", data_dir=tmp_path / "explicit").status == "ready"
    predictor = PhageScout(profile="peptide-deseq2", profile_dir=result.path)
    assert predictor.predict("AAACD").sites[0].score == -2
    monkeypatch.setenv("PHAGESCOUT_HOME", result.path)
    assert fetch("phagescout").manager == "user"


def test_missing_asset_has_actionable_selection_error(small_profiles):
    with pytest.raises(RuntimeError, match="Run mhctools fetch phagescout"):
        PhageScout(profile="peptide-deseq2")
    with pytest.raises(ValueError, match="only to the optional"):
        PhageScout(profile_dir="unused")


def test_full_profile_native_mean_gaps_negative_scores_and_chemistry(small_profiles):
    fetch("phagescout")
    predictor = PhageScout(profile="peptide-deseq2")
    result = predictor.predict(CleavageInput("AAACDE", source_start=7))
    assert [site.score for site in result.sites] == [-2, 1, 1, 1, 4]
    assert [site["source_bond"] for site in result.to_dict()["sites"]] == [8, 9, 10, 11, 12]
    assert "bundled_data_sha256" not in dict(result.conditions)
    assert dict(result.conditions)["normalization"] == "none"
    assert predictor.predict("AAAAA").unsupported_reason
    assert predictor.predict(CleavageInput("AAACD", n_term="acetylated")).unsupported_reason
    assert predictor.predict("AAAA").unsupported_reason
    with pytest.raises(TypeError):
        predictor.weights[(0, 1, 2, 3, 4)]["AAACD"] = (0, 0)
    with pytest.raises(ValueError):
        predictor.predict("AAAXD")


def test_cli_download_inventory_native_scores_and_batch(small_profiles, capsys):
    fetch_main(["phagescout", "--json"])
    assert json.loads(capsys.readouterr().out)["status"] == "ready"
    ls_main(["phagescout", "--json"])
    assert json.loads(capsys.readouterr().out)[0]["version"] == assets.RECORD
    name = "phagescout-elane-peptide-deseq2"
    main(["--sequence", "AAACD", "--model", name])
    report = json.loads(capsys.readouterr().out)
    assert [site["score"] for site in report["results"][0]["sites"]] == [-2] * 4
    batch = predict_cleavage_batch(
        [dict(id="synthetic", scope="construct", sequence="AAACD", n_term="free", c_term="free",
              epitopes=[dict(id="core", start=1, end=4, sequence="AAC")])],
        [dict(id="inflammation", context="extracellular", models=[name])])
    internal = batch["assessments"][0]["overlays"][0]["internal"]
    assert [(site["bond"], site["score"]) for site in internal] == [(2, -2), (3, -2)]


@pytest.mark.parametrize("raw", ["wrong\theader\n", "peptide\tlog2FoldChange\nAAACD\tnan\n",
                                "peptide\tlog2FoldChange\nAAAXD\t2\n"])
def test_parser_rejects_unusable_source_rows(raw):
    with pytest.raises(ValueError):
        _peptide_weights(io.StringIO(raw), 5)


def test_parser_preserves_first_source_row():
    weights = _peptide_weights(io.StringIO(
        "peptide\tlog2FoldChange\nAAACD\t-2\nAAACD\t4\n"), 5)
    assert weights[(0, 1, 2, 3, 4)]["AAACD"] == (0, -2)
