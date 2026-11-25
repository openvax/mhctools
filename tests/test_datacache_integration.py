# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Unit tests for Wave 12 datacache integration and CLI help delegation."""

import hashlib
import json
from pathlib import Path

import datacache
from datacache import OfflineError
import datacache.cli_help as dc_cli_help
import pytest

from mhctools import artifacts, mixtcrpred, phagescout_artifacts as ps_assets
from mhctools.cli import artifacts as cli_artifacts
from mhctools.cli import help as cli_help
from mhctools.cli import integrations as cli_integrations
from mhctools.cli.script import main


def _isolate_env(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg-cache"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg-data"))
    for key in (
        "MHCTOOLS_DATA_DIR",
        "OPENVAX_DATA_CACHE",
        "OPENVAX_OFFLINE",
        "DATACACHE_OFFLINE",
        "PHAGESCOUT_HOME",
        "MIXTCRPRED_HOME",
    ):
        monkeypatch.delenv(key, raising=False)
    return home


def test_data_path_honors_openvax_data_cache_and_precedence(monkeypatch, tmp_path):
    _isolate_env(monkeypatch, tmp_path)
    shared = tmp_path / "shared-openvax"
    pkg_env = tmp_path / "mhctools-env"
    explicit = tmp_path / "explicit-dir"

    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    assert artifacts.data_path() == (shared / "mhctools").resolve()

    # MHCTOOLS_DATA_DIR overrides OPENVAX_DATA_CACHE.
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(pkg_env))
    assert artifacts.data_path() == pkg_env.resolve()

    # Explicit data_dir overrides both environment variables.
    assert artifacts.data_path(explicit) == explicit.resolve()


@pytest.mark.parametrize("legacy_kind", ["user_data", "platform_cache", "dot_mhctools"])
def test_data_path_preserves_populated_legacy_directories(
        monkeypatch, tmp_path, legacy_kind):
    home = _isolate_env(monkeypatch, tmp_path)
    # Point platformdirs/datacache defaults inside tmp_path so real host caches
    # never influence this test.
    default_platform = tmp_path / "default-platform-openvax"
    legacy_user_data = tmp_path / "legacy-user-data"
    legacy_platform_cache = tmp_path / "legacy-platform-cache"
    legacy_dot = home / ".mhctools"

    monkeypatch.setattr(
        datacache.common,
        "platform_cache_dir",
        lambda name: str(
            default_platform if name == "openvax" else legacy_platform_cache
        ),
    )
    monkeypatch.setattr(
        artifacts,
        "user_data_path",
        lambda name, appauthor=False: legacy_user_data,
    )

    chosen = {
        "user_data": legacy_user_data,
        "platform_cache": legacy_platform_cache,
        "dot_mhctools": legacy_dot,
    }[legacy_kind]
    (chosen / "artifacts" / "eramer").mkdir(parents=True)
    (chosen / "artifacts" / "eramer" / "marker.txt").write_text("cached")

    assert artifacts.data_path() == chosen.resolve()

    # Setting OPENVAX_DATA_CACHE takes precedence over populated legacy dirs.
    shared = tmp_path / "shared-cache"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    assert artifacts.data_path() == (shared / "mhctools").resolve()


def test_cli_cache_root_alias_matches_data_dir():
    for parser_fn, argv_prefix in (
        (cli_artifacts.make_ls_parser, []),
        (cli_artifacts.make_fetch_parser, ["eramer"]),
        (cli_integrations.make_parser, []),
    ):
        parser = parser_fn()
        via_data_dir = parser.parse_args(argv_prefix + ["--data-dir", "/tmp/a"])
        via_cache_root = parser.parse_args(argv_prefix + ["--cache-root", "/tmp/a"])
        assert via_data_dir.data_dir == "/tmp/a"
        assert via_cache_root.data_dir == "/tmp/a"


@pytest.mark.parametrize("env_key", ["OPENVAX_OFFLINE", "DATACACHE_OFFLINE"])
def test_phagescout_offline_blocks_download_and_allows_cache_hit(
        monkeypatch, tmp_path, env_key):
    _isolate_env(monkeypatch, tmp_path)
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(tmp_path / "data"))

    payloads = {}
    files = []
    for enzyme, name, _, _ in ps_assets.FILES:
        raw = b"peptide\tlog2FoldChange\nAAACD\t-2\nAACDE\t4\n"
        payloads[name] = raw
        files.append((enzyme, name, len(raw), hashlib.sha256(raw).hexdigest()))
    monkeypatch.setattr(ps_assets, "FILES", tuple(files))

    fetch_calls = []

    def fake_fetch_file(url, *, destination, expected_sha256, expected_size,
                        max_bytes, force, raw, **kwargs):
        assert force is True
        assert raw is True
        assert max_bytes == expected_size
        fetch_calls.append(url)
        name = next(
            fname for fname in payloads
            if ps_assets.asset(next(r[0] for r in files if r[1] == fname))[3] == url
        )
        data = payloads[name]
        assert len(data) == expected_size
        assert hashlib.sha256(data).hexdigest() == expected_sha256
        Path(destination).write_bytes(data)
        return str(destination)

    monkeypatch.setattr(ps_assets, "fetch_file", fake_fetch_file)

    # Offline mode blocks uncached download before touching disk or fetch_file.
    monkeypatch.setenv(env_key, "1")
    with pytest.raises(OfflineError, match=env_key):
        ps_assets.fetch_profiles()
    assert fetch_calls == []
    assert not ps_assets.managed_directory().parent.exists()

    # Online fetch routes through datacache.fetch_file and installs atomically.
    monkeypatch.delenv(env_key)
    status = ps_assets.fetch_profiles()
    assert status.status == "ready"
    assert len(fetch_calls) == 2

    # Once cached, offline mode succeeds without touching the network.
    monkeypatch.setenv(env_key, "1")
    cached = ps_assets.fetch_profiles()
    assert cached.status == "ready"
    assert len(fetch_calls) == 2


def test_phagescout_cli_offline_exits_cleanly(monkeypatch, tmp_path, capsys):
    _isolate_env(monkeypatch, tmp_path)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(tmp_path / "shared"))
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")
    with pytest.raises(SystemExit) as caught:
        main(["fetch", "phagescout"])
    assert caught.value.code == 2
    assert "OPENVAX_OFFLINE" in capsys.readouterr().err


def test_mixtcrpred_delegates_to_datacache_and_enforces_offline(
        monkeypatch, tmp_path):
    _isolate_env(monkeypatch, tmp_path)
    home = tmp_path / "MixTCRpred"
    models_dir = home / "pretrained_models"
    models_dir.mkdir(parents=True)
    (home / "MixTCRpred.py").write_text("# stub\n")
    (home / "LICENSE.md").write_text("license\n")
    src = home / "src"
    src.mkdir()
    for fname in ("dataloaders.py", "models.py", "utils.py"):
        (src / fname).write_text("# stub\n")
    (models_dir / "anchors_perc_rank.pickle").write_bytes(b"anchors")
    (models_dir / "model_A0201_ELAGIGILTV.ckpt").write_bytes(b"shipped")
    (models_dir / "model_A0201_GILGFVFTL.ckpt").write_bytes(b"shipped")
    model = "DPB10401_TFEYVSQPFLMDLE"
    (models_dir / "info_models.csv").write_text(
        "MixTCRpred_model_name,Peptide,MHC,MHC_class,Host_species,Origin,"
        "Number_training_abTCR,AUC_5fold\n"
        "A0201_GILGFVFTL,GILGFVFTL,HLA-A*02:01,MHCI,Human,Influenza,120,0.95\n"
        "%s,TFEYVSQPFLMDLE,HLA-DPA1*01:03;HLA-DPB1*04:01,MHCII,Human,SARS-CoV-2,80,0.90\n"
        % model
    )

    content = b"datacache downloaded checkpoint"
    md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    record = {
        "id": int(mixtcrpred.ZENODO_RECORD),
        "files": [{
            "key": "model_%s.ckpt" % model,
            "size": len(content),
            "checksum": "md5:%s" % md5,
            "links": {"self": "https://zenodo.example/model.ckpt"},
        }],
    }

    bytes_calls = []
    file_calls = []

    def fake_fetch_bytes(url, **kwargs):
        bytes_calls.append(url)
        return json.dumps(record).encode("utf-8")

    def fake_fetch_file(url, *, destination, expected_size, max_bytes, **kwargs):
        assert expected_size == len(content)
        assert max_bytes == len(content)
        file_calls.append(url)
        Path(destination).write_bytes(content)
        return str(destination)

    monkeypatch.setattr(mixtcrpred, "fetch_bytes", fake_fetch_bytes)
    monkeypatch.setattr(mixtcrpred, "fetch_file", fake_fetch_file)

    # Offline mode blocks uncached download before calling datacache transfer.
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")
    with pytest.raises(OfflineError, match="OPENVAX_OFFLINE"):
        mixtcrpred.fetch_models(home, models=[model])
    assert bytes_calls == []
    assert file_calls == []

    # Online download uses datacache.fetch_bytes + datacache.fetch_file.
    monkeypatch.delenv("OPENVAX_OFFLINE")
    fetched = mixtcrpred.fetch_models(home, models=[model])
    assert fetched[0].status == "ready"
    assert Path(fetched[0].path).read_bytes() == content
    assert bytes_calls == [mixtcrpred.ZENODO_API_URL]
    assert file_calls == ["https://zenodo.example/model.ckpt"]

    # Verified cached model works offline; refresh_metadata=True is blocked.
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")
    assert mixtcrpred.fetch_models(home, models=[model])[0].status == "ready"
    with pytest.raises(OfflineError, match="OPENVAX_OFFLINE"):
        mixtcrpred.fetch_models(home, models=[model], refresh_metadata=True)


def test_snapshot_and_mhcflurry_fetch_enforce_offline(monkeypatch, tmp_path):
    _isolate_env(monkeypatch, tmp_path)
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")

    def unexpected_git(*args, **kwargs):
        pytest.fail("git should not be invoked in offline mode")

    monkeypatch.setattr(artifacts, "_run_git", unexpected_git)
    with pytest.raises(OfflineError, match="OPENVAX_OFFLINE"):
        artifacts.fetch("eramer", data_dir=tmp_path / "cache")

    with pytest.raises(OfflineError, match="OPENVAX_OFFLINE"):
        artifacts._fetch_mhcflurry(
            name="mhcflurry",
            download_name="models_class1_presentation",
            relative_path="models",
            version="2.2.0",
        )


def test_cli_help_delegates_to_datacache_cli_help():
    assert cli_help.HelpFormatter is dc_cli_help.HelpFormatter
    assert cli_help.enable_argcomplete is dc_cli_help.enable_argcomplete
    assert cli_help.style_table_header is dc_cli_help.style_table_header
    assert cli_help.style_status is dc_cli_help.style_status
    assert cli_help.classify_status is dc_cli_help.classify_status
    assert cli_help.OPTION_PATTERN is dc_cli_help.OPTION_PATTERN
    assert issubclass(cli_help.HelpArgumentParser, dc_cli_help.HelpArgumentParser)
