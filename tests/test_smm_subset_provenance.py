# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""The hosted SMM subset must stay traceable to one upstream IEDB release.

The subset builder and the fixture recorder each pin the official archive by
SHA-256, and the builder also pins the digest of the subset it derives from
it. The installer pins only that derived digest, which carries no upstream
version of its own -- so the chain from "what CI installs" back to "which
IEDB release" holds only if the installer's digest equals the one the builder
declares. Assert that, or a rebuilt-and-reuploaded asset could change what CI
installs while the recorded fixtures keep describing the old release.
"""

import importlib.util
import hashlib
import io
from pathlib import Path
import re
import sys
import tarfile

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name):
    path = SCRIPTS / ("%s.py" % name)
    if not path.exists():
        pytest.skip("%s not present" % path)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    # These scripts import each other as top-level modules when run from
    # scripts/, so reproduce that import path here.
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


def test_subset_builder_and_fixture_recorder_pin_the_same_release():
    builder = _load("build_iedb_smm_subset")
    recorder = _load("record_smm_fixtures")

    assert builder.ARCHIVE_SHA256 == recorder.ARCHIVE_SHA256


def test_installer_pins_a_subset_checksum():
    installer = _load("setup_test_backends")

    assert re.fullmatch(r"[0-9a-f]{64}", installer.SMM_SUBSET_SHA256)


def test_installer_digest_is_the_one_the_builder_produces():
    """Ties the installed artifact to a build from the pinned IEDB release.

    Without this, uploading a modified asset and updating only
    SMM_SUBSET_SHA256 passes every other check here.
    """
    builder = _load("build_iedb_smm_subset")
    installer = _load("setup_test_backends")

    assert installer.SMM_SUBSET_SHA256 == builder.EXPECTED_SUBSET_SHA256


def test_subset_carries_the_upstream_licenses():
    """NPOSL-3.0 material is redistributed, so its notices must travel."""
    builder = _load("build_iedb_smm_subset")

    for notice in ("LIAI_license.txt", "Copenhagen_license.txt", "README"):
        assert notice in builder.INCLUDED_PREFIXES


def test_installer_downloads_the_release_the_builder_produces():
    """The installed artifact must be the one build_iedb_smm_subset writes."""
    installer = _load("setup_test_backends")

    assert installer.SMM_SUBSET_URL.startswith(
        "https://github.com/openvax/mhctools/releases/download/")
    assert installer.SMM_SUBSET_URL.endswith(
        "IEDB_MHC_I-3.1.7-smm-subset.tar.gz")


def test_subset_allowlist_keeps_the_license_and_percentile_data():
    """Both are load-bearing and easy to drop while pruning for size."""
    builder = _load("build_iedb_smm_subset")

    assert "LIAI_license.txt" in builder.INCLUDED_PREFIXES
    # Percentile ranks read distribution_consensus_bin.cpickle from here.
    assert any(
        prefix.endswith("consensus/")
        for prefix in builder.INCLUDED_PREFIXES)
    # The bulk of the upstream archive must stay out.
    assert not any(
        name in prefix
        for prefix in builder.INCLUDED_PREFIXES
        for name in ("netmhcpan", "netmhccons", "netmhcstabpan", "pickpocket"))


@pytest.fixture
def source_files():
    return {
        "mhc_i/LIAI_license.txt": b"LIAI notice\n",
        "mhc_i/Copenhagen_license.txt": b"Copenhagen notice\n",
        "mhc_i/README": b"Upstream readme\n",
        "mhc_i/src/setupinfo.template": b"home = '%s'\n",
        "mhc_i/method/allele-info/alleles.txt": b"allele metadata\n",
        "mhc_i/method/iedbtools-utilities/util.py": b"# utility\n",
        "mhc_i/data/MHCI_mhcibinding20130222/smm/model": b"SMM matrix\n",
        "mhc_i/data/MHCI_mhcibinding20130222/smmpmbec/model": b"PMBEC matrix\n",
        "mhc_i/data/MHCI_mhcibinding20130222/consensus/distribution": b"ranks\n",
    }


def _write_archive(path, files, extra_members=()):
    with tarfile.open(path, "w:gz") as archive:
        for name, payload in reversed(list(files.items())):
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            info.mode = 0o755 if name.endswith(".py") else 0o664
            info.mtime = 123456789
            info.uid = info.gid = 42
            info.uname = info.gname = "upstream"
            archive.addfile(info, io.BytesIO(payload))
        for member in extra_members:
            archive.addfile(member)


def _build_subset(monkeypatch, archive, output, **pins):
    builder = _load("build_iedb_smm_subset")
    monkeypatch.setattr(builder, "ARCHIVE_SHA256", pins.get(
        "source", hashlib.sha256(archive.read_bytes()).hexdigest()))
    # Synthetic archives exercise the real builder without the 341 MB download.
    monkeypatch.setattr(builder, "EXPECTED_SUBSET_SHA256",
                        pins.get("subset", "PLACEHOLDER"))
    monkeypatch.setattr(sys, "argv", [
        str(SCRIPTS / "build_iedb_smm_subset.py"),
        "--archive", str(archive), "--output", str(output)])
    builder.main()


def test_builder_preserves_payload_and_normalizes_metadata(
        tmp_path, monkeypatch, source_files):
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, dict(source_files, **{
        "mhc_i/README.backup": b"not a notice",
        "mhc_i/src/__pycache__/cached.pyc": b"not source",
        "mhc_i/method/netmhcpan/program": b"not needed",
    }))
    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"
    _build_subset(monkeypatch, archive, first)
    _build_subset(monkeypatch, archive, second,
                  subset=hashlib.sha256(first.read_bytes()).hexdigest())
    assert first.read_bytes() == second.read_bytes()
    with tarfile.open(first) as subset:
        assert subset.getnames() == sorted(source_files)
        for member in subset:
            assert subset.extractfile(member).read() == source_files[member.name]
            assert (member.mtime, member.uid, member.gid) == (0, 0, 0)
            assert (member.uname, member.gname) == ("", "")
            assert member.mode == (0o755 if member.name.endswith(".py") else 0o644)


@pytest.mark.parametrize("empty_directory", [False, True])
def test_builder_requires_files_from_every_prefix(
        tmp_path, monkeypatch, source_files, empty_directory):
    prefix = "mhc_i/data/MHCI_mhcibinding20130222/smm/"
    del source_files[prefix + "model"]
    directory = tarfile.TarInfo(prefix)
    directory.type = tarfile.DIRTYPE
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, source_files, [directory] if empty_directory else [])
    output = tmp_path / "subset.tar.gz"
    with pytest.raises(SystemExit, match="No regular archive file matched: .*smm/"):
        _build_subset(monkeypatch, archive, output)
    assert not output.exists()


@pytest.mark.parametrize("member_type", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE])
def test_builder_rejects_unsupported_members(
        tmp_path, monkeypatch, source_files, member_type):
    member = tarfile.TarInfo("mhc_i/src/extra")
    member.type = member_type
    member.linkname = "setupinfo.template"
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, source_files, [member])
    output = tmp_path / "subset.tar.gz"
    with pytest.raises(SystemExit, match="Unsupported archive member type"):
        _build_subset(monkeypatch, archive, output)
    assert not output.exists()


def test_builder_rejects_wrong_source_checksum(tmp_path, monkeypatch, source_files):
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, source_files)
    output = tmp_path / "subset.tar.gz"
    with pytest.raises(SystemExit, match="Expected the pinned IEDB"):
        _build_subset(monkeypatch, archive, output, source="0" * 64)
    assert not output.exists()


def test_builder_rejects_wrong_subset_checksum(tmp_path, monkeypatch, source_files):
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, source_files)
    with pytest.raises(SystemExit, match="does not match EXPECTED_SUBSET_SHA256"):
        _build_subset(monkeypatch, archive, tmp_path / "subset.tar.gz", subset="0" * 64)


def test_installer_replaces_stale_tree_after_verifying_checksum(
        tmp_path, monkeypatch, source_files):
    archive = tmp_path / "upstream.tar.gz"
    _write_archive(archive, source_files)
    subset = tmp_path / "IEDB_MHC_I-3.1.7-smm-subset.tar.gz"
    _build_subset(monkeypatch, archive, subset)
    installer = _load("setup_test_backends")
    config = {}
    stale = tmp_path / "iedb-3.1.7/mhc_i/data/old-model"
    stale.parent.mkdir(parents=True)
    stale.write_text("old data")
    monkeypatch.setattr(installer, "SMM_SUBSET_SHA256", "0" * 64)
    with pytest.raises(SystemExit, match="Unexpected IEDB SMM subset checksum"):
        installer.smm(tmp_path, sys.executable, config)
    assert stale.read_text() == "old data"
    assert not config

    monkeypatch.setattr(installer, "SMM_SUBSET_SHA256",
                        hashlib.sha256(subset.read_bytes()).hexdigest())
    installer.smm(tmp_path, sys.executable, config)
    assert not stale.exists()
    for name, payload in source_files.items():
        assert (tmp_path / "iedb-3.1.7" / name).read_bytes() == payload
    source = tmp_path / "iedb-3.1.7/mhc_i"
    assert (source / "src/setupinfo.py").read_text() == "home = '%s'\n" % source
    assert Path(config["IEDB_MHCI_EXECUTABLE"]).is_file()
