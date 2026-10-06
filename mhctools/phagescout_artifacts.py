"""Pinned, optional CC BY 4.0 PhageScout DESeq2 peptide tables."""

import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from urllib.parse import quote
from urllib.request import urlopen


RECORD = "21387981"
LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/legalcode"
# enzyme, source filename, bytes, SHA-256: unmodified Zenodo record files.
FILES = (
    ("ELANE", "peptide_profile_elastase.txt", 16321225,
     "43f0eab839fcfed19d6f2d6c696b379167919be6a33d44a3a92e31a50085f681"),
    ("CTSG", "peptide_profile_cathepsin G.txt", 16304718,
     "e8d3a86c7d3fd4e44735f8036bc35f038342ebf8ea95c1f8669187337f7b6b0c"),
)


def asset(enzyme):
    """Return the pinned file identity and immutable-record download URL."""
    _, name, size, digest = next(row for row in FILES if row[0] == enzyme)
    url = "https://zenodo.org/api/records/%s/files/%s/content" % (RECORD, quote(name))
    return name, size, digest, url


def managed_directory(data_dir=None):
    from .artifacts import data_path
    return data_path(data_dir) / "artifacts" / "phagescout" / RECORD


def profile_directory(profile_dir=None, data_dir=None):
    """Explicit profile path, then PHAGESCOUT_HOME, then managed data path."""
    configured = profile_dir
    if configured is None and data_dir is None:
        configured = os.environ.get("PHAGESCOUT_HOME")
    return (Path(configured).expanduser().resolve() if configured is not None else
            managed_directory(data_dir))


def verify_file(path, size, expected):
    """Verify exact source bytes without reading a whole table into memory."""
    if path.stat().st_size != size:
        raise ValueError("PhageScout file size mismatch: %s" % path)
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise ValueError("PhageScout file checksum mismatch: %s" % path)


def _manifest():
    return dict(name="phagescout", version=RECORD,
                source="https://zenodo.org/records/" + RECORD,
                license="CC-BY-4.0", license_url=LICENSE_URL,
                attribution="Enoch Yu, Colin Kretz and Matt Holding; PhageScout",
                files=[dict(name=name, size=size, sha256=digest, url=url)
                       for name, size, digest, url in (asset(row[0]) for row in FILES)])


def _directory_error(path, managed=False):
    try:
        for enzyme, _, _, _ in FILES:
            name, size, digest, _ = asset(enzyme)
            verify_file(path / name, size, digest)
        if managed:
            if json.loads((path / ".mhctools-artifact.json").read_text()) != _manifest():
                return "Invalid PhageScout artifact provenance"
            expected_license = (Path(__file__).parent / "data/PHAGESCOUT_LICENSE.txt").read_bytes()
            if (path / "PHAGESCOUT_LICENSE.txt").read_bytes() != expected_license:
                return "Invalid PhageScout license/attribution file"
    except (OSError, ValueError) as error:
        return str(error)
    return None


def status(data_dir=None):
    """Inventory verifies hashes, without parsing or downloading profiles."""
    from .artifacts import ArtifactStatus
    path = profile_directory(data_dir=data_dir)
    managed = path == managed_directory(data_dir)
    error = _directory_error(path, managed=managed)
    return ArtifactStatus(
        name="phagescout", status="missing" if error else "ready",
        manager="mhctools" if managed else "user", version=RECORD,
        path=str(path), fetchable=True,
        detail=(error + "; run mhctools fetch phagescout" if error else
                "Verified full ELANE/CTSG DESeq2 peptide tables; CC BY 4.0"))


def fetch_profiles(version=None, data_dir=None):
    """Stream, verify and atomically install both tables and provenance."""
    if version is not None and str(version) != RECORD:
        raise ValueError("PhageScout supports only Zenodo record " + RECORD)
    current = status(data_dir)
    if current.status == "ready":
        return current
    target = managed_directory(data_dir)
    if Path(current.path) != target:
        raise RuntimeError("PHAGESCOUT_HOME is incomplete or invalid: %s. Correct it or unset it "
                           "before running mhctools fetch phagescout." % current.path)
    if target.exists():
        raise RuntimeError("PhageScout directory exists but is incomplete or invalid: %s. "
                           "Move it aside and fetch again." % target)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".phagescout-", dir=target.parent))
    try:
        for enzyme, _, _, _ in FILES:
            name, size, expected, url = asset(enzyme)
            digest, downloaded = hashlib.sha256(), 0
            with urlopen(url, timeout=30) as source, (staging / name).open("wb") as output:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    downloaded += len(block)
                    if downloaded > size:
                        raise ValueError("PhageScout download exceeds pinned size: " + name)
                    digest.update(block)
                    output.write(block)
                output.flush()
                os.fsync(output.fileno())
            if downloaded != size or digest.hexdigest() != expected:
                raise ValueError("PhageScout download checksum/size mismatch: " + name)
        (staging / ".mhctools-artifact.json").write_text(
            json.dumps(_manifest(), indent=2, sort_keys=True) + "\n")
        shutil.copyfile(Path(__file__).parent / "data/PHAGESCOUT_LICENSE.txt",
                        staging / "PHAGESCOUT_LICENSE.txt")
        try:
            staging.rename(target)
        except OSError as error:
            if error.errno not in (errno.ENOTEMPTY, errno.EEXIST):
                raise
            # Another complete fetch may win the atomic directory rename.
            if _directory_error(target, managed=True):
                raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return status(data_dir)


def selected_profile(enzyme, profile_dir=None):
    """Resolve one explicit model, verifying its exact source identity."""
    path = profile_directory(profile_dir)
    name, size, digest, _ = asset(enzyme)
    try:
        verify_file(path / name, size, digest)
    except (OSError, ValueError) as error:
        raise RuntimeError("Full PhageScout profile unavailable or invalid: %s. "
                           "Run mhctools fetch phagescout; configure PHAGESCOUT_HOME or "
                           "PhageScout(profile_dir=...) for a custom directory." % error) from error
    return path / name, digest
