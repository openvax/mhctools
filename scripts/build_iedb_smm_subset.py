#!/usr/bin/env python3
"""Derive the SMM/SMM-PMBEC subset of the official IEDB MHC-I 3.1.7 bundle.

The full release is a 341 MB download, 96% of which is method data mhctools
never uses from this bundle: pickpocket, netmhccons, netmhcpan and
netmhcstabpan account for 177 MB, and the netMHC family is wrapped separately
through its own licensed distribution. SMM and SMM-PMBEC need 2.3 MB of
matrices.

CI downloading the full archive from downloads.iedb.org made every merge
depend on a third-party academic host being reachable, which it was not on
2026-09-28 (curl exit 28, connection timeout, four attempts). This script
produces a byte-reproducible subset that can be hosted where CI can always
reach it, pinned by its own checksum.

The subset is a verbatim copy of the allowlisted paths; nothing is rewritten.
LIAI_license.txt is included because the bundle is Non-Profit Open Software
License 3.0 and the terms travel with the material.

Reproducibility: entries are sorted, and mtime, mode, uid/gid and uname/gname
are normalized, so a rebuild from the same official archive yields the same
SHA-256 and the hosted copy can be audited against upstream.

    python scripts/build_iedb_smm_subset.py \
        --archive env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
        --output dist/IEDB_MHC_I-3.1.7-smm-subset.tar.gz
"""

import argparse
import gzip
import hashlib
from pathlib import Path
import tarfile
import tempfile

# The official release this subset is derived from.
ARCHIVE_SHA256 = "1cea64173886cc612d686313d4cb035c986c908e9042dab9cfa9a2bd492d2e31"

# Paths kept, relative to the archive's mhc_i/ directory. An allowlist rather
# than an exclude list: a new method directory upstream must be considered
# deliberately instead of silently inflating the subset.
# These mirror the prefixes setup_test_backends.smm() already installs, with
# data/ narrowed from every method's training data to the two SMM methods.
# method/ as a whole is not included: in the archive it holds the bundled tool
# implementations, which is most of the 341 MB.
INCLUDED_PREFIXES = (
    "LIAI_license.txt",
    "src/",
    "method/allele-info/",
    "method/iedbtools-utilities/",
    "data/MHCI_mhcibinding20130222/smm/",
    "data/MHCI_mhcibinding20130222/smmpmbec/",
    # Percentile ranks come from the consensus score distributions
    # (distribution_consensus_bin.cpickle), so this is required for SMM even
    # though the consensus method itself is not used.
    "data/MHCI_mhcibinding20130222/consensus/",
)

EXCLUDED_SUFFIXES = ("__pycache__", ".pyc")


def _digest(path):
    sha256 = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(block)
    return sha256.hexdigest()


def _is_included(name):
    """Is this archive member part of the subset?

    Members arrive as ``mhc_i/...`` or with a leading release directory; match
    on the portion at and below ``mhc_i/``.
    """
    marker = "mhc_i/"
    index = name.find(marker)
    if index == -1:
        return False
    relative = name[index + len(marker):]
    if not relative:
        return False
    if any(part in relative for part in EXCLUDED_SUFFIXES):
        return False
    return any(
        relative == prefix or relative.startswith(prefix)
        for prefix in INCLUDED_PREFIXES)


def _normalize(info):
    """Strip filesystem noise so the same input yields the same checksum."""
    info.mtime = 0
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.mode = 0o755 if info.isdir() or info.mode & 0o100 else 0o644
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive", type=Path, required=True,
        help="Official IEDB_MHC_I-3.1.7.tar.gz")
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Subset tarball to write")
    parser.add_argument(
        "--skip-archive-checksum", action="store_true",
        help="Skip verifying the source archive (for local iteration only)")
    args = parser.parse_args()

    if not args.skip_archive_checksum:
        digest = _digest(args.archive)
        if digest != ARCHIVE_SHA256:
            raise SystemExit(
                "Expected the pinned IEDB 3.1.7 archive (%s), got %s"
                % (ARCHIVE_SHA256, digest))

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Stage in one streaming pass. Seeking back to each member in a gzip
    # stream re-decompresses from the start, so collecting the member list
    # first and then extracting one by one reads the 341 MB archive once per
    # file; staging to disk reads it once in total.
    with tempfile.TemporaryDirectory() as staging_name:
        staging = Path(staging_name)
        with tarfile.open(args.archive, "r|gz") as source:
            for member in source:
                if not _is_included(member.name):
                    continue
                relative = member.name[member.name.find("mhc_i/"):]
                target = (staging / relative).resolve()
                if not str(target).startswith(str(staging.resolve())):
                    raise SystemExit("Unsafe archive path: %s" % member.name)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif member.isreg():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with source.extractfile(member) as handle:
                        target.write_bytes(handle.read())
                    target.chmod(0o755 if member.mode & 0o100 else 0o644)

        staged = sorted(
            path for path in staging.rglob("*")
            if not any(part in str(path) for part in EXCLUDED_SUFFIXES))
        if not staged:
            raise SystemExit("No members matched the subset allowlist")

        # gzip's mtime and the stored filename are both part of the output
        # bytes. Without filename="" the header records the output path, so
        # building to two different names gives two different checksums.
        with open(args.output, "wb") as raw:
            with gzip.GzipFile(
                    filename="", fileobj=raw, mode="wb",
                    compresslevel=9, mtime=0) as compressed:
                with tarfile.open(
                        fileobj=compressed, mode="w",
                        format=tarfile.GNU_FORMAT) as subset:
                    for path in staged:
                        info = subset.gettarinfo(
                            str(path), arcname=str(path.relative_to(staging)))
                        if info.isreg():
                            with open(path, "rb") as handle:
                                subset.addfile(_normalize(info), handle)
                        else:
                            subset.addfile(_normalize(info))

        members = staged
        total = sum(
            path.stat().st_size for path in staged if path.is_file())

    print("members:  %d" % len(members))
    print("bytes:    %d (%.1f MB uncompressed)" % (total, total / 1e6))
    print("output:   %s (%.1f MB)"
          % (args.output, args.output.stat().st_size / 1e6))
    print("sha256:   %s" % _digest(args.output))


if __name__ == "__main__":
    main()
