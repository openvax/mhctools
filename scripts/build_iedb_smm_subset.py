#!/usr/bin/env python3
"""Derive the SMM/SMM-PMBEC subset of the official IEDB MHC-I 3.1.7 bundle.

The official release unpacks to 1031 MB across 38,236 members, dominated by
the bundled DTU executables under method/ (netmhc-4.0 alone is 210 MB,
netmhc-3.4 192 MB, netmhcpan-4.1 114 MB) and by 192 MB of per-method training
data under data/. mhctools does not run those executables from this bundle;
the netMHC family is wrapped through its own licensed distribution. The paths
SMM and SMM-PMBEC actually need come to 9.4 MB.

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
import io
from pathlib import Path
import tarfile

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
    # The upstream licenses and release notes travel with the material.
    "LIAI_license.txt",
    "Copenhagen_license.txt",
    "README",
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

EXCLUDED_NAMES = ("__pycache__", ".pyc")

# The digest this script produces from the pinned archive. setup_test_backends
# installs an artifact by checksum; asserting that checksum equals this one is
# what ties the hosted copy to a build from the official release rather than to
# whatever was uploaded.
EXPECTED_SUBSET_SHA256 = (
    "eef720a71991a76ceee64ce32fbb52e19684b91d9d44eb078c18c96f395e192e")


def _digest(path):
    sha256 = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha256.update(block)
    return sha256.hexdigest()


MARKER = "mhc_i/"


def _matched_prefix(name):
    """Which allowlist prefix admits this archive member, if any.

    Members arrive as ``mhc_i/...``; match on the portion at and below that.
    A prefix ending in "/" matches a directory subtree, anything else must
    match exactly, so "README" does not also admit "README.backup".
    """
    index = name.find(MARKER)
    if index == -1:
        return None
    relative = name[index + len(MARKER):]
    if not relative:
        return None
    if any(excluded in relative for excluded in EXCLUDED_NAMES):
        return None
    for prefix in INCLUDED_PREFIXES:
        if prefix.endswith("/"):
            if relative == prefix.rstrip("/") or relative.startswith(prefix):
                return prefix
        elif relative == prefix:
            return prefix
    return None


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

    # One streaming pass. Seeking back to a member in a gzip stream
    # re-decompresses from the start, so collecting the member list first and
    # extracting one by one would read the 341 MB archive once per file. The
    # selected payload is ~9 MB, so it is held in memory rather than staged
    # through a temporary directory.
    entries = []
    per_prefix = {prefix: 0 for prefix in INCLUDED_PREFIXES}
    with tarfile.open(args.archive, "r|gz") as source:
        for member in source:
            prefix = _matched_prefix(member.name)
            if prefix is None:
                continue
            per_prefix[prefix] += 1
            arcname = member.name[member.name.find(MARKER):]
            if member.isdir():
                entries.append((arcname, None, 0o755))
            elif member.isreg():
                entries.append((
                    arcname,
                    source.extractfile(member).read(),
                    0o755 if member.mode & 0o100 else 0o644))
            else:
                # Symlinks, hardlinks, devices and fifos. The 3.1.7 archive has
                # 15 symlinks, none under these prefixes. Silently dropping one
                # would ship a subset that installs and then fails at runtime.
                raise SystemExit(
                    "Unsupported archive member type %r for %s; the subset "
                    "builder only copies directories and regular files"
                    % (member.type, member.name))

    empty = [prefix for prefix, count in per_prefix.items() if not count]
    if empty:
        # An upstream rename (the data directory carries a dataset date) would
        # otherwise produce a cheerful build with no model data at all.
        raise SystemExit(
            "No archive member matched: %s" % ", ".join(sorted(empty)))

    # Sort on the archive name itself so ordering does not depend on pathlib's
    # comparison semantics, which have been reworked across CPython releases.
    entries.sort(key=lambda entry: entry[0])

    # gzip's mtime and stored filename are both part of the output bytes.
    # Without filename="" the header records the output path, so building the
    # same content to two paths gives two different checksums.
    with open(args.output, "wb") as raw:
        with gzip.GzipFile(
                filename="", fileobj=raw, mode="wb",
                compresslevel=9, mtime=0) as compressed:
            with tarfile.open(
                    fileobj=compressed, mode="w",
                    format=tarfile.GNU_FORMAT) as subset:
                for arcname, payload, mode in entries:
                    info = tarfile.TarInfo(arcname)
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = mode
                    if payload is None:
                        info.type = tarfile.DIRTYPE
                        subset.addfile(info)
                    else:
                        info.size = len(payload)
                        subset.addfile(info, io.BytesIO(payload))

    files = [entry for entry in entries if entry[1] is not None]
    total = sum(len(entry[1]) for entry in files)
    digest = _digest(args.output)
    print("files:      %d" % len(files))
    print("directories: %d" % (len(entries) - len(files)))
    print("bytes:      %d (%.1f MB uncompressed)" % (total, total / 1e6))
    print("output:     %s (%.1f MB)"
          % (args.output, args.output.stat().st_size / 1e6))
    print("sha256:     %s" % digest)
    if EXPECTED_SUBSET_SHA256 not in ("PLACEHOLDER", digest):
        raise SystemExit(
            "Built subset does not match EXPECTED_SUBSET_SHA256 (%s)"
            % EXPECTED_SUBSET_SHA256)


if __name__ == "__main__":
    main()
