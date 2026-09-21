#!/usr/bin/env python3
"""Record local SMM outputs for the source-linked osteosarc class-I panel."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

from osteosarc_fixtures import tree_identity


ARCHIVE_SHA256 = "1cea64173886cc612d686313d4cb035c986c908e9042dab9cfa9a2bd492d2e31"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--installation", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--program", default=os.environ.get("IEDB_MHCI_EXECUTABLE", "iedb-mhci"))
    args = parser.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("Expected the pinned IEDB 3.1.7 archive")
    panel = Path(__file__).resolve().parents[1] / "tests/data/osteosarc/inputs/class_i.txt"
    peptides = panel.read_text().splitlines()
    args.output.mkdir(parents=True, exist_ok=False)
    captures = []
    for length in sorted(set(map(len, peptides))):
        fasta = args.output / ("peptides-%d.fasta" % length)
        batch = [p for p in peptides if len(p) == length]
        fasta.write_text("".join(">seq%d\n%s\n" % (i + 1, p) for i, p in enumerate(batch)))
        for method in ("smm", "smmpmbec"):
            stem = "%s-%d" % (method, length)
            command = [args.program, method, "HLA-A*01:01,HLA-B*08:01",
                       "%d,%d" % (length, length), str(fasta.resolve())]
            result = subprocess.run(command, check=True, capture_output=True)
            (args.output / (stem + ".tsv")).write_bytes(result.stdout)
            (args.output / (stem + ".stderr")).write_bytes(result.stderr)
            captures.append(dict(method=method, length=length, command=[
                arg.replace(str(args.output.resolve()), "{outputs}") for arg in command]))
    manifest = dict(
        archive_url="https://downloads.iedb.org/tools/mhci/3.1.7/IEDB_MHC_I-3.1.7.tar.gz",
        archive_sha256=ARCHIVE_SHA256, model_version="1.0", bundle_version="3.1.7",
        recorded_at=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
        platform=platform.platform(), captures=captures,
        installation=tree_identity(args.installation),
        source_panel_sha256=hashlib.sha256(panel.read_bytes()).hexdigest(),
        files={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.iterdir()})
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
