#!/usr/bin/env python3
"""Export a pinned osteosarc snapshot and record real predictor CLI outputs.

Regeneration is explicit and writes a NEW directory. Neither osteosarc nor
predictor executables are needed by the offline regression tests.
"""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export_inputs(args):
    from osteosarc import Dataset
    from osteosarc.cache import Cache

    data = Dataset.open(args.snapshot, cache=Cache(args.cache), offline=True)
    if data.id != args.snapshot_id:
        raise ValueError("Snapshot ID differs from the requested immutable snapshot")
    annotations = list(data.annotations)
    peptides = list(data.vaccine_peptides())
    records, missing = [], []
    cursor = 0
    for variant in annotations:
        disclosed = variant.get("vaccine_peptides", [])
        for index, peptide in enumerate(disclosed, 1):
            api_record = peptides[cursor]
            cursor += 1
            assert api_record == dict(
                peptide, variant_id=variant["id"], gene=variant["gene"])
            sequence = api_record["sequence"]
            if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
                raise ValueError("Noncanonical disclosed sequence: %r" % sequence)
            records.append(dict(
                api_record,
                record_id="%s:vaccine-peptide-%d" % (variant["id"], index),
                source_url="https://osteosarc.com/variant/%s/" % variant["id"],
                validation=variant["validation"],
                corrections=variant["corrections"],
            ))
        for vaccine, member in variant["vaccines"].items():
            if member and not any(vaccine in p["in_vaccines"] for p in disclosed):
                missing.append(dict(variant_id=variant["id"], vaccine=vaccine))
    assert cursor == len(peptides)

    # Keep record identity separate from sequence identity. This preserves the
    # EPG5/EXOC4 source conflict instead of silently choosing one variant.
    groups = {name: defaultdict(list) for name in
              ("class_i", "cons_9mer", "class_ii", "cleavage")}
    for record in records:
        if not record["in_vaccines"]:
            continue  # retain associated test peptides in source.json only
        sequence = record["sequence"]
        reference = dict(record_id=record["record_id"], offset=0)
        groups["cleavage"][sequence].append(reference)
        if 8 <= len(sequence) <= 11:
            groups["class_i"][sequence].append(reference)
        if len(sequence) == 9:
            groups["cons_9mer"][sequence].append(reference)
        if len(sequence) >= 15:
            offset = (len(sequence) - 15) // 2
            groups["class_ii"][sequence[offset:offset + 15]].append(
                dict(record_id=record["record_id"], offset=offset))
    panels = {
        name: [dict(id="p%03d" % index, sequence=sequence, sources=refs)
               for index, (sequence, refs) in enumerate(sorted(sequences.items()))]
        for name, sequences in groups.items()
    }
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "source.json", dict(
        schema_version=1,
        osteosarc_version=version("osteosarc"),
        snapshot=data.manifest,
        curation_report=list(data.corrections),
        records=records,
        undisclosed_memberships=missing,
    ))
    write_json(output / "panels.json", panels)
    for name, panel in panels.items():
        (output / (name + ".txt")).write_text(
            "".join(row["sequence"] + "\n" for row in panel))
        (output / (name + ".fasta")).write_text(
            "".join(">%s\n%s\n" % (row["id"], row["sequence"]) for row in panel))


def tree_identity(root):
    """Hash installed code/models without redistributing licensed assets.

    Hash = SHA256 of sorted relative-path + NUL + file-SHA256 + newline.
    Scratch directories, VCS files and macOS metadata are excluded.
    """
    digest = hashlib.sha256()
    count = 0
    size = 0
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part in {".git", "tmp", "__pycache__"} or part.startswith(".")
               for part in relative.parts) or not path.is_file():
            continue
        digest.update((relative.as_posix() + "\0" + sha256(path) + "\n").encode())
        count += 1
        size += path.stat().st_size
    return dict(sha256=digest.hexdigest(), files=count, bytes=size)


def record_outputs(args):
    from mhctools import NetChop

    inputs = Path(args.inputs).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    bundle = Path(os.environ["NETMHC_BUNDLE_HOME"])
    mix = Path(os.environ["MIXMHCPRED_PATH"])
    mix2 = Path(os.environ["MIXMHC2PRED_EXECUTABLE"])
    prime = Path(os.environ["PRIME_EXECUTABLE"])
    class_i = ["HLA-A*01:01", "HLA-B*08:01"]
    class_ii = ["HLA-DRB1*03:01", "HLA-DRB1*08:01"]
    captures = []

    def run(name, panel, command, roots, alleles=(), native_names=(), result=None,
            runtime=None):
        print("Recording " + name, flush=True)
        completed = subprocess.run(command, capture_output=True, timeout=600, check=True)
        (output / (name + ".stdout")).write_bytes(completed.stdout)
        (output / (name + ".stderr")).write_bytes(completed.stderr)
        if result:
            if not (output / result).is_file():
                raise RuntimeError("Predictor produced no result: " + name)
            raw = result
        else:
            raw = name + ".stdout"
        executable = Path(shutil.which(command[0]) or command[0]).resolve()
        captures.append(dict(
            name=name, panel=panel, alleles=list(alleles),
            native_alleles=list(native_names), raw=raw,
            command=[str(arg).replace(str(inputs), "{inputs}").replace(
                str(output), "{outputs}") for arg in command],
            executable_sha256=sha256(executable),
            installations={label: tree_identity(Path(root)) for label, root in roots.items()},
            runtime=runtime,
        ))

    run("netmhcpan", "class_i", [
        "netMHCpan-4.1", "-f", str(inputs / "class_i.fasta"),
        "-l", "9,10", "-a", "HLA-A01:01,HLA-B08:01", "-BA",
    ], {"netMHCpan-4.1": bundle / "netMHCpan-4.1"}, class_i)
    run("netmhciipan", "class_ii", [
        "netMHCIIpan-4.3", "-f", str(inputs / "class_ii.fasta"),
        "-length", "15", "-a", "DRB1_0301,DRB1_0801", "-BA",
    ], {"netMHCIIpan-4.3": bundle / "netMHCIIpan-4.3"}, class_ii)
    legacy_image = os.environ.get("MHCTOOLS_LEGACY_NETMHC_IMAGE", "mhctools-test-netmhc-legacy:1")
    legacy_identity = subprocess.check_output(
        ["docker", "image", "ls", "--no-trunc", "--format", "{{.ID}}", legacy_image],
        text=True).strip()
    if not legacy_identity:
        raise RuntimeError("Preload the legacy predictor image: " + legacy_image)
    for allele, (label, native) in zip(class_i, (
            ("a0101", "HLA-A01:01"), ("b0801", "HLA-B08:01"))):
        run("netmhccons-" + label, "cons_9mer", [
            "netMHCcons", "-f", str(inputs / "cons_9mer.fasta"),
            "-length", "9", "-a", native,
        ], {name: bundle / name for name in (
            "netMHCcons-1.1", "netMHCpan-2.8", "netMHC-3.4", "pickpocket-1.1")},
            [allele], runtime=dict(image=legacy_image, image_id=legacy_identity))
    run("mixmhcpred", "class_i", [
        str(mix), "-i", str(inputs / "class_i.txt"),
        "-o", str(output / "mixmhcpred.tsv"), "-a", ",".join(class_i),
    ], {"MixMHCpred": mix.parent}, class_i, result="mixmhcpred.tsv")
    run("prime", "class_i", [
        str(prime), "-i", str(inputs / "class_i.txt"),
        "-o", str(output / "prime.tsv"), "-a", ",".join(class_i), "-mix", str(mix),
    ], {"PRIME": prime.parent, "MixMHCpred": mix.parent}, class_i, result="prime.tsv")
    native_ii = ["DRB1_03_01", "DRB1_08_01"]
    run("mixmhc2pred", "class_ii", [
        str(mix2), "-i", str(inputs / "class_ii.txt"),
        "-o", str(output / "mixmhc2pred.tsv"), "-a", *native_ii,
        "--no_context", "--extra_out",
    ], {"MixMHC2pred": mix2.parent}, class_ii, native_ii, result="mixmhc2pred.tsv")
    for model, label in ((0, "cterm"), (1, "20s")):
        predictor = NetChop(model_variant=model, execution="container")
        run("netchop-" + label, "cleavage",
            predictor._container_command(inputs, "cleavage.fasta"),
            {"NetChop-3.1": predictor.netchop_dir},
            runtime=dict(image=predictor.container_image, platform="linux/386"))
    write_json(output / "manifest.json", dict(
        schema_version=1, recorded_at=datetime.now(timezone.utc).isoformat(),
        mhctools_version=version("mhctools"), python=platform.python_version(),
        platform=platform.platform(),
        inputs={p.name: sha256(p) for p in sorted(inputs.iterdir()) if p.is_file()},
        captures=captures,
        files={p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()},
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    export = commands.add_parser("inputs")
    export.add_argument("--snapshot", required=True)
    export.add_argument("--snapshot-id", required=True)
    export.add_argument("--cache", required=True)
    export.add_argument("--output", required=True)
    record = commands.add_parser("record")
    record.add_argument("--inputs", required=True)
    record.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "inputs":
        export_inputs(args)
    else:
        record_outputs(args)


if __name__ == "__main__":
    main()
