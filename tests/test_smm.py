"""Offline regressions using recorded upstream SMM vaccine predictions."""

import csv
import hashlib
import json
import math
from pathlib import Path
import socket
import subprocess
import sys
from types import SimpleNamespace

import pytest

from mhctools import SMM, SMMPMBEC, artifact_status
from mhctools.cli.args import make_mhc_arg_parser, predictors_from_args
from mhctools.pred import Kind, value_unit
from mhctools.smm import parse_smm_output

ROOT = Path(__file__).parent / "data/osteosarc"
FIXTURES = ROOT / "smm"
PEPTIDES = (ROOT / "inputs/class_i.txt").read_text().splitlines()
ALLELES = ["HLA-A*01:01", "HLA-B*08:01"]


def native_rows(method):
    return [row for length in (9, 10) for row in csv.DictReader(
        (FIXTURES / ("%s-%d.tsv" % (method, length))).read_text().splitlines(), delimiter="\t")]


@pytest.fixture(autouse=True)
def offline_only(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Offline tests must not execute predictors or contact a service")
    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket.socket, "connect_ex", forbidden)


def recorded_run(command, **kwargs):
    _, method, alleles, lengths, fasta = command
    length = int(lengths.split(",")[0])
    assert alleles == ",".join(ALLELES)
    assert lengths == "%d,%d" % (length, length)
    assert Path(fasta).read_text() == (FIXTURES / ("peptides-%d.fasta" % length)).read_text()
    return SimpleNamespace(returncode=0, stdout=(FIXTURES / ("%s-%d.tsv" % (method, length))).read_text(), stderr="")


def test_capture_integrity_and_source_lineage():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    assert hashlib.sha256((ROOT / "inputs/class_i.txt").read_bytes()).hexdigest() == manifest["source_panel_sha256"]
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest() == digest
    assert sorted(p for length in (9, 10) for p in
                  (FIXTURES / ("peptides-%d.fasta" % length)).read_text().splitlines()[1::2]) == PEPTIDES


@pytest.mark.parametrize("name,method", [
    ("smm", "smm"), ("smm-iedb", "smm"),
    ("smm-pmbec", "smmpmbec"), ("smm-pmbec-iedb", "smmpmbec"),
])
def test_cli_runs_local_vaccine_batch_and_preserves_affinity(name, method, monkeypatch):
    monkeypatch.setattr(subprocess, "run", recorded_run)
    args = make_mhc_arg_parser().parse_args([
        "--mhc-predictor", name, "--mhc-alleles", ",".join(ALLELES),
        "--mhc-predictor-path", sys.executable])
    predictor, = predictors_from_args(args)
    submitted = PEPTIDES + PEPTIDES[:1]
    results = predictor.predict(submitted)
    assert [r.preds[0].peptide for r in results] == submitted
    assert [r.preds[0].source_sequence_name for r in results] == [
        "seq%d" % (i + 1) for i in range(len(submitted))]
    predictions = [prediction for result in results for prediction in result.preds]
    expected = {(r["peptide"], r["allele"]): r for r in native_rows(method)}
    assert len(expected) == 38
    assert len(predictions) == 40
    for prediction in predictions:
        row = expected[prediction.peptide, prediction.allele]
        assert prediction.kind == Kind.pMHC_affinity
        assert value_unit(prediction.kind) == "nM"
        assert prediction.value == float(row["ic50"])
        assert prediction.percentile_rank == float(row["rank"])
        assert prediction.score == pytest.approx(1 - math.log(float(row["ic50"])) / math.log(50000))
    # Preserve the model distinction; these are actual upstream outputs.
    anchor = expected["AAPVATPAL", "HLA-A*01:01"]
    assert float(anchor["ic50"]) == pytest.approx(
        153656.16863502964 if method == "smm" else 176368.08525326516)


@pytest.mark.parametrize("change", [
    lambda rows: rows[1:],  # a missing peptide/allele pair
    lambda rows: rows + rows[:1],  # duplicate row
    lambda rows: [dict(rows[0], seq_num="0")] + rows[1:],
    lambda rows: [dict(rows[0], start="2")] + rows[1:],
    lambda rows: [dict(rows[0], peptide="SIINFEKLA")] + rows[1:],
    lambda rows: [dict(rows[0], ic50="nan")] + rows[1:],
    lambda rows: [dict(rows[0], allele="HLA-A*02:01")] + rows[1:],
])
def test_rejects_corrupt_or_incomplete_output(change, monkeypatch):
    import io
    rows = native_rows("smm")
    rows = [r for r in rows if r["length"] == "9"]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(change(rows))
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: SimpleNamespace(
        returncode=0, stdout=output.getvalue(), stderr="upstream detail"))
    with pytest.raises(ValueError, match="Local IEDB stderr: upstream detail"):
        SMM(alleles=ALLELES, program_name=sys.executable).predict_peptides(
            [p for p in PEPTIDES if len(p) == 9])


@pytest.mark.parametrize("returncode,stdout", [(0, "Unsupported allele/length"), (1, "")])
def test_upstream_errors_are_not_suppressed(returncode, stdout, monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: SimpleNamespace(
        returncode=returncode, stdout=stdout, stderr="missing model file"))
    with pytest.raises((RuntimeError, ValueError), match="missing model file"):
        SMMPMBEC(alleles=ALLELES, program_name=sys.executable).predict_peptides(PEPTIDES)


def test_missing_explicit_executable_never_falls_back_to_path(monkeypatch, tmp_path):
    missing = str(tmp_path / "missing-iedb")
    monkeypatch.setenv("IEDB_MHCI_EXECUTABLE", sys.executable)
    with pytest.raises(FileNotFoundError, match="Local IEDB SMM tools not found"):
        SMM(alleles=ALLELES, program_name=missing)
    assert artifact_status("smm-pmbec").path == sys.executable
    # A stale explicit override must not make inventory find a different
    # launcher on PATH while inference correctly rejects the override.
    launcher = tmp_path / "iedb-mhci"
    launcher.write_text("#!/bin/sh\nexit 0\n")
    launcher.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setenv("IEDB_MHCI_EXECUTABLE", missing)
    assert artifact_status("smm").status == "missing"
    with pytest.raises(FileNotFoundError):
        SMM(alleles=ALLELES)


def test_no_alleles_fails_before_inference():
    with pytest.raises(ValueError, match="requires at least one"):
        SMM(program_name=sys.executable).predict_peptides(PEPTIDES)


def test_parser_uses_native_sequence_numbers_not_score_sort_order():
    peptides = [p for p in PEPTIDES if len(p) == 9]
    parsed = parse_smm_output((FIXTURES / "smm-9.tsv").read_text(), peptides, "smm")
    assert parsed[0].peptide != peptides[0]
    assert all(p.peptide == peptides[int(p.source_sequence_name[3:]) - 1] for p in parsed)
