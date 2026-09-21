"""Real vaccine inputs and recorded outputs, with no network or model runtime."""

from collections import Counter
import csv
import hashlib
import importlib
import json
import math
from pathlib import Path
import socket
import subprocess
import sys

import pytest

from mhctools import MixMHC2pred, NetChop, PRIME
from mhctools.allele_normalization import normalize_allele_name
from mhctools.mixmhc2pred import parse_mixmhc2pred_results
from mhctools.mixmhcpred import parse_mixmhcpred_output
from mhctools.parsing import (
    parse_netmhccons_stdout,
    parse_netmhciipan43_stdout,
    parse_netmhcpan_to_preds,
)
from mhctools.pred import Kind, value_unit
from mhctools.prime import parse_prime_results


ROOT = Path(__file__).parent / "data" / "osteosarc"
PANELS = json.loads((ROOT / "inputs" / "panels.json").read_text())


@pytest.fixture(autouse=True)
def offline_only(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Recorded fixture tests must not execute tools or use the network")

    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket.socket, "connect_ex", forbidden)


def native_table(name):
    with (ROOT / "outputs" / (name + ".tsv")).open() as stream:
        return list(csv.DictReader(
            (line for line in stream if not line.startswith("#")), delimiter="\t"))


def capture(name):
    manifest = json.loads((ROOT / "outputs" / "manifest.json").read_text())
    return next(row for row in manifest["captures"] if row["name"] == name)


def test_source_identity_disclosure_and_experimental_states():
    source = json.loads((ROOT / "inputs" / "source.json").read_text())
    records = source["records"]
    assert len(records) == len({r["record_id"] for r in records}) == 94
    assert sum(bool(r["in_vaccines"]) for r in records) == 78
    assert len(source["undisclosed_memberships"]) == 33
    assert Counter(r["validation"]["elispot_status"] for r in records) == {
        "positive": 46, "negative": 4, "not_tested": 44}
    for record in records:
        if record["validation"]["elispot_status"] == "not_tested":
            assert record["validation"]["elispot_response"] is None
    conflict = [r for r in records if r["sequence"] == "KKSVIRTLSTIDDVEDRENEKGR"]
    assert {r["variant_id"] for r in conflict} == {
        "EPG5-chr18-45925851", "EXOC4-chr7-133274996"}
    # Associated experimental peptides are retained without becoming vaccine members.
    associated = next(r for r in records if r["sequence"] == "EALLSKFSL")
    assert associated["in_vaccines"] == []
    assert associated["experiments"][0]["result"] == "Negative"


def test_panel_lineage_and_duplicate_sequence_provenance():
    records = {r["record_id"]: r for r in json.loads(
        (ROOT / "inputs" / "source.json").read_text())["records"]}
    assert {name: len(rows) for name, rows in PANELS.items()} == {
        "class_i": 19, "cons_9mer": 17, "class_ii": 55, "cleavage": 77}
    for panel in PANELS.values():
        assert len(panel) == len({row["id"] for row in panel})
        for row in panel:
            for reference in row["sources"]:
                parent = records[reference["record_id"]]
                offset = reference["offset"]
                assert parent["in_vaccines"]
                assert row["sequence"] == parent["sequence"][offset:offset + len(row["sequence"])]
    assert {ref["record_id"] for row in PANELS["cleavage"] for ref in row["sources"]} == {
        key for key, row in records.items() if row["in_vaccines"]}


def test_capture_integrity():
    manifest = json.loads((ROOT / "outputs" / "manifest.json").read_text())
    for directory, entries in (("inputs", manifest["inputs"]), ("outputs", manifest["files"])):
        for name, digest in entries.items():
            assert hashlib.sha256((ROOT / directory / name).read_bytes()).hexdigest() == digest
    assert len(manifest["captures"]) == 9
    for row in manifest["captures"]:
        assert row["raw"] in manifest["files"]
        assert row["installations"]
        assert all(asset["files"] > 0 for asset in row["installations"].values())


def test_netmhcpan_both_endpoints_and_window_offsets():
    text = (ROOT / "outputs" / "netmhcpan.stdout").read_text()
    names = {row["id"]: "source:" + row["id"] for row in PANELS["class_i"]}
    predictions = parse_netmhcpan_to_preds(text, sequence_key_mapping=names)
    kinds = (Kind.pMHC_affinity, Kind.pMHC_presentation)
    expected = {
        (names[row["id"]], offset, row["sequence"][offset:offset + length], allele, kind)
        for row in PANELS["class_i"] for length in (9, 10)
        for offset in range(len(row["sequence"]) - length + 1)
        for allele in capture("netmhcpan")["alleles"] for kind in kinds
    }
    assert len(predictions) == len(expected) == 92
    assert {(p.source_sequence_name, p.offset, p.peptide, p.allele, p.kind)
            for p in predictions} == expected
    anchors = {p.kind: p for p in predictions
               if p.peptide == "AAPVATPAL" and p.allele == "HLA-A*01:01"}
    # Independently read from the native BA and EL columns, which differ.
    affinity, presentation = anchors[Kind.pMHC_affinity], anchors[Kind.pMHC_presentation]
    assert affinity.value == 32539.90
    assert affinity.score == pytest.approx(1 - math.log(32539.90) / math.log(50000))
    assert affinity.percentile_rank == 35.431
    assert value_unit(affinity.kind) == "nM"
    assert presentation.value is None
    assert presentation.score == 0.0023290
    assert presentation.percentile_rank == 11.362


@pytest.mark.parametrize("mode,score,rank,affinity", [
    ("elution_score", 0.013925, 22.77, None),
    ("binding_affinity", 0.181729, 51.30, 6998.92),
])
def test_netmhciipan_column_order_and_source_mapping(mode, score, rank, affinity):
    names = {row["id"]: "source:" + row["id"] for row in PANELS["class_ii"]}
    predictions = parse_netmhciipan43_stdout(
        (ROOT / "outputs" / "netmhciipan.stdout").read_text(),
        sequence_key_mapping=names, mode=mode)
    alleles = [normalize_allele_name(a) for a in capture("netmhciipan")["alleles"]]
    expected = {(names[row["id"]], row["sequence"], allele)
                for row in PANELS["class_ii"] for allele in alleles}
    assert len(predictions) == len(expected) == 110
    assert {(p.source_sequence_name, p.peptide, p.allele) for p in predictions} == expected
    assert all(p.offset == 0 for p in predictions)
    anchor = next(p for p in predictions if p.peptide == "AAKAVKPKVVKPKKA" and p.allele == alleles[0])
    assert (anchor.score, anchor.percentile_rank, anchor.affinity) == (score, rank, affinity)


@pytest.mark.parametrize("label,anchor", [
    ("a0101", (0.060, 26123.53, 50.00)),
    ("b0801", (0.099, 17130.58, 32.00)),
])
def test_netmhccons_real_successful_output(label, anchor):
    name = "netmhccons-" + label
    predictions = parse_netmhccons_stdout((ROOT / "outputs" / (name + ".stdout")).read_text())
    assert len(predictions) == 17
    assert {(p.source_sequence_name, p.offset, p.peptide) for p in predictions} == {
        (row["id"], 0, row["sequence"]) for row in PANELS["cons_9mer"]}
    assert {p.allele for p in predictions} == set(capture(name)["alleles"])
    assert all(p.affinity > 0 and 0 <= p.percentile_rank <= 100 for p in predictions)
    first = next(p for p in predictions if p.peptide == "AAPVATPAL")
    assert (first.score, first.affinity, first.percentile_rank) == anchor


@pytest.mark.parametrize("name", ["mixmhcpred", "mixmhc2pred", "prime"])
def test_gfeller_scores_match_native_columns_for_every_peptide_and_allele(name):
    metadata = capture(name)
    path = ROOT / "outputs" / metadata["raw"]
    alleles = metadata["alleles"]
    if name == "mixmhcpred":
        result = parse_mixmhcpred_output(path, alleles=alleles)
        assert result.version == "3.0"
        predictions = {row.preds[0].peptide: row.preds for row in result.predictions}
        native_alleles = [row.native_allele for row in result.allele_info]
    elif name == "mixmhc2pred":
        native_alleles = metadata["native_alleles"]
        predictions = parse_mixmhc2pred_results(path, alleles, native_alleles)
    else:
        native_alleles = ["A0101", "B0801"]
        predictions = parse_prime_results(path, alleles)
    expected_sequences = {row["sequence"] for row in PANELS[metadata["panel"]]}
    rows = native_table(name)
    assert {r["Peptide"] for r in rows} == set(predictions) == expected_sequences
    assert len(rows) == len(expected_sequences)
    for row in rows:
        preds = predictions[row["Peptide"]]
        assert [p.allele for p in preds] == alleles
        for pred, native in zip(preds, native_alleles):
            assert pred.score == float(row["Score_" + native])
            assert pred.percentile_rank == float(row["%Rank_" + native])
            assert pred.value is None
            assert pred.kind == (Kind.immunogenicity if name == "prime" else Kind.pMHC_presentation)


@pytest.mark.parametrize("name,cls", [("prime", PRIME), ("mixmhc2pred", MixMHC2pred)])
def test_recorded_wrapper_replay_preserves_order_and_duplicates(name, cls, monkeypatch):
    metadata = capture(name)
    panel = [row["sequence"] for row in PANELS[metadata["panel"]]]
    peptides = list(reversed(panel)) + [panel[-1], panel[0]]
    predictor = cls(alleles=metadata["alleles"], program_name="recorded-" + name)
    if name == "prime":
        monkeypatch.setattr(predictor, "_validate_mixmhcpred", lambda: ("recorded-mixmhcpred", "3.0"))

    def replay(args, **kwargs):
        assert Path(args[args.index("-i") + 1]).read_text().splitlines() == peptides
        Path(args[args.index("-o") + 1]).write_bytes((ROOT / "outputs" / metadata["raw"]).read_bytes())

    monkeypatch.setattr(importlib.import_module("mhctools." + name), "run_command", replay)
    results = predictor.predict(peptides)
    assert len(results) == len(peptides)
    for peptide, result in zip(peptides, results):
        assert len(result.preds) == 2
        assert {p.peptide for p in result.preds} == {peptide}
    assert results[0] == results[-2]
    assert results[-1] == results[-3]


@pytest.mark.parametrize("model", ["cterm", "20s"])
def test_netchop_recorded_positions_and_wrapper_batch(model, monkeypatch):
    sequences = [row["sequence"] for row in PANELS["cleavage"]]
    raw = (ROOT / "outputs" / ("netchop-" + model + ".stdout")).read_bytes()
    scores = NetChop.parse_netchop(raw)
    assert len(scores) == len(sequences) == 77
    assert [len(row) for row in scores] == [len(s) for s in sequences]
    assert sum(map(len, scores)) == 1705
    assert all(0 <= score <= 1 for row in scores for score in row)
    assert scores[0][:2] == ({
        "cterm": [0.262091, 0.230793], "20s": [0.081664, 0.315131]}[model])

    def replay(args, **kwargs):
        fasta = Path(args[-1]).read_text().splitlines()
        assert fasta[1::2] == sequences
        return subprocess.CompletedProcess(args, 0, raw, b"")

    module = importlib.import_module("mhctools.netchop")
    monkeypatch.setattr(module, "resolve_netchop_dir", lambda **kwargs: None)
    monkeypatch.setattr(module.subprocess, "run", replay)
    predictor = NetChop(program_name=sys.executable, execution="native")
    assert predictor.cleavage_probs_many(sequences + [sequences[0]]) == dict(zip(sequences, scores))
