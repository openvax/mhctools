"""Window identity and source-core topology in additional native predictions."""

import hashlib
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from mhctools import CleavageInput, ITCellCleavage
from mhctools.cleavenet import CleaveNetScore, ENZYMES

SCRIPT = Path(__file__).parents[1] / "analyses/osteosarc_vaccine_cleavage/stability_route.py"
SPEC = importlib.util.spec_from_file_location("stability_route", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
with pytest.MonkeyPatch.context() as patch:
    patch.syspath_prepend(str(SCRIPT.parent))
    SPEC.loader.exec_module(REPORT)


def record(id_, sequence="ACDEFGHIKLMN"):
    return dict(sequence_record_id=id_, gene="fixture", sequence=sequence, length=str(len(sequence)),
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                minimal_epitope="FGHIK", minimal_epitope_offset="4")


def test_overlapping_windows_preserve_duplicate_occurrences_and_do_not_invent_cuts():
    records = [record("first"), record("second")]
    inputs = REPORT.window_inputs(records)
    assert len(inputs) == 6
    assert len({item.occurrence_id for item in inputs}) == 6
    assert [item.sequence for item in inputs[:3]] == ["ACDEFGHIKL", "CDEFGHIKLM", "DEFGHIKLMN"]
    assert [item.source_start for item in inputs] == [0, 1, 2, 0, 1, 2]
    scores = tuple(CleaveNetScore(enzyme, -0.125, 0.345) for enzyme in ENZYMES)
    native = [SimpleNamespace(peptide_input=item, scores=scores, cache_key="fixture") for item in inputs]
    rows = REPORT.window_rows(records, inputs, native)
    assert len(rows) == 108
    assert {r["source_core_overlap_residues"] for r in rows} == {5}
    assert {r["z_score"] for r in rows} == {-0.125}
    assert {r["ensemble_sd"] for r in rows} == {0.345}
    assert all("bond" not in r for r in rows)
    with pytest.raises(ValueError, match="Missing"):
        REPORT.window_rows(records, inputs, native[:-1])
    native[0].peptide_input = inputs[1]
    with pytest.raises(ValueError, match="identity"):
        REPORT.window_rows(records, inputs, native)


def test_native_bonds_preserve_core_boundaries_and_terminal_initial_step():
    parent = record("first")
    peptide = CleavageInput(parent["sequence"], source_id="first")
    rows = REPORT.bond_rows([parent], [ITCellCleavage("S").predict(peptide)])
    lookup = {r["bond"]: r for r in rows}
    assert lookup[4]["inside_source_minimal_epitope"] is False
    assert lookup[9]["inside_source_minimal_epitope"] is False
    assert all(lookup[b]["inside_source_minimal_epitope"] for b in range(5, 9))
    assert lookup[5]["left_aa"] == "F"
    assert lookup[5]["right_aa"] == "G"
    h = REPORT.bond_rows([parent], [ITCellCleavage("H").predict(peptide)])
    assert [row["bond"] for row in h] == [1]
    assert h[0]["inside_source_minimal_epitope"] is False
    unannotated = dict(parent, minimal_epitope="", minimal_epitope_offset="")
    rows = REPORT.bond_rows([unannotated], [ITCellCleavage("S").predict(peptide)])
    assert {r["inside_source_minimal_epitope"] for r in rows} == {""}


def test_frozen_vaccine_inventory_scope_and_all_window_counts():
    source = SCRIPT.parent / "results/2026-10-06T142010-273003-0400"
    REPORT.verified_source(source)
    records = REPORT.vaccine_records(source)
    assert len(records) == 40
    assert len({r["sequence"] for r in records}) == 39
    assert len(REPORT.window_inputs(records)) == 538
    assert max(len(r["sequence"]) for r in records) == 80
