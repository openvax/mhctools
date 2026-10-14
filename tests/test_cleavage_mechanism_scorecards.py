"""Check coordinate, endpoint and missing-data semantics of the frozen report."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "analyses/osteosarc_vaccine_cleavage/mechanism_scorecards.py"
SPEC = importlib.util.spec_from_file_location("mechanism_scorecards", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
with pytest.MonkeyPatch.context() as patch:
    patch.syspath_prepend(str(SCRIPT.parent))
    SPEC.loader.exec_module(REPORT)
SOURCE = ROOT / "analyses/osteosarc_vaccine_cleavage/results/2026-09-18T175125-855754-0400"


@pytest.fixture(scope="module")
def frozen_report():
    records = [row for row in REPORT.read_csv(SOURCE / "tables/vaccine_sequence_inventory.csv")
               if row["sequence_type"] == "synthetic_long_peptide"]
    quantitative = REPORT.read_csv(SOURCE / "tables/slp_quantitative_bond_scores.csv")
    motifs = REPORT.read_csv(SOURCE / "tables/slp_motif_assessments.csv")
    catalog = REPORT.mechanism_catalog(REPORT.read_csv(SOURCE / "tables/model_catalog.csv"), "human")
    rows = REPORT.enriched_assessments(records, quantitative, motifs, catalog)
    return records, quantitative, motifs, catalog, rows


@pytest.mark.parametrize("organism,expected", [
    ("human", REPORT.HUMAN), ("nonhuman", REPORT.MAMMAL),
    ("mixed", REPORT.MAMMAL), ("unknown", REPORT.MAMMAL),
])
def test_organism_preference_keeps_both_recorded_models(organism, expected):
    catalog = REPORT.mechanism_catalog(REPORT.read_csv(SOURCE / "tables/model_catalog.csv"), organism)
    choices = {row["model"]: row["pepsickle_preference"] for row in catalog
               if row["model"] in (REPORT.HUMAN, REPORT.MAMMAL)}
    assert choices[expected] == "preferred"
    assert set(choices.values()) == {"preferred", "comparison"}


def test_all_native_scores_states_and_coordinates_survive_export(frozen_report, tmp_path):
    records, quantitative, motifs, catalog, rows = frozen_report
    assert len(records) == 40
    assert len(rows) == len(quantitative) + len(motifs) == 7248
    path = tmp_path / "scores.csv"
    REPORT.write_csv(path, rows)
    exported = REPORT.read_csv(path)
    lookup = REPORT.site_lookup(exported)
    for original in quantitative:
        bond = str(int(float(original["bond"]))) if original["bond"] else ""
        row = lookup[(original["sequence_record_id"], original["model"], bond)]
        assert row["score"] == original["score"]
        assert row["score_units"] == original["score_units"]
        assert row["unsupported_reason"] == original["unsupported_reason"]
    for original in motifs:
        bond = str(int(float(original["bond"]))) if original["bond"] else ""
        row = lookup[(original["sequence_record_id"], original["model"], bond)]
        assert row["assessment_state"] == original["status"]
        assert row["score"] == row["score_units"] == ""
        assert row["motif_strictness"] == original["motif_strictness"]
    map2 = "MAP2-chr2-209694768:vaccine-peptide-4"
    assert lookup[(map2, "netchop-3.1-cterm-3.0", "13")]["score"] == "0.728793"
    assert lookup[(map2, "netchop-3.1-20s-3.0", "13")]["score"] == "0.039743"
    assert lookup[(map2, "netchop-3.1-20s-3.0", "13")]["core_bond"] == "5"
    dpp4 = lookup[(map2, "dpp4-qpisa", "2")]
    assert dpp4["score"] == "0.4856"
    assert dpp4["inside_source_minimal_epitope"] == "False"
    assert dpp4["topology"] == "n_terminal"
    metadata = {row["model"]: row for row in catalog}
    assert "ligand" in metadata["netchop-3.1-cterm-3.0"]["training_or_assay"]
    assert "in-vitro" in metadata["netchop-3.1-20s-3.0"]["training_or_assay"]
    assert "type agnostic" in metadata[REPORT.HUMAN]["enzyme_attribution"]


def test_unassessed_zero_and_motif_nonmatch_remain_distinct(frozen_report):
    rows = frozen_report[-1]
    unsupported = next(row for row in rows if row["assessment_state"] == "unassessed" and not row["bond"])
    assert REPORT.cell_text(unsupported) == "NA"
    assert REPORT.cell_text(None) == "NA"
    assert REPORT.cell_text({"assessment_state": "scored", "score": "0"}) == "0.000"
    assert REPORT.cell_text({"assessment_state": "matched"}) == "Match"
    assert REPORT.cell_text({"assessment_state": "not_matched"}) == "No match"
    assert REPORT.cell_text({"model": "dpp4-qpisa", "assessment_state": "scored", "score": ".686"}) == "38% predicted loss"
    assert REPORT.cell_text({"model": "dpp4-qpisa", "assessment_state": "scored", "score": "-.5799"}) == "No predicted loss"
    lookup = REPORT.site_lookup(rows)
    assert REPORT.assessment_at(lookup, unsupported["sequence_record_id"], unsupported["model"], 13) == unsupported


def test_core_boundary_bonds_are_not_internal_cuts(frozen_report):
    records, quantitative, _, catalog, _ = frozen_report
    parent = next(row for row in records if row["gene"] == "MAP2")
    selected = [row for row in quantitative if row["sequence_record_id"] == parent["sequence_record_id"]
                and row["model"] == "netchop-3.1-cterm-3.0"]
    rows = REPORT.enriched_assessments([parent], selected, [], catalog)
    states = {row["bond"]: row["inside_source_minimal_epitope"] for row in rows}
    assert states["8"] == states["17"] == "False"
    assert all(states[str(bond)] == "True" for bond in range(9, 17))
    bad = deepcopy(parent)
    bad["minimal_epitope_offset"] = "7"
    with pytest.raises(ValueError, match="Core annotation"):
        REPORT.core_bounds(bad)


def test_terminal_scores_cannot_be_projected_to_embedded_core(frozen_report):
    records, quantitative, _, catalog, _ = frozen_report
    parent = next(row for row in records if row["gene"] == "MAP2")
    row = deepcopy(next(row for row in quantitative if row["sequence_record_id"] == parent["sequence_record_id"]
                        and row["model"] == "dpp4-qpisa"))
    row.update(bond="13", left_residue="F", right_residue="T")
    with pytest.raises(ValueError, match="Terminal model assessed"):
        REPORT.enriched_assessments([parent], [row], [], catalog)


def test_unknown_models_and_duplicate_assessments_are_rejected(frozen_report):
    source_catalog = REPORT.read_csv(SOURCE / "tables/model_catalog.csv")
    bad = deepcopy(source_catalog[0])
    bad["model"] = "future-model-without-mechanism"
    with pytest.raises(ValueError, match="Unannotated source model"):
        REPORT.mechanism_catalog([bad], "human")
    records, quantitative, _, catalog, _ = frozen_report
    row = quantitative[0]
    with pytest.raises(ValueError, match="Duplicate source assessment"):
        REPORT.enriched_assessments(records, [row, row], [], catalog)


def test_complete_source_manifest_is_verified_and_paths_cannot_escape(tmp_path):
    assert len(REPORT.verified_source(SOURCE)) == 107
    payload = tmp_path / "payload.txt"
    payload.write_text("original")
    manifest = tmp_path / "SHA256SUMS.json"
    manifest.write_text(json.dumps({"payload.txt": REPORT.sha256(payload)}))
    payload.write_text("altered")
    with pytest.raises(ValueError, match="checksum verification"):
        REPORT.verified_source(tmp_path)
    manifest.write_text(json.dumps({"../escape.txt": "irrelevant"}))
    with pytest.raises(ValueError, match="checksum verification"):
        REPORT.verified_source(tmp_path)


def test_scorecard_pages_cover_every_internal_bond_exactly_once(frozen_report, tmp_path):
    pytest.importorskip("reportlab")
    records, _, _, catalog, rows = frozen_report
    pdf = tmp_path / "scorecards.pdf"
    index = REPORT.render_pdf(pdf, records, rows, catalog, "human", SOURCE.name)
    assert pdf.read_bytes().startswith(b"%PDF-")
    for record in records:
        segments = [row for row in index if row["sequence_record_id"] == record["sequence_record_id"]]
        covered = [bond for row in segments for bond in range(row["first_bond"], row["last_bond"] + 1)]
        assert covered == list(range(1, len(record["sequence"])))
    assert len({row["pdf_page"] for row in index}) == len(index)
