"""Synthetic exact-bond oracles; native evidence must not become serum rates."""

import pytest

from mhctools import (
    CleavageInput, CleavageModel, CleavageResult, CleavageSite,
    DegradationTarget, annotate_target_cleavage, predict_cleavage,
)


def test_new_termini_and_original_coordinates_distinguish_loss_from_release():
    peptide = CleavageInput("APACDEFG", source_id="synthetic", source_start=4)
    model = CleavageModel(
        "synthetic", "1", "synthetic enzyme", "", "synthetic", ("extracellular",),
        "substrate_reference", ("https://example.org/assay",), "synthetic assay", "no kinetics")
    results = [CleavageResult(peptide, model, tuple(
        CleavageSite(b, "reported", "Synthetic site oracle") for b in (1, 2, 3, 6, 7)))]
    rows = annotate_target_cleavage(DegradationTarget("exact", 6, 10), results)
    assert [r["source_bond"] for r in rows] == [5, 6, 7, 10, 11]
    assert [r["target_effect"] for r in rows] == [
        "flank_trim", "boundary_release", "target_split", "boundary_release", "flank_trim"]
    assert rows[2]["bond_label"] == "A7|C8"
    assert rows[2]["score"] is None
    assert rows[2]["evidence"] == "substrate_reference"


def test_native_dpp4_score_and_motif_have_separate_semantics():
    results = predict_cleavage("APACDEFG", models=("dpp4-qpisa", "fap-dipeptidyl"))
    rows = annotate_target_cleavage(DegradationTarget("target", 2, 6), results)
    scored = next(r for r in rows if r["model"] == "dpp4-qpisa")
    motif = next(r for r in rows if r["model"] == "fap-dipeptidyl")
    assert scored["score"] == results[0].sites[0].score
    assert scored["scored_endpoint"] == "substrate_depletion"
    assert scored["score_units"] == results[0].model.score_units
    assert scored["target_effect"] == motif["target_effect"] == "boundary_release"
    assert motif["score"] is None and motif["status"] == "matched"
    assert "probability" not in scored and "rate" not in scored


def test_missing_scope_and_nonmatch_do_not_become_protection():
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
    results = predict_cleavage(sequence, models=("mme-hydrophobic", "fap-dipeptidyl"))
    rows = annotate_target_cleavage(DegradationTarget("target", 5, 14), results)
    unavailable = next(r for r in rows if r["model"] == "mme-hydrophobic")
    nonmatch = next(r for r in rows if r["model"] == "fap-dipeptidyl")
    assert unavailable["status"] == "unavailable"
    assert unavailable["source_bond"] is None and unavailable["target_effect"] is None
    assert "30" in unavailable["reason"]
    assert nonmatch["status"] == "not_matched"
    assert nonmatch["target_effect"] == "flank_trim"


def test_different_chemistry_occurrences_and_invalid_targets_are_rejected():
    a = CleavageInput("APACDEFG", source_id="one")
    b = CleavageInput("APACDEFG", source_id="two")
    results = [predict_cleavage(x, models=("fap-dipeptidyl",))[0] for x in (a, b)]
    with pytest.raises(ValueError, match="Do not mix"):
        annotate_target_cleavage(DegradationTarget("target", 2, 6), results)
    with pytest.raises(ValueError, match="integer"):
        annotate_target_cleavage(DegradationTarget("bad", True, 6), results[:1])
    with pytest.raises(ValueError):
        annotate_target_cleavage(DegradationTarget("absent", 8, 12), results[:1])
    with pytest.raises(ValueError, match="At least one"):
        annotate_target_cleavage(DegradationTarget("target", 2, 6), [])
