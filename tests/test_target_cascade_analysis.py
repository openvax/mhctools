"""Synthetic source mapping and current-fragment recognition checks."""

import importlib.util
from pathlib import Path


spec = importlib.util.spec_from_file_location(
    "target_cascade", Path(__file__).parents[1] / "analyses" /
    "osteosarc_vaccine_cleavage" / "target_cascade.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_exact_source_highlight_mapping_and_unknown_mutations():
    record = dict(sequence="DEFGHIK", minimal_epitope="FGHI", minimal_epitope_offset=2)
    variant = dict(aa_sequence="ACDEFGHIKLM", aa_hi_start=5, aa_hi_end=6)
    assert module.mutant_region(record, variant)[0] == (3, 4)
    assert module.mutant_region(record, {})[0] is None
    # Altered vaccine flanks can be mapped only inside a unique exact core.
    assert module.mutant_region(dict(record, sequence="QQFGHIR"), variant)[0] == (3, 4)
    assert module.mutant_region(record, dict(variant, aa_hi_start=0, aa_hi_end=1))[0] is None


def test_no_alignment_guessed_for_repeated_or_absent_anchors():
    record = dict(sequence="DEFGHIK", minimal_epitope="FGHI", minimal_epitope_offset=2)
    assert module.mutant_region(record, dict(aa_sequence="FGHIACFGHI", aa_hi_start=1, aa_hi_end=2))[0] is None


def test_new_n_terminus_changes_dpp4_evidence():
    parent = "RPTLSACDEFGHIK"
    # Removing the first dipeptide exposes TLS, which must be rescored.
    old = next(r for r in module.extracellular_evidence(parent) if r.model.name == "dpp4-qpisa")
    new = next(r for r in module.extracellular_evidence(parent[2:]) if r.model.name == "dpp4-qpisa")
    assert old.peptide.sequence != new.peptide.sequence
    assert old.sites[0].score != new.sites[0].score
    assert all(w > 0 for w in module.recognition_weights(parent, 10))


def test_class_ii_tracks_core_and_does_not_guess_mutant_status():
    record = dict(sequence_record_id="s1", variant_id="v1", gene="test", sequence="ACDEFGHIKLMNPQR",
                  minimal_epitope="", minimal_epitope_offset="")
    variant = dict(aa_sequence=record["sequence"], aa_hi_start=6, aa_hi_end=7, peptides={})
    ligand = dict(sequence_record_id="s1", gene="test", mhc_class="II", allele="test allele",
                  start="1", end="15", peptide=record["sequence"], binding_core="EFGHIKLMN",
                  percentile_rank="1", rank_threshold="5", predictor="test", model_version="test")
    targets, candidates, _ = module.target_inventory([record], {"v1": variant}, [ligand])
    assert len(targets) == 1
    assert (targets[0]["start"], targets[0]["end"]) == (3, 12)
    assert targets[0]["target_sequence"] == "EFGHIKLMN"
    assert candidates[0]["mutant_overlap"] is True
    targets, candidates, _ = module.target_inventory([record], {"v1": {}}, [ligand])
    assert targets == []
    assert candidates[0]["mutant_overlap"] is None


def test_multiple_disclosed_source_targets_are_retained():
    record = dict(sequence_record_id="s1", variant_id="v1", gene="test", sequence="ACDEFGHIKLMNPQR",
                  minimal_epitope="EFGHIKLMN", minimal_epitope_offset=3)
    targets, _, _ = module.target_inventory([record], {"v1": {"peptides": {"CeGaT_Class_I": "GHIKLMNPQ"}}}, [])
    assert len(targets) == 2
    assert targets[0]["kind"] == "source minimal epitope"
    assert targets[1]["kind"] == "source class-I candidate"


def test_class_ii_mutant_flank_is_preserved_with_its_binding_core():
    record = dict(sequence_record_id="s1", variant_id="v1", gene="test", sequence="ACDEFGHIKLMNPQR",
                  minimal_epitope="", minimal_epitope_offset="")
    variant = dict(aa_sequence=record["sequence"], aa_hi_start=1, aa_hi_end=2, peptides={})
    ligand = dict(sequence_record_id="s1", gene="test", mhc_class="II", allele="test allele",
                  start="1", end="15", peptide=record["sequence"], binding_core="EFGHIKLMN",
                  percentile_rank="1", rank_threshold="5", predictor="test", model_version="test")
    targets, candidates, _ = module.target_inventory([record], {"v1": variant}, [ligand])
    assert len(targets) == 1
    assert (targets[0]["start"], targets[0]["end"]) == (1, 12)
    assert targets[0]["binding_core_start"] == 3
    assert targets[0]["kind"] == "predicted class-II core + mutant flank"
    assert candidates[0]["mutant_overlap"] is True
    assert candidates[0]["mutant_binding_core_overlap"] is False
