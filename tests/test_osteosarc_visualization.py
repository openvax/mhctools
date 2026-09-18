from datetime import datetime, timezone
import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pytest


pytest.importorskip("matplotlib")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
pytest.importorskip("scipy")


SCRIPT_PATH = (
    Path(__file__).parents[1] / "analyses" / "osteosarc_vaccine_cleavage" / "analyze.py"
)
SPEC = spec_from_file_location("osteosarc_vaccine_cleavage_analysis", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
ANALYSIS = module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def test_sequence_segments_cover_every_internal_bond_once():
    for length in (2, 17, 28, 80):
        segments = ANALYSIS.sequence_segments("A" * length)
        observed = [
            bond
            for start_bond, end_bond in segments
            for bond in range(start_bond, end_bond + 1)
        ]
        assert observed == list(range(1, length))


def test_timestamped_output_directory_is_unique_and_timezone_explicit(tmp_path):
    generated_at = datetime(2026, 9, 16, 3, 45, 12, 123456, tzinfo=timezone.utc)
    assert ANALYSIS.timestamped_output_dir(tmp_path, generated_at) == (
        tmp_path / "2026-09-16T034512-123456+0000"
    )
    with pytest.raises(ValueError, match="timezone"):
        ANALYSIS.timestamped_output_dir(tmp_path, generated_at.replace(tzinfo=None))


def test_segment_ligands_use_half_open_boundaries_without_duplicates():
    ligand_df = pd.DataFrame(
        [
            {
                "sequence_record_id": "record",
                "mhc_class": "I",
                "display_candidate": True,
                "start": 24,
                "end": 32,
                "percentile_rank": 0.2,
                "allele": "HLA-A*01:01",
                "overlaps_disclosed_minimal_epitope": False,
            },
            {
                "sequence_record_id": "record",
                "mhc_class": "I",
                "display_candidate": True,
                "start": 24,
                "end": 32,
                "percentile_rank": 0.4,
                "allele": "HLA-A*01:01",
                "overlaps_disclosed_minimal_epitope": False,
            },
        ]
    )
    first = ANALYSIS.segment_ligand_candidates(
        ligand_df, "record", "I", residue_start=1, residue_end=28
    )
    second = ANALYSIS.segment_ligand_candidates(
        ligand_df, "record", "I", residue_start=28, residue_end=55
    )
    assert first.empty
    assert len(second) == 1
    assert second.iloc[0]["percentile_rank"] == 0.2


def test_ligand_display_prioritizes_intended_window_then_allele_diversity():
    candidates = pd.DataFrame(
        [
            {
                "allele": "HLA-A*01:01",
                "start": 1,
                "end": 9,
                "percentile_rank": 0.01,
                "overlaps_disclosed_minimal_epitope": False,
            },
            {
                "allele": "HLA-A*01:01",
                "start": 12,
                "end": 20,
                "percentile_rank": 0.02,
                "overlaps_disclosed_minimal_epitope": False,
            },
            {
                "allele": "HLA-B*08:01",
                "start": 2,
                "end": 10,
                "percentile_rank": 0.20,
                "overlaps_disclosed_minimal_epitope": False,
            },
            {
                "allele": "HLA-C*07:01",
                "start": 3,
                "end": 11,
                "percentile_rank": 0.30,
                "overlaps_disclosed_minimal_epitope": True,
            },
            {
                "allele": "HLA-C*01:02",
                "start": 13,
                "end": 21,
                "percentile_rank": 0.40,
                "overlaps_disclosed_minimal_epitope": False,
            },
        ]
    )
    selected = ANALYSIS.select_ligand_candidates_for_display(candidates)
    assert len(selected) == 5
    assert selected.iloc[0]["allele"] == "HLA-C*07:01"
    assert selected.iloc[0]["selection_reason"] == "intended_epitope"
    assert selected["allele"].nunique() >= 3
    assert selected[["start", "end"]].duplicated().sum() == 0
    assert set(selected["display_lane"]) <= {0, 1, 2, 3, 4}


def test_ligand_display_never_draws_same_span_twice_across_alleles():
    candidates = pd.DataFrame(
        [
            {
                "allele": allele,
                "start": 4,
                "end": 12,
                "percentile_rank": rank,
                "overlaps_disclosed_minimal_epitope": False,
            }
            for allele, rank in (
                ("HLA-A*01:01", 0.1),
                ("HLA-B*08:01", 0.2),
                ("HLA-C*07:01", 0.3),
            )
        ]
    )
    selected = ANALYSIS.select_ligand_candidates_for_display(candidates)
    assert len(selected) == 1
    assert selected.iloc[0]["allele"] == "HLA-A*01:01"


def test_manuscript_caption_defines_visuals_and_selection_audit(tmp_path):
    selection = pd.DataFrame(
        [
            {
                "panel": "A",
                "gene": "GENE1",
                "protein_change": "p.Arg1Gly",
                "selection_reason": "most MHC-I candidates",
                "intended_epitope_internal_conservative_cuts": 2,
                "intended_epitope_overlapping_mhc_candidates": 7,
                "mhc_i_candidate_count": 4,
                "mhc_ii_candidate_count": 3,
            }
        ]
    )
    path = tmp_path / "caption.md"
    ANALYSIS.write_manuscript_caption(path, selection)
    text = path.read_text()
    assert "**A, GENE1 p.Arg1Gly.**" in text
    assert "Most MHC-I candidates" in text
    assert "not probabilities" in text
    assert "agreement does not establish in-vivo degradation" in text


def test_score_matrix_keeps_zero_distinct_from_unassessed():
    models = ANALYSIS.QUANTITATIVE_SITE_MODELS[:2]
    scores = pd.DataFrame(
        [
            {
                "sequence_record_id": "record",
                "model": models[0],
                "bond": 1,
                "score": 0.0,
                "assessable": True,
            },
            {
                "sequence_record_id": "record",
                "model": models[0],
                "bond": 2,
                "score": 0.8,
                "assessable": False,
            },
        ]
    )
    matrix = ANALYSIS._score_matrix("record", scores, models, [1, 2])
    assert matrix[0, 0] == 0.0
    assert np.isnan(matrix[0, 1])
    assert np.isnan(matrix[1]).all()


def test_continuous_score_runs_preserve_values_and_break_at_missing_bonds():
    runs = ANALYSIS.continuous_score_runs(
        [1, 2, 3, 4, 6], np.array([0.0, 0.25, np.nan, 0.75, 1.0])
    )
    assert len(runs) == 3
    assert runs[0][0].tolist() == [1.5, 2.5]
    assert runs[0][1].tolist() == [0.0, 0.25]
    assert runs[1][0].tolist() == [4.5]
    assert runs[1][1].tolist() == [0.75]
    assert runs[2][0].tolist() == [6.5]
    assert runs[2][1].tolist() == [1.0]


def test_standalone_map_stem_is_stable_safe_and_segments_long_sequences():
    stem = ANALYSIS.standalone_map_stem(7, "MT_ND5-chrM-12994:vaccine-peptide-1", 2, 3)
    assert stem == "07-mt-nd5-chrm-12994-vaccine-peptide-1-segment-2-of-3"


def test_figure_model_set_keeps_distinct_20s_but_omits_redundant_all_mammal():
    assert ANALYSIS.FIGURE_CYTOSOL_MODELS == [
        "netchop-3.1-20s-3.0",
        "netcleave-i-hla",
        "pepsickle-in-vivo-human-only",
        "netchop-3.1-cterm-3.0",
    ]
    assert "pepsickle-in-vivo-all-mammal" not in ANALYSIS.FIGURE_QUANTITATIVE_MODELS


def test_display_uses_top_ten_cap_and_opinionated_slp_enzyme_tracks():
    assert ANALYSIS.MHC_DISPLAY_MAX_WINDOWS == 10
    assert ANALYSIS.MHC_DISPLAY_LANES == 5
    displayed = {
        model
        for _, models, _ in ANALYSIS.FIGURE_SLP_ENZYME_TRACKS
        for model in models
    }
    assert displayed == {
        "mme-hydrophobic",
        "fap-dipeptidyl",
        "fap-endo-gp",
        "anpep-ala",
        "enpep-acidic",
    }
    assert not displayed.intersection({
        "ace-dipeptidyl", "cpb2-basic", "cpn-basic", "app1-xp", "erap2-basic"
    })


def test_red_cut_requires_three_hits_and_all_four_display_tracks_assessed():
    ligand = pd.DataFrame(
        [{"sequence_record_id": "record", "mhc_class": "I", "start": 1, "end": 3}]
    )
    rows = [
        {
            "sequence_record_id": "record",
            "model": model,
            "bond": 2,
            "score": 0.8,
            "assessable": True,
        }
        for model in ANALYSIS.FIGURE_CYTOSOL_MODELS[:3]
    ]
    incomplete = ANALYSIS.annotate_ligand_cleavage_exposure(ligand, pd.DataFrame(rows))
    assert incomplete.iloc[0]["internal_candidate_cleavage_bonds"] == ""

    rows.append(
        {
            "sequence_record_id": "record",
            "model": ANALYSIS.FIGURE_CYTOSOL_MODELS[2],
            "bond": 2,
            "score": 0.2,
            "assessable": True,
        }
    )
    complete = ANALYSIS.annotate_ligand_cleavage_exposure(ligand, pd.DataFrame(rows))
    assert complete.iloc[0]["internal_candidate_cleavage_bonds"] == "2(3/4)"


def test_conservative_cut_support_counts_models_not_duplicate_rows():
    rows = [
        {
            "sequence_record_id": "record",
            "model": model,
            "bond": 2,
            "score": 0.8 if index < 3 else 0.2,
            "assessable": True,
        }
        for index, model in enumerate(ANALYSIS.FIGURE_CYTOSOL_MODELS)
    ]
    rows.append(dict(rows[0]))
    support = ANALYSIS.conservative_cut_support(
        pd.DataFrame(rows),
        "record",
        [1, 2],
        ANALYSIS.FIGURE_CYTOSOL_MODELS,
        ANALYSIS.FIGURE_RED_SUPPORT_REQUIRED,
    )
    assert support == {2: (3, 4)}


def test_merge_residue_spans_forms_a_binary_coverage_overlay():
    assert ANALYSIS.merge_residue_spans([(8, 12), (1, 3), (3, 7), (15, 16)]) == [
        (1, 12),
        (15, 16),
    ]
    with pytest.raises(ValueError, match="invalid residue span"):
        ANALYSIS.merge_residue_spans([(4, 3)])


def test_motif_matrix_keeps_no_match_distinct_from_unsupported():
    models = ANALYSIS.EXTRACELLULAR_MOTIF_MODELS[:2]
    motifs = pd.DataFrame(
        [
            {
                "sequence_record_id": "record",
                "model": models[0],
                "bond": 1,
                "status": "matched",
            },
            {
                "sequence_record_id": "record",
                "model": models[0],
                "bond": 2,
                "status": "not_matched",
            },
            {
                "sequence_record_id": "record",
                "model": models[1],
                "bond": np.nan,
                "status": "unsupported",
            },
        ]
    )
    matrix = ANALYSIS._motif_matrix("record", motifs, models, [1, 2])
    assert matrix[0].tolist() == [1.0, 0.0]
    assert np.isnan(matrix[1]).all()


def test_minimal_epitope_bounds_are_one_based_and_validated():
    record = pd.Series(
        {
            "sequence_record_id": "BMP1",
            "sequence": "RISVTPGEKIILNFTTLDLYRSR",
            "minimal_epitope": "KIILNFTTL",
            "minimal_epitope_offset": 8,
        }
    )
    assert ANALYSIS._minimal_epitope_bounds(record) == (9, 17)
    record["minimal_epitope_offset"] = 7
    with pytest.raises(ValueError, match="offset does not match"):
        ANALYSIS._minimal_epitope_bounds(record)


def test_netmhciipan_parser_uses_el_score_and_rank_columns():
    stdout = """
       1 DRB1_0301 RISVTPGEKIILNFT 1 ISVTPGEKI 0.740 0 Sequence 0.079862 9.55 0.000 0.5 10.0 500.0 <=WB
    """
    assert ANALYSIS._parse_netmhciipan_rows(stdout) == [
        {
            "allele": "DRB1_0301",
            "peptide": "RISVTPGEKIILNFT",
            "binding_core": "ISVTPGEKI",
            "score": 0.079862,
            "percentile_rank": 9.55,
        }
    ]


def test_netmhciipan_cli_alleles_do_not_expose_shell_globs():
    assert ANALYSIS.netmhciipan_cli_allele("HLA-DRB1*03:01") == "DRB1_0301"
    assert (
        ANALYSIS.netmhciipan_cli_allele("HLA-DPA1*01:03-DPB1*04:01")
        == "HLA-DPA10103-DPB10401"
    )


def test_vulnerable_bonds_require_context_separated_support():
    record_id = "record"
    records = pd.DataFrame(
        [
            {
                "sequence_record_id": record_id,
                "gene": "GENE",
                "sequence": "ABCDEFGHI",
                "minimal_epitope": "CDEFG",
                "minimal_epitope_offset": 2,
            }
        ]
    )
    quantitative = pd.DataFrame(
        [
            {
                "sequence_record_id": record_id,
                "model": model,
                "bond": 4,
                "score": 0.8,
                "assessable": True,
            }
            for model in ANALYSIS.MHC_I_CLEAVAGE_MODELS[:3]
        ]
    )
    motifs = pd.DataFrame(
        [
            {
                "sequence_record_id": record_id,
                "model": model,
                "bond": 6,
                "status": "matched",
            }
            for model in ANALYSIS.EXTRACELLULAR_MOTIF_MODELS[:2]
        ]
    )
    observed = ANALYSIS.vulnerable_bond_table(records, quantitative, motifs)
    assert observed["biological_context"].tolist() == [
        "cytosolic/proteasome and MHC-I processing",
        "extracellular/plasma recognition motifs",
    ]
    assert observed["support_count"].tolist() == [3, 2]
    assert observed["in_disclosed_minimal_epitope"].tolist() == [True, True]


def _load_rerender():
    """Import rerender.py, which imports ``analyze`` as a sibling script."""
    import sys

    directory = str(SCRIPT_PATH.parent)
    inserted = directory not in sys.path
    if inserted:
        sys.path.insert(0, directory)
    try:
        spec = spec_from_file_location(
            "osteosarc_vaccine_cleavage_rerender", SCRIPT_PATH.with_name("rerender.py")
        )
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if inserted:
            sys.path.remove(directory)


def _frozen_run(directory, contents):
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, text in contents.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        manifest[name] = ANALYSIS.sha256_file(path)
    (directory / "SHA256SUMS.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def test_rerender_reuses_a_frozen_run_only_when_every_checksum_matches(tmp_path):
    rerender = _load_rerender()
    source = tmp_path / "run"
    manifest = _frozen_run(
        source, {"tables/scores.csv": "bond,score\n1,0.5\n", "REPORT.md": "# report\n"}
    )
    assert rerender._verified_source_manifest(source) == manifest

    # A table edited after the run was frozen must not be silently re-rendered
    # as though it were the checksummed prediction output.
    (source / "tables" / "scores.csv").write_text(
        "bond,score\n1,0.9\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="tables/scores.csv"):
        rerender._verified_source_manifest(source)


def test_rerender_reports_a_missing_frozen_file_instead_of_skipping_it(tmp_path):
    rerender = _load_rerender()
    source = tmp_path / "run"
    _frozen_run(source, {"tables/scores.csv": "bond,score\n1,0.5\n"})
    (source / "tables" / "scores.csv").unlink()
    with pytest.raises(RuntimeError, match="tables/scores.csv"):
        rerender._verified_source_manifest(source)


def test_netcleave_revision_comes_from_the_managed_snapshot_record(tmp_path):
    """A fetched snapshot has no .git, so its pinned revision must be read
    from the artifact record rather than from `git log`."""
    revision = "bc90bf8490dcc5bc7628b1f0c06f860a50e6f338"
    snapshot = tmp_path / revision
    snapshot.mkdir()
    (snapshot / ".mhctools-artifact.json").write_text(
        json.dumps({"name": "netcleave", "revision": revision}), encoding="utf-8"
    )
    directory, resolved = ANALYSIS.resolve_netcleave_dir(snapshot)
    assert directory == snapshot.resolve()
    assert resolved == revision


def test_netcleave_directory_must_exist(tmp_path):
    with pytest.raises(SystemExit, match="does not exist"):
        ANALYSIS.resolve_netcleave_dir(tmp_path / "absent")


def test_netcleave_without_git_or_artifact_record_is_refused(tmp_path):
    """Recording an unidentifiable NetCleave install would silently break the
    run's model provenance, so the run stops instead."""
    plain = tmp_path / "unidentified"
    plain.mkdir()
    with pytest.raises(SystemExit, match="neither a git checkout"):
        ANALYSIS.resolve_netcleave_dir(plain)


def test_a_user_managed_netcleave_install_is_accepted(monkeypatch, tmp_path):
    """Any ready NetCleave counts, not only an mhctools-managed snapshot.

    Requiring manager == "mhctools" rejected a perfectly usable NETCLEAVE_DIR
    with the self-contradictory "NetCleave is not available ... Status: ready",
    aborting the run even though the wrapper resolves that same checkout.
    """
    revision = "bc90bf8490dcc5bc7628b1f0c06f860a50e6f338"
    checkout = tmp_path / "user-netcleave"
    checkout.mkdir()
    (checkout / ".mhctools-artifact.json").write_text(
        json.dumps({"revision": revision}), encoding="utf-8")
    status = SimpleNamespace(
        status="ready", manager="user", path=str(checkout))
    monkeypatch.setattr(ANALYSIS, "artifact_status", lambda name: status)

    directory, resolved = ANALYSIS.resolve_netcleave_dir(None)

    assert directory == checkout.resolve()
    assert resolved == revision


def test_missing_netcleave_still_stops_the_run(monkeypatch):
    status = SimpleNamespace(status="missing", manager="mhctools", path="")
    monkeypatch.setattr(ANALYSIS, "artifact_status", lambda name: status)
    with pytest.raises(SystemExit, match="NetCleave is not available"):
        ANALYSIS.resolve_netcleave_dir(None)


def test_corrupt_artifact_record_is_diagnosed_not_a_traceback(tmp_path):
    """An interrupted manifest write must not surface as JSONDecodeError."""
    checkout = tmp_path / "snapshot"
    checkout.mkdir()
    (checkout / ".mhctools-artifact.json").write_text("{trunc", encoding="utf-8")
    with pytest.raises(SystemExit, match="is not valid JSON"):
        ANALYSIS.resolve_netcleave_dir(checkout)
