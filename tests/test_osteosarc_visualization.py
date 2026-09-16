from datetime import datetime, timezone
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


pytest.importorskip("matplotlib")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
pytest.importorskip("scipy")


SCRIPT_PATH = (
    Path(__file__).parents[1]
    / "analyses"
    / "osteosarc_vaccine_cleavage"
    / "analyze.py"
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
            },
            {
                "sequence_record_id": "record",
                "mhc_class": "I",
                "display_candidate": True,
                "start": 24,
                "end": 32,
                "percentile_rank": 0.4,
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
