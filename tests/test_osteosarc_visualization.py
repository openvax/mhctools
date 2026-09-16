from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd


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
