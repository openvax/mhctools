"""Assay scope, source censoring and reference-clock interpretation."""

import math

import pytest

from mhctools import serum_assay_parent_reference, serum_calibration_evidence
from scripts.serum_calibration_review import measurement_records


def test_reference_point_reproduces_half_life_without_identifying_enzyme_or_target():
    result = serum_assay_parent_reference('yi2015:t2:05', [0, 4, 8, 100])
    assert result['measurement']['matrix'] == 'serum'
    assert result['measurement']['substrate'] == 'G36A'
    assert result['measurement']['replicate_sd_hours'] == .2
    assert [r['parent_remaining'] for r in result['curves'][:3]] == pytest.approx([1, .5, .25])
    assert result['curves'][-1]['beyond_study_maximum_horizon']
    assert result['measurement']['observed_horizon_hours'] is None
    assert result['enzyme_cut_rates'] is None
    assert result['target_survival'] is None
    assert not result['vaccine_transfer_validated']


def test_censored_and_donor_range_references_do_not_invent_midpoints():
    right = serum_assay_parent_reference('yi2015:t2:06', [0, 96])
    assert right['measurement']['perturbation'] == 'proprietary inhibitor cocktail'
    assert right['curves'][0]['parent_remaining_lower'] == 1
    assert right['curves'][1]['parent_remaining_lower'] == pytest.approx(.5)
    assert right['curves'][1]['parent_remaining_upper'] == 1
    assert right['curves'][1]['parent_remaining'] is None
    left = serum_assay_parent_reference('yi2015:t2:24', [0, 24])
    assert left['curves'][0]['parent_remaining_lower'] == 1
    assert left['curves'][1]['parent_remaining_lower'] == 0
    assert left['curves'][1]['parent_remaining_upper'] == pytest.approx(.5)
    donor = serum_assay_parent_reference('yi2015:t2:01', [4])
    assert donor['curves'][0]['parent_remaining'] is None
    assert donor['curves'][0]['parent_remaining_lower'] == pytest.approx(.5)
    assert donor['curves'][0]['parent_remaining_upper'] == pytest.approx(math.exp(-math.log(2) / 6))
    assert donor['measurement']['concentration_um'] is None  # unresolved MS/ELISA row


@pytest.mark.parametrize('times', [[], [-1], [math.inf], [math.nan]])
def test_invalid_reference_times_are_rejected(times):
    with pytest.raises(ValueError):
        serum_assay_parent_reference('yi2015:t2:05', times)


def test_reference_does_not_match_a_sequence_or_another_assay_by_substrate_name():
    with pytest.raises(ValueError, match='Unknown serum assay'):
        serum_assay_parent_reference('HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR', [1])
    evidence = serum_calibration_evidence()
    assert len(evidence['measurements']) == 27
    assert len({r['measurement_id'] for r in evidence['measurements']}) == 27
    assert evidence['measurements'][20]['substrate'] == 'GIP(1–42)'
    assert all(not s['vaccine_rate_calibration'] for s in evidence['studies'].values())
    assert 'PREP' in evidence['studies']['bainbridge2017']['limitations']
    assert 'confound' in evidence['studies']['dufresne2017']['limitations']
    assert 'infusate' in evidence['studies']['torang2016']['limitations']
    evidence['measurements'].clear()
    assert len(serum_calibration_evidence()['measurements']) == 27


def test_transcription_preserves_replicate_sd_donor_range_and_source_censoring():
    rows = [
        ['G36A', 'Serum', 'RT', '4 ± 0.2', 'MS'],
        ['G36A', 'EDTA plasma', 'RT', '4–24', 'MS, ELISA'],
        ['G36A', 'P800 plasma', 'RT', '>96', 'MS'],
    ]
    records = measurement_records(rows)
    assert records[0]['half_life_kind'] == 'mean_with_replicate_sd'
    assert records[0]['replicate_sd_hours'] == .2
    assert records[1]['half_life_kind'] == 'donor_range'
    assert records[1]['concentration_um'] is None
    assert records[2]['half_life_point_hours'] is None
    assert records[2]['half_life_upper_hours'] is None
    assert not records[2]['identifiable_enzyme_hazards']
