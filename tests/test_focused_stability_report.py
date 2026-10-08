"""Synthetic selection and endpoint checks for the focused report exporter."""

import importlib.util
import json
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    "focused_stability_report", Path(__file__).parents[1] /
    "analyses/osteosarc_vaccine_cleavage/focused_stability_report.py")
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def target(mhc_class="I", rank=1, mutant=True, label=None, allele="synthetic"):
    return dict(sequence_record_id="synthetic:one", gene="synthetic",
                mhc_class=mhc_class, percentile_rank=rank, mutant_overlap=mutant,
                target_label=label or "Predicted mutant " + mhc_class,
                start=2, end=6, allele=allele)


def test_best_mutant_is_selected_per_class_without_a_cross_class_winner():
    rows = [target(rank=1), target(rank=.4, allele="better I"),
            target("II", rank=3), target("II", rank=1, allele="better II"),
            target(rank=.01, mutant=False), target(rank=.01, mutant=None),
            target(rank=.01, label="Source target"), target("II", rank=6)]
    selected = report.select_primary_targets(rows)
    assert [r["allele"] for r in selected] == ["better I", "better II"]
    assert report.select_primary_targets([target(rank=3), target(mutant=None)]) == []


@pytest.mark.parametrize("value", [-1, 0, float("nan"), float("inf")])
def test_invalid_native_durations_never_become_zero_or_protection(value):
    with pytest.raises(ValueError, match="positive"):
        report.positive_hours(value)


def test_missing_duration_is_not_a_score_and_long_lifetimes_are_allowed():
    assert report.positive_hours("") is None
    assert report.positive_hours(None) is None
    assert report.positive_hours(72) == 72
    assert report.fmt_hours(72) == "72.0 h"


def test_a_scored_site_is_not_generically_called_a_high_probability_cut():
    assert report.is_candidate(dict(model="dpp4-qpisa", status="scored", score=1))
    assert not report.is_candidate(dict(model="dpp4-qpisa", status="scored", score=-1))
    assert not report.is_candidate(dict(model="other native endpoint", status="scored", score=99))
    assert not report.is_candidate(dict(model="motif", status="not_matched", score=None))


def test_export_cannot_overwrite_or_nest_inside_frozen_inputs(tmp_path):
    root = tmp_path / "frozen"
    for output in (root, root / "child", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            report.export(root, output, "synthetic", overwrite=True)


def test_manifest_tampering_and_path_escape_are_rejected(tmp_path):
    item = tmp_path / "native.csv"
    item.write_text("synthetic native evidence")
    manifest = tmp_path / "SHA256SUMS.json"
    manifest.write_text(json.dumps({"files": {item.name: report.sha256(item)}}))
    report.verify_manifest(tmp_path)
    item.write_text("changed evidence")
    with pytest.raises(ValueError, match="checksum mismatch"):
        report.verify_manifest(tmp_path)
    manifest.write_text(json.dumps({"../outside.csv": "untrusted"}))
    with pytest.raises(ValueError, match="leaves frozen"):
        report.verify_manifest(tmp_path)


def test_censored_retention_is_a_bound_not_a_measured_lifetime():
    assert report.fmt_retention(dict(median_status="beyond_horizon",
                                     retention_median_lower_bound_hours=72)) == "> 72.0 h"
    assert report.fmt_retention(dict(median_status="unassessed_paths",
                                     retention_median_hours=None)) == "Unavailable"


def test_class_two_mutant_flank_is_not_called_a_minimal_epitope():
    assert "minimal class-I" in report.tracked_description(dict(mhc_class="I"))
    minimal = dict(mhc_class="II", binding_core="ACDEFGHIK", target_sequence="ACDEFGHIK")
    assert "minimal class-II" in report.tracked_description(minimal)
    extended = dict(minimal, target_sequence="PACDEFGHIK")
    assert report.tracked_description(extended) == "Class-II binding core plus mutant flank"
