"""Synthetic coordinate and kinetic oracles for exploratory fragment paths."""

import math

import pytest

from mhctools import PeptideInput
from mhctools.serum_degradation import (
    DegradationTarget, degradation_curve, simulate_target_degradation, target_fragments,
)


def run(sequence="ACDEFGHIK", target=None, hours=1.0, weights=None, **kwargs):
    target = target or DegradationTarget("target", 2, 7)
    half_lives = {f.sequence: hours for f in target_fragments(sequence, target)}
    return simulate_target_degradation(
        sequence, target, half_lives, weights or (lambda s: [1.0] * (len(s) - 1)),
        estimator="synthetic", scenario="test", **kwargs)


def test_boundaries_retain_target_and_recalculate_new_termini():
    seen = []

    def weights(sequence):
        seen.append(sequence)
        # Remove N-flank, then C-flank, then split the target itself.
        bond = {"ACDEFGHIK": 2, "DEFGHIK": 5, "DEFGH": 2}[sequence]
        return [int(i == bond) for i in range(1, len(sequence))]

    paths = run(weights=weights, n_paths=1, horizon_hours=100)
    steps = paths[0].steps
    assert seen == ["ACDEFGHIK", "DEFGHIK", "DEFGH"]
    assert [s.cut_bond for s in steps] == [2, 7, 4]
    assert [s.event for s in steps] == ["cut", "cut", "target_destroyed"]
    assert [(s.start, s.end) for s in steps] == [(0, 9), (2, 9), (2, 7)]
    curve = degradation_curve(paths, [0, steps[0].event_hours, steps[-1].event_hours])
    assert curve[0]["parent_remaining"] == 1
    assert curve[1]["parent_remaining"] == 0
    assert curve[1]["target_in_circulation"] == 1
    assert curve[2]["target_destroyed"] == 1


def test_each_fragment_uses_its_own_half_life():
    target = DegradationTarget("target", 2, 7)
    half_lives = {"ACDEFGHIK": 0.001, "DEFGHIK": 1e12}
    paths = simulate_target_degradation(
        "ACDEFGHIK", target, half_lives,
        lambda s: [int(b == 2) for b in range(1, len(s))],
        estimator="synthetic", scenario="long lived product", n_paths=1)
    assert [s.half_life_hours for s in paths[0].steps] == [0.001, 1e12]
    assert [s.event for s in paths[0].steps] == ["cut", "censored"]


def test_unknown_fragment_is_not_falsely_protected():
    target = DegradationTarget("target", 2, 7)
    paths = simulate_target_degradation(
        "ACDEFGHIK", target, {"ACDEFGHIK": 0.001},
        lambda s: [int(b == 2) for b in range(1, len(s))],
        estimator="one", scenario="unknown product", n_paths=20)
    row = degradation_curve(paths, [24])[0]
    assert row["unknown"] == 1
    assert row["target_in_circulation"] == 0
    assert row["target_destroyed"] == 0


def test_exact_epitope_has_analytic_exponential_survival():
    target = DegradationTarget("whole", 0, 9)
    paths = run(target=target, n_paths=15000, seed=3)
    row = degradation_curve(paths, [1])[0]
    assert row["parent_remaining"] == pytest.approx(0.5, abs=0.015)
    assert row["target_in_circulation"] == row["parent_remaining"]
    assert row["target_in_fragment"] == 0


def test_competing_clearance_and_uptake_analytic_limit():
    paths = run(target=DegradationTarget("whole", 0, 9), n_paths=20000,
                clearance_rate_per_hour=math.log(2), uptake_rate_per_hour=math.log(2), seed=8)
    row = degradation_curve(paths, [1])[0]
    assert row["target_in_circulation"] == pytest.approx(0.125, abs=0.012)
    for key in ["cleared", "taken_up", "target_destroyed"]:
        assert row[key] == pytest.approx(0.875 / 3, abs=0.012)
    assert sum(row[k] for k in ["parent_remaining", "target_in_fragment", "cleared",
                                "taken_up", "target_destroyed", "unknown"]) == pytest.approx(1)


def test_uniform_cuts_can_preserve_target_longer_than_parent():
    paths = run(n_paths=8000, seed=29)
    rows = degradation_curve(paths, [0, 0.5, 1, 2, 4, 24])
    assert rows[2]["target_in_circulation"] > rows[2]["parent_remaining"] + 0.1
    assert all(a["target_in_circulation"] >= b["target_in_circulation"] for a, b in zip(rows, rows[1:]))
    for row in rows:
        assert row["target_in_circulation"] >= row["parent_remaining"]


def test_occurrences_and_estimator_identity_are_preserved():
    target = DegradationTarget("second occurrence", 4, 7)
    fragments = target_fragments("ACDEACDEF", target)
    assert any(f.start == 4 and f.sequence == "ACDEF" for f in fragments)
    assert all(f.start <= 4 and f.end >= 7 for f in fragments)
    paths = run(n_paths=2, seed=7)
    assert paths == run(n_paths=2, seed=7)
    other = simulate_target_degradation(
        "ACDEFGHIK", DegradationTarget("target", 2, 7), {}, lambda s: [],
        estimator="other", scenario="test", n_paths=1)
    with pytest.raises(ValueError, match="Do not pool"):
        degradation_curve(paths + other, [0])
    with pytest.raises(ValueError, match="horizon"):
        degradation_curve(paths, [25])


@pytest.mark.parametrize("weights", [[1], [0] * 8, [-1] * 8, [float("nan")] * 8])
def test_invalid_cut_allocation_is_rejected(weights):
    with pytest.raises(ValueError):
        run(weights=lambda s: weights, n_paths=1)


@pytest.mark.parametrize("hours", [-1, float("inf"), float("nan"), True])
def test_invalid_duration_rejected(hours):
    with pytest.raises(ValueError):
        run(hours=hours, n_paths=1)


@pytest.mark.parametrize("hours", [None, 0])
def test_missing_or_zero_half_life_abstains(hours):
    paths = run(hours=hours, n_paths=1)
    assert degradation_curve(paths, [0])[0]["unknown"] == 1


@pytest.mark.parametrize("target", [DegradationTarget("bad", -1, 3), DegradationTarget("bad", 3, 3),
                                  DegradationTarget("bad", 1, 10), DegradationTarget("bad", True, 3)])
def test_invalid_target_rejected(target):
    with pytest.raises(ValueError):
        target_fragments("ACDEFGHIK", target)


def test_modified_parent_is_not_silently_treated_as_free():
    with pytest.raises(ValueError, match="free termini"):
        target_fragments(PeptideInput("ACDEFGHIK", n_term="acetylated"),
                         DegradationTarget("target", 2, 7))


def test_distinct_occurrences_cannot_be_silently_pooled():
    paths = []
    for occurrence in ("first", "second"):
        paths.extend(run(sequence=PeptideInput("ACDEFGHIK", occurrence_id=occurrence), n_paths=1))
    assert paths[0].peptide_input.occurrence_id == "first"
    with pytest.raises(ValueError, match="Do not pool"):
        degradation_curve(paths, [1])


@pytest.mark.parametrize("kwargs", [dict(n_paths=0), dict(n_paths=True), dict(seed=True),
                                   dict(horizon_hours=-1), dict(clearance_rate_per_hour=-1),
                                   dict(uptake_rate_per_hour=float("nan"))])
def test_invalid_scenario_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        run(**kwargs)
