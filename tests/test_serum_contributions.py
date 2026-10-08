"""Kinetic units and independent analytic oracles for enzyme perturbations."""

import math

import pytest

from mhctools import (
    DegradationTarget, EnzymeCutRate, PeptideInput, degradation_curve,
    empirical_terminal_cut_rates, enzyme_removal_effects,
    serum_contribution_evidence, serum_reference_kinetics,
    simulate_enzyme_degradation,
)


def rate(bond, enzyme, value):
    return EnzymeCutRate(bond, enzyme, value, "synthetic kinetic oracle", "test")


def test_serum_reference_units_and_substrate_specific_competition():
    # 19.9 pmol/min/ul serum = 19.9 uM/min per unit serum fraction.
    r = serum_reference_kinetics("oxidized-NPY-1-36", concentration_um=5, serum_fraction=0.1)
    dpp4, app = r["enzymes"]
    assert dpp4["local_rate_per_hour"] == pytest.approx(60 * .1 * 19.9 / 14.4)
    assert app["local_rate_per_hour"] == pytest.approx(60 * .1 * 3 / 40.4)
    assert dpp4["share_of_reported_channels"] == pytest.approx(.949, abs=.001)
    assert [e["enzyme"] for e in r["enzymes"]] == ["DPP4", "aminopeptidase P"]
    product = serum_reference_kinetics("oxidized-NPY-3-36", concentration_um=5, serum_fraction=.1)
    assert product["enzymes"][0]["enzyme"] == "KLKB1"
    assert r["extrapolated_concentration"] is False
    assert serum_reference_kinetics(
        "oxidized-NPY-1-36", concentration_um=.01, serum_fraction=.1)["extrapolated_concentration"]
    with pytest.raises(ValueError, match="exact assay substrate"):
        serum_reference_kinetics("YPSKPDNPGEDAPAEDMARYYSALRHYINLITRQRY", concentration_um=5, serum_fraction=.1)


@pytest.mark.parametrize("kwargs", [dict(concentration_um=0, serum_fraction=.1),
    dict(concentration_um=5, serum_fraction=0), dict(concentration_um=5, serum_fraction=2),
    dict(concentration_um=float("inf"), serum_fraction=.1)])
def test_invalid_reference_conditions(kwargs):
    with pytest.raises(ValueError):
        serum_reference_kinetics("oxidized-NPY-1-36", **kwargs)


def test_library_spectrum_is_not_a_rate_or_an_enzyme_identification():
    evidence = serum_contribution_evidence()["plasma_library"]
    assert sum(evidence["boundary_context_counts_by_bond"].values()) == 113
    cuts = empirical_terminal_cut_rates("ACDEFGHIKLNPQR", total_rate_per_hour=2)
    assert {r.bond for r in cuts} == {1, 2, 12, 13}
    assert math.fsum(r.rate_per_hour for r in cuts) == pytest.approx(2)
    assert [r.enzyme for r in cuts] == ["n_mono", "n_di", "c_mono", "c_di"]
    assert cuts[0].rate_per_hour == pytest.approx(2 * 46 / 90)
    assert empirical_terminal_cut_rates("ACDEFGHIK", total_rate_per_hour=2) is None
    assert empirical_terminal_cut_rates("ACDEFGHIK", total_rate_per_hour=2, allow_transfer=True)
    assert empirical_terminal_cut_rates("ACDE", total_rate_per_hour=2, allow_transfer=True) is None


def test_parallel_enzyme_causes_and_removal_have_analytic_survival():
    def cuts(s):
        return (rate(2, "A", math.log(2)), rate(2, "B", math.log(2)))
    result = enzyme_removal_effects("ACDEF", DegradationTarget("whole", 0, 5), cuts,
        enzymes=("A", "B", "absent"), times_hours=[1], scenario="analytic", n_paths=18000)
    assert result["baseline"][0]["target_in_circulation"] == pytest.approx(.25, abs=.012)
    for row in result["enzyme_removal_effects"][:2]:
        assert row["target_gain_percentage_points"] == pytest.approx(25, abs=1.5)
        assert 0 < row["target_gain_mc_se_upper_bound_pp"] < 1
        assert row["parent_remaining_without_enzyme"] == pytest.approx(.5, abs=.012)
    assert result["enzyme_removal_effects"][2]["target_gain_percentage_points"] == 0
    causes = result["destructive_cut_attribution"][0]["destroyed_fraction_by_enzyme"]
    assert causes["A"] == pytest.approx(.375, abs=.012)
    assert causes["B"] == pytest.approx(.375, abs=.012)


def test_sequential_trim_reassesses_fragment_and_attributes_later_loss():
    target = DegradationTarget("core", 2, 7)
    def cuts(s):
        return (rate(2, "trimmer", 100),) if s == "ACDEFGHIK" else (rate(3, "destroyer", 1),)
    paths = simulate_enzyme_degradation("ACDEFGHIK", target, cuts,
                                       scenario="serial", n_paths=1, horizon_hours=100)
    assert [s.mechanism for s in paths[0].steps] == ["trimmer", "destroyer"]
    assert [s.cut_bond for s in paths[0].steps] == [2, 5]
    assert [s.event for s in paths[0].steps] == ["cut", "target_destroyed"]
    assert "synthetic kinetic oracle" in paths[0].steps[-1].reason
    result = enzyme_removal_effects("ACDEFGHIK", target, cuts, enzymes=[],
        times_hours=[100], scenario="serial", n_paths=10)
    attribution = result["destructive_cut_attribution"][0]
    assert attribution["first_cut_fraction_by_enzyme"] == {"trimmer": 1}
    assert attribution["destroyed_fraction_by_enzyme"] == {"destroyer": 1}


def test_removing_a_protective_route_can_reduce_target_availability():
    target = DegradationTarget("core", 2, 7)
    def cuts(s):
        return (rate(2, "protective_trim", 10), rate(4, "destructive", 1)) if s == "ACDEFGHIK" else ()
    result = enzyme_removal_effects("ACDEFGHIK", target, cuts,
        enzymes=["protective_trim"], times_hours=[1], scenario="competing routes", n_paths=14000)
    # Baseline survival: harmless trimming or no cut. Without trimming: exp(-t).
    expected = 100 * (math.exp(-1) - (10 / 11 + math.exp(-11) / 11))
    assert result["enzyme_removal_effects"][0]["target_gain_percentage_points"] == pytest.approx(expected, abs=1.5)


def test_missing_rates_stay_unknown_and_zero_rates_remain_distinct():
    target = DegradationTarget("core", 2, 7)
    def cuts(s):
        return (rate(2, "trim", 100),) if s == "ACDEFGHIK" else None
    r = enzyme_removal_effects("ACDEFGHIK", target, cuts,
        enzymes=["trim"], times_hours=[1], scenario="missing product", n_paths=10)
    assert r["baseline"][0]["unknown"] == 1
    effect = r["enzyme_removal_effects"][0]
    assert effect["target_gain_percentage_points"] is None
    assert effect["target_gain_bounds_percentage_points"] == (0, 100)
    assert effect["parent_gain_percentage_points"] == 100
    zero = simulate_enzyme_degradation("ACDEFGHIK", target, lambda s: (), scenario="assessed zero")
    assert degradation_curve(zero, [24])[0]["target_in_circulation"] == 1
    assert zero[0].steps[-1].event == "censored"


def test_zero_cleavage_still_allows_competing_clearance():
    paths = simulate_enzyme_degradation("ACDEF", DegradationTarget("whole", 0, 5), lambda s: (),
        scenario="clearance only", clearance_rate_per_hour=math.log(2), n_paths=12000, seed=12)
    row = degradation_curve(paths, [1])[0]
    assert row["target_in_circulation"] == pytest.approx(.5, abs=.015)
    assert row["cleared"] + row["target_in_circulation"] == pytest.approx(1)


def test_chemistry_and_bad_channels_are_not_silently_accepted():
    with pytest.raises(ValueError, match="free termini"):
        simulate_enzyme_degradation(PeptideInput("ACDEF", c_term="amidated"),
            DegradationTarget("whole", 0, 5), lambda s: (), scenario="bad chemistry")
    with pytest.raises(ValueError, match="internal"):
        simulate_enzyme_degradation("ACDEF", DegradationTarget("whole", 0, 5),
                                   lambda s: (rate(5, "A", 1),), scenario="bad bond")
    with pytest.raises(ValueError, match="not a string"):
        simulate_enzyme_degradation("ACDEF", DegradationTarget("whole", 0, 5),
                                   lambda s: (), scenario="bad names", disabled_enzymes="A")
    with pytest.raises(ValueError, match="cut rate"):
        rate(2, "A", -1)
