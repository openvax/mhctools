"""First passage, blood-core accounting and intracellular mRNA route oracles."""

from dataclasses import replace
import math

import pytest

from mhctools import DegradationTarget, simulate_vaccine_trajectory
from mhctools.vaccine_exposure import summarize_core_exposure
from .test_vaccine_trajectory import model, set_rates, no_cuts, assessed, cut, erlang_cdf


def test_slp_node_entry_counts_uptake_before_loading_and_competing_clearance():
    m = set_rates(model(), lymph_entry=1, lymph_transit=1, node_antigen_uptake=1, interstitium_clearance=1)
    out = simulate_vaccine_trajectory(m, [0, 100], no_cuts)
    row = out["curves"][-1]
    assert row["node_apc_entry_fraction_of_dosed_peptide"] == pytest.approx(.5)
    assert row["cumulative_node_apc_entry_by_route"]["node_fluid_uptake"] == pytest.approx(.5)
    assert row["cumulative_loaded_copies"] == 0
    scaled = simulate_vaccine_trajectory(replace(m, initial_copies=6e12), [100], no_cuts)["curves"][0]
    assert scaled["node_apc_entry_fraction_of_dosed_peptide"] == pytest.approx(.5)
    assert scaled["cumulative_node_apc_entry_copies"] == pytest.approx(3e12)


def test_apc_migration_counts_first_arrival_even_if_core_is_later_destroyed():
    m = set_rates(model(), local_antigen_uptake=1, migration_endosome=1)
    def cuts(compartment, fragment):
        return assessed(cut(4)) if compartment == "node_endosome" else no_cuts(compartment, fragment)
    row = simulate_vaccine_trajectory(m, [3], cuts)["curves"][0]
    assert row["cumulative_node_apc_entry_by_route"]["local_apc_migration"] == pytest.approx(erlang_cdf(2, 3))
    assert row["node_apc_target_copies"] == pytest.approx(math.exp(-3) * 3**2 / 2)


def test_node_translation_secretion_and_reuptake_do_not_recount_arrival():
    m = set_rates(model("secreted_mrna"), carrier_drainage=1, node_rna_uptake=1,
                  node_rna_escape=1, node_translation=10, node_rna_decay=1,
                  node_signal_entry=1, node_secretion=1, node_antigen_uptake=1)
    row = simulate_vaccine_trajectory(m, [40], no_cuts)["curves"][0]
    assert row["antigen_copies_supplied"] == pytest.approx(10)
    assert row["cumulative_node_apc_entry_copies"] == pytest.approx(10)
    assert row["cumulative_node_apc_entry_by_route"]["node_fluid_uptake"] == 0
    assert row["node_apc_entry_fraction_of_dosed_peptide"] is None
    assert row["node_apc_entry_copies_per_input"] == pytest.approx(10)
    assert row["node_apc_entry_fraction_of_supplied_antigen"] == pytest.approx(1)


def test_cytosolic_mrna_can_be_processed_and_presented_without_secretion():
    m = set_rates(model("cytosolic_mrna", target=DegradationTarget("whole", 0, 9)),
                  carrier_drainage=1, node_rna_uptake=1, node_rna_escape=1, node_translation=2,
                  node_rna_decay=1, node_tap=1, node_er_loading=1, node_surface_export=1)
    out = simulate_vaccine_trajectory(m, [40], no_cuts)
    assert not any("secretion" in s["name"] for s in out["steps"])
    assert out["curves"][0]["surface_pmhc_by_apc_location"]["node"] == pytest.approx(2)


def test_producer_cross_dressing_preserves_antigen_and_does_not_recount_apc_loading():
    m = set_rates(model("cytosolic_mrna", target=DegradationTarget("whole", 0, 9)),
                  producer_rna_uptake=1, producer_rna_escape=1, producer_translation=3,
                  producer_rna_decay=1, producer_tap=1, producer_er_loading=1,
                  producer_surface_export=1, cross_dressing_node=1)
    out = simulate_vaccine_trajectory(m, [40], no_cuts)
    row = out["curves"][0]
    assert row["cumulative_loaded_copies"] == 0
    assert row["cumulative_producer_loaded_copies"] == pytest.approx(3)
    assert row["cumulative_node_apc_entry_by_route"]["cross_dressing"] == pytest.approx(3)
    assert row["surface_pmhc_by_apc_location"]["node"] == pytest.approx(3)
    assert row["antigen_balance_error"] == pytest.approx(0, abs=1e-10)


def test_blood_tracks_core_across_flank_cuts_and_reports_absorption_curve():
    m = set_rates(model(), vascular_absorption=1, blood_clearance=1)
    def cuts(compartment, fragment):
        if compartment == "blood":
            if fragment.start == 0:
                return assessed(cut(2, 10))
            if fragment.end == 9:
                return assessed(cut(5, 10))
        return no_cuts(compartment, fragment)
    out = simulate_vaccine_trajectory(m, [0, 1, 2, 3, 4], cuts)
    for row in out["curves"]:
        t = row["time_hours"]
        assert row["blood_target_copies"] == pytest.approx(t * math.exp(-t))
        assert row["blood_target_copies"] == pytest.approx(row["blood_free_core_copies"] + row["blood_precursor_copies"])
    blood = summarize_core_exposure(out)["blood"]
    assert blood["sampled_peak_time_hours"] == 1
    assert blood["post_peak_half_decline"]["hours_after_sampled_peak_interval"] == [1, 2]
    assert blood["intrinsic_circulation_half_life_hours"] is None


def test_missing_kinetics_and_absent_blood_do_not_acquire_half_lives():
    missing = simulate_vaccine_trajectory(replace(model(), rates={}), [24])
    assert summarize_core_exposure(missing)["node_apc_entry"] is None
    out = simulate_vaccine_trajectory(model(), [0, 24], no_cuts)
    assert summarize_core_exposure(out)["blood"]["post_peak_half_decline"]["status"] == "no_blood_target_observed"
