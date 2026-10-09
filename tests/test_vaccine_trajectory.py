"""Independent kinetic oracles and mechanistic guards for vaccine trajectories."""

from dataclasses import replace
import json
import math

import pytest

from mhctools import (
    DegradationTarget, EnzymeCutRate, PeptideInput, VaccineCleavageRates,
    VaccineRate, VaccineTrajectoryInput, simulate_vaccine_trajectory, vaccine_route_steps,
)
from mhctools.cli.vaccine_trajectory import main


def model(delivery="slp", mhc_class="I", **kwargs):
    options = dict(peptide="ACDEFGHIK", target=DegradationTarget("chosen", 2, 7),
                   allele="HLA-A*02:01" if mhc_class == "I" else "HLA-DRB1*04:01",
                   mhc_class=mhc_class, delivery=delivery, scenario="synthetic analytic test",
                   rates={s.name: VaccineRate(0, "disabled", "Excluded in this synthetic oracle")
                          for s in vaccine_route_steps(delivery, mhc_class)},
                   tap_max_length=9 if mhc_class == "I" else None,
                   class_ii_max_length=20 if mhc_class == "II" else None,
                   signal_end=2 if delivery == "secreted_mrna" else None)
    options.update(kwargs)
    return VaccineTrajectoryInput(**options)


def set_rates(m, **values):
    return replace(m, rates={**m.rates, **{name: VaccineRate(value, "assumed", "Synthetic oracle")
                                          for name, value in values.items()}})


def no_cuts(compartment, fragment):
    return VaccineCleavageRates((), "assumed", "No cleavage in this synthetic oracle")


def cut(bond, value=1):
    return EnzymeCutRate(bond, "synthetic activity", value, "assumed", "Analytic test")


def assessed(*channels):
    return VaccineCleavageRates(channels, "assumed", "Synthetic fragment-specific kinetic oracle")


def amount(output, row, compartment, start=None, end=None):
    return sum(value for state, value in zip(output["states"], output["curves"][row]["state_copies"])
               if state["compartment"] == compartment and
               (start is None or state["start"] == start) and (end is None or state["end"] == end))


def erlang_cdf(n, t):
    return 1 - math.exp(-t) * sum(t**k / math.factorial(k) for k in range(n))


def test_lymph_and_vascular_absorption_are_parallel_not_serial():
    m = set_rates(model(), lymph_entry=2, vascular_absorption=3)
    output = simulate_vaccine_trajectory(m, [1, 0, 1], no_cuts)
    assert output["status"] == "conditional_scenario"
    assert amount(output, 0, "interstitium") == pytest.approx(math.exp(-5))
    assert amount(output, 0, "afferent_lymph") == pytest.approx(2/5 * (1-math.exp(-5)))
    assert amount(output, 0, "blood") == pytest.approx(3/5 * (1-math.exp(-5)))
    assert output["curves"][0] == output["curves"][2]
    assert output["curves"][1]["extracellular_target_copies"] == 1
    assert abs(output["curves"][0]["antigen_balance_error"]) < 1e-11


def test_sequential_flank_then_internal_cut_uses_current_coordinates_and_kinetics():
    def cuts(compartment, fragment):
        return assessed(cut(2)) if fragment.start == 0 else assessed(cut(3, 2))
    output = simulate_vaccine_trajectory(model(), [1], cuts)
    row = output["curves"][0]
    assert row["extracellular_parent_copies"] == pytest.approx(math.exp(-1))
    assert row["target_retained_copies"] == pytest.approx(2*math.exp(-1)-math.exp(-2))
    assert row["target_destroyed_copies"] == pytest.approx(1-row["target_retained_copies"])
    assert any(s["start"] == 2 and s["end"] == 9 for s in output["states"])
    assert output["cuts"][1]["channels"][0]["rate_per_hour"] == 2


def test_missing_product_or_lymph_kinetics_abstains_instead_of_assuming_protection():
    def cuts(compartment, fragment):
        return assessed(cut(2)) if fragment.start == 0 and compartment == "interstitium" else None
    output = simulate_vaccine_trajectory(set_rates(model(), lymph_entry=1), [1], cuts)
    assert output["status"] == "unassessed"
    assert output["curves"] == []
    assert "cuts:interstitium:2:9" in output["missing_kinetics"]
    assert "cuts:afferent_lymph:0:9" in output["missing_kinetics"]


def test_missing_transport_rates_and_downstream_requirements_remain_visible():
    output = simulate_vaccine_trajectory(replace(model(), rates={}), [1])
    assert output["status"] == "unassessed" and not output["curves"]
    assert "step:lymph_entry" in output["missing_kinetics"]
    assert "step:node_antigen_uptake" in output["missing_kinetics"]


def test_mhc_i_needs_c_terminal_release_tap_n_trimming_and_exact_loading():
    m = model(peptide="ACDEFGHIKLM", target=DegradationTarget("chosen", 2, 9))
    m = set_rates(m, local_antigen_uptake=1, local_cross_escape=1, local_tap=1,
                  local_er_loading=1, local_surface_export=1)
    def cuts(compartment, fragment):
        if compartment == "local_cytosol" and fragment.end == 11:
            return assessed(cut(9))
        if compartment == "local_er" and fragment.start < 2:
            return assessed(cut(1))
        return no_cuts(compartment, fragment)
    output = simulate_vaccine_trajectory(m, [3], cuts)
    assert output["curves"][0]["cumulative_loaded_copies"] == pytest.approx(erlang_cdf(7, 3))
    assert output["curves"][0]["surface_pmhc_copies"] == pytest.approx(erlang_cdf(8, 3))
    assert all(s["end"] == 9 for s in output["states"] if s["compartment"] == "local_er")
    blocked = simulate_vaccine_trajectory(m, [3], no_cuts)
    assert blocked["curves"][0]["cumulative_loaded_copies"] == 0


def test_vacuolar_loading_requires_exact_class_i_ligand():
    m = set_rates(model(), local_antigen_uptake=1, local_vacuolar_loading=1)
    assert simulate_vaccine_trajectory(m, [2], no_cuts)["curves"][0]["cumulative_loaded_copies"] == 0
    m = replace(m, target=DegradationTarget("whole", 0, 9))
    assert simulate_vaccine_trajectory(m, [2], no_cuts)["curves"][0]["cumulative_loaded_copies"] == pytest.approx(erlang_cdf(2, 2))


def test_mhc_ii_can_load_a_longer_ligand_with_target_intact_and_competes_with_cuts():
    m = set_rates(model(mhc_class="II"), local_antigen_uptake=1, local_ii_loading=1)
    def cuts(compartment, fragment):
        return assessed(cut(4)) if compartment == "local_endosome" else no_cuts(compartment, fragment)
    output = simulate_vaccine_trajectory(m, [30], cuts)
    assert output["curves"][0]["cumulative_loaded_copies"] == pytest.approx(.5)
    assert output["curves"][0]["target_destroyed_copies"] == pytest.approx(.5)
    blocked = simulate_vaccine_trajectory(replace(m, class_ii_max_length=5), [30], no_cuts)
    assert blocked["curves"][0]["cumulative_loaded_copies"] == 0


def test_apc_migration_does_not_recount_loading_or_surface_export():
    m = model(target=DegradationTarget("whole", 0, 9))
    m = set_rates(m, local_antigen_uptake=1, local_vacuolar_loading=1,
                  local_surface_export=1, migration_surface=1)
    output = simulate_vaccine_trajectory(m, [3], no_cuts)
    row = output["curves"][0]
    assert row["cumulative_loaded_copies"] == pytest.approx(erlang_cdf(2, 3))
    assert row["cumulative_surface_export_copies"] == pytest.approx(erlang_cdf(3, 3))
    assert row["surface_pmhc_by_apc_location"]["node"] == pytest.approx(erlang_cdf(4, 3))
    assert row["surface_pmhc_copies"] == pytest.approx(erlang_cdf(3, 3))


def test_pmhc_surface_turnover_is_separate_from_free_target_cleavage():
    m = model(target=DegradationTarget("whole", 0, 9))
    m = set_rates(m, local_antigen_uptake=1, local_vacuolar_loading=1,
                  local_surface_export=1, local_surface_loss=1)
    output = simulate_vaccine_trajectory(m, [3], no_cuts)
    row = output["curves"][0]
    assert row["surface_pmhc_copies"] == pytest.approx(math.exp(-3)*3**3/6)
    assert row["pmhc_lost_copies"] == pytest.approx(erlang_cdf(4, 3))
    assert row["target_destroyed_copies"] == 0
    assert row["antigen_balance_error"] == pytest.approx(0, abs=1e-11)


@pytest.mark.parametrize("decay", [0, 1])
def test_rna_translation_creates_multiple_antigens_without_consuming_transcript(decay):
    m = set_rates(model("secreted_mrna"), producer_rna_uptake=1, producer_rna_escape=1,
                  producer_rna_decay=decay, producer_translation=3)
    output = simulate_vaccine_trajectory(m, [3], no_cuts)
    expected = 3*erlang_cdf(3, 3) if decay else 3*(3-2+(3+2)*math.exp(-3))
    assert output["curves"][0]["antigen_copies_supplied"] == pytest.approx(expected)
    expected_rna = math.exp(-3)*3**2/2 if decay else erlang_cdf(2, 3)
    assert amount(output, 0, "producer_rna") == pytest.approx(expected_rna)
    assert output["source_units"] == "input carrier-bound RNA copies"
    assert output["curves"][0]["antigen_balance_error"] == pytest.approx(0, abs=1e-10)


def test_secreted_rna_removes_signal_and_feeds_interstitium_not_mandatory_blood():
    m = set_rates(model("secreted_mrna"), producer_rna_uptake=2, producer_rna_escape=2,
                  producer_rna_decay=1, producer_translation=4,
                  producer_signal_entry=2, producer_secretion=2, lymph_entry=2, lymph_transit=2,
                  node_antigen_uptake=2, node_cross_escape=2)
    output = simulate_vaccine_trajectory(m, [12], no_cuts)
    assert amount(output, 0, "interstitium", 2, 9) > 0
    assert amount(output, 0, "node_cytosol", 2, 9) > 0
    assert amount(output, 0, "blood") == 0
    assert output["curves"][0]["antigen_balance_error"] == pytest.approx(0, abs=1e-10)


def test_direct_rna_apc_processing_does_not_require_secretion_or_blood():
    m = set_rates(model("secreted_mrna"), local_rna_uptake=2, local_rna_escape=2,
                  local_rna_decay=1, local_translation=3, local_failed_entry=2, local_tap=2,
                  local_er_loading=2, local_surface_export=2)
    m = replace(m, target=DegradationTarget("mature target", 2, 9))
    def cuts(compartment, fragment):
        return assessed(cut(1, 2)) if compartment == "local_er" and fragment.start < m.target.start else no_cuts(compartment, fragment)
    output = simulate_vaccine_trajectory(m, [12], cuts)
    assert output["curves"][0]["surface_pmhc_copies"] > 2.9
    assert output["curves"][0]["extracellular_target_copies"] == 0


def test_secretory_cuts_need_their_own_assessment_and_can_split_target():
    m = set_rates(model("secreted_mrna"), producer_rna_uptake=2, producer_rna_escape=2,
                  producer_rna_decay=1, producer_translation=3, producer_signal_entry=2)
    def cuts(compartment, fragment):
        return assessed(cut(3, 2)) if compartment == "producer_secretory_er" else no_cuts(compartment, fragment)
    output = simulate_vaccine_trajectory(m, [30], cuts)
    assert output["curves"][0]["target_destroyed_copies"] == pytest.approx(3, abs=1e-8)
    unassessed = simulate_vaccine_trajectory(m, [1], lambda c, f: None if c == "producer_secretory_er" else no_cuts(c, f))
    assert "cuts:producer_secretory_er:2:9" in unassessed["missing_kinetics"]
    assert unassessed["curves"] == []


@pytest.mark.parametrize("kwargs", [dict(signal_end=None), dict(signal_end=9),
    dict(signal_end=True), dict(tap_max_length=None), dict(allele="HLA-DRB1*04:01"),
    dict(initial_copies=-1), dict(initial_copies=True), dict(scenario=""),
    dict(target=DegradationTarget("inside signal", 0, 2)),
    dict(target=DegradationTarget("crossing signal boundary", 1, 5))])
def test_invalid_mrna_inputs_are_rejected(kwargs):
    with pytest.raises(ValueError):
        simulate_vaccine_trajectory(model("secreted_mrna", **kwargs), [1], no_cuts)


def test_chemistry_unknown_rate_keys_illegal_er_cuts_and_size_guard():
    with pytest.raises(ValueError, match="free termini"):
        simulate_vaccine_trajectory(model(peptide=PeptideInput("ACDEFGHIK", n_term="acetylated")), [1], no_cuts)
    with pytest.raises(ValueError, match="Unknown kinetic"):
        simulate_vaccine_trajectory(replace(model(), rates={"typo": VaccineRate(1, "assumed", "Test")}), [1], no_cuts)
    with pytest.raises(ValueError, match="max_states"):
        simulate_vaccine_trajectory(set_rates(model(), lymph_entry=1), [1], no_cuts, max_states=1)
    m = set_rates(model(target=DegradationTarget("whole", 0, 9)), local_antigen_uptake=1,
                  local_cross_escape=1, local_tap=1)
    with pytest.raises(ValueError, match="N-terminal"):
        simulate_vaccine_trajectory(m, [1], lambda c, f: assessed(cut(2)) if c == "local_er" else no_cuts(c, f))


@pytest.mark.parametrize("value,basis,source", [(1, "unassessed", "Test"), (None, "assumed", "Test"),
    (-1, "assumed", "Test"), (1, "disabled", "Test"), (float("nan"), "measured", "Test"), (0, "assumed", "")])
def test_rate_evidence_cannot_contradict_value(value, basis, source):
    with pytest.raises(ValueError):
        VaccineRate(value, basis, source)


def test_cli_json_roundtrip_preserves_cut_assumptions_and_refuses_overwrite(tmp_path):
    payload = dict(schema_version=1, peptide="ACDEFGHIK", target=dict(label="chosen", start=2, end=7),
                   allele="HLA-A*02:01", mhc_class="I", delivery="slp", scenario="CLI synthetic",
                   tap_max_length=9, times_hours=[0, 1], rates={s.name: dict(value_per_hour=0,
                       basis="disabled", source="Synthetic disabled route") for s in vaccine_route_steps("slp", "I")},
                   cleavage_assessments=[dict(compartment="interstitium", start=0, end=9,
                       sequence="ACDEFGHIK", channels=[], basis="assumed", source="Explicit no-cut example")])
    source, output = tmp_path / "input.json", tmp_path / "output.json"
    source.write_text(json.dumps(payload))
    main(["--input", str(source), "--output", str(output)])
    result = json.loads(output.read_text())
    assert result["status"] == "conditional_scenario"
    assert result["cuts"][0]["source"] == "Explicit no-cut example"
    assert result["curves"][1]["target_retained_copies"] == pytest.approx(1)
    with pytest.raises(SystemExit):
        main(["--input", str(source), "--output", str(output)])
