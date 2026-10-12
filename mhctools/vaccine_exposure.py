"""Core-bearing blood exposure and first nodal-APC arrival summaries."""


def summarize_core_exposure(result):
    """Describe the observed blood curve and cumulative nodal APC entry.

    The half-decline interval is measured from the sampled maximum and brackets
    the first subsequent downward half-maximum crossing. It is not a fitted
    terminal elimination half-life: injection-site absorption, secretion and
    precursor processing can all supply blood antigen throughout observation.
    No interpolation, extrapolation or physiological rate fitting is performed.
    """
    if "node_entry_definition" not in result:
        raise ValueError("Rerun the scenario with mhctools 3.52.0 or newer for core delivery accounting")
    summary = {"status": result["status"], "node_apc_entry": None, "blood": None}
    if result["status"] == "unassessed":
        return summary
    rows = sorted({r["time_hours"]: r for r in result["curves"]}.values(), key=lambda r: r["time_hours"])
    final = rows[-1]
    peak_index = max(range(len(rows)), key=lambda i: rows[i]["blood_target_copies"])
    peak = rows[peak_index]
    maximum = peak["blood_target_copies"]
    decline = {"status": "no_blood_target_observed" if maximum == 0 else "not_observed_within_horizon",
               "hours_after_sampled_peak_interval": None,
               "definition": "First sampled downward half-maximum crossing of all blood fragments retaining this core; not intrinsic elimination half-life"}
    if maximum > 0:
        for previous, current in zip(rows[peak_index:], rows[peak_index + 1:]):
            if current["blood_target_copies"] <= maximum / 2:
                decline.update(status="bracketed_on_observation_grid", hours_after_sampled_peak_interval=[
                    previous["time_hours"] - peak["time_hours"], current["time_hours"] - peak["time_hours"]])
                break
    summary["blood"] = {
        "core_bearing_copies_at_horizon": final["blood_target_copies"],
        "exact_free_core_copies_at_horizon": final["blood_free_core_copies"],
        "longer_precursor_copies_at_horizon": final["blood_precursor_copies"],
        "auc_copy_hours": final["extracellular_target_auc_by_compartment"]["blood"],
        "sampled_peak_copies": maximum, "sampled_peak_time_hours": peak["time_hours"] if maximum > 0 else None,
        "post_peak_half_decline": decline,
        "intrinsic_circulation_half_life_hours": None,
        "intrinsic_half_life_status": "not_identified_from_an_IM_or_SC_input_curve",
    }
    summary["node_apc_entry"] = {
        "cumulative_core_bearing_copies": final["cumulative_node_apc_entry_copies"],
        "by_route": final["cumulative_node_apc_entry_by_route"],
        "fraction_of_dosed_peptide_copies": final["node_apc_entry_fraction_of_dosed_peptide"],
        "copies_per_input_molecule": final["node_apc_entry_copies_per_input"],
        "fraction_of_antigen_supplied_by_horizon": final["node_apc_entry_fraction_of_supplied_antigen"],
        "current_core_bearing_copies": final["node_apc_target_copies"],
        "definition": result["node_entry_definition"],
    }
    return summary
