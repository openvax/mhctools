"""Readable summaries of explicit vaccine scenarios, without efficacy scores."""

from collections import Counter
import math


def summarize_vaccine_trajectory(result):
    """Summarize the last observation and finite-horizon presentation exposure.

    The result must come from ``simulate_vaccine_trajectory``. Missing kinetics
    leave numerical summaries absent. Means at a zero-hour horizon are undefined.
    Source normalization is per peptide or per transcript, never per APC.
    """
    if not result.get("times_hours") or "dose_basis" not in result:
        raise ValueError("Rerun the original scenario with mhctools 3.50.0 or newer to report dose and integrated exposure")
    summary = {"status": result["status"], "horizon_hours": max(result["times_hours"]),
               "initial_copies": result["initial_copies"], "source_units": result["source_units"],
               "presentation": None, "antigen_fates": None}
    if result["status"] == "unassessed":
        return summary
    row = max(result["curves"], key=lambda r: r["time_hours"])
    horizon, initial = row["time_hours"], result["initial_copies"]
    summary["presentation"] = {
        location: {"surface_copies_at_horizon": row["surface_pmhc_by_apc_location"][location],
                   "auc_copy_hours": auc, "auc_copy_hours_per_input_copy": auc / initial,
                   "mean_surface_copies": auc / horizon if horizon > 0 else None}
        for location, auc in row["surface_pmhc_auc_by_apc_location"].items()}
    summary["cumulative_loaded_copies"] = row["cumulative_loaded_copies"]
    summary["extracellular_target_auc_copy_hours"] = row["extracellular_target_auc_copy_hours"]
    summary["antigen_fates"] = {key: row[key] for key in (
        "antigen_copies_supplied", "target_retained_copies", "target_destroyed_copies",
        "antigen_cleared_copies", "pmhc_lost_copies", "antigen_balance_error")}
    return summary


def _text(value):
    return str(value).replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def _number(value):
    return "undefined" if value is None else format(value, ".6g")


def format_vaccine_trajectory_report(result):
    """Return a Markdown report with dose, endpoints, missing inputs and sources."""
    summary = summarize_vaccine_trajectory(result)
    target = result["target"]
    lines = ["# Vaccine trajectory", "", _text(result["scenario"]), "",
             "Status: **%s**. All numerical results are conditional on supplied rates." % result["status"], "",
             "- Delivery: %s; local administration: %s; organism: %s." % (
                 "free SLP" if result["delivery"] == "slp" else "secretion-tagged mRNA-LNP",
                 _text(result["administration"]), _text(result.get("organism") or "unspecified")),
             "- Target: %s (%s), [%s, %s) in the input; %s, MHC-%s." % (
                 _text(target["label"]), target["sequence"], target["start"], target["end"],
                 _text(result["allele"]), result["mhc_class"]),
             "- Starting amount: %s %s (%s)." % (
                 _number(result["initial_copies"]), result["source_units"], result["dose_basis"]),
             "- Observation horizon: 0 to %s hours after injection." % _number(summary["horizon_hours"])]
    dose = result.get("dose")
    if dose:
        lines += ["- Supplied amount: %s %s; %s; source: %s." % (
            _number(dose["amount"]), dose["unit"], dose["basis"], _text(dose["source"]))]
        if dose["molecular_weight_g_per_mol"] is not None:
            lines += ["- Starting-material molecular weight: %s g/mol." % _number(dose["molecular_weight_g_per_mol"])]
    else:
        lines += ["- No administered mass or molar amount is documented."]
    lines += ["", "## Presentation", ""]
    if summary["presentation"] is None:
        lines += ["No numerical prediction: reachable kinetics are missing.", "",
                  "## Required inputs", ""]
        steps = {step["name"]: step for step in result["steps"]}
        for missing in result["missing_kinetics"]:
            if missing.startswith("step:"):
                step = steps[missing[5:]]
                lines += ["- %s (`%s`): %s → %s." % (
                    step["mechanism"], step["name"], step["origin"], step["destination"])]
            else:
                _, compartment, start, end = missing.split(":")
                lines += ["- Current-fragment cleavage in %s, input interval [%s, %s)." % (compartment, start, end)]
    else:
        lines += ["| APC location | Surface copies at horizon | Exposure (complex-hours) | Mean surface copies | Exposure per input copy |",
                  "|---|---:|---:|---:|---:|"]
        labels = {"local": "Injection site", "node": "Draining node", "systemic": "Systemic"}
        for location, values in summary["presentation"].items():
            lines += ["| %s | %s | %s | %s | %s |" % (labels[location], *(
                _number(values[key]) for key in ("surface_copies_at_horizon", "auc_copy_hours",
                                                "mean_surface_copies", "auc_copy_hours_per_input_copy")))]
        lines += ["", "Exposure integrates display from injection to the horizon. Mean copies are totals across the modeled APC compartment, not copies per cell.",
                  "Per-input exposure uses peptide copies for SLP and carrier-bound transcript copies for mRNA.", "",
                  "Cumulative loading: %s complexes. Extracellular target availability: %s copy-hours." % (
                      _number(summary["cumulative_loaded_copies"]), _number(summary["extracellular_target_auc_copy_hours"])),
                  "", "## Antigen accounting at the horizon", "",
                  "| Quantity | Copies |", "|---|---:|"]
        labels = ("Supplied or translated", "Target retained, including loaded/displayed complexes",
                  "Target destroyed by internal cuts", "Antigen cleared", "Loaded/displayed pMHC lost", "Balance error")
        for label, value in zip(labels, summary["antigen_fates"].values()):
            lines += ["| %s | %s |" % (label, _number(value))]
    rate_counts = Counter(step["rate"]["basis"] if step["rate"] else "missing" for step in result["steps"])
    cut_counts = Counter(cut["basis"] for cut in result["cuts"])
    channel_counts = Counter(channel["basis"] for cut in result["cuts"] for channel in cut["channels"])
    lines += ["", "## Evidence and assumptions", "",
              "Rate inventory (includes inactive routes): " + ", ".join("%s %s" % (n, basis) for basis, n in sorted(rate_counts.items())) + ".",
              "Current-fragment assessments: " + (", ".join("%s %s" % (n, basis) for basis, n in sorted(cut_counts.items())) or "none") + ".",
              "Cleavage channels: " + (", ".join("%s %s" % (n, basis) for basis, n in sorted(channel_counts.items())) or "none") + ".",
              "These labels record input provenance; they do not certify biological calibration.", ""]
    if any(step["rate"] and step["rate"]["value_per_hour"] for step in result["steps"]):
        lines += ["| Nonzero rate | Value | Basis | Source |", "|---|---:|---|---|"]
    else:
        lines += ["No nonzero transport, expression or loading rates were supplied."]
    for step in result["steps"]:
        rate = step["rate"]
        if rate and rate["value_per_hour"]:
            unit = "antigens/transcript/hour" if step["kind"] == "birth" else "per hour"
            lines += ["| %s (`%s`) | %s %s | %s | %s |" % (
                step["mechanism"], step["name"], _number(rate["value_per_hour"]), unit, rate["basis"], _text(rate["source"]))]
    cut_sources = sorted({cut["source"] for cut in result["cuts"]} |
                         {channel["source"] for cut in result["cuts"] for channel in cut["channels"]})
    if cut_sources:
        lines += ["", "Cleavage assessment/channel sources:", ""] + ["- " + _text(source) for source in cut_sources]
    lines += ["", "## Interpretation", "",
              "Presentation exposure is an amount-time descriptor, not a validated effective immune dose. Equal AUC can hide different timing; retain the time courses.",
              "Parent half-life, target retention, RNA decay and pMHC turnover describe different quantities. Ongoing mRNA expression prevents a single decay clock from summarizing the trajectory.",
              "Dose scaling assumes the same dilute, linear kinetics. It does not establish equal-dose efficacy or a clinical dose-response relationship.",
              "The JSON retains all rates, fragment assessments, curves and mechanism evidence.", "",
              "Model exclusions:", ""] + ["- " + item for item in result["exclusions"]]
    return "\n".join(lines) + "\n"


def compare_vaccine_trajectories(results):
    """Compare matched endpoints without ranking efficacy or equating RNA/SLP.

    Require the same exact target sequence, allele, MHC class and horizon.
    Known organism differences are rejected. Route differences and unspecified
    context stay visible. Unassessed scenarios remain rows with absent outcomes.
    Results keep their input order; there is no winner or combined score.
    """
    results = list(results)
    if len(results) < 2:
        raise ValueError("A comparison requires at least two trajectory results")
    for result in results:
        if result.get("schema_version") != 1 or result.get("status") not in ("conditional_scenario", "unassessed"):
            raise ValueError("Expected vaccine trajectory result JSON")
        summarize_vaccine_trajectory(result)
    first = results[0]
    def identity(result):
        return (result["target"]["sequence"], result["allele"], result["mhc_class"], max(result["times_hours"]))
    if any(identity(result) != identity(first) for result in results[1:]):
        raise ValueError("Compare the same target sequence, allele, MHC class and observation horizon")
    organisms = {r["organism"] for r in results if r.get("organism") is not None}
    if len(organisms) > 1:
        raise ValueError("Known organism contexts differ; do not pool cross-species exposures")
    notes = ["Conditional scenarios; presentation exposure is not a calibrated immune dose.",
             "Equal area can hide different timing; retain the original time courses.",
             "One peptide and one carrier-bound transcript are different source units."]
    if any(not r.get("organism") or r["administration"] == "unspecified" for r in results):
        notes.append("Organism or administration context is unspecified in at least one scenario.")
    if len({r["administration"] for r in results}) > 1:
        notes.append("Administration routes differ; route-specific rates and evidence require separate review.")
    if any(r["dose_basis"] == "one_input_normalization" for r in results):
        basis = "Includes one-input normalization; no equal administered-dose comparison is established"
    elif all(r.get("dose") and r["dose"]["unit"] == "ug" for r in results) and all(
            math.isclose(r["dose"]["amount"], first["dose"]["amount"], rel_tol=1e-12) for r in results):
        basis = "Equal documented starting-material mass; formulation delivery and efficacy remain conditional"
    elif all(math.isclose(r["initial_copies"], first["initial_copies"], rel_tol=1e-12) for r in results):
        basis = "Equal starting molecule counts; peptide and transcript identities remain distinct"
    else:
        basis = "Unequal specified starting amounts; absolute exposure includes this dose difference"
    scenarios = []
    for result in results:
        scenarios.append({"scenario": result["scenario"], "delivery": result["delivery"],
                          "administration": result["administration"], "organism": result.get("organism"),
                          "dose": result.get("dose"), "dose_basis": result["dose_basis"],
                          "mhctools_version": result["mhctools_version"],
                          "summary": summarize_vaccine_trajectory(result),
                          "missing_kinetics": result["missing_kinetics"],
                          "rates": {s["name"]: s["rate"] for s in result["steps"]},
                          "cleavage_sources": sorted({c["source"] for c in result["cuts"]} |
                                                     {v["source"] for c in result["cuts"] for v in c["channels"]})})
    return {"schema_version": 1, "target_sequence": first["target"]["sequence"],
            "allele": first["allele"], "mhc_class": first["mhc_class"],
            "horizon_hours": max(first["times_hours"]), "dose_comparison": basis,
            "scenarios": scenarios, "interpretation": notes,
            "calibration_requirements": first["mechanism_evidence"].get("calibration_requirements", [])}


def format_vaccine_comparison_report(comparison):
    """Render the matched endpoint comparison in supplied scenario order."""
    lines = ["# Vaccine scenario comparison", "",
             "Target %s; %s, MHC-%s; horizon 0–%s hours." % (
                 comparison["target_sequence"], comparison["allele"], comparison["mhc_class"],
                 _number(comparison["horizon_hours"])), "", comparison["dose_comparison"] + ".", "",
             "| Scenario | Delivery / route | Starting copies | Node exposure (complex-hours) | Node exposure per input copy | Status |",
             "|---|---|---:|---:|---:|---|"]
    for scenario in comparison["scenarios"]:
        summary = scenario["summary"]
        node = summary["presentation"]["node"] if summary["presentation"] is not None else None
        lines += ["| %s | %s / %s | %s | %s | %s | %s |" % (
            _text(scenario["scenario"]), scenario["delivery"], scenario["administration"],
            _number(summary["initial_copies"]), _number(node["auc_copy_hours"]) if node else "unassessed",
            _number(node["auc_copy_hours_per_input_copy"]) if node else "unassessed", summary["status"])]
    lines += ["", "## Interpretation", ""] + ["- " + item for item in comparison["interpretation"]]
    lines += ["", "## Scenario inputs", ""]
    for index, scenario in enumerate(comparison["scenarios"], 1):
        lines += ["%s. %s" % (index, _text(scenario["scenario"])),
                  "   Organism: %s; amount basis: %s; package: %s." % (
                      _text(scenario["organism"] or "unspecified"), scenario["dose_basis"], scenario["mhctools_version"])]
        dose = scenario["dose"]
        if dose:
            lines += ["   Starting amount: %s %s (%s); source: %s." % (
                _number(dose["amount"]), dose["unit"], dose["basis"], _text(dose["source"]))]
        if scenario["missing_kinetics"]:
            lines += ["   Missing kinetics: " + ", ".join(scenario["missing_kinetics"]) + "."]
        lines.append("")
    lines += ["The comparison JSON preserves each scenario's rate values and sources. Individual trajectory reports explain all fragment assessments and exclusions.", ""]
    if comparison["calibration_requirements"]:
        lines += ["## Evidence needed to calibrate these quantities", ""]
        for item in comparison["calibration_requirements"]:
            lines += ["- **%s:** %s" % (item["stage"], item["measurement"])]
    return "\n".join(lines) + "\n"
