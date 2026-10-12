"""Readable summaries of explicit vaccine scenarios, without efficacy scores."""

from collections import Counter


def summarize_vaccine_trajectory(result):
    """Summarize the last observation and finite-horizon presentation exposure.

    The result must come from ``simulate_vaccine_trajectory``. Missing kinetics
    leave numerical summaries absent. Means at a zero-hour horizon are undefined.
    Source normalization is per peptide or per transcript, never per APC.
    """
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
