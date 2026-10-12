"""Dose conversions, endpoint interpretation and readable missing-input reports."""

from dataclasses import replace
import json

import pytest

from mhctools import VaccineDose, VaccineTrajectoryInput, simulate_vaccine_trajectory
from mhctools.cli.vaccine_trajectory import main
from mhctools.vaccine_trajectory_report import format_vaccine_trajectory_report, summarize_vaccine_trajectory
from .test_vaccine_trajectory import model, no_cuts, set_rates


def test_dose_conversion_and_legacy_normalization():
    assert VaccineDose(1, "pmol", "measured", "Assay record").copies == pytest.approx(6.02214076e11)
    dose = VaccineDose(1, "ug", "assumed", "Synthetic material", 1000)
    assert dose.copies == pytest.approx(6.02214076e14)
    m = replace(model(), dose=dose, administration="SC", organism="Synthetic example")
    result = simulate_vaccine_trajectory(m, [0, 2], no_cuts)
    assert result["initial_copies"] == pytest.approx(dose.copies)
    assert result["curves"][-1]["extracellular_target_auc_copy_hours"] == pytest.approx(2 * dose.copies)
    assert result["dose"]["source"] == "Synthetic material"
    normalized = simulate_vaccine_trajectory(model(), [0], no_cuts)
    assert normalized["initial_copies"] == 1
    assert normalized["dose_basis"] == "one_input_normalization"
    assert simulate_vaccine_trajectory(model(initial_copies=4), [0], no_cuts)["dose_basis"] == "explicit_copies"


@pytest.mark.parametrize("kwargs", [dict(amount=0), dict(amount=True), dict(amount=float("nan")),
    dict(unit="mg"), dict(unit="ug"), dict(basis="calibrated"), dict(source=""),
    dict(molecular_weight_g_per_mol=1000), dict(unit="ug", molecular_weight_g_per_mol=0),
    dict(unit="ug", molecular_weight_g_per_mol=True)])
def test_invalid_or_ambiguous_doses_fail(kwargs):
    values = dict(amount=1, unit="copies", basis="assumed", source="Synthetic")
    with pytest.raises(ValueError):
        VaccineDose(**{**values, **kwargs})


def test_dose_cannot_conflict_with_copy_count_or_guess_rna_molecular_weight():
    with pytest.raises(ValueError, match="not both"):
        simulate_vaccine_trajectory(model(initial_copies=3, dose=VaccineDose(1, "copies", "assumed", "Toy")), [1])
    with pytest.raises(ValueError, match="molecular_weight"):
        VaccineDose(1, "ug", "assumed", "RNA amount without transcript mass")
    with pytest.raises(ValueError, match="administration"):
        simulate_vaccine_trajectory(model(administration="IV"), [1])


def test_summary_horizon_is_not_input_order_and_zero_time_mean_is_undefined():
    m = set_rates(model(), interstitium_clearance=1)
    result = simulate_vaccine_trajectory(m, [2, 0, 1], no_cuts)
    summary = summarize_vaccine_trajectory(result)
    assert summary["horizon_hours"] == 2
    assert summary["antigen_fates"]["antigen_cleared_copies"] > .8
    zero = summarize_vaccine_trajectory(simulate_vaccine_trajectory(m, [0], no_cuts))
    assert zero["presentation"]["node"]["mean_surface_copies"] is None
    report = format_vaccine_trajectory_report(result)
    assert "Draining node" in report and "one_input_normalization" in report
    assert "not copies per cell" in report and "not a validated effective immune dose" in report


def test_missing_report_explains_steps_and_fragments_without_numerical_exposure():
    result = simulate_vaccine_trajectory(replace(model(), rates={}), [12])
    summary = summarize_vaccine_trajectory(result)
    assert summary["presentation"] is None
    report = format_vaccine_trajectory_report(result)
    assert "No numerical prediction" in report
    assert "Entry into afferent lymph" in report
    assert "Current-fragment cleavage in interstitium" in report
    assert "| Draining node | 0" not in report


def test_cli_dose_report_roundtrip_and_no_overwrite(tmp_path):
    value = dict(schema_version=1, peptide="ACDEFGHIK", target=dict(label="chosen", start=2, end=7),
                 allele="HLA-A*02:01", mhc_class="I", delivery="slp", scenario="Unassessed example",
                 tap_max_length=9, times_hours=[0, 24], administration="IM", organism="Homo sapiens",
                 dose=dict(amount=2, unit="pmol", basis="assumed", source="Illustrative amount"))
    parsed = VaccineTrajectoryInput.from_dict(value)
    assert parsed.source_copies() == pytest.approx(2 * 6.02214076e11)
    source, output, report = (tmp_path / name for name in ("input.json", "out.json", "report.md"))
    source.write_text(json.dumps(value))
    main(["--input", str(source), "--output", str(output), "--report", str(report)])
    assert "2 pmol" in report.read_text() and "Homo sapiens" in report.read_text()
    assert json.loads(output.read_text())["administration"] == "IM"
    fresh = tmp_path / "fresh.json"
    with pytest.raises(SystemExit):
        main(["--input", str(source), "--output", str(fresh), "--report", str(report)])
    assert not fresh.exists()
