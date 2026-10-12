"""Matched endpoints and explicitly unequal source/dose bases in comparisons."""

from copy import deepcopy
from dataclasses import replace
import hashlib
import json

import pytest

from mhctools import (
    VaccineDose, compare_vaccine_trajectories, format_vaccine_comparison_report,
    simulate_vaccine_trajectory,
)
from mhctools.cli.vaccine_compare import main
from .test_vaccine_trajectory import model, no_cuts


def result(delivery="slp", **kwargs):
    return simulate_vaccine_trajectory(model(delivery, **kwargs), [0, 12], no_cuts)


def test_normalized_slp_mrna_comparison_keeps_distinct_units_and_no_ranking():
    inputs = [result(), result("secreted_mrna")]
    comparison = compare_vaccine_trajectories(inputs)
    assert "no equal administered-dose" in comparison["dose_comparison"]
    assert [r["delivery"] for r in comparison["scenarios"]] == ["slp", "secreted_mrna"]
    assert "transcript" in comparison["interpretation"][2]
    report = format_vaccine_comparison_report(comparison)
    assert "Node exposure per input copy" in report
    assert "Starting material" in report


@pytest.mark.parametrize("field,value", [("allele", "HLA-A*03:01"), ("mhc_class", "II"),
                                        ("times_hours", [0, 24]), ("target", {"sequence": "OTHER"})])
def test_incompatible_endpoints_are_rejected(field, value):
    a, b = result(), result()
    b[field] = value
    with pytest.raises(ValueError, match="same target"):
        compare_vaccine_trajectories([a, b])


def test_mass_comparison_requires_documented_mass_and_preserves_conversion():
    a = result(dose=VaccineDose(1, "ug", "assumed", "Synthetic peptide MW", 1000))
    b = result("secreted_mrna", dose=VaccineDose(1, "ug", "assumed", "Synthetic full RNA MW", 1000000))
    comparison = compare_vaccine_trajectories([a, b])
    assert comparison["dose_comparison"].startswith("Equal documented starting-material mass")
    assert comparison["scenarios"][0]["summary"]["initial_copies"] == pytest.approx(1000 * comparison["scenarios"][1]["summary"]["initial_copies"])
    unequal = deepcopy(b)
    unequal["dose"]["amount"] = 2
    unequal["initial_copies"] *= 2
    assert compare_vaccine_trajectories([a, unequal])["dose_comparison"].startswith("Unequal")


def test_unknowns_stay_absent_and_species_context_cannot_be_pooled():
    a = result(organism="Homo sapiens", administration="IM")
    b = simulate_vaccine_trajectory(replace(model(), rates={}), [0, 12])
    comparison = compare_vaccine_trajectories([a, b])
    assert comparison["scenarios"][1]["summary"]["presentation"] is None
    assert "unassessed" in format_vaccine_comparison_report(comparison)
    assert any("unspecified" in n for n in comparison["interpretation"])
    b["organism"] = "Mus musculus"
    with pytest.raises(ValueError, match="organism"):
        compare_vaccine_trajectories([a, b])


def test_cli_comparison_preserves_input_hashes_and_outputs(tmp_path):
    paths = [tmp_path / "slp.json", tmp_path / "rna.json"]
    for path, value in zip(paths, [result(), result("secreted_mrna")]):
        path.write_text(json.dumps(value))
    out, report = tmp_path / "comparison.json", tmp_path / "comparison.md"
    args = ["--input", *map(str, paths), "--output", str(out), "--report", str(report)]
    main(args)
    comparison = json.loads(out.read_text())
    assert comparison["input_files"][0]["sha256"] == hashlib.sha256(paths[0].read_bytes()).hexdigest()
    before = report.read_text()
    with pytest.raises(SystemExit):
        main(args)
    assert report.read_text() == before
    with pytest.raises(ValueError, match="at least two"):
        compare_vaccine_trajectories([result()])


def test_old_results_request_recomputation_instead_of_inventing_exposure():
    old = result()
    del old["times_hours"]
    with pytest.raises(ValueError, match="Rerun the original scenario"):
        compare_vaccine_trajectories([old, result()])
