from datetime import datetime, timezone
import json

import pytest

from mhctools.vaccine_report import (
    VaccineReportInput,
    generate_vaccine_report,
    placement_assessments,
    processing_route_policy,
    select_mhc_windows,
)


def manifest(delivery="synthetic_long_peptide", routing=None):
    record = {
        "id": "construct-1",
        "sequence": "ACDEFGHIK",
        "delivery": delivery,
        "intended_epitopes": [{"label": "target", "start": 3, "end": 7}],
        "mhc_windows": [
            {"mhc_class": "I", "allele": "HLA-A*01:01", "start": 1,
             "end": 2, "percentile_rank": 0.01, "predictor": "p"},
            {"mhc_class": "I", "allele": "HLA-B*08:01", "start": 3,
             "end": 7, "percentile_rank": 0.2, "predictor": "p"},
            {"mhc_class": "I", "allele": "HLA-C*07:01", "start": 2,
             "end": 6, "percentile_rank": 0.3, "predictor": "p"},
        ],
        "cleavage_tracks": [
            {"name": "proteasome", "context": "cytosolic_proteasome",
             "evidence_type": "quantitative_model",
             "scores": [0.1, 0.6, None, 0.8, 0.2, 0.9, 0.3, 0.4],
             "threshold": 0.5, "units": "native"},
            {"name": "serum", "context": "circulation",
             "evidence_type": "motif_rule",
             "scores": [0, 0, 0, 1, 0, 0, 0, 0], "threshold": 1},
        ],
    }
    if routing is not None:
        record["routing"] = routing
    return {"schema_version": 1, "constructs": [record]}


def test_delivery_route_policy_distinguishes_slp_and_cytosolic_rna():
    slp = VaccineReportInput.from_dict(manifest()).constructs[0]
    rna = VaccineReportInput.from_dict(
        manifest("rna_encoded", "cytosolic")
    ).constructs[0]
    slp_policy = {item.context: item.relevance for item in processing_route_policy(slp)}
    rna_policy = {item.context: item.relevance for item in processing_route_policy(rna)}
    assert slp_policy["endolysosomal"] == "primary"
    assert slp_policy["cytosolic_proteasome"] == "conditional"
    assert slp_policy["circulation"] == "not_applicable"
    assert rna_policy["cytosolic_proteasome"] == "primary"
    assert rna_policy["endolysosomal"] == "conditional"
    assert rna_policy["extracellular_interstitial"] == "not_applicable"


def test_rna_defaults_to_cytosolic_but_can_declare_other_routing():
    record = VaccineReportInput.from_dict(manifest("rna_encoded")).constructs[0]
    assert record.routing == "cytosolic"
    with pytest.raises(ValueError, match="routing"):
        VaccineReportInput.from_dict(manifest("rna_encoded", "blood"))


def test_cleavage_arrays_are_internal_bonds_only():
    value = manifest()
    value["constructs"][0]["cleavage_tracks"][0]["scores"].append(0.0)
    with pytest.raises(ValueError, match="exactly 8 internal-bond"):
        VaccineReportInput.from_dict(value)


def test_mhc_selection_keeps_intended_overlap_then_allele_diversity():
    record = VaccineReportInput.from_dict(manifest()).constructs[0]
    selected = select_mhc_windows(record, "I", maximum=2)
    assert selected[0].allele == "HLA-B*08:01"
    assert len(selected) == 2
    assert len({item.allele for item in selected}) == 2


def test_placement_keeps_internal_boundary_and_route_separate():
    record = VaccineReportInput.from_dict(manifest()).constructs[0]
    rows = placement_assessments(record)
    proteasome = next(row for row in rows if row["track"] == "proteasome")
    assert proteasome["internal_assessed_bonds"] == [4, 5, 6]
    assert proteasome["internal_supported_bonds"] == [4, 6]
    assert proteasome["boundary_assessed_bonds"] == [2, 7]
    assert proteasome["boundary_supported_bonds"] == [2]
    assert proteasome["route_relevance"] == "conditional"


def test_generate_report_uses_timestamped_directory_and_checksums(tmp_path):
    pytest.importorskip("matplotlib")
    generated_at = datetime(2026, 9, 17, 12, 34, 56, 123456, tzinfo=timezone.utc)
    output = generate_vaccine_report(
        manifest(), tmp_path, generated_at=generated_at, maximum_mhc_windows=2
    )
    assert output.name == "2026-09-17T123456-123456+0000"
    assert (output / "vaccine-processing-report.pdf").stat().st_size > 1000
    assert json.loads((output / "SHA256SUMS.json").read_text())["vaccine-processing-report.pdf"]


def test_figure_lanes_keep_every_selected_window_and_its_audit_rank():
    from mhctools.vaccine_report import _lane_assign

    value = manifest()
    # Nine mutually overlapping windows cannot share a lane, so a fixed lane
    # budget would have to drop some of them.
    value["constructs"][0]["mhc_windows"] = [
        {"mhc_class": "I", "allele": "HLA-A*%02d:01" % index, "start": index,
         "end": 9, "percentile_rank": 0.1 * index, "predictor": "p"}
        for index in range(1, 10)
    ]
    record = VaccineReportInput.from_dict(value).constructs[0]
    selected = select_mhc_windows(record, "I", maximum=9)
    assigned = _lane_assign(selected)
    assert len(assigned) == len(selected) == 9
    assert len({lane for _, _, lane in assigned}) == 9
    # The rank drawn on each bar is the row's display_rank in the audit CSV.
    for rank, window, _ in assigned:
        assert selected[rank - 1] is window
