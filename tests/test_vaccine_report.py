from datetime import datetime, timezone
import json
from pathlib import Path

import pytest

from mhctools.vaccine_report import (
    VaccineReportInput,
    generate_vaccine_report,
    placement_assessments,
    processing_route_policy,
    select_mhc_windows,
    compact_target_summary,
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
    generated_at = datetime(2026, 9, 17, 12, 34, 56, 123456, tzinfo=timezone.utc)
    output = generate_vaccine_report(
        manifest(), tmp_path, generated_at=generated_at, maximum_mhc_windows=2
    )
    assert output.name == "2026-09-17T123456-123456+0000"
    assert (output / "vaccine-processing-report.pdf").stat().st_size > 1000
    assert json.loads((output / "SHA256SUMS.json").read_text())["vaccine-processing-report.pdf"]
    assert json.loads((output / "target-summary.json").read_text())[0]['serum']['status'] == 'not_applicable'


def test_documented_simple_manifest_uses_actual_tool_output(tmp_path, capsys):
    from mhctools.cli.cleavage import main as cleavage_main
    from mhctools.cli.vaccine_report import main as report_main

    docs = Path(__file__).parents[1] / "docs"
    guide = (docs / "vaccine-reports.md").read_text()
    inline = guide.split("cat > vaccine-report.json <<'JSON'\n", 1)[1].split("\nJSON\n", 1)[0]
    source = docs / "vaccine-report-simple.json"
    assert json.loads(inline) == json.loads(source.read_text())
    manifest_path = tmp_path / "vaccine-report.json"
    manifest_path.write_text(source.read_text())
    cleavage_main(["--sequence", "AAAKIILNFTTLAAA", "--model", "mme-hydrophobic",
                   "--out", str(tmp_path / "cleavage.json")])
    inline_output = guide.split("## Inspect `cleavage.json`\n", 1)[1].split(
        "```json\n", 1)[1].split("\n```", 1)[0]
    assert inline_output + "\n" == (tmp_path / "cleavage.json").read_text()
    report_main(["--input", str(manifest_path), "--output-dir", str(tmp_path)])
    output = Path(capsys.readouterr().out.strip())
    assert output.parent == tmp_path
    assert (output / "vaccine-processing-report.pdf").read_bytes().startswith(b"%PDF-")
    summary, = json.loads((output / "placement-assessments.json").read_text())
    assert summary["internal_supported_bonds"] == [4, 5, 6, 8, 11]
    assert summary["boundary_supported_bonds"] == []
    assert summary["route_relevance"] == "primary"
    assert "categorical_sites" not in summary
    normalized = json.loads((output / "normalized-input.json").read_text())
    serialized_track = normalized["constructs"][0]["cleavage_tracks"][0]
    assert "scores" not in serialized_track
    assert set(serialized_track["evidence"]["sites"][0]) == {
        "bond", "match", "source_bond"}
    saved = VaccineReportInput.from_dict(normalized)
    track, = saved.constructs[0].cleavage_tracks
    assert track.name == "mme-hydrophobic"
    assert track.evidence_type == "motif_rule"
    assert track.evidence.model.name == "mme-hydrophobic"
    assert track.evidence.model.references
    assert track.scores == (None,) * 14


def test_vaccine_manifest_reads_legacy_cleavage_json(tmp_path):
    from dataclasses import asdict
    from mhctools import get_cleavage_model

    result = get_cleavage_model("mme-hydrophobic").predict("AAAKIILNFTTLAAA")
    old = asdict(result)
    model = old.pop("model")
    model.pop("rule")
    for site in old["sites"]:
        site["source_bond"] = site["bond"]
    old["model"] = model["name"]
    (tmp_path / "cleavage.json").write_text(json.dumps({
        "schema_version": 2, "models": {model["name"]: model}, "results": [old]}))
    manifest = tmp_path / "vaccine-report.json"
    manifest.write_text((Path(__file__).parents[1] / "docs" /
                         "vaccine-report-simple.json").read_text())
    report = VaccineReportInput.from_json(manifest)
    track, = report.constructs[0].cleavage_tracks
    assert track.categorical_sites[4] == "matched"


def test_cleavage_json_reference_requires_a_unique_matching_result(tmp_path):
    from mhctools.cli.cleavage import main as cleavage_main

    manifest_path = tmp_path / "report.json"
    value = json.loads((Path(__file__).parents[1] / "docs" /
                        "vaccine-report-simple.json").read_text())
    manifest_path.write_text(json.dumps(value))
    output_path = tmp_path / "cleavage.json"
    cleavage_main(["--sequence", "AAAKIILNFTTLAAA", "--model", "mme-hydrophobic",
                   "--out", str(output_path)])
    output = json.loads(output_path.read_text())
    output["results"].append(output["results"][0])
    output_path.write_text(json.dumps(output))
    with pytest.raises(ValueError, match="exactly one result"):
        VaccineReportInput.from_json(manifest_path)
    output["results"] = output["results"][:1]
    output["results"][0]["peptide"]["sequence"] = "AAAAAAAAAAAAAAA"
    output_path.write_text(json.dumps(output))
    with pytest.raises(ValueError, match="exactly one result"):
        VaccineReportInput.from_json(manifest_path)


def test_cleavage_json_model_selection_and_metadata_are_checked(tmp_path):
    from mhctools.cli.cleavage import main as cleavage_main

    value = json.loads((Path(__file__).parents[1] / "docs" /
                        "vaccine-report-simple.json").read_text())
    manifest_path = tmp_path / "report.json"
    output_path = tmp_path / "cleavage.json"
    cleavage_main(["--sequence", "AAAKIILNFTTLAAA", "--model", "mme-hydrophobic",
                   "--model", "fap-endo-gp", "--out", str(output_path)])
    manifest_path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="exactly one result"):
        VaccineReportInput.from_json(manifest_path)
    track = value["constructs"][0]["cleavage_tracks"][0]
    track["model"] = "mme-hydrophobic"
    manifest_path.write_text(json.dumps(value))
    selected = VaccineReportInput.from_json(manifest_path)
    assert selected.constructs[0].cleavage_tracks[0].name == "mme-hydrophobic"
    track["evidence_type"] = "quantitative_model"
    manifest_path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="unsupported fields"):
        VaccineReportInput.from_json(manifest_path)
    track.pop("evidence_type")
    output = json.loads(output_path.read_text())
    output["models"]["mme-hydrophobic"]["name"] = "different-model"
    output_path.write_text(json.dumps(output))
    manifest_path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="metadata contradicts"):
        VaccineReportInput.from_json(manifest_path)


def test_compact_summary_keeps_internal_and_release_cuts_and_conditional_routes_separate():
    construct = VaccineReportInput.from_dict(manifest()).constructs[0]
    summary, = compact_target_summary(construct)
    assert summary['target_sequence'] == 'DEFGH'
    assert summary['mrna_expression']['status'] == 'not_applicable'
    assert summary['antigen_processing']['status'] == 'candidate_internal_cuts'
    track, = summary['antigen_processing']['evidence']
    assert track['relevance'] == 'conditional'
    assert track['internal_cut_flags'] == [4, 6]
    assert track['boundary_cut_flags'] == [2]
    assert track['internal_unassessed_bonds'] == [3]
    assert track['threshold'] == .5
    assert summary['antigen_processing']['mhc_loading'] == 'unassessed'
    assert summary['serum']['target_survival'] is None
    assert summary['matching_mhc_windows'][0]['allele'] == 'HLA-B*08:01'


def test_compact_summary_never_calls_missing_or_below_threshold_scores_protection():
    value = manifest('rna_encoded')
    track = value['constructs'][0]['cleavage_tracks'][0]
    track['scores'] = [.1] * 8
    summary, = compact_target_summary(VaccineReportInput.from_dict(value).constructs[0])
    assert summary['mrna_expression']['status'] == 'unassessed'
    assert summary['antigen_processing']['status'] == 'site_evidence_only'
    assert summary['serum']['status'] == 'not_applicable'
    track['scores'] = [None] * 8
    summary, = compact_target_summary(VaccineReportInput.from_dict(value).constructs[0])
    assert summary['antigen_processing']['status'] == 'unassessed'
    assert summary['antigen_processing']['evidence'][0]['internal_unassessed_bonds'] == [3, 4, 5, 6]


def test_compact_circulation_summary_requires_rates_even_when_motif_flag_is_present():
    value = manifest()
    value['constructs'][0]['exposures'] = ['circulation']
    summary, = compact_target_summary(VaccineReportInput.from_dict(value).constructs[0])
    assert summary['serum']['relevance'] == 'conditional'
    assert summary['serum']['status'] == 'unassessed'
    assert summary['serum']['evidence'][0]['internal_cut_flags'] == [4]
    assert summary['serum']['target_survival'] is None


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
def test_canonical_categorical_tracks_survive_report_round_trip(tmp_path):
    from mhctools import PeptideInput, get_cleavage_model
    from mhctools.vaccine_report import (
        CleavageTrack, VaccineConstruct, VaccineReportInput, placement_assessments,
    )
    from dataclasses import asdict

    result = get_cleavage_model("cpn-basic").predict(PeptideInput("RPPGFSPFR"))
    track = CleavageTrack.from_result(result, "circulation", "AARPPGFSPFRGG", start=2,
                                     conditional_on="If boundary cleavage releases bradykinin")
    value = dict(id="peptide", sequence="AARPPGFSPFRGG", delivery="synthetic_long_peptide",
                 exposures=["circulation"], intended_epitopes=[dict(label="target", start=3, end=11)],
                 cleavage_tracks=[asdict(track)])
    construct = VaccineConstruct.from_dict(value)
    assert construct.cleavage_tracks[0].evidence == result
    row = placement_assessments(construct)[0]
    assert row["internal_supported_bonds"] == [10]
    assert "canonical_evidence" not in row
    assert row["conditional_on"]
    assert not any(score is not None for score in track.scores)
    report = VaccineReportInput.from_dict(dict(schema_version=1, constructs=[value]))
    restored = VaccineReportInput.from_dict(report.to_dict())
    assert restored.constructs[0].cleavage_tracks[0] == track
    output = generate_vaccine_report(report, tmp_path / "canonical-report")
    assert (output / "vaccine-processing-report.pdf").stat().st_size > 0
    saved = json.loads((output / "normalized-input.json").read_text())
    evidence = saved["constructs"][0]["cleavage_tracks"][0]["evidence"]
    assert evidence["sites"][0]["match"] == 1
    assert "score" not in evidence["sites"][0]
    assert VaccineReportInput.from_dict(saved).constructs[0].cleavage_tracks[0] == track


def test_canonical_track_keeps_unsupported_and_substrate_only_observations():
    import pytest
    from mhctools import get_cleavage_model
    from mhctools.vaccine_report import CleavageTrack

    source = get_cleavage_model("nln-observed").predict("YGGFLRRIR")
    track = CleavageTrack.from_result(source, "cytosolic_proteasome", "YGGFLRRIR")
    assert track.evidence.substrate_observation == "cleavage_reported"
    assert track.scores == (None,) * 8 and track.categorical_sites == {}
    unsupported = get_cleavage_model("lnpep-observed").predict("AAAAAAAAA")
    assert CleavageTrack.from_result(unsupported, "endolysosomal", "AAAAAAAAA").evidence.unsupported_reason
    with pytest.raises(ValueError, match="threshold"):
        CleavageTrack.from_result(source, "cytosolic_proteasome", "YGGFLRRIR", threshold=0.5)
