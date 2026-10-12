"""Dose normalization, multi-core preservation and composed batch interface."""

import copy
import hashlib
import json
import math

from Bio.SeqUtils import molecular_weight
import pytest

from mhctools import assess_vaccine_windows, vaccine_window_rows, format_vaccine_window_report, vaccine_route_steps
from mhctools.cli.script import main


def profile(sequence, delivery="slp", changes=None, cut_bond=None):
    rates = {step.name: {"value_per_hour": 0, "basis": "disabled", "source": "Synthetic test exclusion"}
             for step in vaccine_route_steps(delivery, "I")}
    for name, value in (changes or {}).items():
        rates[name] = {"value_per_hour": value, "basis": "assumed", "source": "Synthetic test rate, not human kinetics"}
    compartments = ("interstitium", "afferent_lymph", "node_fluid", "blood", "producer_cytosol", "producer_er") + tuple(
        region + "_" + phase for region in ("local", "node", "systemic") for phase in ("endosome", "cytosol", "er"))
    cuts = []
    for start in ([0, cut_bond] if cut_bond is not None else [0]):
        for compartment in compartments:
            channels = ([{"bond": cut_bond, "enzyme": "synthetic protease", "rate_per_hour": 1,
                          "basis": "assumed", "source": "Analytical competition oracle"}]
                        if cut_bond is not None and start == 0 and compartment == "interstitium" else [])
            cuts.append({"compartment": compartment, "start": start, "end": len(sequence), "sequence": sequence[start:],
                         "channels": channels, "basis": "assumed", "source": "Synthetic cleavage oracle"})
    return {"I": {"rates": rates, "tap_max_length": len(sequence), "cleavage_assessments": cuts}}


def request(sequence="ACDEFGHIKLMN"):
    return {"schema_version": 1, "scenario": "Synthetic validation only", "times_hours": [0, 1, 2, 4, 60],
            "administration": "IM", "candidates": [{"id": "window1", "peptide": sequence, "delivery": "slp",
                "cores": [{"label": "core1", "start": 2, "end": 5, "allele": "HLA-A*02:01", "mhc_class": "I", "score_contribution": 1},
                          {"label": "core2", "start": 7, "end": 10, "allele": "HLA-A*02:01", "mhc_class": "I", "score_contribution": 3}],
                "kinetics": profile(sequence, changes={"lymph_entry": 1, "lymph_transit": 1, "node_antigen_uptake": 1})}]}


def test_one_microgram_uses_actual_window_mass_and_marginal_fractions():
    data = request()
    other = copy.deepcopy(data["candidates"][0])
    other.update(id="longer", peptide="ACDEFGHIKLMNPQRST", kinetics=profile(
        "ACDEFGHIKLMNPQRST", changes={"lymph_entry": 1, "lymph_transit": 1, "node_antigen_uptake": 1}))
    data["candidates"].append(other)
    out = assess_vaccine_windows(data)
    rows = vaccine_window_rows(out)
    assert len(rows) == 4
    for window in out["candidates"]:
        assert window["dose"]["amount"] == 1
        expected = 6.02214076e17 / molecular_weight(window["peptide_input"]["sequence"], seq_type="protein")
        assert window["initial_copies"] == pytest.approx(expected)
        for core in window["cores"]:
            assert core["exposure"]["node_apc_entry"]["fraction_of_dosed_peptide_copies"] == pytest.approx(1)
    assert rows[0]["node_apc_entry_copies"] > rows[2]["node_apc_entry_copies"]
    assert rows[0]["intrinsic_circulation_half_life_hours"] is None


def test_flank_cut_preserves_one_core_but_destroys_another_without_union_double_count():
    data = request()
    data["candidates"][0]["kinetics"] = profile("ACDEFGHIKLMN", cut_bond=4,
        changes={"lymph_entry": 1, "lymph_transit": 1, "node_antigen_uptake": 1})
    window = assess_vaccine_windows(data)["candidates"][0]
    assert [c["exposure"]["node_apc_entry"]["copies_per_input_molecule"] for c in window["cores"]] == pytest.approx([.5, 1])
    assert window["aggregate"]["score_weighted_mean_node_entry_per_input"] == pytest.approx(.875)
    assert window["aggregate"]["any_core_molecule_fraction"] is None
    assert window["cores"][1]["sequence"] == "IKL"


def test_rna_needs_its_own_dose_and_can_produce_multiple_antigen_copies():
    data = request()
    candidate = data["candidates"][0]
    candidate.update(delivery="cytosolic_mrna", kinetics=profile("ACDEFGHIKLMN", "cytosolic_mrna", {
        "carrier_drainage": 1, "node_rna_uptake": 1, "node_rna_escape": 1, "node_translation": 5, "node_rna_decay": 1}))
    missing = assess_vaccine_windows(data)["candidates"][0]
    assert missing["dose"] is None
    assert missing["status"] == "unassessed"
    candidate["dose"] = {"amount": 1, "unit": "ug", "basis": "assumed", "source": "Synthetic transcript mass", "molecular_weight_g_per_mol": 1e6}
    rows = vaccine_window_rows(assess_vaccine_windows(data))
    assert rows[0]["node_apc_entries_per_input"] == pytest.approx(5)
    assert rows[0]["node_apc_fraction_of_dosed_peptide"] is None
    assert rows[0]["initial_copies"] == pytest.approx(6.02214076e11)


def test_missing_core_kinetics_never_become_zero_or_partial_weighted_score():
    data = request()
    data["candidates"][0]["kinetics"]["I"]["rates"].pop("lymph_entry")
    out = assess_vaccine_windows(data)
    window = out["candidates"][0]
    assert window["status"] == "unassessed"
    assert window["aggregate"]["score_weighted_mean_node_entry_per_input"] is None
    assert all(r["node_apc_entry_copies"] is None for r in vaccine_window_rows(out))
    assert "lymph_entry" in format_vaccine_window_report(out)
    del data["candidates"][0]["kinetics"]
    assert "tap_max_length" in assess_vaccine_windows(data)["candidates"][0]["cores"][0]["missing_inputs"][0]


def test_repeated_sequence_occurrences_and_alleles_are_not_collapsed():
    data = request("ACDACDACDACD")
    cores = [{"label": str(i), "start": i, "end": i+3, "allele": "HLA-A*02:01", "mhc_class": "I"} for i in (0, 3)]
    cores.append({**cores[0], "label": "other allele", "allele": "HLA-A*03:01"})
    data["candidates"][0]["cores"] = cores
    result = assess_vaccine_windows(data)
    assert len(vaccine_window_rows(result)) == 3
    assert result["candidates"][0]["aggregate"]["status"] == "unassessed"
    cores.append(dict(cores[0]))
    with pytest.raises(ValueError, match="Duplicate core"):
        assess_vaccine_windows(data)


def test_modified_input_is_retained_and_unassessed_not_silently_unmodified():
    data = request()
    data["candidates"][0]["peptide"] = {"sequence": "ACDEFGHIKLMN", "c_term": "amidated"}
    window = assess_vaccine_windows(data)["candidates"][0]
    assert window["peptide_input"]["c_term"] == "amidated"
    assert window["dose"] is None
    assert window["cores"][0]["missing_inputs"][0] == "unsupported C-terminal chemistry: amidated"


@pytest.mark.parametrize("value", [math.nan, math.inf, True])
def test_invalid_score_rejected(value):
    data = request()
    data["candidates"][0]["cores"][0]["score_contribution"] = value
    with pytest.raises(ValueError, match="score_contribution"):
        assess_vaccine_windows(data)


def test_cli_composes_audit_csv_and_report_and_preserves_outputs(tmp_path):
    source, output, table, report = [tmp_path / name for name in ("input.json", "out.json", "out.csv", "out.md")]
    source.write_text(json.dumps(request()))
    args = ["vaccine-windows", "--input", str(source), "--output", str(output), "--table", str(table), "--report", str(report)]
    main(args)
    out = json.loads(output.read_text())
    assert out["input_file"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert len(table.read_text().splitlines()) == 3
    assert "conditional_scenario" in report.read_text()
    before = output.read_bytes()
    with pytest.raises(SystemExit):
        main(args)
    assert output.read_bytes() == before
    output.unlink()
    with pytest.raises(SystemExit):
        main(args[:-1] + [str(table)])
    assert not output.exists()


def test_core_sequence_must_match_explicit_occurrence():
    data = request()
    data["candidates"][0]["cores"][0]["sequence"] = "AAA"
    with pytest.raises(ValueError, match="Core sequence"):
        assess_vaccine_windows(data)
