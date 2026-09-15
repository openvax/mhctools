"""Exact peptide form, context, provenance, and cache identity contracts."""

import json

import pytest

from mhctools.optional_backend import (
    BackendSpec,
    backend_inventory,
    inspect_artifact,
    prediction_cache_key,
)
from mhctools.peptide_input import (
    PeptideContext,
    PeptideInput,
    coerce_peptide_inputs,
    sequence_only_chemistry_error,
)
from mhctools.pred import Kind, PeptideResult, Prediction


_SPEC = BackendSpec(
    name="fixture",
    endpoint="cellular_uptake",
    developed_against="fixture-v1",
    license="MIT",
    serialization="JSON",
    entry_point="prediction_only",
    supported_platforms=("linux",),
    supported_interpreters=("Python 3.10",),
)


def _inventory(tmp_path, setting="one"):
    artifact = tmp_path / f"model-{setting}.bin"
    artifact.write_bytes(setting.encode("utf-8"))
    return backend_inventory(
        _SPEC,
        [inspect_artifact("weights", "model_weights", artifact)],
        settings={"mode": setting},
    )


def _context(**updates):
    values = {
        "administration_route": "intravenous",
        "formulation": "saline",
        "cargo": "none",
        "administered_material": "intact peptide",
        "released_material": "free peptide",
        "measured_analyte": "free parent peptide",
        "study_id": "study-7",
        "assay_species": "Homo sapiens",
        "matrix": "plasma",
        "cell_type": "dendritic cell",
        "cell_subtype": "conventional type 1",
        "maturation_state": "mature",
        "readout": "LC-MS concentration",
        "timepoint": 1.0,
        "time_unit": "hours",
        "conditions": {"temperature": "37 C", "pH": "7.4"},
    }
    values.update(updates)
    return PeptideContext(**values)


def test_context_unknowns_remain_none_not_defaults():
    context = PeptideContext()
    assert context.administration_route is None
    assert context.formulation is None
    assert context.assay_species is None
    assert context.matrix is None
    assert context.cell_type is None
    assert context.timepoint is None
    assert PeptideInput("SIINFEKL").inference_identity_sha256 == \
        PeptideInput("SIINFEKL", context=context).inference_identity_sha256


def test_exact_input_json_round_trip_preserves_all_fields():
    peptide_input = PeptideInput(
        sequence="SIINFEKL",
        n_term="acetylated",
        c_term="amidated",
        attachments=(("residue:4", "FITC isomer I"),),
        occurrence_id="sample-1:occurrence-2",
        source_sequence_name="protein-9",
        source_start=12,
        source_gene="GENE1",
        source_species="Homo sapiens",
        context=_context(),
    )

    restored = PeptideInput.from_dict(
        json.loads(json.dumps(peptide_input.to_dict())))

    assert restored == peptide_input
    assert restored.schema_version == 1
    assert restored.context.schema_version == 1
    assert restored.context.conditions == (("pH", "7.4"),
                                            ("temperature", "37 C"))


def test_attachment_mapping_and_order_have_one_canonical_identity():
    mapping = PeptideInput(
        "SIINFEKL", attachments={"residue:4": "FITC", "N-term": "PEG2"})
    pairs = PeptideInput(
        "SIINFEKL",
        attachments=(("N-term", "PEG2"), ("residue:4", "FITC")))
    assert mapping == pairs
    assert mapping.chemical_identity_sha256 == pairs.chemical_identity_sha256


@pytest.mark.parametrize("changed", [
    PeptideInput("SIINFEKL", c_term="amidated"),
    PeptideInput("SIINFEKL", n_term="acetylated"),
    PeptideInput("SIINFEKL", attachments=(("residue:1", "FITC"),)),
    PeptideInput("SIINFEKL", context=_context(
        administration_route="subcutaneous")),
    PeptideInput("SIINFEKL", context=_context(matrix="serum")),
    PeptideInput("SIINFEKL", context=_context(cell_type="HeLa")),
    PeptideInput("SIINFEKL", context=_context(timepoint=4.0)),
])
def test_same_sequence_different_chemistry_or_context_has_distinct_identity(
        changed):
    baseline = PeptideInput("SIINFEKL", context=_context())
    assert changed.inference_identity_sha256 != baseline.inference_identity_sha256


def test_occurrence_and_source_provenance_do_not_poison_inference_cache(
        tmp_path):
    first = PeptideInput(
        "SIINFEKL", occurrence_id="first", source_sequence_name="protein-a",
        source_start=3, source_gene="GENE-A", source_species="mouse",
        context=_context())
    second = PeptideInput(
        "SIINFEKL", occurrence_id="second", source_sequence_name="protein-b",
        source_start=8, source_gene="GENE-B", source_species="human",
        context=_context())
    inventory = _inventory(tmp_path)

    assert first.record_identity_sha256 != second.record_identity_sha256
    assert first.inference_identity_sha256 == second.inference_identity_sha256
    assert prediction_cache_key(first, inventory) == prediction_cache_key(
        second, inventory)


def test_cache_key_includes_context_assets_and_settings(tmp_path):
    baseline = PeptideInput("SIINFEKL", context=_context())
    changed_context = PeptideInput(
        "SIINFEKL", context=_context(matrix="serum"))
    first_inventory = _inventory(tmp_path, "one")
    second_inventory = _inventory(tmp_path, "two")

    keys = {
        prediction_cache_key(baseline, first_inventory),
        prediction_cache_key(changed_context, first_inventory),
        prediction_cache_key(baseline, second_inventory),
    }
    assert len(keys) == 3


def test_prediction_round_trip_and_dataframe_preserve_exact_input_and_cache():
    peptide_input = PeptideInput(
        "SIINFEKL", occurrence_id="two", source_sequence_name="protein",
        source_start=7, context=_context())
    prediction = Prediction(
        kind=Kind.serum_half_life,
        score=2.0,
        value=2.0,
        peptide_input=peptide_input,
        cache_key="cache-identity",
    )
    result = PeptideResult(preds=(prediction, prediction))

    restored = PeptideResult.from_dict(
        json.loads(json.dumps(result.to_dict())))
    frame = restored.to_dataframe()

    assert len(restored.preds) == 2
    assert restored.preds[0].peptide == "SIINFEKL"
    assert restored.preds[0].source_sequence_name == "protein"
    assert restored.preds[0].offset == 7
    assert frame["cache_key"].tolist() == ["cache-identity"] * 2
    assert [item["occurrence_id"] for item in frame["peptide_input"]] == [
        "two", "two"]


def test_prediction_rejects_inconsistent_or_untyped_identity_fields():
    peptide_input = PeptideInput("SIINFEKL")
    with pytest.raises(ValueError, match="differs from peptide_input"):
        Prediction(
            kind=Kind.serum_half_life, score=1.0, peptide="GILGFVFTL",
            peptide_input=peptide_input)
    with pytest.raises(ValueError, match="nonempty string"):
        Prediction(
            kind=Kind.serum_half_life, score=1.0,
            peptide_input=peptide_input, cache_key=0)


def test_reordered_batch_keeps_each_stable_identity():
    inputs = [
        PeptideInput("SIINFEKL", context=_context(matrix="serum")),
        PeptideInput("SIINFEKL", context=_context(matrix="plasma")),
    ]
    assert [item.inference_identity_sha256 for item in inputs] == list(
        reversed([
            item.inference_identity_sha256 for item in reversed(inputs)]))


def test_string_shorthand_is_explicit_natural_free_peptide():
    peptide_input = coerce_peptide_inputs(" siinfekl ")[0]
    assert peptide_input == PeptideInput(
        "SIINFEKL", n_term="free", c_term="free")
    assert sequence_only_chemistry_error(peptide_input) is None


@pytest.mark.parametrize("peptide_input, message", [
    (PeptideInput("SIINFEKL", n_term="unknown"), "N-terminal"),
    (PeptideInput("SIINFEKL", c_term="amidated"), "C-terminal"),
    (PeptideInput(
        "SIINFEKL", attachments=(("residue:4", "FITC"),)), "attachments"),
])
def test_sequence_only_compatibility_rejects_unsupported_chemistry(
        peptide_input, message):
    assert message in sequence_only_chemistry_error(peptide_input)


def test_noncanonical_or_ambiguous_sequence_is_rejected_not_stripped():
    with pytest.raises(ValueError, match="non-standard residues"):
        PeptideInput("SIINFEKLB")
    with pytest.raises(ValueError, match="uppercase"):
        PeptideInput("siinfekl")


def test_timepoint_requires_units_and_context_is_descriptive_only():
    with pytest.raises(ValueError, match="provided together"):
        PeptideContext(timepoint=1.0)
    assert not hasattr(PeptideContext(), "recommended_route")
    assert not hasattr(PeptideContext(), "patient_parameters")
