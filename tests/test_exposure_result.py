"""Contracts for PK, uptake, and tissue-exposure result semantics."""

import json

import pytest

from mhctools.pred import (
    CONTEXT_DEPENDENT_KINDS,
    Kind,
    MeasurementContext,
    PeptideResult,
    Prediction,
    best_direction,
)


def _context(**updates):
    values = {
        "estimate_type": "ml_predicted",
        "analyte": "free parent peptide",
        "compartment": "plasma",
        "matrix": "plasma",
        "unit": "hours",
        "transform": "linear",
        "pk_scope": "systemic",
    }
    values.update(updates)
    return MeasurementContext(**values)


def test_endpoint_kinds_are_distinct_and_no_ambiguous_delivery_score_exists():
    kinds = {
        Kind.systemic_elimination_half_life,
        Kind.systemic_clearance,
        Kind.distribution_volume,
        Kind.systemic_exposure,
        Kind.cpp_classification,
        Kind.cellular_uptake,
        Kind.tissue_concentration,
        Kind.serum_half_life,
        Kind.plasma_half_life,
        Kind.blood_half_life,
        Kind.pMHC_stability,
    }
    assert len(kinds) == 11
    assert not hasattr(Kind, "delivery_score")


def test_measurement_context_json_and_prediction_round_trip():
    context = _context()
    prediction = Prediction(
        kind=Kind.systemic_elimination_half_life,
        score=0.72,
        value=4.5,
        peptide="SIINFEKL",
        predictor_name="fixture",
        predictor_version="model@sha256:abc",
        measurement_context=context,
    )

    encoded = json.loads(json.dumps(prediction.to_dict()))
    restored = Prediction.from_dict(encoded)

    assert restored == prediction
    assert restored.measurement_context.schema_version == 1
    assert restored.measurement_context.pk_scope == "systemic"
    assert restored.measurement_context.unit == "hours"


def test_time_series_identity_and_dataframe_preserve_every_occurrence():
    first = MeasurementContext(
        estimate_type="observed",
        analyte="total parent peptide",
        compartment="liver",
        matrix="tissue homogenate",
        unit="ng/g",
        transform="linear",
        concentration_basis="total",
        timepoint=1.0,
        time_unit="hours",
        time_origin="dose administration",
        series_id="study-7:subject-2:liver:parent",
    )
    second = MeasurementContext.from_dict({
        **first.to_dict(),
        "timepoint": 4.0,
    })
    result = PeptideResult(preds=(
        Prediction(
            kind=Kind.tissue_concentration,
            score=12.0,
            value=12.0,
            peptide="SIINFEKL",
            offset=3,
            source_sequence_name="occurrence-a",
            measurement_context=first,
        ),
        Prediction(
            kind=Kind.tissue_concentration,
            score=8.0,
            value=8.0,
            peptide="SIINFEKL",
            offset=3,
            source_sequence_name="occurrence-a",
            measurement_context=second,
        ),
    ))

    restored = PeptideResult.from_dict(
        json.loads(json.dumps(result.to_dict())))
    frame = restored.to_dataframe(sample_name="sample")

    assert len(restored.preds) == 2
    assert len(frame) == 2
    assert frame["source_sequence_name"].tolist() == [
        "occurrence-a", "occurrence-a"]
    contexts = frame["measurement_context"].tolist()
    assert [c["timepoint"] for c in contexts] == [1.0, 4.0]
    assert {c["series_id"] for c in contexts} == {
        "study-7:subject-2:liver:parent"}
    assert all(pred.allele == "" for pred in restored.preds)


@pytest.mark.parametrize("status", [
    "unsupported", "missing", "out_of_domain", "failed",
])
def test_unavailable_results_are_explicit_and_carry_no_number(status):
    context = _context(
        status=status,
        unit=None,
        transform=None,
        detail="fixture reason",
    )
    prediction = Prediction(
        kind=Kind.systemic_clearance,
        score=None,
        peptide="SIINFEKL",
        measurement_context=context,
    )

    assert prediction.measurement_context.status == status
    assert prediction.score is None
    assert prediction.value is None


def test_unavailable_result_rejects_stale_numeric_output():
    with pytest.raises(ValueError, match="cannot carry score or value"):
        Prediction(
            kind=Kind.cellular_uptake,
            score=0.8,
            measurement_context=_context(
                status="out_of_domain", unit=None, transform=None),
        )


def test_value_requires_explicit_linear_unit_and_transform():
    with pytest.raises(ValueError, match="unit is required"):
        Prediction(
            kind=Kind.systemic_clearance,
            score=1.0,
            value=1.0,
            measurement_context=_context(unit=None),
        )
    with pytest.raises(ValueError, match="linear transform"):
        Prediction(
            kind=Kind.systemic_clearance,
            score=1.0,
            value=1.0,
            measurement_context=_context(transform="log10"),
        )


def test_partial_time_series_identity_is_rejected():
    with pytest.raises(ValueError, match="must be provided together"):
        _context(timepoint=1.0, time_unit="hours")


def test_mhc_independent_endpoint_rejects_allele_duplication():
    with pytest.raises(ValueError, match="MHC-independent"):
        Prediction(
            kind=Kind.cellular_uptake,
            score=0.4,
            value=0.4,
            allele="HLA-A*02:01",
            measurement_context=_context(
                unit="relative fluorescence",
                score_semantics="model-native uptake score",
            ),
        )


def test_context_dependent_endpoints_have_no_implicit_best_direction():
    for kind in CONTEXT_DEPENDENT_KINDS:
        for field in ("score", "value", "percentile_rank"):
            with pytest.raises(ValueError, match="context-dependent"):
                best_direction(kind, field)


def test_mixed_kind_aggregation_keeps_legacy_policy_but_refuses_pk_policy():
    result = PeptideResult(preds=(
        Prediction(
            kind=Kind.pMHC_affinity,
            score=0.8,
            value=25.0,
            allele="HLA-A*02:01",
        ),
        Prediction(
            kind=Kind.systemic_exposure,
            score=30.0,
            value=30.0,
            measurement_context=_context(unit="ng*h/mL"),
        ),
    ))

    assert result.best_by_value(Kind.pMHC_affinity).value == 25.0
    with pytest.raises(ValueError, match="context-dependent"):
        result.best_by_value(Kind.systemic_exposure)


def test_cpp_confidence_is_not_a_physical_value():
    prediction = Prediction(
        kind=Kind.cpp_classification,
        score=0.91,
        value=None,
        peptide="RKKRRQRRR",
        measurement_context=MeasurementContext(
            estimate_type="ml_predicted",
            score_semantics="confidence for class_label",
            class_label="cell_penetrating",
        ),
    )

    assert prediction.score == 0.91
    assert prediction.value is None
    assert prediction.measurement_context.class_label == "cell_penetrating"

    with pytest.raises(ValueError, match="belongs in score"):
        Prediction(
            kind=Kind.cpp_classification,
            score=0.91,
            value=91.0,
            measurement_context=MeasurementContext(
                estimate_type="ml_predicted",
                unit="percent",
                transform="linear",
                score_semantics="confidence for class_label",
                class_label="cell_penetrating",
            ),
        )


def test_quantitative_endpoints_cannot_be_reduced_to_native_scores():
    with pytest.raises(ValueError, match="requires a quantitative value"):
        Prediction(
            kind=Kind.cellular_uptake,
            score=0.83,
            measurement_context=MeasurementContext(
                estimate_type="ml_predicted",
                score_semantics="model-native regression output",
            ),
        )


def test_cpp_classification_requires_label_and_confidence_semantics():
    with pytest.raises(ValueError, match="class_label and score_semantics"):
        Prediction(
            kind=Kind.cpp_classification,
            score=0.83,
            measurement_context=MeasurementContext(
                estimate_type="ml_predicted",
            ),
        )
