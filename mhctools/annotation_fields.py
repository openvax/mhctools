"""Prediction-table field metadata shared with the lightweight CLI parser."""

from .pred import Kind


# Output-field token -> (kind, Prediction attribute). ``kind`` is ``None`` when
# the field's "best" direction is kind-independent (``score`` is always
# higher-better, ``percentile_rank`` always lower-better); ``affinity`` maps to
# the ``value`` of ``pMHC_affinity`` predictions (IC50 nM, lower-better).
OUTPUT_FIELDS = {
    "affinity": (Kind.pMHC_affinity, "value"),
    "score": (None, "score"),
    "percentile_rank": (None, "percentile_rank"),
    # convenience aliases
    "rank": (None, "percentile_rank"),
    "value": (Kind.pMHC_affinity, "value"),
    "presentation": (Kind.pMHC_presentation, "score"),
    "processing": (Kind.antigen_processing, "score"),
    # NetMHCstabpan's Thalf(h) lands in both `score` and `value`; read it from
    # `score`, which is where this field has always read it.
    "stability": (Kind.pMHC_stability, "score"),
    "immunogenicity": (Kind.immunogenicity, "score"),
    "foreignness": (Kind.foreignness, "score"),
    "tap_transport": (Kind.tap_transport, "score"),
    "erap_trimming": (Kind.erap_trimming, "score"),
    # Parent-peptide half-life. Matrix-specific names remain CLI conveniences;
    # the prediction's MeasurementContext carries the actual matrix.
    "peptide_half_life": (Kind.peptide_half_life, "value"),
    "serum_half_life": (Kind.peptide_half_life, "value"),
    "blood_half_life": (Kind.peptide_half_life, "value"),
}


# Units of the serialized field selected by each CLI token. An empty value
# means dimensionless or predictor-specific, not unknown canonical metadata.
OUTPUT_UNITS = {
    "affinity": "nM",
    "value": "nM",
    "percentile_rank": "percentile",
    "rank": "percentile",
    "stability": "hours",
    "peptide_half_life": "hours",
    "serum_half_life": "hours",
    "blood_half_life": "hours",
}


def output_field_tokens():
    """Sorted list of accepted output-field tokens (for CLI help / errors)."""
    return sorted(OUTPUT_FIELDS)

