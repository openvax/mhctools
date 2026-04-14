# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the MHCflurry wrapper's presentation-score lookup.

Verifies that when the presentation predictor's output allele string differs
from the affinity predictor's output allele string (e.g. different
normalization), the wrapper fails loudly instead of silently returning
score=0.0 for all presentations.
"""

import types

import pandas as pd
import pytest

from mhctools import MHCflurry


def _make_fake_predictor(aff_allele_str, pres_allele_str, supported):
    """Build a fake mhcflurry Class1PresentationPredictor with configurable
    allele string in the affinity and presentation outputs."""
    affinity_predictor = types.SimpleNamespace(
        predict_to_dataframe=lambda peptides, alleles: pd.DataFrame({
            "peptide": peptides,
            "allele": [aff_allele_str] * len(peptides),
            "prediction": [500.0] * len(peptides),
            "prediction_percentile": [1.5] * len(peptides),
        }),
        supported_alleles=supported,
    )
    def predict(peptides, alleles, include_affinity_percentile=False, verbose=0):
        return pd.DataFrame({
            "peptide": list(peptides),
            "allele": [pres_allele_str] * len(peptides),
            "presentation_score": [0.75] * len(peptides),
            "presentation_percentile": [2.0] * len(peptides),
        })
    return types.SimpleNamespace(
        affinity_predictor=affinity_predictor,
        predict=predict,
        supported_alleles=supported,
    )


def test_consistent_allele_strings_produce_correct_scores():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",
        pres_allele_str="HLA-A*02:01",
        supported=["HLA-A*02:01"])
    p = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)
    results = p.predict(["SIINFEKLA"])
    assert len(results) == 1
    r = results[0]
    assert r.presentation.score == 0.75
    assert r.presentation.allele == "HLA-A*02:01"


def test_inconsistent_allele_strings_raise_instead_of_silently_returning_zero():
    """Regression: previously the wrapper silently returned (0.0, None) for
    presentation when the allele string in aff_df didn't match the key in
    pres_by_pep_allele. Now it raises."""
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",    # affinity output
        pres_allele_str="HLA-A0201",     # presentation output (different!)
        supported=["HLA-A*02:01"])
    p = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)
    with pytest.raises(ValueError, match="missing presentation score"):
        p.predict(["SIINFEKLA"])
