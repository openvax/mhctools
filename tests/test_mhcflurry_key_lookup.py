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

from mhctools import MHCflurry, MHCflurry_Affinity


def _make_fake_predictor(
        aff_allele_str,
        pres_allele_str,
        supported,
        percent_rank_transforms=None,
        allele_to_sequence=None):
    """Build a fake mhcflurry Class1PresentationPredictor with configurable
    allele string in the affinity and presentation outputs."""
    predict_calls = []

    def predict_to_dataframe(
            peptides, alleles, include_percentile_ranks=True):
        data = {
            "peptide": peptides,
            "allele": [aff_allele_str] * len(peptides),
            "prediction": [500.0] * len(peptides),
        }
        if include_percentile_ranks:
            data["prediction_percentile"] = [1.5] * len(peptides)
        return pd.DataFrame(data)

    affinity_predictor = types.SimpleNamespace(
        predict_to_dataframe=predict_to_dataframe,
        supported_alleles=supported,
    )
    if percent_rank_transforms is not None:
        affinity_predictor.allele_to_percent_rank_transform = \
            percent_rank_transforms
    if allele_to_sequence is not None:
        affinity_predictor.allele_to_sequence = allele_to_sequence
        affinity_predictor.canonicalize_allele_name = lambda allele: allele

    def predict(
            peptides,
            alleles,
            sample_names=None,
            n_flanks=None,
            c_flanks=None,
            include_affinity_percentile=False,
            verbose=0,
            throw=True,
            affinity_model_kwargs=None,
            processing_batch_size="auto"):
        predict_calls.append({
            "peptides": list(peptides),
            "alleles": list(alleles),
            "n_flanks": None if n_flanks is None else list(n_flanks),
            "c_flanks": None if c_flanks is None else list(c_flanks),
        })
        return pd.DataFrame({
            "peptide": list(peptides),
            "allele": [pres_allele_str] * len(peptides),
            "presentation_score": [0.75] * len(peptides),
            "presentation_percentile": [2.0] * len(peptides),
        })
    return types.SimpleNamespace(
        affinity_predictor=affinity_predictor,
        predict=predict,
        predict_calls=predict_calls,
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


def test_accepts_affinity_percentile_calibration_from_same_pseudosequence():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-C*15:05",
        pres_allele_str="HLA-C*15:05",
        supported=["HLA-C*15:05"],
        percent_rank_transforms={"HLA-C*15:99": object()},
        allele_to_sequence={
            "HLA-C*15:05": "PSEUDOSEQ",
            "HLA-C*15:99": "PSEUDOSEQ",
        })
    predictor = MHCflurry(alleles=["HLA-C*15:05"], predictor=fake)
    results = predictor.predict(["SIINFEKLA"])
    assert results[0].affinity.percentile_rank == 1.5


def test_missing_affinity_percentile_calibration_raises_early():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-C*15:05",
        pres_allele_str="HLA-C*15:05",
        supported=["HLA-C*15:05"],
        percent_rank_transforms={},
        allele_to_sequence={"HLA-C*15:05": "PSEUDOSEQ"})
    with pytest.raises(ValueError, match="affinity percentile ranks"):
        MHCflurry(alleles=["HLA-C*15:05"], predictor=fake)


def test_can_disable_missing_affinity_percentile_ranks():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-C*15:05",
        pres_allele_str="HLA-C*15:05",
        supported=["HLA-C*15:05"],
        percent_rank_transforms={},
        allele_to_sequence={"HLA-C*15:05": "PSEUDOSEQ"})
    predictor = MHCflurry(
        alleles=["HLA-C*15:05"],
        predictor=fake,
        include_affinity_percentile_ranks=False)
    results = predictor.predict(["SIINFEKLA"])
    assert results[0].affinity.percentile_rank is None


def test_mhcflurry_forwards_flanks_to_presentation_predictor():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",
        pres_allele_str="HLA-A*02:01",
        supported=["HLA-A*02:01"])
    predictor = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)

    results = predictor.predict(
        ["SIINFEKLA", "SIINFEKLA"],
        n_flanks=["NN", "XX"],
        c_flanks=["CC", "YY"])

    assert fake.predict_calls[0]["n_flanks"] == ["NN", "XX"]
    assert fake.predict_calls[0]["c_flanks"] == ["CC", "YY"]
    assert len(results) == 2
    assert results[0].presentation.n_flank == "NN"
    assert results[0].presentation.c_flank == "CC"
    assert results[1].presentation.n_flank == "XX"
    assert results[1].presentation.c_flank == "YY"


def test_mhcflurry_predict_proteins_threads_flanking_context():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",
        pres_allele_str="HLA-A*02:01",
        supported=["HLA-A*02:01"])
    predictor = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)

    result = predictor.predict_proteins({"protein": "MSIINFEKLAC"})

    assert fake.predict_calls[0]["peptides"] == [
        "MSIINFEKL",
        "SIINFEKLA",
        "IINFEKLAC",
    ]
    assert fake.predict_calls[0]["n_flanks"] == ["", "M", "MS"]
    assert fake.predict_calls[0]["c_flanks"] == ["AC", "C", ""]
    middle = result["protein"][1]
    assert middle.peptide == "SIINFEKLA"
    assert middle.presentation.n_flank == "M"
    assert middle.presentation.c_flank == "C"


def test_affinity_only_disabled_percentile_ranks_convert_to_none():
    def predict_to_dataframe(
            peptides, alleles, include_percentile_ranks=True):
        data = {
            "peptide": list(peptides),
            "allele": list(alleles),
            "prediction": [500.0] * len(peptides),
        }
        if include_percentile_ranks:
            data["prediction_percentile"] = [1.5] * len(peptides)
        return pd.DataFrame(data)

    fake = types.SimpleNamespace(
        predict_to_dataframe=predict_to_dataframe,
        supported_alleles=["HLA-A*02:232"],
    )
    predictor = MHCflurry_Affinity(
        alleles=["HLA-A*02:232"],
        predictor=fake,
        include_affinity_percentile_ranks=False)

    result = predictor.predict(["SIINFEKLA"])[0]

    assert result.affinity.percentile_rank is None
    assert result.best_affinity_by_rank is None
