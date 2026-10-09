# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the MHCflurry wrapper's presentation and affinity bookkeeping."""

import types

import pandas as pd
import pytest

from mhctools import (
    MHC_CLASS_VALUES,
    MHC_DEPENDENCE_VALUES,
    MHCflurry,
    MHCflurry_Affinity,
)
from mhctools.pred import Kind


def _make_fake_predictor(
        aff_allele_str,
        pres_allele_str,
        supported,
        percent_rank_transforms=None,
        allele_to_sequence=None):
    """Build a fake mhcflurry Class1PresentationPredictor with configurable
    allele string in the affinity and presentation outputs."""
    predict_calls = []
    affinity_calls = []

    def predict_to_dataframe(
            peptides, alleles, include_percentile_ranks=True):
        affinity_calls.append({
            "peptides": list(peptides), "alleles": list(alleles),
            "include_percentile_ranks": include_percentile_ranks,
        })
        output_alleles = (
            [aff_allele_str] * len(peptides)
            if aff_allele_str is not None else list(alleles))
        data = {
            "peptide": peptides,
            "allele": output_alleles,
            "prediction": [500.0] * len(peptides),
        }
        if include_percentile_ranks:
            data["prediction_percentile"] = [1.5] * len(peptides)
        return pd.DataFrame(data)

    affinity_predictor = types.SimpleNamespace(
        predict_to_dataframe=predict_to_dataframe,
        supported_alleles=supported,
        affinity_calls=affinity_calls,
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
        allele_arg = {
            k: list(v)
            for (k, v) in alleles.items()
        } if isinstance(alleles, dict) else list(alleles)
        predict_calls.append({
            "peptides": list(peptides),
            "alleles": allele_arg,
            "n_flanks": None if n_flanks is None else list(n_flanks),
            "c_flanks": None if c_flanks is None else list(c_flanks),
        })
        rows = []
        if isinstance(alleles, dict):
            for sample_name in alleles:
                for i, peptide in enumerate(peptides):
                    rows.append({
                        "peptide": peptide,
                        "peptide_num": i,
                        "sample_name": sample_name,
                        "best_allele": sample_name,
                        "presentation_score": 0.75,
                        "presentation_percentile": 2.0,
                    })
        else:
            for i, peptide in enumerate(peptides):
                rows.append({
                    "peptide": peptide,
                    "peptide_num": i,
                    "best_allele": pres_allele_str,
                    "presentation_score": 0.75,
                    "presentation_percentile": 2.0,
                })
        return pd.DataFrame(rows)
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


def test_presentation_best_allele_is_attribution_not_affinity_lookup_key():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",    # affinity output
        pres_allele_str="HLA-A0201",     # presentation attribution
        supported=["HLA-A*02:01"])
    p = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)
    result = p.predict(["SIINFEKLA"])[0]
    assert result.presentation.allele == "HLA-A0201"
    assert result.presentation.score == 0.75


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


def test_mhcflurry_kind_support_marks_presentation_as_haplotype_level():
    fake = _make_fake_predictor(
        aff_allele_str="HLA-A*02:01",
        pres_allele_str="HLA-A*02:01",
        supported=["HLA-A*02:01"])
    predictor = MHCflurry(alleles=["HLA-A*02:01"], predictor=fake)

    support = predictor.kind_support()

    assert support[Kind.pMHC_affinity]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_presentation]["mhc_dependence"] == "haplotype"
    assert support[Kind.pMHC_presentation]["mhc_class"] == "I"
    assert support[Kind.pMHC_presentation]["mhc_dependence"] in (
        MHC_DEPENDENCE_VALUES)
    assert support[Kind.pMHC_presentation]["mhc_class"] in MHC_CLASS_VALUES


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


def test_mhcflurry_presentation_uses_full_haplotype_once():
    fake = _make_fake_predictor(
        aff_allele_str=None,
        pres_allele_str="HLA-B*07:02",
        supported=["HLA-A*02:01", "HLA-B*07:02"])
    predictor = MHCflurry(
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        predictor=fake)

    result = predictor.predict(["SIINFEKLA"])[0]

    assert len(fake.predict_calls) == 1
    assert set(fake.predict_calls[0]["alleles"]) == {
        "HLA-A*02:01",
        "HLA-B*07:02",
    }
    assert len(result.filter(kind=Kind.pMHC_affinity)) == 2
    assert len(result.filter(kind=Kind.pMHC_presentation)) == 1
    assert result.presentation.allele == "HLA-B*07:02"


def test_mhcflurry_large_panel_uses_one_sample_per_allele():
    alleles = [
        "HLA-A*01:01",
        "HLA-A*02:01",
        "HLA-A*03:01",
        "HLA-B*07:02",
        "HLA-B*08:01",
        "HLA-C*07:01",
        "HLA-C*07:02",
    ]
    fake = _make_fake_predictor(
        aff_allele_str=None,
        pres_allele_str="unused",
        supported=alleles)
    predictor = MHCflurry(alleles=alleles, predictor=fake)

    result = predictor.predict(["SIINFEKLA"])[0]

    assert predictor.presentation_allele_mode == "per_allele"
    assert predictor.kind_support()[Kind.pMHC_presentation][
        "mhc_dependence"] == "single_allele"
    assert len(fake.predict_calls) == 1
    presentation_alleles = fake.predict_calls[0]["alleles"]
    assert set(presentation_alleles) == set(alleles)
    assert all(v == [k] for k, v in presentation_alleles.items())
    assert len(result.filter(kind=Kind.pMHC_affinity)) == len(alleles)
    presentations = result.filter(kind=Kind.pMHC_presentation)
    assert len(presentations) == len(alleles)
    assert {p.allele for p in presentations} == set(alleles)


def test_mhcflurry_rejects_large_explicit_haplotype():
    alleles = [
        "HLA-A*01:01",
        "HLA-A*02:01",
        "HLA-A*03:01",
        "HLA-B*07:02",
        "HLA-B*08:01",
        "HLA-C*07:01",
        "HLA-C*07:02",
    ]
    fake = _make_fake_predictor(
        aff_allele_str=None,
        pres_allele_str="unused",
        supported=alleles)

    with pytest.raises(ValueError, match="haplotype mode accepts at most"):
        MHCflurry(
            alleles=alleles,
            predictor=fake,
            presentation_allele_mode="haplotype")


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


@pytest.fixture(params=[MHCflurry, MHCflurry_Affinity])
def allele_key_wrapper(request):
    return request.param


def _key_predictor(wrapper, supported, **kwargs):
    fake = _make_fake_predictor(
        aff_allele_str=None, pres_allele_str=supported[0],
        supported=supported, **kwargs)
    return fake if wrapper is MHCflurry else fake.affinity_predictor


@pytest.mark.parametrize('requested,key', [
    ('DLA-88*001:01', 'DLA-88*001:01'),
    ('DLA-88*01:01', 'DLA-88*001:01'),
    ('dla-88*001:01', 'DLA-88*001:01'),
    ('DLA-88*003:02', 'DLA-88*003:02'),
    ('DLA-88*03:02', 'DLA-88*003:02'),
    ('DLA-88*001:01', 'DLA-88*01:01'),
    ('HLA-A0201', 'HLA-A*02:01'),
    ('hla-a*02:01', 'HLA-A*02:01'),
    ('HLA-A*02:01', 'HLA-A0201'),
    ('H2-Kb', 'H-2-Kb'),
    ('H-2-Kb', 'H2-Kb'),
    ('predictor-native-key', 'predictor-native-key'),
])
def test_resolves_actual_supported_keys(allele_key_wrapper, requested, key):
    fake = _key_predictor(allele_key_wrapper, [key])
    model = allele_key_wrapper(
        [requested], predictor=fake, include_affinity_percentile_ranks=False)
    assert model.alleles == [key]


@pytest.mark.parametrize('requested', [
    'DLA-88*999:99', 'HLA-A*99:99', 'foo', 'A2', 'HLA-A*02',
    'HLA-A*02:01:01', 'HLA-A*02:01N', 'HLA-DQA1*01:01',
])
def test_rejects_unknown_or_non_equivalent_identities(allele_key_wrapper, requested):
    from mhctools.unsupported_allele import UnsupportedAllele
    fake = _key_predictor(allele_key_wrapper, ['DLA-88*001:01', 'HLA-A*02:01'])
    with pytest.raises(UnsupportedAllele):
        allele_key_wrapper(
            [requested], predictor=fake, include_affinity_percentile_ranks=False)


def test_exact_key_wins_but_ambiguous_non_exact_identity_is_rejected(allele_key_wrapper):
    from mhctools.unsupported_allele import UnsupportedAllele
    supported = ['DLA-88*001:01', 'DLA-88*01:01']
    fake = _key_predictor(allele_key_wrapper, supported)
    for exact in supported:
        model = allele_key_wrapper(
            [exact], predictor=fake, include_affinity_percentile_ranks=False)
        assert model.alleles == [exact]
    with pytest.raises(UnsupportedAllele, match='Ambiguous'):
        allele_key_wrapper(
            ['dla-88*001:01'], predictor=fake,
            include_affinity_percentile_ranks=False)


def test_preserves_requested_order_and_forwards_resolved_keys(allele_key_wrapper):
    expected = ['DLA-88*003:02', 'H-2-Kb', 'HLA-A*02:01', 'DLA-88*001:01']
    fake = _key_predictor(allele_key_wrapper, expected)
    requested = 'DLA-88*03:02,H2-Kb,A0201,DLA-88*01:01,DLA-88*001:01'
    model = allele_key_wrapper(
        requested, predictor=fake, include_affinity_percentile_ranks=False)
    peptides = ['SIINFEKLA', 'SIINFEKLL']
    results = model.predict(peptides)
    assert model.alleles == expected
    affinity = fake.affinity_predictor if allele_key_wrapper is MHCflurry else fake
    assert affinity.affinity_calls[-1]['alleles'] == [a for a in expected for _ in peptides]
    assert affinity.affinity_calls[-1]['peptides'] == peptides * len(expected)
    assert all({p.allele for p in r.filter(kind=Kind.pMHC_affinity)} == set(expected)
               for r in results)
    if allele_key_wrapper is MHCflurry:
        assert fake.predict_calls[-1]['alleles'] == expected


def test_percent_rank_support_uses_actual_dla_key(allele_key_wrapper):
    key = 'DLA-88*001:01'
    fake = _key_predictor(
        allele_key_wrapper, [key], percent_rank_transforms={key: object()})
    model = allele_key_wrapper(['DLA-88*01:01'], predictor=fake)
    assert model.alleles == [key]
