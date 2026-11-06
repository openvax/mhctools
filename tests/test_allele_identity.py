"""The shared, independently specified allele identity contract for #631."""

import json
from pathlib import Path

import pytest
from mhcgnomes import AlleleWithoutGene

from mhctools.allele_normalization import (
    AlleleIdentityError, allele_identity_provenance, native_allele_identity,
    normalize_allele_name, normalize_allele_name_or_raw, parse_allele,
)
from mhctools.allele_formats import format_allele, netmhc_allele_name
from mhctools.base_predictor import BasePredictor
from mhctools.cli.args import make_mhc_arg_parser, mhc_alleles_from_args

CONTRACT = json.loads((Path(__file__).parent / 'data' /
                       'canonical_allele_identities.json').read_text())


@pytest.mark.parametrize('raw,canonical', CONTRACT['accepted'])
def test_lossless_canonical_round_trip(raw, canonical):
    for normalizer in (normalize_allele_name, normalize_allele_name_or_raw):
        assert normalizer(raw) == canonical
        assert normalizer(canonical) == canonical
    assert parse_allele(raw).to_string() == canonical


@pytest.mark.parametrize('raw', CONTRACT['rejected'])
def test_ambiguity_and_invalid_pairs_cannot_escape_validation(raw):
    for normalizer in (normalize_allele_name, normalize_allele_name_or_raw):
        with pytest.raises(AlleleIdentityError):
            normalizer(raw)
    with pytest.raises(AlleleIdentityError):
        BasePredictor._check_hla_alleles([raw], keep_unparseable=True)
    with pytest.raises(AlleleIdentityError):
        format_allele(raw, netmhc_allele_name, fallback=lambda n: n)
    args = make_mhc_arg_parser().parse_args(['--mhc-predictor', 'netmhcpan', '--mhc-alleles', raw])
    with pytest.raises(AlleleIdentityError):
        mhc_alleles_from_args(args)


def test_missing_genes_do_not_invent_identity():
    assert isinstance(parse_allele('BoLA-T2b'), AlleleWithoutGene)
    assert isinstance(parse_allele('BoLA-T2c'), AlleleWithoutGene)


@pytest.mark.parametrize('keyword', ['infer_class2_pair', 'omit_dra1'])
def test_legacy_identity_projection_requires_explicit_migration(keyword):
    assert normalize_allele_name('HLA-DRB1*01:01', **{keyword: False}) == 'HLA-DRB1*01:01'
    with pytest.raises(AlleleIdentityError, match='explicit'):
        normalize_allele_name('HLA-DRB1*01:01', **{keyword: True})


def test_curated_aliases_are_scoped_to_native_catalogs():
    assert native_allele_identity('Mamu-A01') == 'Mamu-A1*001:01'
    assert native_allele_identity('HLA-C*04:22') == 'HLA-C*04:22'
    for name in ['A2', 'HLA-A*02:01P', 'H-2-Qa1', 'Mamu-B17']:
        with pytest.raises(AlleleIdentityError):
            native_allele_identity(name)


def test_unparseable_raw_fallback_preserves_case_and_star_matching():
    # A synthetic tool-catalog entry, not a parsed gene or group.
    assert normalize_allele_name_or_raw('  native-X*foo  ') == 'native-Xfoo'
    assert normalize_allele_name_or_raw('native-xfoo') == 'native-xfoo'


def test_cache_provenance_names_policy_parser_and_options():
    provenance = allele_identity_provenance()
    assert provenance['policy'] == CONTRACT['policy']
    assert provenance['mhcgnomes']
    assert provenance['parser_options']['infer_class2_pairing'] is False
    assert provenance['parser_options']['use_allele_aliases'] is False



def test_native_catalog_cache_cannot_mix_catalogs_or_parser_versions(monkeypatch):
    import mhctools.base_commandline_predictor as module
    cls = module.BaseCommandlinePredictor
    first = cls._normalized_supported_alleles('catalog-boundary', '-list', {'HLA-A0201'})
    second = cls._normalized_supported_alleles('catalog-boundary', '-list', {'HLA-B0702'})
    assert first == {'HLA-A*02:01': 'HLA-A0201'}
    assert second == {'HLA-B*07:02': 'HLA-B0702'}
    original = module.allele_identity_provenance
    monkeypatch.setattr(module, 'allele_identity_provenance',
                        lambda: dict(original(), mhcgnomes='future-parser'))
    third = cls._normalized_supported_alleles('catalog-boundary', '-list', {'HLA-A0201'})
    assert third == first
    assert third is not first


def test_native_catalog_alias_and_excluded_entries_do_not_escape_to_users(caplog):
    from mhctools.base_commandline_predictor import BaseCommandlinePredictor
    mapping = BaseCommandlinePredictor._normalized_supported_alleles(
        'mixed-catalog', '-list', {'Mamu-A01', 'Mamu-B17', 'H-2-Qa1', 'HLA-A0201', 'BoLA-T2b'})
    assert mapping == {'Mamu-A1*001:01': 'Mamu-A01', 'HLA-A*02:01': 'HLA-A0201',
                       'BoLA-T2b': 'BoLA-T2b'}
    assert 'excluded 2' in caplog.text
