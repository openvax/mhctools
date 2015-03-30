from mhctools.common import (normalize_hla_allele_name,
                              compact_hla_allele_name,
                              mhc_class_from_normalized_allele_name)

from nose.tools import assert_equals, raises

hla_alleles = [
	'HLA-A*02:01',
	'HLA-A*0201',
	'A*02:01',
	'A*0201',
	'HLA-A02:01',
	'A0201',
	'HLA-A0201',
	'A2',
	'A2:01',
	'HLA-A2',
    'A0201',
]

def test_long_names():
	expected = 'HLA-A*02:01'
	for name in hla_alleles:
		result = normalize_hla_allele_name(name)
		assert expected == result, result

def test_short_names():
	expected = 'A0201'
	for name in hla_alleles:
		result = compact_hla_allele_name(name)
		assert expected == result, result

@raises(AssertionError)
def test_mhc_from_not_normalized():
	mhc_class = mhc_class_from_normalized_allele_name('HLAA*03:02')

def test_mhc_from_class_1():
	assert_equals(mhc_class_from_normalized_allele_name('HLA-A*03:02'), 1)
	assert_equals(mhc_class_from_normalized_allele_name('HLA-K*03:02'), 1)

@raises(ValueError)
def test_mhc_class_1_error():
	mhc_class = mhc_class_from_normalized_allele_name('HLA-AA*03:02')

def test_mhc_from_class_2():
	assert_equals(mhc_class_from_normalized_allele_name('HLA-DM*03:02'), 2)
	assert_equals(mhc_class_from_normalized_allele_name('HLA-DQ*03:02'), 2)

@raises(ValueError)
def test_mhc_class_2_error():
	mhc_class = mhc_class_from_normalized_allele_name('HLA-DS*03:02')
