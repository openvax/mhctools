from nose.tools import eq_
from mhctools.alleles import (
    normalize_allele_name,
    compact_allele_name
)

hla_alleles = [
    "HLA-A*02:01",
    "HLA-A*0201",
    "A*02:01",
    "A*0201",
    "HLA-A02:01",
    "A0201",
    "HLA-A0201",
    "A2",
    "A2:01",
    "HLA-A2",
    "A0201",
]

def test_hla_long_names():
    expected = "HLA-A*02:01"
    for name in hla_alleles:
        result = normalize_allele_name(name)
        eq_(expected, result)

def test_hla_short_names():
    expected = "A0201"
    for name in hla_alleles:
        result = compact_allele_name(name)
        eq_(expected, result)

def test_macaque_alleles():
    allele_name = "Mamu-B*082:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*82:02")
    eq_(compact_allele_name(allele_name), "B8202")

    # expect 3rd zero in the family "007" to be trimmed in the normalized form
    # of this allele
    allele_name = "Mamu-B*007:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*07:02")
    eq_(compact_allele_name(allele_name), "B0702")

def test_swine_alleles():
    pass
    """
    allele_name = "Mamu-B*082:02"
    eq_(normalize_allele_name(allele_name), allele_name)
    eq_(compact_allele_name(allele_name), "Mamu-B8202")
    """