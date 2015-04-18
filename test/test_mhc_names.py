from mhctools.common import (
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

def test_long_names():
    expected = "HLA-A*02:01"
    for name in hla_alleles:
        result = normalize_allele_name(name)
        assert expected == result, result

def test_short_names():
    expected = "A0201"
    for name in hla_alleles:
        result = compact_allele_name(name)
        assert expected == result, result
