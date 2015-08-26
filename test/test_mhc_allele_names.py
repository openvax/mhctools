from nose.tools import eq_
from mhctools.alleles import (
    normalize_allele_name,
    parse_allele_name,
    compact_allele_name,
    AlleleName,
)


# TODO: test swine and rat alleles since those
# can be significantly different from other species we've implementeds

# Rat alleles:
# - RT1-Bb*u
# - RT1-Db1*a
# - RT1-DMa*a
# - RT1-9.5*f
# - RT1-M3-1*av1

hla_02_01_names = [
    "HLA-A*02:01",
    "HLA-A*0201",
    "A*02:01",
    "A*0201",
    "HLA-A02:01",
    # no punctuation
    "A0201",
    "HLA-A0201",
    "A0201",
    "A2",
    "A2:01",
    "HLA-A2",
    # lower case
    "hla-a*0201",
    "a*0201",
    "a*02:01",
    "a0201"
]

def test_hla_long_names():
    expected = "HLA-A*02:01"
    for name in hla_02_01_names:
        result = normalize_allele_name(name)
        eq_(result, expected)

def test_hla_short_names():
    expected = "A0201"
    for name in hla_02_01_names:
        result = compact_allele_name(name)
        eq_(result, expected)

def test_macaque_alleles():
    allele_name = "Mamu-B*082:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*82:02")
    eq_(compact_allele_name(allele_name), "B8202")

    # expect 3rd zero in the family "007" to be trimmed in the normalized form
    # of this allele
    allele_name = "Mamu-B*007:02"
    eq_(normalize_allele_name(allele_name), "Mamu-B*07:02")
    eq_(compact_allele_name(allele_name), "B0702")

def test_dog_class2_allele():
    eq_(parse_allele_name("DLA-DQA1*00101"),
        AlleleName("DLA", "DQA1", "01", "01"))

def test_sheep_class1_allele():
    eq_(parse_allele_name("Ovar-N*50001"),
        AlleleName("Ovar", "N", "500", "01"))

def test_sheep_class2_allele():
    eq_(parse_allele_name("Ovar-DRB1*0804"),
        AlleleName("Ovar", "DRB1", "08", "04"))

def test_mouse_class1_alleles_H2_Kk():
    # H2-Kk
    eq_(parse_allele_name("H2-Kk"),
        AlleleName("H-2", "K", "", "k"))
    eq_(normalize_allele_name("H2-Kk"), "H-2-Kk")
    eq_(compact_allele_name("H-2-Kk"), "Kk")

    # with a hyphen in "H-2"
    eq_(parse_allele_name("H-2-Kk"),
        AlleleName("H-2", "K", "", "k"))
    eq_(normalize_allele_name("H-2-Kk"), "H-2-Kk")
    eq_(compact_allele_name("H-2-Kk"), "Kk")

def test_mouse_class1_alleles_H2_Db():
    # H2-Db
    eq_(parse_allele_name("H2-Db"),
        AlleleName("H-2", "D", "", "b"))
    eq_(normalize_allele_name("H2-Db"), "H-2-Db")
    eq_(compact_allele_name("H2-Db"), "Db")

    # with hyphen in "H-2"
    eq_(parse_allele_name("H-2-Db"),
        AlleleName("H-2", "D", "", "b"))
    eq_(normalize_allele_name("H-2-Db"), "H-2-Db")
    eq_(compact_allele_name("H-2-Db"), "Db")

def test_human_class2():
    expected = "HLA-DRB1*01:02"
    expected_compact = "DRB10102"
    for name in ["DRB1_0102",
                 "DRB101:02",
                 "HLA-DRB1_0102",
                 "DRB10102",
                 "DRB1*0102",
                 "HLA-DRB1*0102",
                 "HLA-DRB1*01:02"]:
        eq_(normalize_allele_name(name), expected)
        eq_(compact_allele_name(name), expected_compact)

def test_human_class2_alpha_beta():
    expected = "HLA-DPA1*01:05-DPB1*100:01"
    expected_compact = "DPA10105-DPB110001"
    for name in ["DPA10105-DPB110001",
                 "HLA-DPA1*01:05-DPB1*100:01",
                 "hla-dpa1*0105-dpb1*10001",
                 "dpa1*0105-dpb1*10001",
                 "HLA-DPA1*01:05/DPB1*100:01"]:
        eq_(normalize_allele_name(name), expected)
        eq_(compact_allele_name(name), expected_compact)

def test_mouse_class2_alleles():
    # H2-IAb
    eq_(parse_allele_name("H2-IAb"),
        AlleleName("H-2", "IA", "", "b"))
    eq_(normalize_allele_name("H2-IAb"), "H-2-IAb")
    eq_(compact_allele_name("H2-IAb"), "IAb")

    # with hyphen in "H-2"
    eq_(parse_allele_name("H-2-IAb"),
        AlleleName("H-2", "IA", "", "b"))
    eq_(normalize_allele_name("H-2-IAb"), "H-2-IAb")
    eq_(compact_allele_name("H-2-IAb"), "IAb")
