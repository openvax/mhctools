from nose.tools import eq_

from mhctools import NetMHCIIpan
from mhctools.alleles import normalize_allele_name

def test_netmhcii_pan_DRB():
    alleles = [normalize_allele_name("HLA-DRB1*01:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles,
        epitope_lengths=[15, 16])
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    epitope_collection = ii_pan_predictor.predict(
        fasta_dictionary=fasta_dictionary)

    unique_lengths = {x.length for x in epitope_collection}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in epitope_collection}
    eq_(unique_alleles, {"HLA-DRB1*01:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(epitope_collection) == 52, \
        "Expected 52 epitopes from %s" % (epitope_collection,)

def test_netmhcii_pan_alpha_beta():
    alleles = [normalize_allele_name("HLA-DPA1*01:05-DPB1*100:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles,
        epitope_lengths=[15, 16])
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }

    epitope_collection = ii_pan_predictor.predict(
        fasta_dictionary=fasta_dictionary)
    unique_lengths = {x.length for x in epitope_collection}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in epitope_collection}
    eq_(unique_alleles, {"HLA-DPA1*01:05-DPB1*100:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(epitope_collection) == 52, \
        "Expected 52 epitopes from %s" % (epitope_collection,)

def test_netmhcii_pan_multiple_alleles():
    alleles = [
        normalize_allele_name("HLA-DPA1*01:05-DPB1*100:01"),
        normalize_allele_name("HLA-DQA1*05:11-DQB1*03:02"),
        normalize_allele_name("HLA-DRB1*01:01")
    ]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles,
        epitope_lengths=[15, 16])
    fasta_dictionary = {
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    epitope_collection = ii_pan_predictor.predict(
        fasta_dictionary=fasta_dictionary)

    unique_lengths = {x.length for x in epitope_collection}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in epitope_collection}
    eq_(unique_alleles, {
        "HLA-DPA1*01:05-DPB1*100:01",
        "HLA-DQA1*05:11-DQB1*03:02",
        "HLA-DRB1*01:01"
    })

    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect 3 * ((21-15+1) + (21-16+1)) = 39 entries
    assert len(epitope_collection) == 39, \
        "Expected 39 epitopes from %s" % (epitope_collection,)

def test_netmhcii_pan_mouse():
    alleles = [normalize_allele_name("H2-IAb")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles,
        epitope_lengths=[15, 16])
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    epitope_collection = ii_pan_predictor.predict(
        fasta_dictionary=fasta_dictionary)

    unique_lengths = {x.length for x in epitope_collection}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in epitope_collection}
    eq_(unique_alleles, {"H-2-IAb"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(epitope_collection) == 52, \
        "Expected 52 epitopes from %s" % (epitope_collection,)
