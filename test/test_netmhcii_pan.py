from mhctools import NetMHCIIpan
from mhctools.alleles import normalize_allele_name

def test_netmhcii_pan():
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

    assert len(epitope_collection) == 27, \
        "Expected 27 epitopes from %s" % (epitope_collection,)

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

    assert len(epitope_collection) == 27, \
        "Expected 27 epitopes from %s" % (epitope_collection,)
