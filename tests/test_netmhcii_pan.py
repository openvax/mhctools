from nose.tools import eq_

from mhctools import NetMHCIIpan
from mhcnames import normalize_allele_name

def test_netmhcii_pan_DRB():
    alleles = [normalize_allele_name("HLA-DRB1*01:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DRA1*01:01-DRB1*01:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_alpha_beta():
    alleles = [normalize_allele_name("HLA-DPA1*01:05-DPB1*100:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }

    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])
    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DPA1*01:05-DPB1*100:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_multiple_alleles():
    alleles = [
        normalize_allele_name("HLA-DPA1*01:05-DPB1*100:01"),
        normalize_allele_name("HLA-DQA1*05:11-DQB1*03:02"),
        normalize_allele_name("HLA-DRA1*01:01-DRB1*01:01")
    ]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {
        "HLA-DPA1*01:05-DPB1*100:01",
        "HLA-DQA1*05:11-DQB1*03:02",
        "HLA-DRA1*01:01-DRB1*01:01"
    })

    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect 3 * ((21-15+1) + (21-16+1)) = 39 entries
    assert len(binding_predictions) == 39, \
        "Expected 39 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_mouse():
    alleles = [normalize_allele_name("H2-IAb")]
    ii_pan_predictor = NetMHCIIpan(alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"H-2-IAb"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

# Tests will also be run using the following program names if they are installed.
# In any case, a program called "netMHCIIpan" MUST be installed and working for
# this test suite to succeeed.
OPTIONAL_NETMHCIIPAN_PROGRAM_NAMES = [
    "netMHCIIpan-4.0",
]

DEFAULT_ALLELE = "HLA-DPA1*01:05-DPB1*100:01"

protein_sequence_dict = {
    "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
    "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
}

# follows similar model to test_netmhc_pan.py
def test_netmhcii_pan():
    yield check_netmhcii_pan, "netMHCIIpan", True  # required

    for program_name in OPTIONAL_NETMHCIIPAN_PROGRAM_NAMES:
        yield check_netmhcii_pan, program_name, False  # optional

# TODO: add more actual tests to this
def check_netmhcii_pan(program_name, fail_if_no_such_program=True):
    try:
        predictor = NetMHCIIpan(
            alleles=[DEFAULT_ALLELE], program_name=program_name)
    except FileNotFoundError:
        if fail_if_no_such_program:
            raise
        print("Skipping because no such program: %s" % program_name)
        return