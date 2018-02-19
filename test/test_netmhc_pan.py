from nose.tools import eq_

from mhctools import NetMHCpan

# Defining FileNotFoundError for Python 2.x
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


DEFAULT_ALLELE = 'HLA-A*02:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

# Tests will also be run using the following program names if they are installed.
# In any case, a program called "netMHCpan" MUST be installed and working for
# this test suite to succeeed.
OPTIONAL_NETMHCPAN_PROGRAM_NAMES = [
    "netMHCpan-2.8",
    "netMHCpan-3.0",
    "netMHCpan-4.0",
]


def test_netmhc_pan():
    yield check_netmhc_pan, "netMHCpan", True  # required

    for program_name in OPTIONAL_NETMHCPAN_PROGRAM_NAMES:
        yield check_netmhc_pan, program_name, False  # optional


def check_netmhc_pan(program_name, fail_if_no_such_program=True):
    try:
        predictor = NetMHCpan(
            alleles=[DEFAULT_ALLELE], program_name=program_name)
    except FileNotFoundError:
        if fail_if_no_such_program:
            raise
        print("Skipping because no such program: %s" % program_name)
        return

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)
    for x in binding_predictions:
        # recompute the peptide from the offset and starting sequence, and make sure it matches.
        # this is currently wrong in netMHCpan-3.0 and we want to test our wrapper fix to that
        offset = x.offset
        length = x.length
        seq_name = x.source_sequence_name
        expected_peptide = protein_sequence_dict[seq_name][offset:offset + length]
        eq_(expected_peptide, x.peptide,
            "Peptide mismatch: expected %s but got %s in binding prediction '%s'" % (
                expected_peptide, x.peptide, x,))
