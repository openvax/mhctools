from nose.tools import eq_

from mhctools import NetMHCpan


DEFAULT_ALLELE = 'HLA-A*02:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

# This test assumes you have versions of netMHCpan installed in your PATH
# with the following names:
NETMHCPAN_BINARIES = [
    "netMHCpan",
    "netMHCpan-2.8",
    "netMHCpan-3.0",
    "netMHCpan-4.0",
]


def test_netmhc_pan():
    for program_name in NETMHCPAN_BINARIES:
        yield check_netmhc_pan, program_name


def check_netmhc_pan(program_name):
    predictor = NetMHCpan(alleles=[DEFAULT_ALLELE], program_name=program_name)
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
