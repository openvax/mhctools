from nose.tools import eq_

from mhctools import NetMHCpan
from mhcnames import normalize_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

def run_class_with_executable(mhcpan_class, mhcpan_executable):
    alleles = [normalize_allele_name("HLA-A*02:01")]
    predictor = mhcpan_class(
        alleles=alleles,
        program_name=mhcpan_executable)
    return predictor.predict_subsequences(
        sequence_dict=protein_sequence_dict,
        peptide_lengths=[9])

def test_netmhc_pan():
    binding_predictions = run_class_with_executable(NetMHCpan, "netMHCpan")
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
