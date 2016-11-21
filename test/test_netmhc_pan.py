from nose.tools import raises, eq_

from mhctools import NetMHCpan, NetMHCpan28, NetMHCpan3
from mhctools.alleles import normalize_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'


def run_class_with_executable(mhcpan_class, mhcpan_executable):
    alleles = [normalize_allele_name("HLA-A*02:01")]
    predictor = mhcpan_class(
        alleles=alleles,
        epitope_lengths=[9],
        program_name=mhcpan_executable)
    fasta_dictionary = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    return predictor.predict(fasta_dictionary=fasta_dictionary)

def test_netmhc_pan():
    epitope_collection = run_class_with_executable(NetMHCpan, "netMHCpan")
    assert len(epitope_collection) == 4, \
            "Expected 4 epitopes from %s" % (epitope_collection,)
    for epitope in epitope_collection:
        # recompute the peptide from the offset and starting sequence, and make sure it matches.
        # this is currently wrong in netMHCpan-3.0 and we want to test our wrapper fix to that
        offset = epitope.offset
        length = epitope.length
        expected_peptide = epitope.source_sequence[offset:offset+length]
        eq_(expected_peptide, epitope.peptide,
            "Peptide mismatch: expected %s but got %s in binding prediction '%s'" % (
                expected_peptide, epitope.peptide, epitope,))
