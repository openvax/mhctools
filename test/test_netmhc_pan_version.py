from nose.tools import raises, eq_

from mhctools import NetMHCpan, NetMHCpan2, NetMHCpan3
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

@raises(SystemError)
def test_executable_mismatch_2_3():
    run_class_with_executable(NetMHCpan2, "netMHCpan3")

@raises(SystemError)
def test_executable_mismatch_3_2():
    run_class_with_executable(NetMHCpan3, "netMHCpan2")

def test_netmhc_pan():
    epitope_collection = run_class_with_executable(NetMHCpan, "netMHCpan")
    assert len(epitope_collection) == 4, \
            "Expected 4 epitopes from %s" % (epitope_collection,)
