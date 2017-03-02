from nose.tools import raises, eq_

from mhctools import NetMHC, NetMHC3, NetMHC4
from mhcnames import normalize_allele_name


def run_class_with_executable(mhc_class, mhc_executable):
    alleles = [normalize_allele_name("HLA-A*02:01")]
    predictor = mhc_class(
        alleles=alleles,
        program_name=mhc_executable)
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    return predictor.predict_subsequences(
        sequence_dict=sequence_dict,
        peptide_lengths=[9])

@raises(SystemError)
def test_executable_mismatch_3_4():
    run_class_with_executable(NetMHC3, "netMHC")

@raises(SystemError)
def test_executable_mismatch_4_3():
    run_class_with_executable(NetMHC4, "netMHC-3.4")

def test_wrapper_function():
    alleles = [normalize_allele_name("HLA-A*02:01")]
    wrapped_4 = NetMHC(
        alleles=alleles,
        default_peptide_lengths=[9],
        program_name="netMHC")
    eq_(type(wrapped_4), NetMHC4)
    wrapped_3 = NetMHC(
        alleles=alleles,
        default_peptide_lengths=[9],
        program_name="netMHC-3.4")
    eq_(type(wrapped_3), NetMHC3)

@raises(SystemError, OSError)
def test_wrapper_failure():
    alleles = [normalize_allele_name("HLA-A*02:01")]
    NetMHC(alleles=alleles,
           default_peptide_lengths=[9],
           program_name="netMHC-none")
