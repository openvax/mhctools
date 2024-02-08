from nose.tools import assert_raises

from mhctools import IedbNetMHCpan


DEFAULT_ALLELE = 'HLA-A*02:01'
UNSUPPORTED_ALLELE = 'HLA-A*24:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}


def test_netmhcpan_iedb():
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[8,9])
    assert len(binding_predictions) == 10, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)


def test_netmhcpan_iedb_unsupported_allele():
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE, UNSUPPORTED_ALLELE], raise_on_error=False)
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)

    # check that the error is raised when raise_on_error is left at default (True)
    predictor = IedbNetMHCpan(alleles=[DEFAULT_ALLELE, UNSUPPORTED_ALLELE])
    with assert_raises(ValueError):
        binding_predictions = predictor.predict_subsequences(
            protein_sequence_dict,
            peptide_lengths=[9])
