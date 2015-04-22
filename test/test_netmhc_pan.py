from mhctools import NetMHCpan
from mhctools.alleles import normalize_allele_name


DEFAULT_ALLELE = 'HLA-A*02:01'

def test_netmhc_pan():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    pan_predictor = NetMHCpan(
        alleles=alleles,
        epitope_lengths=[9])
    fasta_dictionary = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    epitope_collection = pan_predictor.predict(
        fasta_dictionary=fasta_dictionary)

    assert len(epitope_collection) == 4, \
        "Expected 4 epitopes from %s" % (epitope_collection,)