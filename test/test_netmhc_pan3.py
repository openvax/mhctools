from mhctools import NetMHCpan3
from mhctools.alleles import normalize_allele_name


DEFAULT_ALLELE = 'HLA-B*18:01'

def test_netmhc_pan3():
    alleles = [normalize_allele_name(DEFAULT_ALLELE)]
    pan_predictor = NetMHCpan3(
        alleles=alleles,
        epitope_lengths=[8])
    fasta_dictionary = {
    	'sequence_0': 'MFCQLAKTY',
    }
    epitope_collection = pan_predictor.predict(
    	fasta_dictionary=fasta_dictionary)

    assert len(epitope_collection) == 2, \
    	'Expected 2 epitopes from %s' % (epitope_collection,)
