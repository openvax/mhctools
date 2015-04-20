
from nose.tools import eq_
from mhctools import RandomBindingPredictor

alleles = [
    "HLA-A*02:01",
    "HLA-B*05:02"
]
predictor = RandomBindingPredictor(alleles)

fasta_dict = {
    "seq0": "SIINKFKEEL",
    "seq1": "AAAPPPLLLTTT",
}
def test_random_mhc_binding_predictions():
    epitopes = predictor.predict(fasta_dict)
    df = epitopes.dataframe()
    # make sure we have predictions for all the alleles
    eq_(set(df.allele), set(alleles))
    # make sure all entries from the fasta dict are present
    eq_(set(df.source_sequence), set(fasta_dict.values()))
    eq_(set(df.source_sequence_key), set(fasta_dict.keys()))
