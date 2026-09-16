
from mhctools import RandomBindingPredictor
from .common import eq_


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
    df = predictor.predict_subsequences_dataframe(fasta_dict)
    print(df)
    # make sure we have predictions for all the alleles
    eq_(set(df.allele), set(alleles))
    # make sure all entries from the fasta dict are present
    eq_(set(df.source_sequence_name), set(fasta_dict.keys()))


def test_default_peptide_lengths_are_not_shared_between_instances():
    first = RandomBindingPredictor()
    second = RandomBindingPredictor()

    first.default_peptide_lengths.append(10)

    assert first.default_peptide_lengths == [9, 10]
    assert second.default_peptide_lengths == [9]
    assert RandomBindingPredictor().default_peptide_lengths == [9]
