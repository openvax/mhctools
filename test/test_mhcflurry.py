from nose.tools import eq_
from numpy import testing

from mhcflurry import Class1AffinityPredictor
from mhctools import MHCflurry

DEFAULT_ALLELE = "HLA-A*02:01"

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

def test_mhcflurry():
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    eq_(4, len(binding_predictions),
        "Expected 4 binding predictions from %s" % (binding_predictions,))

    prediction_scores = {
        (x.peptide, x.allele): x.affinity for x in binding_predictions
    }

    predictor = Class1AffinityPredictor.load()
    # test one prediction at a time to make sure there's no peptide/allele mixup
    for (peptide, allele), affinity in prediction_scores.items():
        prediction = predictor.predict([peptide], allele=allele)
        assert len(prediction) == 1
        # we've seen results differ a bit so doing an approximate check, not an error condition
        testing.assert_almost_equal(prediction[0], affinity, decimal=0)
