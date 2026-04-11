# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



from .common import eq_
from numpy import testing

from mhcflurry import Class1AffinityPredictor
from mhctools import MHCflurry, MHCflurry_Affinity
from mhctools.pred import Kind

DEFAULT_ALLELE = "HLA-A*02:01"

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}


def test_mhcflurry_presentation_predict_subsequences():
    """Legacy predict_subsequences still works with the presentation-based MHCflurry."""
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    eq_(4, len(binding_predictions),
        "Expected 4 binding predictions from %s" % (binding_predictions,))

    # Verify affinity values are reasonable (positive nM values)
    for bp in binding_predictions:
        assert bp.affinity > 0, "Affinity should be positive, got %s" % bp.affinity
        assert bp.allele == DEFAULT_ALLELE


def test_mhcflurry_presentation_predict():
    """New predict() API returns PeptideResults with both affinity and presentation."""
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    peptides = ["SIINFEKL", "ASILLLVFY"]
    results = predictor.predict(peptides)

    eq_(len(peptides), len(results),
        "Expected one PeptideResult per peptide")

    for r in results:
        assert r.peptide in peptides
        assert r.affinity is not None, "Should have affinity prediction"
        assert r.presentation is not None, "Should have presentation prediction"

        # Check affinity prediction
        assert r.affinity.kind == Kind.pMHC_affinity
        assert r.affinity.value > 0, "Affinity value should be positive nM"
        assert r.affinity.score >= 0
        assert r.affinity.allele == DEFAULT_ALLELE
        assert r.affinity.predictor_name == "mhcflurry"

        # Check presentation prediction
        assert r.presentation.kind == Kind.pMHC_presentation
        assert r.presentation.score >= 0
        assert r.presentation.allele == DEFAULT_ALLELE
        assert r.presentation.predictor_name == "mhcflurry"

        # Check that percentile ranks are present
        assert r.affinity.percentile_rank is not None
        assert r.presentation.percentile_rank is not None

    # Verify kinds
    for r in results:
        assert Kind.pMHC_affinity in r.kinds
        assert Kind.pMHC_presentation in r.kinds


def test_mhcflurry_presentation_affinity_matches_old_api():
    """Affinity values from the presentation predictor should be close to the
    old Class1AffinityPredictor values."""
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_peptides(["SIINFEKL"])

    old_predictor = Class1AffinityPredictor.load()
    old_prediction = old_predictor.predict(["SIINFEKL"], allele=DEFAULT_ALLELE)

    assert len(binding_predictions) == 1
    # Approximate check -- presentation predictor uses the same affinity model
    # but values can differ slightly
    testing.assert_almost_equal(
        binding_predictions[0].affinity, old_prediction[0], decimal=0)


def test_mhcflurry_affinity_only():
    """MHCflurry_Affinity subclass only produces affinity predictions."""
    predictor = MHCflurry_Affinity(alleles=[DEFAULT_ALLELE])
    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    eq_(4, len(binding_predictions),
        "Expected 4 binding predictions from %s" % (binding_predictions,))

    prediction_scores = {
        (x.peptide, x.allele): x.affinity for x in binding_predictions
    }

    old_predictor = Class1AffinityPredictor.load()
    # test one prediction at a time to make sure there's no peptide/allele mixup
    for (peptide, allele), affinity in prediction_scores.items():
        prediction = old_predictor.predict([peptide], allele=allele)
        assert len(prediction) == 1
        testing.assert_almost_equal(prediction[0], affinity, decimal=0)


def test_mhcflurry_multiple_alleles():
    """MHCflurry with multiple alleles produces predictions for each allele."""
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    predictor = MHCflurry(alleles=alleles)
    results = predictor.predict(["SIINFEKL"])

    eq_(1, len(results), "Expected one PeptideResult")
    r = results[0]

    # Should have 2 alleles x 2 kinds = 4 predictions total
    eq_(4, len(r.preds), "Expected 4 predictions (2 alleles x 2 kinds)")

    # Both alleles should be present
    assert r.alleles == set(alleles)

    # Both kinds should be present
    assert r.kinds == {Kind.pMHC_affinity, Kind.pMHC_presentation}
