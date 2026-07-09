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


def test_mhcflurry_predict_proteins_matches_direct_flanked_prediction():
    """Protein scanning should forward the same flanks as direct MHCflurry."""
    predictor = MHCflurry(alleles=[DEFAULT_ALLELE])
    peptide = "SIINFEKL"
    protein = "MDSKG%sGSRLL" % peptide
    offset = protein.index(peptide)
    n_flank_length, c_flank_length = predictor._predict_protein_flank_lengths()
    n_flank = protein[max(0, offset - n_flank_length):offset]
    c_flank = protein[offset + len(peptide):
                      offset + len(peptide) + c_flank_length]

    protein_results = predictor.predict_proteins(
        {"protein": protein},
        peptide_lengths=[len(peptide)])
    protein_result = [
        pp for pp in protein_results["protein"] if pp.offset == offset
    ][0]

    direct = predictor.predictor.predict(
        peptides=[peptide],
        alleles=[DEFAULT_ALLELE],
        n_flanks=[n_flank],
        c_flanks=[c_flank],
        include_affinity_percentile=False,
        verbose=0)

    assert protein_result.presentation.n_flank == n_flank
    assert protein_result.presentation.c_flank == c_flank
    testing.assert_allclose(
        protein_result.presentation.score,
        direct.presentation_score.iloc[0],
        rtol=1e-6,
        atol=1e-6)
    testing.assert_allclose(
        protein_result.presentation.percentile_rank,
        direct.presentation_percentile.iloc[0],
        rtol=1e-6,
        atol=1e-6)


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
    """MHCflurry affinity is per-allele; presentation is haplotype-level."""
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    predictor = MHCflurry(alleles=alleles)
    results = predictor.predict(["SIINFEKL"])

    eq_(1, len(results), "Expected one PeptideResult")
    r = results[0]

    # Two per-allele affinity predictions, one haplotype-level presentation
    # prediction with best_allele attribution, and one allele-independent
    # antigen-processing prediction.
    eq_(4, len(r.preds), "Expected 4 predictions")
    eq_(2, len(r.filter(kind=Kind.pMHC_affinity)))
    eq_(1, len(r.filter(kind=Kind.pMHC_presentation)))
    eq_(1, len(r.filter(kind=Kind.antigen_processing)))

    # Both alleles should be present through affinity predictions
    # (the processing prediction is allele-less, so it doesn't add one).
    assert r.alleles == set(alleles)

    # All three kinds should be present
    assert r.kinds == {
        Kind.pMHC_affinity, Kind.pMHC_presentation, Kind.antigen_processing}


def test_mhcflurry_processing_score():
    """MHCflurry surfaces its antigen-processing (cleavage) score.

    The score is allele-independent (one per peptide, no allele) and matches
    the ``processing_score`` MHCflurry's presentation predictor computes.
    """
    from mhcflurry import Class1PresentationPredictor

    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    peptides = ["SIINFEKL", "GILGFVFTL"]
    n_flanks = ["AAA", "CCC"]
    c_flanks = ["KKK", "DDD"]

    predictor = MHCflurry(alleles=alleles)
    results = predictor.predict(peptides, n_flanks=n_flanks, c_flanks=c_flanks)

    # Ground truth straight from MHCflurry.
    raw = Class1PresentationPredictor.load().predict(
        peptides=peptides,
        alleles={a: [a] for a in alleles},
        n_flanks=n_flanks,
        c_flanks=c_flanks,
        include_affinity_percentile=False,
        verbose=0)
    expected = dict(zip(raw["peptide"], raw["processing_score"]))

    for r in results:
        processing = r.filter(kind=Kind.antigen_processing)
        eq_(1, len(processing),
            "Expected exactly one processing prediction per peptide")
        pred = processing[0]
        assert pred.allele == "", "Processing score is allele-independent"
        assert pred.predictor_name == "mhcflurry"
        # flanks are carried through as provenance
        idx = peptides.index(r.peptide)
        assert pred.n_flank == n_flanks[idx]
        assert pred.c_flank == c_flanks[idx]
        testing.assert_allclose(pred.score, expected[r.peptide], rtol=1e-5)
        # the convenience accessor resolves to the same prediction
        assert r.processing is not None
        assert r.processing.kind == Kind.antigen_processing

    # antigen_processing is advertised as an allele-independent supported kind
    support = predictor.kind_support()[Kind.antigen_processing]
    eq_("none", support["mhc_dependence"])
    eq_("none", support["mhc_class"])


def test_mhcflurry_processing_score_without_flanks():
    """The processing score is still emitted when no flanks are provided."""
    predictor = MHCflurry(alleles=["HLA-A*02:01"])
    r = predictor.predict(["SIINFEKL"])[0]
    processing = r.filter(kind=Kind.antigen_processing)
    eq_(1, len(processing), "Expected one processing prediction")
    pred = processing[0]
    assert pred.allele == ""
    assert pred.n_flank == "" and pred.c_flank == ""
    assert pred.score is not None


def test_mhcflurry_legacy_predict_peptides_unchanged():
    """Surfacing processing on predict() must not change the legacy path.

    The CLI (--output-csv) and predict_peptides go through the affinity-only
    BindingPredictionCollection path, which should stay one affinity record
    per (peptide, allele) with no antigen_processing rows.
    """
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    peptides = ["SIINFEKL", "GILGFVFTL"]
    predictor = MHCflurry(alleles=alleles)
    collection = predictor.predict_peptides(peptides)

    eq_(len(peptides) * len(alleles), len(collection),
        "Legacy path should emit one affinity record per (peptide, allele)")
    for bp in collection:
        assert bp.affinity is not None
        assert bp.prediction_method_name == "mhcflurry"
