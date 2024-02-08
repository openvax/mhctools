from mhctools.binding_prediction import BindingPrediction
from nose.tools import eq_

def test_binding_prediction_fields():
    bp = BindingPrediction(
            source_sequence_name="seq",
            offset=0,
            peptide="SIINFEKL",
            allele="H-2-K-d",
            affinity=200.0,
            percentile_rank=0.3)
    eq_(bp.source_sequence_name, "seq")
    eq_(bp.offset, 0)
    eq_(bp.peptide, "SIINFEKL")
    eq_(bp.allele, "H-2-K-d")
    eq_(bp.affinity, 200.0)
    eq_(bp.percentile_rank, 0.3)

def test_binding_prediction_str_repr():
    bp = BindingPrediction(
            source_sequence_name="seq",
            offset=0,
            peptide="SIINFEKL",
            allele="H-2-K-d",
            affinity=200.0,
            percentile_rank=0.3)
    eq_(str(bp), repr(bp))
    assert "SIINFEKL" in str(bp)
    assert "200.0" in str(bp)

def test_binding_predicton_eq():
    bp1 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKL",
        allele="H-2-K-d",
        affinity=200.0,
        percentile_rank=0.3)
    bp2 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKLY",
        allele="H-2-K-d",
        affinity=5000.0,
        percentile_rank=9.3)
    eq_(bp1, bp1)
    eq_(bp2, bp2)
    assert bp1 != bp2

def test_binding_predicton_hash():
    bp1 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKL",
        allele="H-2-K-d",
        affinity=200.0,
        percentile_rank=0.3)
    bp2 = BindingPrediction(
        source_sequence_name="seq",
        offset=0,
        peptide="SIINFEKLY",
        allele="H-2-K-d",
        affinity=5000.0,
        percentile_rank=9.3)
    eq_(hash(bp1), hash(bp1))
    assert hash(bp1) != hash(bp2)
