from nose.tools import eq_
from mhctools import BindingPrediction, BindingPredictionCollection

def test_collection_to_dataframe():
    bp = BindingPrediction(
        peptide="SIINFEKL",
        allele="A0201",
        affinity=1.5,
        percentile_rank=0.1)
    collection = BindingPredictionCollection([bp])
    df = collection.to_dataframe()
    eq_(df.peptide.iloc[0], "SIINFEKL")
    eq_(df.affinity.iloc[0], 1.5)
    eq_(df.allele.iloc[0], "A0201")
    eq_(df.percentile_rank.iloc[0], 0.1)
