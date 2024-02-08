
"""
Class I predictors typically allow 8mer or longer peptides, whereas
NetMHCIIpan allows 9mer or longer. So far only MHCflurry has a max
length (of 15mer).
"""

from mhctools import NetMHCIIpan, NetMHC
from nose.tools import eq_, assert_raises


def test_class2_9mer_success():
    ii_pan_predictor = NetMHCIIpan(alleles=["HLA-DRB1*01:01"])
    predictions = ii_pan_predictor.predict_peptides(["A" * 9])
    eq_(len(predictions), 1)

def test_class2_8mer_fails():
    ii_pan_predictor = NetMHCIIpan(alleles=["HLA-DRB1*01:01"])
    with assert_raises(ValueError):
        ii_pan_predictor.predict_peptides(["A" * 8])

def test_class1_8mer_success():
    netmhc = NetMHC(alleles=["HLA-A0201"])
    predictions = netmhc.predict_peptides(["A" * 8])
    eq_(len(predictions), 1)

def test_class1_7mer_failure():
    netmhc = NetMHC(alleles=["HLA-A0201"])
    with assert_raises(ValueError):
        netmhc.predict_peptides(["A" * 7])
