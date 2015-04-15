from collections import namedtuple

"""
Given a longer amino acid sequence such as MEEPQSDPSV
result of a binding predictor can be represented as a dictionary
from alleles to a list of ShortPeptidePrediction objects, such as:
    {
        "HLA-A*02:01" : [
            ShortPeptidePrediction(
                allele="HLA-A*02:01",
                peptide="MEEPQSDPS",
                length=9,
                start=0,
                end=9,
                value=0.9,
                prediction_method_name="NetMHC",
                measure=ic50_nM,
                prediction_method_name="NetMHC"),
            ShortPeptidePrediction(
                allele="HLA-A*02:01",
                peptide="EEPQSDPSV",
                length=9,
                start=1,
                end=10,
                value=20.9,
                percentile_rank=
                measure=ic50_nM,
                prediction_method_name="NetMHC"),
        ]
    }
"""

BindingPrediction = namedtuple("ShortPeptidePrediction",
    [
        # HLA allele, e.g. "HLA-A*02:01"
        "allele",
        # peptide sequence, e.g. "SIINFKELL"
        "peptide",
        # length of peptide
        "length",
        # base 0 half-open start position of this peptide in a larger sequence
        "start",
        # base 0 half-open end position of this peptide in a larger sequence
        "end",
        # predicted binding value
        "value",
        # what is the predicted value measured (a BindingMeasure object)
        "measure",
        # percentile rank of predicted value, if available (lower is better)
        "percentile_rank",
        # name of predictor e.g. "NetMHC"
        "prediction_method_name",
    ])
