from collections import namedtuple

"""
Given the following single sequence FASTA file:
    >seq0
    MEEPQSDPSV
the result of a binding predictor for HLA-A*02:01 will be
a collection of the following BindingPrediction objects:
    [
        BindingPrediction(
            allele="HLA-A*02:01",
            peptide="MEEPQSDPS",
            length=9,
            base0_start=0,
            base0_end=9,
            value=0.9,
            percentile_rank=1.3,
            measure=ic50_nM,
            prediction_method_name="NetMHC",
            source_sequence="MEEPQSDPSV",
            source_sequence_key="seq0",
        ),
        BindingPrediction(
            allele="HLA-A*02:01",
            peptide="EEPQSDPSV",
            length=9,
            base0_start=1,
            base0_end=10,
            value=20.9,
            percentile_rank=39.9,
            measure=ic50_nM,
            prediction_method_name="NetMHC",
            source_sequence="MEEPQSDPSV",
            source_sequence_key="seq0",
        ),
    ]
"""

BindingPrediction = namedtuple("BindingPrediction",
    [
        # HLA allele, e.g. "HLA-A*02:01"
        "allele",
        # peptide sequence, e.g. "SIINFKELL"
        "peptide",
        # length of peptide
        "length",
        # start position of this peptide in larger sequence
        "base0_start",
        # end position of this peptide in larger sequence
        "base0_end",
        # predicted binding value
        "value",
        # what is the predicted value measured (a BindingMeasure object)
        "measure",
        # percentile rank of predicted value, if available (lower is better)
        "percentile_rank",
        # name of predictor e.g. "NetMHC"
        "prediction_method_name",
        # longer amino acid sequence from which peptide originated
        "source_sequence",
        # "key" of source sequence is often the ID from a FASTA file
        "source_sequence_key",
    ])
