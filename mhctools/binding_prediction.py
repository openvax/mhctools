# Copyright (c) 2014. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
from collections import namedtuple

"""
Given the following single sequence FASTA file:
    >seq0
    MEEPQSDPSV
the result of a binding predictor for HLA-A*02:01 will be
a collection of the following BindingPrediction objects:
    [
        BindingPrediction(
            source_sequence_key="seq0",
            source_sequence="MEEPQSDPSV",
            offset=0,
            allele="HLA-A*02:01",
            peptide="MEEPQSDPS",
            length=9,
            value=0.9,
            percentile_rank=1.3,
            measure=ic50_nM,
            prediction_method_name="NetMHC",

        ),
        BindingPrediction(
            source_sequence_key="seq0",
            source_sequence="MEEPQSDPSV",
            offset=1,
            allele="HLA-A*02:01",
            peptide="EEPQSDPSV",
            length=9,
            value=20.9,
            percentile_rank=39.9,
            measure=ic50_nM,
            prediction_method_name="NetMHC",
        ),
    ]
"""

BindingPrediction = namedtuple("BindingPrediction",
    [
        # "key" of source sequence is often the ID from a FASTA file
        "source_sequence_key",
        # longer amino acid sequence from which peptide originated
        "source_sequence",
        # base-0 start position of this peptide in larger sequence
        "offset",
        # HLA allele, e.g. "HLA-A*02:01"
        "allele",
        # peptide sequence, e.g. "SIINFKELL"
        "peptide",
        # length of peptide
        "length",
        # predicted binding value
        "value",
        # what is the predicted value measured (a BindingMeasure object)
        "measure",
        # percentile rank of predicted value, if available (lower is better)
        "percentile_rank",
        # name of predictor e.g. "NetMHC"
        "prediction_method_name",

    ])
