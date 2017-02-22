# Copyright (c) 2014-2017. Mount Sinai School of Medicine
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

import numpy as np
import pandas as pd

class BindingPrediction(object):
    """
    Given the following single sequence FASTA file:
        >seq0
        MEEPQSDPSV
    the result of a binding predictor for HLA-A*02:01 will be
    a collection of the following BindingPrediction objects:
    [
        BindingPrediction(
            source_sequence_key="seq0",
            offset=0,
            allele="HLA-A*02:01",
            peptide="MEEPQSDPS",
            value=0.9,
            percentile_rank=1.3,
            measure=ic50_nM,
            prediction_method_name="NetMHC",

        ),
        BindingPrediction(
            source_sequence_key="seq0",
            offset=1,
            allele="HLA-A*02:01",
            peptide="EEPQSDPSV",
            value=20.9,
            percentile_rank=39.9,
            measure=ic50_nM,
            prediction_method_name="NetMHC",
        ),
    ]
    """
    def __init__(
            self,
            source_sequence_name,
            offset,
            peptide,
            allele,
            affinity,
            percentile_rank,
            log_affinity=None,
            prediction_method_name=""):
        """
        Parameters
        ----------
        source_sequence_name : str
            Name of sequence from which peptide was extracted

        offset : int
            Base0 starting position in source sequence that all epitopes were
            extracted from

        peptide : str
            Short amino acid sequence

        allele : str
            HLA allele, e.g. HLA-A*02:01

        affinity : float
            Predicted binding affinity

        percentile_rank : float
            Percentile rank of the binding affinity for that allele

        log_affinity : float, optional
            NetMHC sometimes gives invalid IC50 values but we can still
            reconstruct the value from its (1.0 - log_50000(IC50)) score.

        prediction_method_name : str, optional
            Name of predictor used to generate this prediction.
        """
        # if we have a bad IC50 score we might still get a salvageable
        # log of the score. Strangely, this is necessary sometimes!
        if invalid_binding_score(affinity) and np.isfinite(log_affinity):
            affinity = 50000 ** (-log_affinity + 1)

        # if IC50 is still NaN or otherwise invalid, abort
        if invalid_binding_score(affinity):
            raise ValueError(
                "Invalid IC50 value %0.4f for %s w/ allele %s" % (
                    affinity,
                    peptide,
                    allele))

        if invalid_binding_score(percentile_rank) or percentile_rank > 100:
            raise ValueError(
                "Invalid percentile rank %s for %s w/ allele %s" % (
                    percentile_rank, peptide, allele))

        self.source_sequence_name = source_sequence_name
        self.offset = offset
        self.allele = allele
        self.peptide = peptide
        self.affinity = affinity
        self.percentile_rank = percentile_rank
        self.prediction_method_name = prediction_method_name,

    def __str__(self):
        format_string = (
            "BindingPrediction("
            "source_sequence_name='%s', "
            "peptide='%s', "
            "allele='%s', "
            "affinity=%0.4f, "
            "percentile_rank=%0.4f, "
            "prediction_method_name='%s')")
        return format_string % (
                self.source_sequence_name,
                self.peptide,
                self.allele,
                self.affinity,
                self.percentile_rank,
                self.prediction_method_name)

    def __repr__(self):
        return str(self)

    @property
    def length(self):
        """Length of peptide, preserved for backwards compatibility"""
        return len(self.peptide)

    @property
    def value(self):
        """Alias for affinity preserved for backwards compatibility"""
        return self.affinity

    fields = (
        "source_sequence_name",
        "offset",
        "peptide",
        "allele",
        "affinity",
        "percentile_rank",
        "prediction_method_name"
    )

    def to_tuple(self):
        return (
            self.source_sequence_name,
            self.offset,
            self.peptide,
            self.allele,
            self.affinity,
            self.percentile_rank,
            self.prediction_method_name)

    def to_dict(self):
        return {k: v for (k, v) in zip(self.fields, self.to_tuple())}

    def __eq__(self, other):
        return (
            other.__class__ is BindingPrediction and
            self.to_tuple() == other.to_tuple())

    def __hash__(self):
        return hash(self.to_tuple())

def invalid_binding_score(x):
    return x < 0 or np.isnan(x) or np.isinf(x)

def binding_predictions_to_dataframe(
        binding_predictions,
        columns=BindingPrediction.fields + ("length",)):
    """
    Converts collection of BindingPrediction objects to DataFrame
    """
    return pd.DataFrame.from_records(
        [tuple([getattr(x, name) for name in columns]) for x in binding_predictions],
        columns=columns)

def update_binding_prediction_fields(binding_predictions, **kwargs):
    """
    Changes fields corresponding to names of keyword arguments in
    each BindingPrediction.
    """
    field_dicts = [x.to_dict() for x in binding_predictions]
    for field_name, values in kwargs.items():
        for i, value in enumerate(values):
            field_dicts[i][field_name] = value
    return [BindingPrediction(**d) for d in field_dicts]
