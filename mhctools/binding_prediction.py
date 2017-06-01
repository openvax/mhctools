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

class BindingPrediction(object):
    def __init__(
            self,
            peptide,
            allele,
            affinity,
            percentile_rank,
            source_sequence_name=None,
            offset=0,
            log_affinity=None,
            prediction_method_name=""):
        """
        Parameters
        ----------
        peptide : str
            Short amino acid sequence

        allele : str
            HLA allele, e.g. HLA-A*02:01

        affinity : float
            Predicted binding affinity

        percentile_rank : float
            Percentile rank of the binding affinity for that allele

        source_sequence_name : str
            Name of sequence from which peptide was extracted

        offset : int
            Base0 starting position in source sequence that all epitopes were
            extracted from

        log_affinity : float, optional
            NetMHC sometimes gives invalid IC50 values but we can still
            reconstruct the value from its (1.0 - log_50000(IC50)) score.

        prediction_method_name : str, optional
            Name of predictor used to generate this prediction.
        """
        # if we have a bad IC50 score we might still get a salvageable
        # log of the score. Strangely, this is necessary sometimes!
        if invalid_affinity(affinity) and np.isfinite(log_affinity):
            # pylint: disable=invalid-unary-operand-type
            affinity = 50000 ** (-log_affinity + 1)

        # if IC50 is still NaN or otherwise invalid, abort
        if invalid_affinity(affinity):
            raise ValueError(
                "Invalid IC50 value %0.4f for %s w/ allele %s" % (
                    affinity,
                    peptide,
                    allele))

        if invalid_percentile_rank(percentile_rank):
            raise ValueError(
                "Invalid percentile rank %s for %s w/ allele %s" % (
                    percentile_rank, peptide, allele))

        self.source_sequence_name = source_sequence_name
        self.offset = offset
        self.allele = allele
        self.peptide = peptide
        self.affinity = affinity
        self.percentile_rank = percentile_rank
        self.prediction_method_name = prediction_method_name

    def __str__(self):
        format_string = (
            "BindingPrediction("
            "peptide='%s', "
            "allele='%s', "
            "affinity=%0.4f, "
            "percentile_rank=%s, "
            "source_sequence_name=%s, "
            "offset=%d, "
            "prediction_method_name='%s')")
        return format_string % (
                self.peptide,
                self.allele,
                self.affinity,
                ('%0.4f' % self.percentile_rank
                    if self.percentile_rank
                    else None),
                ('%s' % self.source_sequence_name
                    if self.source_sequence_name
                    else None),
                self.offset,
                self.prediction_method_name)

    def clone_with_updates(self, **kwargs):
        """Returns new BindingPrediction with updated fields"""
        fields_dict = self.to_dict()
        fields_dict.update(kwargs)
        return BindingPrediction(**fields_dict)

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

    def __lt__(self, other):
        return self.value < other.value

def invalid_affinity(x):
    return x < 0 or np.isnan(x) or np.isinf(x)

def invalid_percentile_rank(x):
    # for now, we accept a null percentile rank - not all predictors generate a value
    if x is None:
        return False
    return x < 0 or x > 100
