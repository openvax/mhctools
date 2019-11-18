# Copyright (c) 2014-2019. Mount Sinai School of Medicine
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
from serializable import Serializable

class BindingPrediction(Serializable):
    def __init__(
            self,
            peptide,
            allele,
            score=None,
            percentile_rank=None,
            affinity=None,
            source_sequence_name=None,
            offset=0,
            prediction_method_name=""):
        """
        Parameters
        ----------
        peptide : str
            Short amino acid sequence

        allele : str
            HLA allele, e.g. HLA-A*02:01

        score : float
            Continuous prediction of peptide-MHC binding where larger values
            indicate either higher affinity or higher probability. For affinity
            predictors this can be 1-log(IC50)/log(max_IC50) For mass spec
            predictors this can be the probability of detection.

        percentile_rank : float
            Percentile rank of the score

        affinity : float
            Predicted binding affinity IC50

        source_sequence_name : str
            Name of sequence from which peptide was extracted

        offset : int
            Base0 starting position in source sequence that all epitopes were
            extracted from

        prediction_method_name : str, optional
            Name of predictor used to generate this prediction.
        """
        self.source_sequence_name = source_sequence_name
        self.offset = offset
        self.allele = allele
        self.peptide = peptide

        if score is None and affinity is not None:
            # make an ascending score by taking 1-log_50k (IC50)
            score = 1.0 - (np.log(affinity) / np.log(50000))

        self.score = score
        self.percentile_rank = percentile_rank
        self.affinity = affinity
        self.prediction_method_name = prediction_method_name

    def __str__(self):
        format_string = (
            "BindingPrediction("
            "peptide='%s', "
            "allele='%s', "
            "score=%s, "
            "percentile_rank=%s, "
            "affinity=%s, "
            "source_sequence_name=%s, "
            "offset=%d, "
            "prediction_method_name='%s')")
        return format_string % (
                self.peptide,
                self.allele,
                ("None" if self.score is None else '%0.3f' % self.score),
                ("None" if self.percentile_rank is None else '%0.3f' % self.percentile_rank),
                ("None" if self.affinity is None else '%0.3f' % self.affinity),
                ("None" if self.source_sequence_name is None else "'%s'" % self.source_sequence_name),
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

    @property
    def elution_score(self):
        """
        Deprecated alias of `score` from when we only considered
        predictors of peptide-MHC binding affinity.

        Returns
        -------
        float
        """
        return self.score


    fields = (
        "source_sequence_name",
        "offset",
        "peptide",
        "allele",
        "score",
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
            self.score,
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

