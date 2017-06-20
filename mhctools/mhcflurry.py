# Copyright (c) 2017. Mount Sinai School of Medicine
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
import logging

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection

from mhcflurry import Class1AffinityPredictor
from mhcflurry.encodable_sequences import EncodableSequences
from numpy import nan

logger = logging.getLogger(__name__)


class MHCflurry(BasePredictor):
    """
    Wrapper around MHCflurry. Users will need to download MHCflurry models
    first.
    See https://github.com/hammerlab/mhcflurry
    """

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            predictor=None):
        """
        Parameters
        -----------
        predictor : mhcflurry.Class1AffinityPredictor (optional)
            MHCflurry predictor to use

        """
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=15)
        if predictor is None:
            predictor = Class1AffinityPredictor.load()
        self.predictor = predictor

    def predict_peptides(self, peptides):
        binding_predictions = []
        encodable_sequences = EncodableSequences.create(peptides)
        for allele in self.alleles:
            predictions = self.predictor.predict(
                encodable_sequences, allele=allele)
            for (i, peptide) in enumerate(peptides):
                binding_prediction = BindingPrediction(
                    allele=allele,
                    peptide=peptide,
                    affinity=predictions[i],

                    # TODO: include percentile rank when MHCflurry supports it
                    percentile_rank=None,
                    prediction_method_name="mhcflurry"
                )
                logger.info(binding_prediction)
                binding_predictions.append(binding_prediction)
        return BindingPredictionCollection(binding_predictions)
