# Copyright (c) 2016-2017. Mount Sinai School of Medicine
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

logger = logging.getLogger(__name__)


class MHCFlurry(BasePredictor):
    """
    Wrapper around MHCFlurry. Assumes you've installed it and downloaded datasets and trained  MHCFlurry: see
    https://github.com/hammerlab/mhcflurry for instructions.
    """

    def __init__(
            self,
            alleles,
            default_peptide_lengths=[9],
            models_dir=None):
        """
        Standard mhctools constructor, with additional ability to pass in a directory containing
        saved MHCFlurry models.
        """
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths)
        self.models_dir = models_dir

    def predict_peptides(self, peptides):
        binding_predictions = []
        predictor = Class1AffinityPredictor.load(models_dir=self.models_dir)
        encodable_sequences = EncodableSequences.create(peptides)
        for allele in self.alleles:
            predictions = predictor.predict(encodable_sequences, allele=allele)
            for (i, peptide) in enumerate(peptides):
                binding_prediction = BindingPrediction(
                    allele=allele,
                    peptide=peptide,
                    affinity=predictions[i],
                    percentile_rank=0.0,  # TODO: grab percentile rank when it's available
                    prediction_method_name="mhcflurry"
                )
                binding_predictions.append(binding_prediction)
        return BindingPredictionCollection(binding_predictions)
