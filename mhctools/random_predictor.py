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
import random

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection

class RandomBindingPredictor(BasePredictor):
    def __init__(
            self,
            alleles=['HLA-A*02:01'],
            default_peptide_lengths=[9]):
        BasePredictor.__init__(
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths)

    def predict_peptides(self, peptides):
        binding_predictions = []
        for allele in self.alleles:
            for p in peptides:
                binding_predictions.append(
                    BindingPrediction(
                        source_sequence_name=None,
                        offset=0,
                        allele=allele,
                        peptide=p,
                        ic50=random.random() * 10000.0,
                        rank=random.randint(0, 99)))
        return BindingPredictionCollection(binding_predictions)
