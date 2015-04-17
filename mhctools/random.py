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

import random

from .binding_prediction import BindingPrediction
from .epitope_collection import EpitopeCollection
from .binding_measure import ic50_nM

class RandomPredictor(object):

    def __init__(
            self,
            alleles=['HLA-A*02:01'],
            epitope_lengths=[9]):
        self.alleles = alleles
        self.epitope_lengths = epitope_lengths

    def predict(self, fasta_dict):
        binding_predictions = []
        # if wer'e not running the MHC prediction then we have to manually
        # extract 9mer substrings
        for key, sequence in fasta_dict.items():
            for epitope_length in self.epitope_lengths:
                for i in xrange(len(sequence) - epitope_length + 1):
                    for allele in self.alleles:
                        binding_predictions.append(
                            BindingPrediction(
                                allele=allele,
                                peptide=sequence[i:i + epitope_length],
                                length=epitope_length,
                                base0_start=i,
                                base0_end=i + epitope_length - 1,
                                value=random.random() * 10000.0,
                                percentile_rank=random.randint(0, 99),
                                prediction_method_name="random",
                                binding_measure=ic50_nM,
                                source_sequence=sequence,
                                source_sequence_key=key))
        return EpitopeCollection(binding_predictions)
