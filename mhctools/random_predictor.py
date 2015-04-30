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
import random

from .epitope_collection_builder import EpitopeCollectionBuilder
from .common import check_sequence_dictionary
from .binding_measure import ic50_nM

class RandomBindingPredictor(object):
    def __init__(
            self,
            alleles=['HLA-A*02:01'],
            epitope_lengths=[9]):
        self.alleles = alleles
        self.epitope_lengths = epitope_lengths

    def predict(self, fasta_dictionary):
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        builder = EpitopeCollectionBuilder(
            fasta_dictionary=fasta_dictionary,
            prediction_method_name="random",
            binding_measure=ic50_nM)
        # if wer'e not running the MHC prediction then we have to manually
        # extract 9mer substrings
        for key, sequence in fasta_dictionary.items():
            for epitope_length in self.epitope_lengths:
                for i in range(len(sequence) - epitope_length + 1):
                    for allele in self.alleles:
                        builder.add_binding_prediction(
                            source_sequence_key=key,
                            offset=i,
                            allele=allele,
                            peptide=sequence[i:i + epitope_length],
                            ic50=random.random() * 10000.0,
                            rank=random.randint(0, 99))
        return builder.get_collection()
