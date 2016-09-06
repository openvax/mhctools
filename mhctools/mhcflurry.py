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

import collections

from mhcflurry import Class1BindingPredictor

from .epitope_collection_builder import EpitopeCollectionBuilder
from .common import check_sequence_dictionary
from .binding_measure import ic50_nM
from .base_predictor import BasePredictor

class MHCFlurry(BasePredictor):
    def __init__(self, alleles, epitope_lengths=[9], valid_alleles=None):
        BasePredictor.__init__(self, alleles, epitope_lengths, valid_alleles)

        def get_predictor(allele):
            return Class1BindingPredictor.from_allele_name(allele)

        self._predictors = dict(
            (allele, get_predictor(allele)) for allele in self.alleles)

    def predict(self, fasta_dictionary):
        fasta_dictionary = check_sequence_dictionary(fasta_dictionary)
        # dict of peptide sequence -> list of (fasta key, offset into sequence)
        # where the peptide came from.
        peptides = collections.defaultdict(list)  
        for (key, sequence) in fasta_dictionary.items():
            for epitope_length in self.epitope_lengths:
                for i in range(len(sequence) - epitope_length + 1):
                    peptide = sequence[i:i + epitope_length]
                    peptides[peptide].append((key, i))
        peptides_list = list(peptides)

        builder = EpitopeCollectionBuilder(
            fasta_dictionary=fasta_dictionary,
            prediction_method_name="mhcflurry",
            binding_measure=ic50_nM)

        for (allele, predictor) in self._predictors.items():
            predictions_df = predictor.predict(peptides_list)
            for (i, record) in predictions_df.iterrows():
                for (key, offset) in peptides[record.Peptide]:
                    # TODO: determine percentile rank somehow.
                    # Currently we just set it to 50.
                    builder.add_binding_prediction(
                            source_sequence_key=key,
                            offset=offset,
                            allele=allele,
                            peptide=record.Peptide,
                            ic50=record.Prediction,
                            rank=50)

        return builder.get_collection()

