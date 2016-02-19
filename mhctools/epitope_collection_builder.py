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

import numpy as np

from .alleles import normalize_allele_name
from .binding_measure import ic50_nM
from .binding_prediction import BindingPrediction
from .epitope_collection import EpitopeCollection


def invalid_binding_score(x):
    return x < 0 or np.isnan(x) or np.isinf(x)

class EpitopeCollectionBuilder(object):
    """Creates BindingPrediction objects from the fields contained in each
    row of a NetMHC output file.

    Simplifies some of the hassle of mapping epitopes back onto their
    source protein sequences by hanging on to the originating FASTA dictionary
    and using sequence keys to look up the original sequences.
    """
    def __init__(
            self,
            fasta_dictionary,
            prediction_method_name,
            binding_measure=ic50_nM):
        self.fasta_dictionary = fasta_dictionary
        self.prediction_method_name = prediction_method_name
        self.binding_measure = binding_measure
        self.binding_predictions = []

    def add_binding_prediction(
            self,
            source_sequence_key,
            offset,
            peptide,
            allele,
            ic50,
            rank,
            log_ic50=None):
        """
        Parameters
        ----------
        source_sequence_key : str
            Unique identifier for source sequence

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
        """
        # if we have a bad IC50 score we might still get a salvageable
        # log of the score. Strangely, this is necessary sometimes!
        if invalid_binding_score(ic50) and log_ic50 is not None:
            ic50 = 50000 ** (-log_ic50 + 1)
        # if IC50 is still NaN or otherwise invalid, abort
        if invalid_binding_score(ic50):
            raise ValueError(
                "Invalid IC50 value %0.4f for %s w/ allele %s" % (
                    ic50,
                    peptide,
                    allele))

        if invalid_binding_score(rank) or rank > 100:
            raise ValueError(
                "Invalid percentile rank %s for %s w/ allele %s" % (
                    rank, peptide, allele))

        if source_sequence_key not in self.fasta_dictionary:
            raise ValueError("Sequence identifier %s not found in %s" % (
                source_sequence_key, self.fasta_dictionary))

        source_sequence = str(self.fasta_dictionary[source_sequence_key])

        binding_prediction = BindingPrediction(
            source_sequence_key=source_sequence_key,
            source_sequence=source_sequence,
            offset=offset,
            allele=normalize_allele_name(allele),
            peptide=peptide,
            length=len(peptide),
            value=ic50,
            percentile_rank=rank,
            prediction_method_name=self.prediction_method_name,
            measure=self.binding_measure)
        self.binding_predictions.append(binding_prediction)

    def get_collection(self):
        return EpitopeCollection(self.binding_predictions)

    def get_dataframe(self):
        return self.get_collection().dataframe()
