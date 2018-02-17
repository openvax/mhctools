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

from numpy import nan

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .unsupported_allele import UnsupportedAllele

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
            predictor=None,
            models_path=None):
        """
        Parameters
        -----------
        predictor : mhcflurry.Class1AffinityPredictor (optional)
            MHCflurry predictor to use

        models_path : string
            Models dir to use if predictor argument is None

        """
        # moving import here since the mhcflurry package imports
        # Keras and its backend (either Theano or TF) which end up
        # slowing down responsive for any CLI application using MHCtools
        from mhcflurry import Class1AffinityPredictor
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=15)
        if predictor:
            self.predictor = predictor
        elif models_path:
            logging.info("Loading MHCflurry models from %s" % models_path)
            self.predictor = Class1AffinityPredictor.load(models_path)
        else:
            self.predictor = Class1AffinityPredictor.load()

        # relying on BasePredictor and MHCflurry to both normalize
        # allele names the same way using mhcnames
        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)

    def predict_peptides(self, peptides):
        """
        Predict MHC affinity for peptides.
        """

        # importing locally to avoid slowing down CLI applications which
        # don't use MHCflurry
        from mhcflurry.encodable_sequences import EncodableSequences

        binding_predictions = []
        encodable_sequences = EncodableSequences.create(peptides)
        for allele in self.alleles:
            predictions_df = self.predictor.predict_to_dataframe(
                encodable_sequences, allele=allele)
            for (_, row) in predictions_df.iterrows():
                binding_prediction = BindingPrediction(
                    allele=allele,
                    peptide=row.peptide,
                    affinity=row.prediction,
                    percentile_rank=(
                        row.prediction_percentile
                        if 'prediction_percentile' in row else nan),
                    prediction_method_name="mhcflurry"
                )
                binding_predictions.append(binding_prediction)
        return BindingPredictionCollection(binding_predictions)
