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

import logging

from numpy import nan

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Prediction, Kind
from .unsupported_allele import UnsupportedAllele

logger = logging.getLogger(__name__)


class MHCflurry(BasePredictor):
    """
    MHCflurry predictor using the modern Class1PresentationPredictor API.

    Produces both ``pMHC_affinity`` and ``pMHC_presentation`` predictions
    per peptide-allele pair. The legacy ``predict_peptides`` method returns
    BindingPrediction objects based on affinity values for backward compat.

    See https://github.com/openvax/mhcflurry
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
        alleles : list of str

        default_peptide_lengths : list of int

        predictor : mhcflurry.Class1PresentationPredictor (optional)
            MHCflurry presentation predictor to use

        models_path : string
            Models dir to use if predictor argument is None
        """
        # Lazy import to avoid importing Keras/TF at module load time
        from mhcflurry import Class1PresentationPredictor
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
            self.predictor = Class1PresentationPredictor.load(models_path)
        else:
            self.predictor = Class1PresentationPredictor.load()

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)

    def predict_peptides(self, peptides):
        """
        Predict MHC binding affinity and presentation for peptides.

        Returns a BindingPredictionCollection (legacy API) using affinity
        values for backward compatibility.
        """
        binding_predictions = []
        for allele in self.alleles:
            df = self.predictor.predict(
                peptides=list(peptides),
                alleles=[allele],
                include_affinity_percentile=True,
                verbose=0,
            )
            for row in df.itertuples(index=False):
                binding_predictions.append(BindingPrediction(
                    allele=allele,
                    peptide=row.peptide,
                    affinity=row.affinity,
                    percentile_rank=(
                        row.affinity_percentile
                        if hasattr(row, 'affinity_percentile') else nan),
                    prediction_method_name="mhcflurry",
                ))
        return BindingPredictionCollection(binding_predictions)

    def predict(self, peptides):
        """
        Predict for a list of peptide sequences.

        Returns a list of PeptideResult, each containing both
        pMHC_affinity and pMHC_presentation Prediction objects per allele.
        """
        from collections import defaultdict
        from .pred import PeptideResult

        groups = defaultdict(list)
        for allele in self.alleles:
            df = self.predictor.predict(
                peptides=list(peptides),
                alleles=[allele],
                include_affinity_percentile=True,
                verbose=0,
            )
            for row in df.itertuples(index=False):
                pep = row.peptide
                affinity_nM = row.affinity
                affinity_pct = (
                    row.affinity_percentile
                    if hasattr(row, 'affinity_percentile') else None)
                presentation_score = row.presentation_score
                presentation_pct = row.presentation_percentile

                # Higher score = better binding. Convert IC50 nM to 0-1 score:
                # score = 1 - log(IC50)/log(50000), clamped to [0,1]
                import math
                aff_score = max(0.0, min(1.0,
                    1.0 - math.log(max(affinity_nM, 1e-6)) / math.log(50000)))

                groups[pep].append(Prediction(
                    kind=Kind.pMHC_affinity,
                    score=aff_score,
                    peptide=pep,
                    allele=allele,
                    value=affinity_nM,
                    percentile_rank=affinity_pct,
                    predictor_name="mhcflurry",
                ))
                groups[pep].append(Prediction(
                    kind=Kind.pMHC_presentation,
                    score=presentation_score,
                    peptide=pep,
                    allele=allele,
                    percentile_rank=presentation_pct,
                    predictor_name="mhcflurry",
                ))

        return [PeptideResult(preds=tuple(preds)) for preds in groups.values()]

    def _default_pred_kind(self):
        return Kind.pMHC_affinity


class MHCflurry_Affinity(BasePredictor):
    """
    MHCflurry predictor using the older Class1AffinityPredictor API.

    Only produces ``pMHC_affinity`` predictions. Use this if you only need
    affinity scores and don't want the presentation model overhead.

    See https://github.com/openvax/mhcflurry
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
        alleles : list of str

        default_peptide_lengths : list of int

        predictor : mhcflurry.Class1AffinityPredictor (optional)
            MHCflurry affinity predictor to use

        models_path : string
            Models dir to use if predictor argument is None
        """
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

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)

    def predict_peptides(self, peptides):
        """
        Predict MHC affinity for peptides.
        """
        from mhcflurry.encodable_sequences import EncodableSequences

        binding_predictions = []
        encodable_sequences = EncodableSequences.create(peptides)
        for allele in self.alleles:
            predictions_df = self.predictor.predict_to_dataframe(
                encodable_sequences, allele=allele)
            for row in predictions_df.itertuples(index=False):
                binding_predictions.append(BindingPrediction(
                    allele=allele,
                    peptide=row.peptide,
                    affinity=row.prediction,
                    percentile_rank=(
                        row.prediction_percentile
                        if hasattr(row, 'prediction_percentile') else nan),
                    prediction_method_name="mhcflurry",
                ))
        return BindingPredictionCollection(binding_predictions)
