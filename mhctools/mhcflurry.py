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
import math
import os

from numpy import nan

from .base_predictor import BasePredictor
from .binding_prediction import BindingPrediction
from .binding_prediction_collection import BindingPredictionCollection
from .pred import Prediction, Kind
from .unsupported_allele import UnsupportedAllele

logger = logging.getLogger(__name__)

# Module-level cache for loaded models. Keyed by (kind, normalized_path).
_model_cache = {}


def _normalize_models_path(models_path):
    """Normalize a models directory path for cache-key deduplication.

    Different strings pointing at the same directory (relative vs.
    absolute, with ``~``, via a symlink) should share a single cache
    entry. ``None`` stays ``None`` — mhcflurry resolves that to the
    package-default path internally, and we treat it as its own key.
    """
    if models_path is None:
        return None
    return os.path.realpath(os.path.expanduser(models_path))


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
        from mhcflurry import Class1PresentationPredictor
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=8,
            max_peptide_length=15)
        if predictor:
            self.predictor = predictor
        else:
            cache_key = ("presentation", _normalize_models_path(models_path))
            if cache_key not in _model_cache:
                if models_path:
                    logger.info(
                        "Loading MHCflurry models from %s", models_path)
                    _model_cache[cache_key] = \
                        Class1PresentationPredictor.load(models_path)
                else:
                    _model_cache[cache_key] = \
                        Class1PresentationPredictor.load()
            self.predictor = _model_cache[cache_key]

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)

    def predict_peptides(self, peptides):
        """
        Predict MHC binding affinity and presentation for peptides.

        Returns a BindingPredictionCollection (legacy API) using affinity
        values for backward compatibility.
        """
        peptide_list = list(peptides)
        allele_list = list(self.alleles)

        # Build cross product for batch prediction
        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]

        df = self.predictor.affinity_predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
        )
        binding_predictions = []
        for row in df.itertuples(index=False):
            binding_predictions.append(BindingPrediction(
                allele=row.allele,
                peptide=row.peptide,
                affinity=row.prediction,
                percentile_rank=(
                    row.prediction_percentile
                    if hasattr(row, 'prediction_percentile') else nan),
                prediction_method_name="mhcflurry",
            ))
        return BindingPredictionCollection(binding_predictions)

    def predict(self, peptides):
        """
        Predict for a list of peptide sequences.

        Returns a list of PeptideResult, each containing both
        pMHC_affinity and pMHC_presentation Prediction objects per allele.

        Uses batch prediction across all alleles in a single call for
        both affinity and presentation scores.
        """
        from collections import defaultdict
        from .pred import PeptideResult

        peptide_list = list(peptides)
        allele_list = list(self.alleles)

        # Build cross product
        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]

        # Single batched call for affinity
        aff_df = self.predictor.affinity_predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
        )

        # Per-allele presentation calls (presentation predictor does
        # deconvolution across alleles, so we call per-allele to get
        # per-allele presentation scores). Key by mhcflurry's output
        # allele string so lookups with aff_df.allele always match.
        pres_by_pep_allele = {}
        for input_allele in allele_list:
            df = self.predictor.predict(
                peptides=peptide_list,
                alleles=[input_allele],
                include_affinity_percentile=False,
                verbose=0,
            )
            for row in df.itertuples(index=False):
                output_allele = getattr(row, 'allele', input_allele)
                pres_by_pep_allele[(row.peptide, output_allele)] = (
                    row.presentation_score,
                    row.presentation_percentile,
                )

        groups = defaultdict(list)
        for row in aff_df.itertuples(index=False):
            pep = row.peptide
            allele = row.allele
            affinity_nM = row.prediction
            affinity_pct = (
                row.prediction_percentile
                if hasattr(row, 'prediction_percentile') else None)

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

            key = (pep, allele)
            if key not in pres_by_pep_allele:
                raise ValueError(
                    "MHCflurry: missing presentation score for "
                    "peptide='%s' allele='%s' (this indicates an allele or "
                    "peptide string mismatch between the affinity and "
                    "presentation predictor outputs)" % (pep, allele))
            pres_score, pres_pct = pres_by_pep_allele[key]
            groups[pep].append(Prediction(
                kind=Kind.pMHC_presentation,
                score=pres_score,
                peptide=pep,
                allele=allele,
                percentile_rank=pres_pct,
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
        else:
            cache_key = ("affinity", _normalize_models_path(models_path))
            if cache_key not in _model_cache:
                if models_path:
                    logger.info(
                        "Loading MHCflurry models from %s", models_path)
                    _model_cache[cache_key] = \
                        Class1AffinityPredictor.load(models_path)
                else:
                    _model_cache[cache_key] = \
                        Class1AffinityPredictor.load()
            self.predictor = _model_cache[cache_key]

        for allele in self.alleles:
            if allele not in self.predictor.supported_alleles:
                raise UnsupportedAllele(allele)

    def predict_peptides(self, peptides):
        """
        Predict MHC affinity for peptides.
        """
        peptide_list = list(peptides)
        allele_list = list(self.alleles)

        batch_peptides = peptide_list * len(allele_list)
        batch_alleles = [a for a in allele_list for _ in peptide_list]

        df = self.predictor.predict_to_dataframe(
            peptides=batch_peptides,
            alleles=batch_alleles,
        )
        binding_predictions = []
        for row in df.itertuples(index=False):
            binding_predictions.append(BindingPrediction(
                allele=row.allele,
                peptide=row.peptide,
                affinity=row.prediction,
                percentile_rank=(
                    row.prediction_percentile
                    if hasattr(row, 'prediction_percentile') else nan),
                prediction_method_name="mhcflurry",
            ))
        return BindingPredictionCollection(binding_predictions)
