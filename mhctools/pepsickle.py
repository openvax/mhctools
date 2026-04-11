# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from .proteasome_predictor import ProteasomePredictor

# Module-level cache for loaded pepsickle models. Keyed by human_only.
_model_cache = {}


class Pepsickle(ProteasomePredictor):
    """
    Proteasomal cleavage predictor using pepsickle's epitope model.

    Uses the in-vivo epitope model from Weeder et al. (Bioinformatics
    2021), which the paper shows outperforms the in-vitro alternatives
    and NetChop for neoantigen identification.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : callable, optional
        See :class:`ProcessingPredictor`.  Default:
        ``score_cterm_anti_max_internal``.

    threshold : float
        Cleavage probability threshold used by pepsickle internally
        (default 0.5).

    human_only : bool
        If True, use human-only trained model instead of all-mammal.
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            threshold=0.5,
            human_only=False):
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.threshold = threshold
        self.human_only = human_only
        self._model = None

    def __str__(self):
        return "%s(scoring=%s)" % (
            self.__class__.__name__,
            getattr(self.scoring, "__name__", repr(self.scoring)))

    def _predictor_name(self):
        return "pepsickle"

    def _load_model(self):
        if self._model is None:
            cache_key = self.human_only
            if cache_key not in _model_cache:
                from pepsickle.model_functions import initialize_epitope_model
                _model_cache[cache_key] = initialize_epitope_model(
                    human_only=self.human_only)
            self._model = _model_cache[cache_key]
        return self._model

    def cleavage_probs(self, sequence):
        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()
        preds_raw = predict_protein_cleavage_locations(
            sequence,
            model,
            mod_type="epitope",
            proteasome_type="C",
            threshold=self.threshold,
        )
        return [entry[2] for entry in preds_raw]
