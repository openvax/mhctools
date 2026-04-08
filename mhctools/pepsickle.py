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


class Pepsickle(ProteasomePredictor):
    """
    Wrapper around the pepsickle proteasomal cleavage predictor.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : callable, optional
        See :class:`ProcessingPredictor`.  Default:
        ``score_cterm_anti_max_internal``.

    model_type : str
        Pepsickle model to use. One of ``"epitope"`` (default),
        ``"in-vitro"`` (gradient-boosted), or ``"in-vitro-2"`` (neural net).

    proteasome_type : str
        ``"C"`` for constitutive (default) or ``"I"`` for immunoproteasome.
        Only used by in-vitro models; ignored by the epitope model.

    threshold : float
        Cleavage probability threshold used by pepsickle internally
        (default 0.5).

    human_only : bool
        If True, use human-only trained models instead of all-mammal.
    """

    VALID_MODEL_TYPES = ("epitope", "in-vitro", "in-vitro-2")
    VALID_PROTEASOME_TYPES = ("C", "I")

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            model_type="epitope",
            proteasome_type="C",
            threshold=0.5,
            human_only=False):
        if model_type not in self.VALID_MODEL_TYPES:
            raise ValueError(
                "model_type must be one of %s, got %r" % (
                    self.VALID_MODEL_TYPES, model_type))
        if proteasome_type not in self.VALID_PROTEASOME_TYPES:
            raise ValueError(
                "proteasome_type must be 'C' or 'I', got %r" % proteasome_type)
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.model_type = model_type
        self.proteasome_type = proteasome_type
        self.threshold = threshold
        self.human_only = human_only
        self._model = None

    def __str__(self):
        return "%s(model_type=%r, proteasome_type=%r, scoring=%s)" % (
            self.__class__.__name__,
            self.model_type,
            self.proteasome_type,
            getattr(self.scoring, "__name__", repr(self.scoring)))

    def _predictor_name(self):
        return "pepsickle"

    def _load_model(self):
        if self._model is None:
            from pepsickle.model_functions import (
                initialize_epitope_model,
                initialize_digestion_model,
                initialize_digestion_gb_model,
            )
            if self.model_type == "epitope":
                self._model = initialize_epitope_model(
                    human_only=self.human_only)
            elif self.model_type == "in-vitro":
                self._model = initialize_digestion_gb_model()
            elif self.model_type == "in-vitro-2":
                self._model = initialize_digestion_model(
                    human_only=self.human_only)
        return self._model

    def cleavage_probs(self, sequence):
        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()
        preds_raw = predict_protein_cleavage_locations(
            sequence,
            model,
            mod_type=self.model_type,
            proteasome_type=self.proteasome_type,
            threshold=self.threshold,
        )
        return [entry[2] for entry in preds_raw]
