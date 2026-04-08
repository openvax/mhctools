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

from collections import defaultdict

from .base_predictor import BasePredictor
from .pred import Pred, Kind, PeptidePreds


class Pepsickle(BasePredictor):
    """
    Wrapper around the pepsickle proteasomal cleavage predictor.

    Pepsickle predicts proteasomal cleavage probabilities for each position
    in a protein sequence. This wrapper exposes those predictions through
    the standard mhctools predictor API.

    For ``predict(peptides)``, returns the C-terminal cleavage probability
    for each peptide (i.e. the probability that the proteasome cleaves after
    the last residue, which is the relevant cleavage event for generating
    that peptide).

    For ``predict_proteins()``, runs pepsickle on the full protein so that
    flanking-sequence context is available, then extracts per-peptide
    C-terminal cleavage scores.

    Parameters
    ----------
    alleles : list of str, optional
        Not used by cleavage predictors but accepted for API compatibility.

    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    model_type : str
        Pepsickle model to use. One of ``"epitope"`` (default),
        ``"in-vitro"`` (gradient-boosted), or ``"in-vitro-2"`` (neural net).

    proteasome_type : str
        ``"C"`` for constitutive (default) or ``"I"`` for immunoproteasome.
        Only used by in-vitro models; ignored by the epitope model.

    threshold : float
        Cleavage probability threshold (default 0.5).

    human_only : bool
        If True, use human-only trained models instead of all-mammal models.
    """

    VALID_MODEL_TYPES = ("epitope", "in-vitro", "in-vitro-2")
    VALID_PROTEASOME_TYPES = ("C", "I")

    def __init__(
            self,
            alleles=None,
            default_peptide_lengths=None,
            model_type="epitope",
            proteasome_type="C",
            threshold=0.5,
            human_only=False):
        if alleles is None:
            alleles = []
        if default_peptide_lengths is None:
            default_peptide_lengths = [9]
        if model_type not in self.VALID_MODEL_TYPES:
            raise ValueError(
                "model_type must be one of %s, got %r" % (
                    self.VALID_MODEL_TYPES, model_type))
        if proteasome_type not in self.VALID_PROTEASOME_TYPES:
            raise ValueError(
                "proteasome_type must be 'C' or 'I', got %r" % proteasome_type)
        BasePredictor.__init__(
            self,
            alleles=alleles,
            default_peptide_lengths=default_peptide_lengths,
            min_peptide_length=None,
            max_peptide_length=None,
            allow_X_in_peptides=True,
        )
        self.model_type = model_type
        self.proteasome_type = proteasome_type
        self.threshold = threshold
        self.human_only = human_only
        self._model = None

    def __str__(self):
        return "%s(model_type=%r, proteasome_type=%r)" % (
            self.__class__.__name__,
            self.model_type,
            self.proteasome_type)

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

    def _default_pred_kind(self):
        return Kind.proteasome_cleavage

    def predict(self, peptides):
        """
        Predict C-terminal cleavage probability for each peptide.

        Parameters
        ----------
        peptides : list of str

        Returns
        -------
        list of PeptidePreds
        """
        self._check_peptide_inputs(peptides)
        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()
        results = []
        for peptide in peptides:
            preds_raw = predict_protein_cleavage_locations(
                peptide,
                model,
                mod_type=self.model_type,
                proteasome_type=self.proteasome_type,
                threshold=self.threshold,
            )
            # C-terminal cleavage score is at the last position
            if preds_raw:
                score = preds_raw[-1][2]
            else:
                score = 0.0
            pred = Pred(
                kind=Kind.proteasome_cleavage,
                score=score,
                peptide=peptide,
                predictor_name="pepsickle",
            )
            results.append(PeptidePreds(preds=(pred,)))
        return results

    def predict_proteins(self, sequence_dict, peptide_lengths=None):
        """
        Predict cleavage for peptides derived from full protein sequences.

        Runs pepsickle on each full protein (preserving flanking context)
        and extracts the C-terminal cleavage score for every sub-peptide.

        Parameters
        ----------
        sequence_dict : dict or str
            Mapping of sequence names to amino acid strings, or a single
            sequence string.

        peptide_lengths : list of int, optional

        Returns
        -------
        dict mapping sequence_name -> list of PeptidePreds
        """
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}
        elif isinstance(sequence_dict, (list, tuple)):
            sequence_dict = {seq: seq for seq in sequence_dict}

        peptide_lengths = self._check_peptide_lengths(peptide_lengths)

        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()

        results = defaultdict(list)
        for name, sequence in sequence_dict.items():
            preds_raw = predict_protein_cleavage_locations(
                sequence,
                model,
                protein_id=name,
                mod_type=self.model_type,
                proteasome_type=self.proteasome_type,
                threshold=self.threshold,
            )
            # Build 1-based position -> cleavage probability map
            cleavage_scores = {entry[0]: entry[2] for entry in preds_raw}

            for peptide_length in peptide_lengths:
                for i in range(len(sequence) - peptide_length + 1):
                    peptide = sequence[i:i + peptide_length]
                    # C-terminal residue is at 1-based position (i + peptide_length)
                    c_term_pos = i + peptide_length
                    score = cleavage_scores.get(c_term_pos, 0.0)
                    pred = Pred(
                        kind=Kind.proteasome_cleavage,
                        score=score,
                        peptide=peptide,
                        source_sequence_name=name,
                        offset=i,
                        predictor_name="pepsickle",
                    )
                    results[name].append(PeptidePreds(preds=(pred,)))
        return dict(results)

    def predict_cleavage_sites(self, sequence_dict):
        """
        Return per-position cleavage probabilities for full protein sequences.

        This is a pepsickle-specific method that exposes the full cleavage
        profile rather than summarising it per peptide.

        Parameters
        ----------
        sequence_dict : dict or str
            Mapping of sequence names to amino acid strings, or a single
            sequence string.

        Returns
        -------
        dict mapping sequence_name -> list of (position, residue, probability)
            position is 1-based.
        """
        if isinstance(sequence_dict, str):
            sequence_dict = {"seq": sequence_dict}

        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()

        results = {}
        for name, sequence in sequence_dict.items():
            preds_raw = predict_protein_cleavage_locations(
                sequence,
                model,
                protein_id=name,
                mod_type=self.model_type,
                proteasome_type=self.proteasome_type,
                threshold=self.threshold,
            )
            results[name] = [
                (entry[0], entry[1], entry[2]) for entry in preds_raw
            ]
        return results
