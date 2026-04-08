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

from .pred import Kind
from .processing_predictor import ProcessingPredictor, resolve_scoring


class ProteasomePredictor(ProcessingPredictor):
    """
    Base class for proteasome-cleavage predictors.

    Defaults *scoring* to ``"cterm_max_internal"``
    (``c_term * (1 - max(internal))``) and emits
    :attr:`Kind.proteasome_cleavage` predictions.

    Subclasses must implement :meth:`cleavage_probs`.
    """

    def __init__(self, default_peptide_lengths=None, scoring=None, **kwargs):
        scoring = resolve_scoring(scoring)
        if scoring is None:
            scoring = "cterm_max_internal"
        ProcessingPredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
            **kwargs,
        )

    def _pred_kind(self):
        return Kind.proteasome_cleavage
