# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

"""Shared helpers for predictors that emit the new ``PeptideResult`` model
directly, rather than going through ``BasePredictor``/``BindingPrediction``.

``NewModelPredictorMixin`` holds the members that were otherwise copy-pasted
across the standalone wrappers (``predict_dataframe``, ``supported_kinds``,
``__repr__``, peptide normalization). ``AlleleFreePredictor`` adds the
allele-independent ``kind_support`` shape used by Calis, DeepTAP, and ERAMER.
"""

import pandas as pd

from .pred import COLUMNS


class NewModelPredictorMixin:
    """Members shared by every new-model wrapper (allele-free or not)."""

    def __repr__(self):
        return str(self)

    @property
    def supported_kinds(self):
        """Prediction kind strings this predictor can emit."""
        return tuple(self.kind_support())

    @staticmethod
    def _normalize_peptides(peptides):
        """Accept a single string or an iterable; strip + upper-case each."""
        if isinstance(peptides, str):
            peptides = [peptides]
        return [str(p).strip().upper() for p in peptides]

    def predict_dataframe(self, peptides, sample_name="", n_flanks=None,
                          c_flanks=None):
        """``predict()`` flattened to a DataFrame.

        ``n_flanks``/``c_flanks`` are accepted for a uniform flank-aware API;
        allele-free wrappers do not use flanking context and ignore them.
        """
        dfs = [pp.to_dataframe(sample_name) for pp in self.predict(peptides)]
        if not dfs:
            return pd.DataFrame(columns=COLUMNS)
        return pd.concat(dfs, ignore_index=True)


class AlleleFreePredictor(NewModelPredictorMixin):
    """Base for allele-independent new-model predictors.

    Subclasses set ``mhc_class`` and implement ``_default_pred_kind()`` plus the
    scoring in ``predict()``; peptide length/character validation stays with the
    subclass since each tool has its own bounds and messages.
    """

    mhc_class = "none"

    def kind_support(self):
        return {
            self._default_pred_kind(): {
                "mhc_dependence": "none",
                "mhc_class": self.mhc_class,
            },
        }
