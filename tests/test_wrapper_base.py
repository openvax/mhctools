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

"""Tests for the shared new-model wrapper base classes."""

from mhctools.pred import COLUMNS, Kind, Prediction, PeptideResult
from mhctools.wrapper_base import AlleleFreePredictor, NewModelPredictorMixin


def test_normalize_peptides_single_string_and_strip_upper():
    norm = NewModelPredictorMixin._normalize_peptides
    assert norm("siinfekl") == ["SIINFEKL"]          # single string -> list
    assert norm([" gilgfvftl ", "NlV"]) == ["GILGFVFTL", "NLV"]


class _FakeAlleleFree(AlleleFreePredictor):
    mhc_class = "I"

    def _default_pred_kind(self):
        return Kind.immunogenicity

    def predict(self, peptides):
        return [
            PeptideResult(preds=(Prediction(
                kind=Kind.immunogenicity, score=1.0, peptide=p,
                predictor_name="fake"),))
            for p in self._normalize_peptides(peptides)]


def test_allele_free_kind_support_and_supported_kinds():
    p = _FakeAlleleFree()
    support = p.kind_support()[Kind.immunogenicity]
    assert support == {"mhc_dependence": "none", "mhc_class": "I"}
    assert p.supported_kinds == (Kind.immunogenicity,)


def test_inherited_predict_dataframe():
    p = _FakeAlleleFree()
    df = p.predict_dataframe(["siinfekl"], sample_name="s")
    assert list(df["peptide"]) == ["SIINFEKL"]
    assert df.iloc[0]["kind"] == Kind.immunogenicity
    # empty input -> empty frame with the canonical columns
    empty = p.predict_dataframe([])
    assert list(empty.columns) == list(COLUMNS)
    assert len(empty) == 0


def test_wrappers_use_the_shared_base():
    # The standalone wrappers inherit the shared helpers rather than
    # re-implementing them.
    from mhctools.calis import Calis
    from mhctools.bigmhc import BigMHC
    assert issubclass(Calis, AlleleFreePredictor)
    assert issubclass(BigMHC, NewModelPredictorMixin)
    # BigMHC standardizes on _default_pred_kind (not the old _pred_kind).
    assert not hasattr(BigMHC, "_pred_kind")
    assert hasattr(BigMHC, "_default_pred_kind")
