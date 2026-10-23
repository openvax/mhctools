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

import pytest

from mhctools.pred import COLUMNS, Kind, Prediction, PeptideResult
from mhctools.wrapper_base import AlleleFreePredictor, NewModelPredictorMixin
from mhctools.bigmhc import BigMHC, BigMHC_EL, BigMHC_IM


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


@pytest.mark.parametrize("factory,kind", [
    (BigMHC_EL, Kind.pMHC_presentation),
    (BigMHC_IM, Kind.immunogenicity),
])
def test_bigmhc_fixed_mode_supported_kinds(factory, kind, tmp_path):
    # Construction checks the source path; model loading stays lazy.
    source = tmp_path / "src"
    source.mkdir()
    (source / "bigmhc.py").touch()
    predictor = factory(alleles=["HLA-A*02:01"], bigmhc_path=str(tmp_path))
    assert predictor.supported_kinds == (kind,)
    assert predictor.kind_support() == {
        kind: {"mhc_dependence": "single_allele", "mhc_class": "I"},
    }


def test_bigmhc_legacy_kind_hook_stays_removed():
    assert not hasattr(BigMHC, "_pred_kind")
