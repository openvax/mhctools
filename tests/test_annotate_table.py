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

"""Binary-free tests for mhctools.annotate.annotate_table.

Uses a deterministic in-memory predictor so best-allele selection and
direction handling can be asserted exactly (RandomBindingPredictor returns
fresh random scores on each call, so it can only be used for a structural
smoke test).
"""

import math

import pandas as pd
import pytest

from mhctools import (
    annotate_table,
    AnnotationSpec,
    RandomBindingPredictor,
)
from mhctools.annotate import parse_annotation_spec, output_field_tokens
from mhctools.pred import Kind, PeptideResult, Prediction


# A fixed (peptide, allele) -> (affinity nM, score, percentile_rank) table.
# Chosen so the best allele differs by field, to catch direction bugs:
#   SIINFEKL: A*02:01 has lower IC50 (better affinity) but LOWER score;
#             B*07:02 has higher IC50 but HIGHER score.
_FIXTURE = {
    ("SIINFEKL", "HLA-A*02:01"): (100.0, 0.20, 5.0),
    ("SIINFEKL", "HLA-B*07:02"): (500.0, 0.90, 1.0),
    ("GILGFVFTL", "HLA-A*02:01"): (50.0, 0.80, 0.5),
    ("GILGFVFTL", "HLA-B*07:02"): (9000.0, 0.05, 40.0),
}


class _FixturePredictor:
    """Deterministic predictor returning predictions from ``_FIXTURE``.

    Only implements the ``predict`` method that ``annotate_table`` calls.
    Emits pMHC_affinity predictions carrying value (IC50), score, and rank.
    """

    def __init__(self, alleles):
        self.alleles = list(alleles) if alleles else []

    def predict(self, peptides):
        results = []
        for peptide in peptides:
            preds = []
            for allele in self.alleles:
                key = (peptide, allele)
                if key not in _FIXTURE:
                    continue
                affinity, score, rank = _FIXTURE[key]
                preds.append(Prediction(
                    kind=Kind.pMHC_affinity,
                    peptide=peptide,
                    allele=allele,
                    score=score,
                    value=affinity,
                    percentile_rank=rank,
                    predictor_name="fixture"))
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def _factory(alleles):
    return _FixturePredictor(alleles)


# A fixed peptide -> processing score table for the allele-free path.
_PROCESSING = {"SIINFEKL": 0.80, "GILGFVFTL": 0.30}


class _ProcessingFixturePredictor:
    """Deterministic allele-free predictor (like a processing predictor).

    Emits antigen_processing predictions with no allele, so annotate_table
    routes them through the by-peptide lookup rather than the (peptide,
    allele) lookup.
    """

    def __init__(self, alleles=None):
        pass

    def predict(self, peptides):
        results = []
        for peptide in peptides:
            preds = []
            if peptide in _PROCESSING:
                preds.append(Prediction(
                    kind=Kind.antigen_processing,
                    peptide=peptide,
                    allele="",
                    score=_PROCESSING[peptide],
                    predictor_name="processing-fixture"))
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def _table():
    return pd.DataFrame({
        "sample_id": ["s1", "s2"],
        "hit": [1, 0],
        "peptide": ["SIINFEKL", "GILGFVFTL"],
        "hla": ["HLA-A*02:01 HLA-B*07:02", "A0201,B0702"],
    })


def test_affinity_picks_lowest_ic50():
    out = annotate_table(
        _table(),
        [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    row = out[out.peptide == "SIINFEKL"].iloc[0]
    # A*02:01 has the lower IC50 (100 < 500) -> better affinity
    assert row["aff"] == 100.0
    assert row["aff_best_allele"] == "HLA-A*02:01"
    row2 = out[out.peptide == "GILGFVFTL"].iloc[0]
    assert row2["aff"] == 50.0
    assert row2["aff_best_allele"] == "HLA-A*02:01"


def test_score_picks_highest_and_differs_from_affinity():
    out = annotate_table(
        _table(),
        [AnnotationSpec(_factory, "sc", field="score")],
        allele_column="hla")
    row = out[out.peptide == "SIINFEKL"].iloc[0]
    # B*07:02 has the higher score (0.90 > 0.20) even though worse affinity
    assert row["sc"] == 0.90
    assert row["sc_best_allele"] == "HLA-B*07:02"


def test_percentile_rank_picks_lowest():
    out = annotate_table(
        _table(),
        [AnnotationSpec(_factory, "pr", field="percentile_rank")],
        allele_column="hla")
    row = out[out.peptide == "SIINFEKL"].iloc[0]
    # B*07:02 has the lower rank (1.0 < 5.0) -> better
    assert row["pr"] == 1.0
    assert row["pr_best_allele"] == "HLA-B*07:02"


def test_preserves_input_columns_and_order():
    df = _table()
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    # original columns come first, unchanged, then the two appended columns
    assert list(out.columns) == list(df.columns) + ["aff", "aff_best_allele"]
    pd.testing.assert_frame_equal(out[df.columns], df)


def test_does_not_mutate_input():
    df = _table()
    before = df.copy()
    annotate_table(df, [AnnotationSpec(_factory, "aff")], allele_column="hla")
    pd.testing.assert_frame_equal(df, before)


def test_multiple_specs_appended_independently():
    out = annotate_table(
        _table(),
        [AnnotationSpec(_factory, "aff", field="affinity"),
         AnnotationSpec(_factory, "sc", field="score", add_best_allele=False)],
        allele_column="hla")
    assert "aff" in out.columns and "aff_best_allele" in out.columns
    assert "sc" in out.columns
    # add_best_allele=False suppresses the provenance column
    assert "sc_best_allele" not in out.columns


def test_missing_allele_yields_nan():
    df = pd.DataFrame({
        "peptide": ["SIINFEKL"],
        "hla": ["HLA-C*07:01"],  # not in the fixture
    })
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    assert math.isnan(out.iloc[0]["aff"])
    assert out.iloc[0]["aff_best_allele"] is None


def test_unknown_peptide_yields_nan():
    df = pd.DataFrame({
        "peptide": ["WWWWWWWWW"],  # not in the fixture
        "hla": ["HLA-A*02:01"],
    })
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    assert math.isnan(out.iloc[0]["aff"])


def test_empty_allele_cell_yields_nan():
    df = pd.DataFrame({
        "peptide": ["SIINFEKL", "GILGFVFTL"],
        "hla": ["HLA-A*02:01", None],
    })
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    assert out.iloc[0]["aff"] == 100.0
    assert math.isnan(out.iloc[1]["aff"])
    assert out.iloc[1]["aff_best_allele"] is None


def test_column_collision_raises():
    df = _table()
    df["aff"] = 0.0
    with pytest.raises(ValueError, match="already exists"):
        annotate_table(
            df, [AnnotationSpec(_factory, "aff", field="affinity")],
            allele_column="hla")


def test_column_collision_overwrite_allowed():
    df = _table()
    df["aff"] = 0.0
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla", overwrite=True)
    assert out.iloc[0]["aff"] == 100.0


def test_missing_peptide_column_raises():
    with pytest.raises(KeyError, match="peptide"):
        annotate_table(
            pd.DataFrame({"seq": ["SIINFEKL"]}),
            [AnnotationSpec(_factory, "aff")])


def test_missing_allele_column_raises():
    with pytest.raises(KeyError, match="nope"):
        annotate_table(
            _table(), [AnnotationSpec(_factory, "aff")],
            allele_column="nope")


def test_non_dataframe_raises():
    with pytest.raises(TypeError, match="DataFrame"):
        annotate_table("not a df", [AnnotationSpec(_factory, "aff")])


def test_custom_best_allele_column_name():
    out = annotate_table(
        _table(),
        [AnnotationSpec(_factory, "aff", field="affinity",
                        best_allele_column="winner")],
        allele_column="hla")
    assert "winner" in out.columns
    assert "aff_best_allele" not in out.columns


def test_prebuilt_predictor_instance_used_directly():
    # Passing a built predictor (not a factory) should work: it is used as-is.
    predictor = _FixturePredictor(["HLA-A*02:01", "HLA-B*07:02"])
    out = annotate_table(
        _table(), [AnnotationSpec(predictor, "aff", field="affinity")],
        allele_column="hla")
    assert out[out.peptide == "SIINFEKL"].iloc[0]["aff"] == 100.0


# --- duplicate output-column guard (cross-spec) -----------------------------

def test_duplicate_output_column_across_specs_raises():
    with pytest.raises(ValueError, match="more than one spec"):
        annotate_table(
            _table(),
            [AnnotationSpec(_factory, "dup", field="affinity"),
             AnnotationSpec(_factory, "dup", field="score")],
            allele_column="hla")


def test_duplicate_output_column_across_specs_raises_even_with_overwrite():
    # overwrite lets you replace an *existing* column, but two specs writing
    # the same new column would still clobber each other -> always an error.
    with pytest.raises(ValueError, match="more than one spec"):
        annotate_table(
            _table(),
            [AnnotationSpec(_factory, "dup", field="affinity"),
             AnnotationSpec(_factory, "dup", field="score")],
            allele_column="hla", overwrite=True)


def test_duplicate_best_allele_column_across_specs_raises():
    # distinct score columns, but the same explicit provenance column
    with pytest.raises(ValueError, match="more than one spec"):
        annotate_table(
            _table(),
            [AnnotationSpec(_factory, "a", field="affinity",
                            best_allele_column="prov"),
             AnnotationSpec(_factory, "b", field="score",
                            best_allele_column="prov")],
            allele_column="hla")


def test_output_column_collides_with_other_spec_best_allele_raises():
    # spec A's score column name equals spec B's provenance column name
    with pytest.raises(ValueError, match="more than one spec"):
        annotate_table(
            _table(),
            [AnnotationSpec(_factory, "b_best_allele", field="affinity"),
             AnnotationSpec(_factory, "b", field="score")],
            allele_column="hla")


# --- peptide normalization (case / whitespace insensitive matching) ---------

def test_peptide_matched_case_insensitively():
    df = pd.DataFrame({"peptide": ["siinfekl"], "hla": ["HLA-A*02:01"]})
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    # matched the fixture (keyed uppercase) despite the lowercase cell
    assert out.iloc[0]["aff"] == 100.0
    # the input peptide column is preserved verbatim, not uppercased
    assert out.iloc[0]["peptide"] == "siinfekl"


def test_peptide_matched_ignoring_surrounding_whitespace():
    df = pd.DataFrame({"peptide": [" SIINFEKL "], "hla": ["HLA-A*02:01"]})
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    assert out.iloc[0]["aff"] == 100.0
    assert out.iloc[0]["peptide"] == " SIINFEKL "


# --- allele-free predictor path (by-peptide lookup) -------------------------

def test_processing_field_token_pins_antigen_processing_kind():
    # the 'processing' field selects Kind.antigen_processing (higher-better)
    spec = AnnotationSpec(_factory, "col", field="processing")
    assert spec.kind == Kind.antigen_processing
    assert spec.prediction_field == "score"
    df = pd.DataFrame({"peptide": ["SIINFEKL", "GILGFVFTL", "WWWWWWWWW"]})
    out = annotate_table(
        df,
        [AnnotationSpec(lambda alleles: _ProcessingFixturePredictor(),
                        "proc", field="processing")])
    assert out.iloc[0]["proc"] == 0.80
    assert out.iloc[1]["proc"] == 0.30
    assert math.isnan(out.iloc[2]["proc"])


def test_allele_free_predictor_uses_by_peptide_lookup():
    df = pd.DataFrame({"peptide": ["SIINFEKL", "GILGFVFTL", "WWWWWWWWW"]})
    out = annotate_table(
        df,
        [AnnotationSpec(lambda alleles: _ProcessingFixturePredictor(),
                        "proc", field="score")])
    # no allele column -> no provenance column
    assert "proc" in out.columns
    assert "proc_best_allele" not in out.columns
    assert out.iloc[0]["proc"] == 0.80
    assert out.iloc[1]["proc"] == 0.30
    # peptide absent from the fixture -> NaN
    assert math.isnan(out.iloc[2]["proc"])


def test_allele_free_predictor_works_even_with_allele_column_present():
    # An allele-free predictor emits allele-less predictions; when a (peptide,
    # allele) lookup finds nothing, we fall back to the by-peptide prediction,
    # so the score is still filled in rather than silently NaN.
    df = pd.DataFrame({"peptide": ["SIINFEKL"], "hla": ["HLA-A*02:01"]})
    out = annotate_table(
        df,
        [AnnotationSpec(lambda alleles: _ProcessingFixturePredictor(),
                        "proc", field="score")],
        allele_column="hla")
    assert out.iloc[0]["proc"] == 0.80
    # the winning prediction had no allele -> provenance is None
    assert out.iloc[0]["proc_best_allele"] is None


def test_binding_predictor_unsupported_allele_still_nan():
    # Regression guard: the by-peptide fallback must NOT rescue a binding
    # predictor's unsupported-allele miss (it never populates by_peptide).
    df = pd.DataFrame({"peptide": ["SIINFEKL"], "hla": ["HLA-C*07:01"]})
    out = annotate_table(
        df, [AnnotationSpec(_factory, "aff", field="affinity")],
        allele_column="hla")
    assert math.isnan(out.iloc[0]["aff"])
    assert out.iloc[0]["aff_best_allele"] is None


# --- spec parsing -----------------------------------------------------------

def test_parse_spec_full():
    spec = parse_annotation_spec("netmhcpan42-ba:mycol:affinity")
    assert spec.output_column == "mycol"
    assert spec.field == "affinity"
    assert spec.kind == Kind.pMHC_affinity
    assert spec.prediction_field == "value"


def test_parse_spec_defaults_column_and_field():
    spec = parse_annotation_spec("mixmhcpred")
    assert spec.field == "affinity"
    assert spec.output_column == "mixmhcpred_affinity"


def test_parse_spec_default_field_only():
    spec = parse_annotation_spec("netmhcpan42-el:elcol")
    assert spec.output_column == "elcol"
    assert spec.field == "affinity"


def test_parse_spec_unknown_predictor():
    with pytest.raises(ValueError, match="Unknown predictor"):
        parse_annotation_spec("does-not-exist:c:affinity")


def test_parse_spec_unknown_field():
    with pytest.raises(ValueError, match="Unknown output field"):
        parse_annotation_spec("mixmhcpred:c:bogus")


def test_bad_field_in_spec_object_raises():
    with pytest.raises(ValueError, match="Unknown output field"):
        AnnotationSpec(_factory, "col", field="nonsense")


def test_output_field_tokens_include_primary_three():
    tokens = output_field_tokens()
    for expected in ("affinity", "score", "percentile_rank"):
        assert expected in tokens


# --- stability token reads `score` (NetMHCstabpan puts Thalf there) ---------

class _StabilityFixturePredictor:
    """Allele-bearing predictor emitting pMHC_stability with half-life in
    ``score`` and no ``value`` — mirroring parse_netmhcstabpan."""

    def __init__(self, alleles):
        self.alleles = list(alleles) if alleles else []

    def predict(self, peptides):
        results = []
        for peptide in peptides:
            preds = [
                Prediction(
                    kind=Kind.pMHC_stability,
                    peptide=peptide,
                    allele=allele,
                    score=12.5,          # half-life (hours), higher = better
                    value=None,          # NetMHCstabpan sets no ic50
                    predictor_name="stab-fixture")
                for allele in self.alleles]
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def test_stability_token_reads_score_not_value():
    # Regression: the `stability` token used to map to `value`, which
    # NetMHCstabpan never populates, so the column was always NaN.
    spec = parse_annotation_spec("netmhcstabpan:stab:stability")
    assert (spec.kind, spec.prediction_field) == (Kind.pMHC_stability, "score")

    out = annotate_table(
        _table(),
        [AnnotationSpec(lambda a: _StabilityFixturePredictor(a), "stab",
                        field="stability")],
        allele_column="hla")
    assert list(out["stab"]) == [12.5, 12.5]  # not NaN


# --- generic score/rank token is rejected on multi-kind predictors ----------

class _MultiKindFixturePredictor:
    """Emits BOTH affinity and presentation for the same (peptide, allele)."""

    def __init__(self, alleles):
        self.alleles = list(alleles) if alleles else []

    def predict(self, peptides):
        results = []
        for peptide in peptides:
            preds = []
            for allele in self.alleles:
                preds.append(Prediction(
                    kind=Kind.pMHC_affinity, peptide=peptide, allele=allele,
                    score=0.2, value=100.0, predictor_name="multi"))
                preds.append(Prediction(
                    kind=Kind.pMHC_presentation, peptide=peptide, allele=allele,
                    score=0.9, predictor_name="multi"))
            results.append(PeptideResult(preds=tuple(preds)))
        return results


def test_generic_score_token_ambiguous_multikind_raises():
    with pytest.raises(ValueError, match="[Aa]mbiguous"):
        annotate_table(
            _table(),
            [AnnotationSpec(lambda a: _MultiKindFixturePredictor(a), "s",
                            field="score")],
            allele_column="hla")


def test_kind_specific_token_unambiguous_on_multikind():
    # The same multi-kind predictor is fine with a kind-specific token.
    out = annotate_table(
        _table(),
        [AnnotationSpec(lambda a: _MultiKindFixturePredictor(a), "pres",
                        field="presentation")],
        allele_column="hla")
    assert list(out["pres"]) == [0.9, 0.9]


# --- integration-flavored smoke test with the real random predictor ---------

def test_random_predictor_smoke():
    df = _table()
    out = annotate_table(
        df,
        [AnnotationSpec(lambda a: RandomBindingPredictor(alleles=a),
                        "rand", field="affinity")],
        allele_column="hla")
    assert list(out.columns) == list(df.columns) + ["rand", "rand_best_allele"]
    # every row got a numeric prediction and a best allele from its own set
    for _, row in out.iterrows():
        assert not math.isnan(row["rand"])
        assert row["rand_best_allele"] in {"HLA-A*02:01", "HLA-B*07:02"}
