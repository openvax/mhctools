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

import pytest

from mhctools.pred import (
    COLUMNS,
    FIELD_BEST_DIRECTIONS,
    Kind,
    MHC_CLASS_VALUES,
    MHC_DEPENDENCE_VALUES,
    PeptideResult,
    Prediction,
    VALUE_BEST_DIRECTIONS,
    best_direction,
    preds_from_rows,
)
from mhctools.sample import MultiSample
from mhctools.base_predictor import BasePredictor
from mhctools.binding_prediction import BindingPrediction
from mhctools.binding_prediction_collection import BindingPredictionCollection
from mhctools.netmhc_pan41 import NetMHCpan41
from mhctools.netmhcstabpan import NetMHCstabpan
from mhctools.netmhcii_pan import NetMHCIIpan4, NetMHCIIpan43
from mhctools.random_predictor import RandomBindingPredictor


# -- Prediction --

def test_pred_basic():
    p = Prediction(
        kind=Kind.pMHC_affinity,
        score=0.85,
        peptide="SIINFEKL",
        allele="HLA-A*02:01",
        value=120.5,
        percentile_rank=0.8,
    )
    assert p.score == 0.85
    assert p.value == 120.5
    assert p.peptide == "SIINFEKL"
    assert p.allele == "HLA-A*02:01"
    assert p.kind == Kind.pMHC_affinity


def test_pred_frozen():
    p = Prediction(kind=Kind.pMHC_affinity, score=0.5)
    try:
        p.score = 0.9
        assert False, "should be frozen"
    except AttributeError:
        pass


def test_pred_to_row():
    p = Prediction(
        kind=Kind.pMHC_affinity,
        score=0.85,
        peptide="SIINFEKL",
        allele="HLA-A*02:01",
    )
    row = p.to_row(sample_name="pat001")
    assert row["sample_name"] == "pat001"
    assert row["peptide"] == "SIINFEKL"
    assert row["kind"] == "pMHC_affinity"
    assert set(row.keys()) == set(COLUMNS)


def test_pred_defaults():
    p = Prediction(kind=Kind.antigen_processing, score=0.7)
    assert p.peptide == ""
    assert p.allele == ""
    assert p.n_flank == ""
    assert p.c_flank == ""
    assert p.value is None
    assert p.percentile_rank is None


def test_pred_repr():
    p = Prediction(kind=Kind.pMHC_affinity, score=0.85, peptide="SIINFEKL",
             allele="HLA-A*02:01", value=120.5, percentile_rank=0.8,
             predictor_name="netMHCpan")
    r = repr(p)
    assert "SIINFEKL" in r
    assert "HLA-A*02:01" in r
    assert "pMHC_affinity" in r
    assert "score=0.85" in r
    assert "value=120.5" in r
    assert "rank=0.80%" in r
    assert "netMHCpan" in r
    assert str(p) == r


def test_pred_repr_minimal():
    p = Prediction(kind=Kind.proteasome_cleavage, score=0.5)
    r = repr(p)
    assert "score=0.5" in r
    assert "value" not in r
    assert "rank" not in r


def test_pred_to_dict_round_trip():
    p = Prediction(kind=Kind.pMHC_affinity, score=0.85, peptide="SIINFEKL",
             allele="HLA-A*02:01", value=120.5, percentile_rank=0.8)
    d = p.to_dict()
    assert d["kind"] == "pMHC_affinity"
    assert d["score"] == 0.85
    p2 = Prediction.from_dict(d)
    assert p == p2


def test_pred_to_dict_json_serializable():
    import json
    p = Prediction(kind=Kind.pMHC_affinity, score=0.85, peptide="SIINFEKL")
    s = json.dumps(p.to_dict())
    p2 = Prediction.from_dict(json.loads(s))
    assert p == p2


def test_pred_eq():
    p1 = Prediction(kind=Kind.pMHC_affinity, score=0.85, peptide="SIINFEKL")
    p2 = Prediction(kind=Kind.pMHC_affinity, score=0.85, peptide="SIINFEKL")
    p3 = Prediction(kind=Kind.pMHC_affinity, score=0.50, peptide="SIINFEKL")
    assert p1 == p2
    assert p1 != p3


# -- PeptideResult --

def _make_pred_set():
    return preds_from_rows(
        [
            dict(kind=Kind.pMHC_affinity, allele="HLA-A*02:01",
                 score=0.85, value=120.5, percentile_rank=0.8),
            dict(kind=Kind.pMHC_affinity, allele="HLA-B*07:02",
                 score=0.42, value=5000.0, percentile_rank=15.0),
            dict(kind=Kind.pMHC_presentation, allele="HLA-A*02:01",
                 score=0.92, percentile_rank=0.3),
            dict(kind=Kind.pMHC_presentation, allele="HLA-B*07:02",
                 score=0.15, percentile_rank=12.0),
            dict(kind=Kind.antigen_processing, score=0.85),
        ],
        peptide="SIINFEKL",
        predictor_name="mhcflurry",
        predictor_version="2.1",
    )


def test_preds_from_rows():
    ps = _make_pred_set()
    assert len(ps.preds) == 5
    for p in ps.preds:
        assert p.peptide == "SIINFEKL"
        assert p.predictor_name == "mhcflurry"


def test_peptide_result_repr():
    ps = _make_pred_set()
    r = repr(ps)
    assert "SIINFEKL" in r
    assert "5 preds" in r
    assert "pMHC_affinity" in r
    assert str(ps) == r


def test_peptide_result_repr_empty():
    assert "empty" in repr(PeptideResult())


def test_peptide_result_to_dict_round_trip():
    ps = _make_pred_set()
    d = ps.to_dict()
    assert len(d["preds"]) == 5
    ps2 = PeptideResult.from_dict(d)
    assert ps == ps2


def test_peptide_result_to_dict_json_serializable():
    import json
    ps = _make_pred_set()
    s = json.dumps(ps.to_dict())
    ps2 = PeptideResult.from_dict(json.loads(s))
    assert ps == ps2


def test_peptide_result_eq():
    ps1 = _make_pred_set()
    ps2 = _make_pred_set()
    assert ps1 == ps2


def test_best_affinity():
    ps = _make_pred_set()
    best = ps.best_affinity
    assert best is not None
    assert best.allele == "HLA-A*02:01"
    assert best.score == 0.85


def test_best_affinity_by_rank():
    ps = _make_pred_set()
    best = ps.best_affinity_by_rank
    assert best is not None
    assert best.allele == "HLA-A*02:01"
    assert best.percentile_rank == 0.8


def test_best_presentation():
    ps = _make_pred_set()
    best = ps.best_presentation
    assert best is not None
    assert best.allele == "HLA-A*02:01"
    assert best.score == 0.92


def test_best_presentation_by_rank():
    ps = _make_pred_set()
    best = ps.best_presentation_by_rank
    assert best is not None
    assert best.allele == "HLA-A*02:01"
    assert best.percentile_rank == 0.3


def test_best_stability_empty():
    ps = _make_pred_set()
    assert ps.best_stability is None
    assert ps.best_stability_by_rank is None


# -- shared fields --

def test_peptide_result_peptide():
    ps = _make_pred_set()
    assert ps.peptide == "SIINFEKL"

def test_peptide_result_offset():
    ps = _make_pred_set()
    assert ps.offset == 0

def test_peptide_result_source_sequence_name():
    ps = _make_pred_set()
    assert ps.source_sequence_name is None

def test_peptide_result_kinds():
    ps = _make_pred_set()
    assert ps.kinds == {
        Kind.pMHC_affinity,
        Kind.pMHC_presentation,
        Kind.antigen_processing,
    }

def test_peptide_result_alleles():
    ps = _make_pred_set()
    assert ps.alleles == {"HLA-A*02:01", "HLA-B*07:02"}

def test_empty_peptide_result_shared_fields():
    ps = PeptideResult()
    assert ps.peptide == ""
    assert ps.offset == 0
    assert ps.source_sequence_name is None
    assert ps.kinds == set()
    assert ps.alleles == set()


# -- kind accessors --

def test_affinity_accessor():
    ps = _make_pred_set()
    assert ps.affinity is not None
    assert ps.affinity.value == 120.5
    assert ps.affinity.score == 0.85
    assert ps.affinity.percentile_rank == 0.8
    assert ps.affinity.allele == "HLA-A*02:01"

def test_presentation_accessor():
    ps = _make_pred_set()
    assert ps.presentation is not None
    assert ps.presentation.score == 0.92
    assert ps.presentation.allele == "HLA-A*02:01"

def test_missing_kind_is_none():
    ps = _make_pred_set()
    assert ps.stability is None
    assert ps.immunogenicity is None

def test_kind_accessor_on_empty_result():
    ps = PeptideResult()
    assert ps.affinity is None


def test_filter():
    ps = _make_pred_set()
    affinity_preds = ps.filter(kind=Kind.pMHC_affinity)
    assert len(affinity_preds) == 2

    a1_preds = ps.filter(allele="HLA-A*02:01")
    assert len(a1_preds) == 2

    processing = ps.filter(kind=Kind.antigen_processing)
    assert len(processing) == 1
    assert processing[0].allele == ""


def test_to_dataframe():
    ps = _make_pred_set()
    df = ps.to_dataframe(sample_name="pat001")
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 5
    assert (df["sample_name"] == "pat001").all()
    assert (df["peptide"] == "SIINFEKL").all()


def test_empty_pred_set():
    ps = PeptideResult()
    assert ps.best_affinity is None
    df = ps.to_dataframe()
    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 0


# -- predict() on a predictor --

def test_predict_returns_peptide_preds():
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        default_peptide_lengths=[9])
    results = predictor.predict(["SIINFEKLL", "GILGFVFTL"])
    assert isinstance(results, list)
    assert all(isinstance(pp, PeptideResult) for pp in results)
    assert len(results) == 2
    for pp in results:
        assert pp.best_affinity is not None


def test_predict_dataframe():
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01"],
        default_peptide_lengths=[9])
    df = predictor.predict_dataframe(["SIINFEKLL"], sample_name="pat001")
    assert list(df.columns) == list(COLUMNS)
    assert (df["sample_name"] == "pat001").all()


def test_predict_proteins():
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01"],
        default_peptide_lengths=[9])
    result = predictor.predict_proteins({"TP53": "SIINFEKLLAA"})
    assert "TP53" in result
    assert isinstance(result["TP53"], list)
    assert all(isinstance(pp, PeptideResult) for pp in result["TP53"])
    for pp in result["TP53"]:
        for pred in pp.preds:
            assert pred.source_sequence_name == "TP53"


def test_binding_predictor_kind_support_defaults_to_single_allele_class_i():
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01"],
        default_peptide_lengths=[9])
    support = predictor.kind_support()
    assert predictor.supported_kinds == (Kind.pMHC_affinity,)
    assert support[Kind.pMHC_affinity]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_affinity]["mhc_class"] == "I"
    assert support[Kind.pMHC_affinity]["mhc_dependence"] in (
        MHC_DEPENDENCE_VALUES)
    assert support[Kind.pMHC_affinity]["mhc_class"] in MHC_CLASS_VALUES


def test_netmhcpan41_kind_support_includes_affinity_and_presentation():
    predictor = NetMHCpan41.__new__(NetMHCpan41)
    predictor.mode = "binding_affinity"

    support = predictor.kind_support()

    assert set(support) == {Kind.pMHC_affinity, Kind.pMHC_presentation}
    assert support[Kind.pMHC_affinity]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_presentation]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_affinity]["mhc_class"] == "I"


def test_netmhciipan4_el_kind_support_is_class_ii_presentation():
    predictor = NetMHCIIpan4.__new__(NetMHCIIpan4)
    predictor.mode = "elution_score"

    support = predictor.kind_support()

    assert set(support) == {Kind.pMHC_presentation}
    assert support[Kind.pMHC_presentation]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_presentation]["mhc_class"] == "II"


def test_netmhciipan43_ba_kind_support_is_class_ii_affinity():
    predictor = NetMHCIIpan43.__new__(NetMHCIIpan43)
    predictor.mode = "binding_affinity"

    support = predictor.kind_support()

    assert set(support) == {Kind.pMHC_affinity}
    assert support[Kind.pMHC_affinity]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_affinity]["mhc_class"] == "II"


def test_netmhcstabpan_kind_support_is_stability():
    predictor = NetMHCstabpan.__new__(NetMHCstabpan)

    support = predictor.kind_support()

    assert set(support) == {Kind.pMHC_stability}
    assert support[Kind.pMHC_stability]["mhc_dependence"] == "single_allele"
    assert support[Kind.pMHC_stability]["mhc_class"] == "I"


class FlankEchoPredictor(BasePredictor):
    uses_flanking_sequences = True
    flank_length = 2

    def __init__(self):
        BasePredictor.__init__(
            self,
            alleles=["HLA-A*02:01"],
            default_peptide_lengths=[3],
            min_peptide_length=3)
        self.calls = []

    def predict(self, peptides, n_flanks=None, c_flanks=None):
        peptide_list, n_flank_list, c_flank_list = self._check_flank_inputs(
            peptides, n_flanks, c_flanks)
        self.calls.append((peptide_list, n_flank_list, c_flank_list))
        results = []
        for i, peptide in enumerate(peptide_list):
            n_flank = n_flank_list[i] if n_flank_list is not None else ""
            c_flank = c_flank_list[i] if c_flank_list is not None else ""
            results.append(PeptideResult(preds=(Prediction(
                kind=Kind.pMHC_presentation,
                score=float(i),
                peptide=peptide,
                allele="HLA-A*02:01",
                n_flank=n_flank,
                c_flank=c_flank,
                predictor_name="flank_echo",
            ),)))
        return results

    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        return self.predict(peptides, n_flanks=n_flanks, c_flanks=c_flanks)


def test_predict_proteins_preserves_distinct_flanking_contexts():
    predictor = FlankEchoPredictor()
    result = predictor.predict_proteins({"protein": "XABCYABCZ"})

    peptides, n_flanks, c_flanks = predictor.calls[0]
    assert peptides.count("ABC") == 2

    abc_results = [pp for pp in result["protein"] if pp.peptide == "ABC"]
    assert len(abc_results) == 2
    assert [(pp.offset, pp.presentation.n_flank, pp.presentation.c_flank)
            for pp in abc_results] == [
                (1, "X", "YA"),
                (5, "CY", "Z"),
            ]
    assert n_flanks[1] == "X"
    assert c_flanks[1] == "YA"
    assert n_flanks[5] == "CY"
    assert c_flanks[5] == "Z"


def test_predict_proteins_dataframe():
    predictor = RandomBindingPredictor(
        alleles=["HLA-A*02:01"],
        default_peptide_lengths=[9])
    df = predictor.predict_proteins_dataframe(
        {"TP53": "SIINFEKLLAA"}, sample_name="pat001")
    assert list(df.columns) == list(COLUMNS)
    assert (df["source_sequence_name"] == "TP53").all()


# -- MultiSample --

def test_multi_sample_predict():
    ms = MultiSample(
        samples={
            "pat001": ["HLA-A*02:01"],
            "pat002": ["HLA-B*07:02"],
        },
        predictor_class=RandomBindingPredictor,
        default_peptide_lengths=[9],
    )
    results = ms.predict(["SIINFEKLL"])
    assert "pat001" in results
    assert "pat002" in results
    assert all(isinstance(pp, PeptideResult) for pp in results["pat001"])


def test_multi_sample_predict_dataframe():
    ms = MultiSample(
        samples={
            "pat001": ["HLA-A*02:01"],
            "pat002": ["HLA-B*07:02"],
        },
        predictor_class=RandomBindingPredictor,
        default_peptide_lengths=[9],
    )
    df = ms.predict_dataframe(["SIINFEKLL"])
    assert set(df["sample_name"]) == {"pat001", "pat002"}


# -- Compat layer --

def test_binding_prediction_to_pred():
    bp = BindingPrediction(
        peptide="SIINFEKL",
        allele="HLA-A*02:01",
        affinity=200.0,
        percentile_rank=0.3,
        source_sequence_name="seq",
        offset=5,
        prediction_method_name="netMHCpan",
    )
    pred = bp.to_pred()
    assert pred.kind == Kind.pMHC_affinity
    assert pred.peptide == "SIINFEKL"
    assert pred.allele == "HLA-A*02:01"
    assert pred.value == 200.0
    assert pred.percentile_rank == 0.3
    assert pred.source_sequence_name == "seq"
    assert pred.offset == 5
    assert pred.predictor_name == "netMHCpan"


def test_binding_prediction_from_pred():
    pred = Prediction(
        kind=Kind.pMHC_affinity,
        score=0.85,
        peptide="SIINFEKL",
        allele="HLA-A*02:01",
        value=120.5,
        percentile_rank=0.8,
        predictor_name="mhcflurry",
    )
    bp = BindingPrediction.from_pred(pred)
    assert bp.peptide == "SIINFEKL"
    assert bp.allele == "HLA-A*02:01"
    assert bp.affinity == 120.5
    assert bp.percentile_rank == 0.8
    assert bp.score == 0.85
    assert bp.prediction_method_name == "mhcflurry"


def test_collection_to_preds():
    bps = BindingPredictionCollection([
        BindingPrediction(peptide="SIINFEKL", allele="HLA-A*02:01", affinity=200.0),
        BindingPrediction(peptide="SIINFEKL", allele="HLA-B*07:02", affinity=5000.0),
    ])
    preds = bps.to_preds()
    assert len(preds) == 2
    assert all(isinstance(p, Prediction) for p in preds)


def test_collection_to_peptide_preds():
    bps = BindingPredictionCollection([
        BindingPrediction(peptide="SIINFEKL", allele="HLA-A*02:01", affinity=200.0, offset=0),
        BindingPrediction(peptide="SIINFEKL", allele="HLA-B*07:02", affinity=5000.0, offset=0),
        BindingPrediction(peptide="GILGFVFTL", allele="HLA-A*02:01", affinity=50.0, offset=10),
    ])
    pp_list = bps.to_peptide_preds()
    assert len(pp_list) == 2  # two distinct peptide positions
    sizes = sorted(len(pp.preds) for pp in pp_list)
    assert sizes == [1, 2]


# -- best_direction --


def test_field_best_directions_constants():
    """score is max-better, percentile_rank is min-better — uniform across kinds."""
    assert FIELD_BEST_DIRECTIONS["score"] == "max"
    assert FIELD_BEST_DIRECTIONS["percentile_rank"] == "min"


def test_value_best_directions_kind_specific():
    """value direction is kind-dependent; affinity uses IC50 (min-better)."""
    assert VALUE_BEST_DIRECTIONS[Kind.pMHC_affinity] == "min"
    assert VALUE_BEST_DIRECTIONS[Kind.pMHC_stability] == "max"


def test_best_direction_field_defaults():
    # score is max for any kind that uses it
    assert best_direction(Kind.pMHC_affinity, "score") == "max"
    assert best_direction(Kind.pMHC_presentation, "score") == "max"
    assert best_direction(Kind.proteasome_cleavage, "score") == "max"
    # percentile_rank is min for any kind that uses it
    assert best_direction(Kind.pMHC_affinity, "percentile_rank") == "min"
    assert best_direction(Kind.pMHC_presentation, "percentile_rank") == "min"


def test_best_direction_value_kind_specific():
    assert best_direction(Kind.pMHC_affinity, "value") == "min"   # IC50
    assert best_direction(Kind.pMHC_stability, "value") == "max"  # half-life


def test_best_direction_value_unregistered_raises():
    """`value` for a kind without a registered direction must raise —
    we refuse to guess the unit semantics."""
    with pytest.raises(ValueError, match="best_direction undefined"):
        best_direction(Kind.pMHC_presentation, "value")
    with pytest.raises(ValueError, match="best_direction undefined"):
        best_direction(Kind.proteasome_cleavage, "value")


def test_best_direction_unknown_field_raises():
    with pytest.raises(ValueError, match="best_direction undefined for field"):
        best_direction(Kind.pMHC_affinity, "made_up_field")


def test_best_direction_accepts_string_kind():
    # Kind constants are plain strings, so callers can pass either form.
    assert best_direction("pMHC_affinity", "value") == "min"
    assert best_direction("pMHC_affinity", "score") == "max"


# -- best_by (public, direction-aware) --


def test_best_by_score_picks_max():
    ps = _make_pred_set()
    best = ps.best_by_score(Kind.pMHC_affinity)
    assert best.allele == "HLA-A*02:01"
    assert best.score == 0.85


def test_best_by_rank_picks_min():
    ps = _make_pred_set()
    best = ps.best_by_rank(Kind.pMHC_affinity)
    assert best.allele == "HLA-A*02:01"
    assert best.percentile_rank == 0.8


def test_best_by_value_affinity_picks_min_ic50():
    ps = _make_pred_set()
    best = ps.best_by_value(Kind.pMHC_affinity)
    assert best.allele == "HLA-A*02:01"
    assert best.value == 120.5


def test_best_by_value_stability_picks_max_half_life():
    ps = preds_from_rows(
        [
            dict(kind=Kind.pMHC_stability, allele="HLA-A*02:01",
                 score=0.6, value=4.0),
            dict(kind=Kind.pMHC_stability, allele="HLA-B*07:02",
                 score=0.9, value=12.0),
        ],
        peptide="SIINFEKL",
    )
    best = ps.best_by_value(Kind.pMHC_stability)
    assert best.allele == "HLA-B*07:02"
    assert best.value == 12.0


def test_best_by_value_unregistered_kind_raises():
    ps = _make_pred_set()
    with pytest.raises(ValueError, match="best_direction undefined"):
        ps.best_by_value(Kind.pMHC_presentation)


def test_best_by_skips_none_field():
    # antigen_processing pred has score but None for percentile_rank
    ps = _make_pred_set()
    assert ps.best_by_rank(Kind.antigen_processing) is None
    assert ps.best_by_score(Kind.antigen_processing) is not None


def test_best_by_falls_back_to_allele_less():
    ps = preds_from_rows(
        [dict(kind=Kind.proteasome_cleavage, score=0.7),
         dict(kind=Kind.proteasome_cleavage, score=0.9)],
        peptide="SIINFEKL",
    )
    best = ps.best_by_score(Kind.proteasome_cleavage)
    assert best is not None
    assert best.score == 0.9
    assert best.allele == ""


def test_best_by_missing_kind_returns_none():
    ps = _make_pred_set()
    assert ps.best_by_score(Kind.immunogenicity) is None
    assert ps.best_by_rank(Kind.pMHC_stability) is None


def test_best_by_accepts_string_kind():
    # Kind constants are plain strings, so callers can pass either form.
    ps = _make_pred_set()
    by_kind = ps.best_by_score(Kind.pMHC_affinity)
    by_str = ps.best_by_score("pMHC_affinity")
    assert by_kind == by_str
    assert ps.best_by("pMHC_affinity", "value") == ps.best_by_value(Kind.pMHC_affinity)


def test_best_by_rank_falls_back_to_allele_less():
    # Behavior change vs the old _best_by_rank: rank predictions without an
    # allele are now considered when no allele-bearing rank pred exists. Pin
    # the new contract so it can't silently regress.
    ps = preds_from_rows(
        [
            dict(kind=Kind.pMHC_presentation, score=0.5, percentile_rank=10.0),
            dict(kind=Kind.pMHC_presentation, score=0.8, percentile_rank=2.0),
        ],
        peptide="SIINFEKL",
    )
    best = ps.best_presentation_by_rank
    assert best is not None
    assert best.allele == ""
    assert best.percentile_rank == 2.0
