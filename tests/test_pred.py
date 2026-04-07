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

from mhctools.pred import Pred, PeptidePreds, Kind, preds_from_rows, COLUMNS
from mhctools.sample import MultiSample
from mhctools.binding_prediction import BindingPrediction
from mhctools.binding_prediction_collection import BindingPredictionCollection
from mhctools.random_predictor import RandomBindingPredictor


# -- Pred --

def test_pred_basic():
    p = Pred(
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
    p = Pred(kind=Kind.pMHC_affinity, score=0.5)
    try:
        p.score = 0.9
        assert False, "should be frozen"
    except AttributeError:
        pass


def test_pred_to_row():
    p = Pred(
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
    p = Pred(kind=Kind.antigen_processing, score=0.7)
    assert p.peptide == ""
    assert p.allele == ""
    assert p.n_flank == ""
    assert p.c_flank == ""
    assert p.value is None
    assert p.percentile_rank is None


# -- PeptidePreds --

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
    ps = PeptidePreds()
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
    assert all(isinstance(pp, PeptidePreds) for pp in results)
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
    assert all(isinstance(pp, PeptidePreds) for pp in result["TP53"])
    for pp in result["TP53"]:
        for pred in pp.preds:
            assert pred.source_sequence_name == "TP53"


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
    assert all(isinstance(pp, PeptidePreds) for pp in results["pat001"])


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
    pred = Pred(
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
    assert all(isinstance(p, Pred) for p in preds)


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
