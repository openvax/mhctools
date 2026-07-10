import pytest

from mhctools.bigmhc import BigMHC
from mhctools.pred import COLUMNS, Kind, PeptideResult


class _FakeBigMHC(BigMHC):
    def __init__(self, alleles=("HLA-A*02:01",), mode="im"):
        self.alleles = list(alleles)
        self.mode = mode
        self.calls = []

    def _predict_raw(self, peptides, alleles):
        self.calls.append((list(peptides), list(alleles)))
        return [i / 10.0 for i in range(len(peptides))]


def test_predict_pairs_scores_parallel_peptide_allele_rows():
    predictor = _FakeBigMHC(mode="im")

    results = predictor.predict_pairs(
        ["siinfekl", "GILGFVFTL"],
        ["HLA-B*07:02", "HLA-A*02:01"])

    assert predictor.calls == [(
        ["SIINFEKL", "GILGFVFTL"],
        ["HLA-B*07:02", "HLA-A*02:01"],
    )]
    assert len(results) == 2
    assert all(isinstance(result, PeptideResult) for result in results)
    assert [result.preds[0].peptide for result in results] == [
        "SIINFEKL", "GILGFVFTL"]
    assert [result.preds[0].allele for result in results] == [
        "HLA-B*07:02", "HLA-A*02:01"]
    assert [result.preds[0].score for result in results] == [0.0, 0.1]
    assert results[0].preds[0].kind == Kind.immunogenicity
    assert results[0].preds[0].predictor_name == "bigmhc_im"


def test_predict_pairs_accepts_pair_iterable():
    predictor = _FakeBigMHC(mode="el")

    results = predictor.predict_pairs([
        ("SIINFEKL", "HLA-B*07:02"),
        ("GILGFVFTL", "HLA-A*02:01"),
    ])

    assert predictor.calls == [(
        ["SIINFEKL", "GILGFVFTL"],
        ["HLA-B*07:02", "HLA-A*02:01"],
    )]
    assert results[0].preds[0].kind == Kind.pMHC_presentation
    assert results[0].preds[0].predictor_name == "bigmhc_el"


def test_predict_pairs_rejects_length_mismatch():
    predictor = _FakeBigMHC()

    with pytest.raises(ValueError) as e:
        predictor.predict_pairs(["SIINFEKL"], ["HLA-A*02:01", "HLA-B*07:02"])

    assert "peptides length 1 != alleles length 2" in str(e.value)
    assert predictor.calls == []


def test_predict_uses_pairs_for_existing_cross_product_api():
    predictor = _FakeBigMHC(
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        mode="el")

    results = predictor.predict(["siinfekl", "GILGFVFTL"])

    assert predictor.calls == [(
        ["SIINFEKL", "SIINFEKL", "GILGFVFTL", "GILGFVFTL"],
        ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*02:01", "HLA-B*07:02"],
    )]
    assert len(results) == 2
    assert [pred.allele for pred in results[0].preds] == [
        "HLA-A*02:01", "HLA-B*07:02"]
    assert [pred.score for pred in results[0].preds] == [0.0, 0.1]
    assert [pred.allele for pred in results[1].preds] == [
        "HLA-A*02:01", "HLA-B*07:02"]
    assert [pred.score for pred in results[1].preds] == [0.2, 0.3]


def test_predict_pairs_dataframe_is_one_row_per_pair():
    predictor = _FakeBigMHC(mode="im")

    df = predictor.predict_pairs_dataframe(
        ["SIINFEKL", "GILGFVFTL"],
        ["HLA-A*02:01", "HLA-B*07:02"],
        sample_name="sample-1")

    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 2
    assert df["sample_name"].tolist() == ["sample-1", "sample-1"]
    assert df["peptide"].tolist() == ["SIINFEKL", "GILGFVFTL"]
    assert df["allele"].tolist() == ["HLA-A*02:01", "HLA-B*07:02"]
