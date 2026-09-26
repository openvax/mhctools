"""Prediction results stay aligned even when backends group or reorder rows."""

import json

import pytest

from mhctools import BindingPrediction, BindingPredictionCollection, RandomBindingPredictor
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.base_predictor import BasePredictor
from mhctools.pred import COLUMNS, Kind, Prediction


ALLELES = ["HLA-A*02:01", "HLA-B*07:02"]
PEPTIDES = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "GILGFVFTL", "SIINFEKL"]
VALUES = {"GILGFVFTL": 10.0, "SIINFEKL": 100.0, "NLVPMVATV": 1000.0}


def _binding_rows(peptides, alleles):
    # Deliberately return allele-major, length-sorted output.
    return [
        BindingPrediction(
            peptide=peptide, allele=allele, affinity=VALUES[peptide],
            prediction_method_name="fixture")
        for allele in alleles
        for peptide in sorted(peptides, key=len)
    ]


class _LegacyPredictor(BasePredictor):
    def predict_peptides(self, peptides):
        return BindingPredictionCollection(_binding_rows(peptides, self.alleles))


@pytest.fixture(params=["base", "commandline-legacy", "commandline-native"])
def predictor(request, monkeypatch):
    if request.param == "base":
        return _LegacyPredictor(alleles=ALLELES, default_peptide_lengths=[9])

    monkeypatch.setattr(
        BaseCommandlinePredictor, "_determine_supported_alleles",
        staticmethod(lambda *args: set(ALLELES)))

    def run_commands(commands, **kwargs):
        for output, command in commands.items():
            with open(command[command.index("-f") + 1]) as source:
                peptides = source.read().splitlines()
            alleles = command[command.index("-a") + 1].split(",")
            rows = []
            for bp in _binding_rows(peptides, alleles):
                rows.append(bp.to_pred().to_dict())
                if request.param == "commandline-native":
                    rows.append(Prediction(
                        peptide=bp.peptide, allele=bp.allele,
                        kind=Kind.pMHC_presentation, score=0.75,
                        predictor_name="fixture", predictor_version="1.2").to_dict())
            output.write(json.dumps(rows))
            output.flush()

    monkeypatch.setattr(
        "mhctools.base_commandline_predictor.run_multiple_commands_redirect_stdout",
        run_commands)

    def parse_legacy(stdout, **kwargs):
        return BindingPredictionCollection([
            BindingPrediction.from_pred(Prediction.from_dict(row))
            for row in json.loads(stdout)])

    def parse_native(stdout, **kwargs):
        return [Prediction.from_dict(row) for row in json.loads(stdout)]

    result = BaseCommandlinePredictor(
        program_name="fixture", alleles=ALLELES,
        parse_output_fn=parse_legacy,
        parse_to_preds_fn=(
            parse_native if request.param == "commandline-native" else None),
        supported_alleles_flag="-list", input_file_flag="-f",
        length_flag="-l", allele_flag="-a",
        group_peptides_by_length=True, max_peptides_per_file=1,
        max_alleles_per_command=1)
    # Keep the fake command's spellings canonical, like its output.
    result._allele_cli_names = {allele: allele for allele in ALLELES}
    return result


def test_predict_preserves_order_duplicates_and_all_scores(predictor):
    results = predictor.predict(iter(PEPTIDES))

    assert [result.peptide for result in results] == PEPTIDES
    for peptide, result in zip(PEPTIDES, results):
        assert result.alleles == set(ALLELES)
        assert len(result.preds) == len(ALLELES) * len(result.kinds)
        assert all(pred.peptide == peptide for pred in result.preds)
        assert all(pred.value == VALUES[peptide]
                   for pred in result.preds if pred.kind == Kind.pMHC_affinity)
        assert all(pred.score == 0.75 and pred.predictor_version == "1.2"
                   for pred in result.preds if pred.kind == Kind.pMHC_presentation)
    assert results[0] == results[3]
    # PeptideResult is mutable: each occurrence must have its own container.
    results[0].preds = ()
    assert results[3].preds


def test_predict_dataframe_preserves_occurrence_blocks(predictor):
    frame = predictor.predict_dataframe(PEPTIDES, sample_name="sample")
    rows_per_input = len(ALLELES) * frame.kind.nunique()
    assert frame.peptide.tolist() == [
        peptide for peptide in PEPTIDES for _ in range(rows_per_input)]
    assert frame.sample_name.tolist() == ["sample"] * len(frame)
    assert list(frame.columns) == list(COLUMNS)


def test_predict_empty_input(predictor):
    assert predictor.predict([]) == []
    frame = predictor.predict_dataframe([])
    assert frame.empty
    assert list(frame.columns) == list(COLUMNS)


@pytest.mark.parametrize("missing", ["peptide", "allele", "all"])
def test_predict_rejects_missing_predictions(predictor, monkeypatch, missing):
    original = _binding_rows

    def incomplete_rows(peptides, alleles):
        return [bp for bp in original(peptides, alleles)
                if missing != "all"
                and not (missing == "peptide" and bp.peptide == PEPTIDES[1])
                and not (missing == "allele" and bp.allele == ALLELES[1])]

    monkeypatch.setattr(__name__ + "._binding_rows", incomplete_rows)
    with pytest.raises(ValueError, match="Missing predictions|No parseable predictions"):
        predictor.predict(PEPTIDES)


def test_random_predictor_preserves_occurrences_without_duplicating_allele_rows():
    results = RandomBindingPredictor(alleles=ALLELES).predict(PEPTIDES)
    assert [result.peptide for result in results] == PEPTIDES
    assert all(len(result.preds) == len(ALLELES) for result in results)
