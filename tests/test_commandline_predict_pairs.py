import os

import pytest

from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.binding_prediction import BindingPrediction
from mhctools.binding_prediction_collection import BindingPredictionCollection
from mhctools.pred import COLUMNS, Kind
from mhctools.unsupported_allele import UnsupportedAllele


class _PairStubPredictor(BaseCommandlinePredictor):
    def __init__(self, supported_alleles=None):
        if supported_alleles is None:
            supported_alleles = ["HLA-A*02:01", "HLA-B*07:02"]
        self.program_name = "fakepan"
        self.peptide_mode_flags = ["-p"]
        self.allele_flag = "-a"
        self.length_flag = "-l"
        self.input_file_flag = "-f"
        self.supported_alleles_flag = "-listMHC"
        self.tempdir_flag = None
        self.extra_flags = []
        self.alleles = ["HLA-A*02:01"]
        self._supported_allele_names = set(supported_alleles)
        self._allele_cli_names = {"HLA-A*02:01": "HLA-A*02:01"}
        self.max_alleles_per_command = 1
        self.max_peptides_per_file = 10 ** 4
        self.process_limit = -1
        self.group_peptides_by_length = False
        self.default_peptide_lengths = [9]
        self.allow_X_in_peptides = False
        self.allow_lowercase_in_peptides = False
        self.min_peptide_length = 8
        self.max_peptide_length = None
        self.parse_output_fn = None
        self.parse_to_preds_fn = None
        self.calls = []

    def prepare_allele_name(self, allele_name):
        return allele_name

    def _run_commands_and_collect_predictions(
            self, commands, input_filenames, temp_dir_list,
            sequence_key_mapping=None):
        predictions = []
        for output_file, command in commands.items():
            allele_arg = command[command.index("-a") + 1]
            alleles = allele_arg.split(",")
            input_filename = command[-1]
            with open(input_filename) as input_file:
                peptides = [
                    line.strip()
                    for line in input_file
                    if line.strip()]
            self.calls.append((alleles, peptides))
            for allele in alleles:
                for peptide in peptides:
                    predictions.append(BindingPrediction(
                        peptide=peptide,
                        allele=allele,
                        affinity=float(len(peptide) + len(allele)),
                        prediction_method_name=self.program_name))
            output_file.close()
            _silent_remove(output_file.name)
        for input_filename in input_filenames:
            _silent_remove(input_filename)
        return BindingPredictionCollection(predictions)


def _silent_remove(path):
    try:
        os.remove(path)
    except OSError:
        pass


def test_predict_pairs_groups_by_allele_and_preserves_input_order():
    predictor = _PairStubPredictor()

    results = predictor.predict_pairs(
        ["SIINFEKLL", "GILGFVFTL", "NLVPMVATV", "SIINFEKLL"],
        ["HLA-B*07:02", "HLA-A*02:01", "HLA-B*07:02", "HLA-A*02:01"])

    assert predictor.calls == [
        (["HLA-B*07:02"], ["SIINFEKLL", "NLVPMVATV"]),
        (["HLA-A*02:01"], ["GILGFVFTL", "SIINFEKLL"]),
    ]
    assert len(results) == 4
    assert [result.preds[0].peptide for result in results] == [
        "SIINFEKLL", "GILGFVFTL", "NLVPMVATV", "SIINFEKLL"]
    assert [result.preds[0].allele for result in results] == [
        "HLA-B*07:02", "HLA-A*02:01", "HLA-B*07:02", "HLA-A*02:01"]
    assert all(result.preds[0].kind == Kind.pMHC_affinity
               for result in results)
    assert all(result.preds[0].predictor_name == "fakepan"
               for result in results)


def test_predict_pairs_accepts_pair_iterable():
    predictor = _PairStubPredictor()

    results = predictor.predict_pairs([
        ("SIINFEKLL", "HLA-B*07:02"),
        ("GILGFVFTL", "HLA-A*02:01"),
    ])

    assert [result.preds[0].peptide for result in results] == [
        "SIINFEKLL", "GILGFVFTL"]
    assert [result.preds[0].allele for result in results] == [
        "HLA-B*07:02", "HLA-A*02:01"]


def test_predict_pairs_rejects_length_mismatch():
    predictor = _PairStubPredictor()

    with pytest.raises(ValueError) as e:
        predictor.predict_pairs(["SIINFEKLL"], [
            "HLA-A*02:01", "HLA-B*07:02"])

    assert "peptides length 1 != alleles length 2" in str(e.value)
    assert predictor.calls == []


def test_predict_pairs_rejects_malformed_pair_iterable():
    predictor = _PairStubPredictor()

    with pytest.raises(ValueError, match="Expected .* pairs"):
        predictor.predict_pairs(["SIINFEKLL"])

    assert predictor.calls == []


def test_predict_pairs_dataframe_is_canonical_schema():
    predictor = _PairStubPredictor()

    df = predictor.predict_pairs_dataframe(
        ["SIINFEKLL", "GILGFVFTL"],
        ["HLA-B*07:02", "HLA-A*02:01"],
        sample_name="sample-1")

    assert list(df.columns) == list(COLUMNS)
    assert len(df) == 2
    assert df["sample_name"].tolist() == ["sample-1", "sample-1"]
    assert df["peptide"].tolist() == ["SIINFEKLL", "GILGFVFTL"]
    assert df["allele"].tolist() == ["HLA-B*07:02", "HLA-A*02:01"]


def test_predict_pairs_raises_unsupported_allele_before_running_commands():
    predictor = _PairStubPredictor(supported_alleles=["HLA-A*02:01"])

    with pytest.raises(UnsupportedAllele, match="HLA-B"):
        predictor.predict_pairs(["SIINFEKLL"], ["HLA-B*07:02"])

    assert predictor.calls == []


def test_predict_pairs_dataframe_empty_result_has_canonical_columns():
    predictor = _PairStubPredictor()

    df = predictor.predict_pairs_dataframe([], [])

    assert list(df.columns) == list(COLUMNS)
    assert df.empty
