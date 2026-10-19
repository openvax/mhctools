"""Runtime input validation must preserve one-shot iterables."""

import subprocess
import sys

import numpy as np
import pytest

from mhctools import RandomBindingPredictor
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.base_predictor import BasePredictor, _check_flank_inputs
from mhctools.mhcflurry import _resolve_supported_alleles


ALLELES = ["HLA-A*02:01", "HLA-B*07:02"]
PEPTIDES = ["SIINFEKLQ", "SIINFEKLA"]


@pytest.fixture
def commandline_predictor(monkeypatch):
    monkeypatch.setattr(
        BaseCommandlinePredictor, "_determine_supported_alleles",
        staticmethod(lambda *args: set()))

    def make(**updates):
        kwargs = dict(
            program_name="fixture", alleles=[], parse_output_fn=lambda _: [],
            supported_alleles_flag="-list", input_file_flag="-f",
            length_flag="-l", allele_flag="-a")
        kwargs.update(updates)
        return BaseCommandlinePredictor(**kwargs)

    return make


def test_constructor_preserves_generator_alleles_and_lengths():
    predictor = RandomBindingPredictor(
        alleles=iter(ALLELES), default_peptide_lengths=iter([8, 9]))
    assert set(predictor.alleles) == set(ALLELES)
    assert predictor.default_peptide_lengths == [8, 9]
    results = predictor.predict_subsequences({"seq": "SIINFEKLQ"})
    assert {len(result.peptide) for result in results} == {8, 9}
    assert len(results) == 6


def test_protein_scan_preserves_explicit_generator_lengths():
    predictor = RandomBindingPredictor(alleles=ALLELES)
    results = predictor.predict_subsequences(
        {"seq": "SIINFEKLQ"}, peptide_lengths=iter([8, 9]))
    assert len(results) == 6


def test_prediction_preserves_generator_peptides_and_flanks():
    predictor = RandomBindingPredictor(alleles=ALLELES)
    results = predictor.predict(
        iter(PEPTIDES), n_flanks=iter(["A", "B"]), c_flanks=iter(["C", "D"]))
    assert [result.preds[0].peptide for result in results] == PEPTIDES
    assert len(predictor.predict_peptides(iter(PEPTIDES))) == 4


@pytest.mark.parametrize("field", ["peptide_mode_flags", "extra_flags"])
def test_constructor_preserves_generator_flags(commandline_predictor, field):
    predictor = commandline_predictor(**{field: iter(["-p", "--BA"])})
    assert getattr(predictor, field) == ["-p", "--BA"]


@pytest.mark.parametrize("field, values", [
    ("alleles", [ALLELES[0], 7]),
    ("default_peptide_lengths", [8, "9"]),
    ("peptide_mode_flags", ["-p", 7]),
    ("extra_flags", ["--BA", 7]),
])
def test_invalid_generator_element_is_reported(commandline_predictor, field, values):
    with pytest.raises(TypeError, match=r"index 1"):
        commandline_predictor(**{field: iter(values)})


@pytest.mark.parametrize("field", [
    "program_name", "supported_alleles_flag", "input_file_flag",
    "length_flag", "allele_flag", "tempdir_flag",
])
def test_non_string_flags_fail_before_running_backend(commandline_predictor, field):
    with pytest.raises(TypeError):
        commandline_predictor(**{field: 7})


@pytest.mark.parametrize("field", [
    "max_peptides_per_file", "max_alleles_per_command", "process_limit",
])
@pytest.mark.parametrize("value", ["3", 3.0, np.int64(3)])
def test_integer_fields_keep_python_int_validation(commandline_predictor, field, value):
    with pytest.raises(TypeError):
        commandline_predictor(**{field: value})


def test_bool_keeps_existing_integer_acceptance(commandline_predictor):
    predictor = commandline_predictor(process_limit=True)
    assert predictor.process_limit is True


def test_invalid_peptide_and_flank_elements_have_type_errors():
    with pytest.raises(TypeError, match=r"index 1"):
        BasePredictor()._check_peptide_inputs([PEPTIDES[0], 7])
    with pytest.raises(TypeError, match=r"n_flanks.*index 1"):
        _check_flank_inputs(PEPTIDES, n_flanks=iter(["A", 7]))


def test_mhcflurry_allele_validation_preserves_generators():
    assert _resolve_supported_alleles(iter(ALLELES), ALLELES) == ALLELES
    with pytest.raises(TypeError, match=r"index 1"):
        _resolve_supported_alleles(iter([ALLELES[0], 7]), ALLELES)


def test_import_and_prediction_without_typechecks():
    result = subprocess.run(
        [sys.executable, "-c", """
import sys
sys.modules["typechecks"] = None
from mhctools import RandomBindingPredictor
predictor = RandomBindingPredictor(alleles=["HLA-A*02:01"])
assert len(predictor.predict(iter(["SIINFEKLQ"]))) == 1
"""], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
