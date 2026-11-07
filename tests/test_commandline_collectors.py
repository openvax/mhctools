"""The public command collectors keep catalog-aware allele normalization."""

import sys
import tempfile

import pytest

from mhctools.base_commandline_predictor import BaseCommandlinePredictor


def stub_predictor():
    """Return a command-line predictor with a one-entry native catalog."""
    predictor = BaseCommandlinePredictor.__new__(BaseCommandlinePredictor)
    predictor.program_name = "fixture-predictor"
    predictor.process_limit = 1
    predictor._supported_allele_names = {"HLA-A30:14L"}
    return predictor


def echo_command():
    """Return one output file mapped to a command that prints a line."""
    output = tempfile.NamedTemporaryFile("w+", delete=False)
    return {output: [sys.executable, "-c", "print('HLA-A30:14L')"]}


@pytest.mark.parametrize("collector, parser_attribute, name_argument", [
    ("run_commands_and_collect_predictions", "parse_output_fn",
     "prediction_method_name"),
    ("run_commands_and_collect_preds", "parse_to_preds_fn", "predictor_name"),
])
def test_collectors_pass_catalog_normalizer(
        collector, parser_attribute, name_argument):
    predictor = stub_predictor()
    calls = []

    def parser(stdout, sequence_key_mapping, allele_normalizer, **names):
        calls.append((names, allele_normalizer))
        return [allele_normalizer(stdout.strip())]

    setattr(predictor, parser_attribute, parser)
    results = getattr(predictor, collector)(echo_command(), [], [])
    assert list(results) == ["HLA-A*30:14L"]
    [(names, normalizer)] = calls
    assert names == {name_argument: "fixture-predictor"}
    assert normalizer == predictor._normalize_native_output


def test_collectors_keep_parsers_without_normalizer_argument_working():
    predictor = stub_predictor()

    def parser(stdout, sequence_key_mapping, prediction_method_name):
        return [stdout.strip()]

    predictor.parse_output_fn = parser
    results = predictor.run_commands_and_collect_predictions(
        echo_command(), [], [])
    assert list(results) == ["HLA-A30:14L"]
