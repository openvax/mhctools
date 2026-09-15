"""Tests that wrapped executables preserve diagnostics on failure."""

import sys
import tempfile

import pytest

from mhctools.base_commandline_predictor import BaseCommandlinePredictor


def _output_file():
    return tempfile.NamedTemporaryFile("w+", delete=False)


def _predictor(parse_output_fn):
    predictor = BaseCommandlinePredictor.__new__(BaseCommandlinePredictor)
    predictor.program_name = "fixture-predictor"
    predictor.process_limit = 1
    predictor.parse_output_fn = parse_output_fn
    return predictor


def test_nonzero_command_includes_combined_output_tail():
    predictor = _predictor(lambda **kwargs: [])
    command = [
        sys.executable,
        "-c",
        "import sys; print('stdout detail'); "
        "print('stderr detail', file=sys.stderr); raise SystemExit(3)",
    ]
    with pytest.raises(RuntimeError) as raised:
        predictor._run_commands_and_collect_predictions(
            {_output_file(): command}, [], [])
    message = str(raised.value)
    assert "stdout detail" in message
    assert "stderr detail" in message


def test_unparseable_success_includes_output_tail():
    predictor = _predictor(lambda **kwargs: [])
    command = [sys.executable, "-c", "print('model diagnostic')"]
    with pytest.raises(ValueError) as raised:
        predictor._run_commands_and_collect_predictions(
            {_output_file(): command}, [], [])
    assert "No parseable predictions" in str(raised.value)
    assert "model diagnostic" in str(raised.value)

