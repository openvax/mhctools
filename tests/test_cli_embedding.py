"""Embedding the shared options in cached-report and fresh-prediction workflows."""

from argparse import ArgumentParser

import pytest

from mhctools.cli import (
    add_mhc_args,
    make_mhc_arg_parser,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
    predictors_from_args,
)
from mhctools.cli.args import _LazyPredictor
from mhctools.cli.script import main


def _optional_parser():
    parser = ArgumentParser()
    add_mhc_args(parser, predictor_required=False)
    return parser


@pytest.mark.parametrize("tokens, expected", [
    (["HLA-A*02:01", "HLA-B*07:02"], ["HLA-A*02:01", "HLA-B*07:02"]),
    (["HLA-A*02:01,HLA-B*07:02"], ["HLA-A*02:01", "HLA-B*07:02"]),
    (["HLA-A*02:01 HLA-B*07:02"], ["HLA-A*02:01", "HLA-B*07:02"]),
    (["HLA-A*02:01,HLA-B*07:02", "HLA-C*07:02"],
     ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]),
    (["H2-Kb", "H2-Db"], ["H-2-Kb", "H-2-Db"]),
])
def test_allele_arguments_preserve_the_string_contract(tokens, expected):
    args = _optional_parser().parse_args(["--mhc-alleles", *tokens])
    assert isinstance(args.mhc_alleles, str)
    assert mhc_alleles_from_args(args) == expected


def test_multiple_allele_arguments_combine_with_file_and_predict(tmp_path):
    allele_file = tmp_path / "alleles.txt"
    allele_file.write_text("HLA-C*07:02\n")
    args = _optional_parser().parse_args([
        "--mhc-alleles", "HLA-A*02:01", "HLA-B*07:02",
        "--mhc-alleles-file", str(allele_file), "--mhc-predictor", "random",
        "--mhc-peptide-lengths", "8",
    ])
    predictor, = predictors_from_args(args)
    predictions = predictor.predict(["SIINFEKL"])
    assert {entry.allele for peptide in predictions for entry in peptide.preds} == {
        "HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"}


def test_repeated_allele_option_keeps_the_last_occurrence():
    args = _optional_parser().parse_args([
        "--mhc-alleles", "HLA-A*02:01", "HLA-B*07:02",
        "--mhc-alleles", "HLA-C*07:02",
    ])
    assert mhc_alleles_from_args(args) == ["HLA-C*07:02"]


def test_allele_option_requires_a_value():
    with pytest.raises(SystemExit) as exc:
        _optional_parser().parse_args(["--mhc-alleles", "--mhc-predictor", "random"])
    assert exc.value.code == 2


def test_cached_report_parses_options_without_loading_models(monkeypatch, tmp_path):
    def unexpected_resolution(self):
        pytest.fail("Parsing cached-report options resolved a predictor")

    monkeypatch.setattr(_LazyPredictor, "_resolve", unexpected_resolution)
    parser = _optional_parser()
    parser.add_argument("--replay", action="store_true")
    missing_alleles = tmp_path / "not-loaded.txt"
    args = parser.parse_args([
        "--replay", "--mhc-alleles-file", str(missing_alleles),
        "--mhc-alleles", "HLA-A*02:01", "--mhc-peptide-lengths", "8-9",
        "--mhc-predictor-models-path", "/unused/models",
    ])

    assert args.replay is True
    assert args.mhc_predictor is None
    assert args.mhc_alleles_file == str(missing_alleles)
    assert args.mhc_alleles == "HLA-A*02:01"
    assert args.mhc_peptide_lengths == [8, 9]
    assert args.mhc_predictor_models_path == "/unused/models"


def test_optional_parser_can_run_an_explicit_predictor():
    args = _optional_parser().parse_args([
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--mhc-peptide-lengths", "8",
    ])

    predictor, = predictors_from_args(args)
    result, = predictor.predict(["SIINFEKL"])

    assert result.peptide == "SIINFEKL"
    assert result.affinity.allele == "HLA-A*02:01"
    assert predictor.default_peptide_lengths == [8]


@pytest.mark.parametrize("factory", [predictors_from_args, mhc_binding_predictor_from_args])
@pytest.mark.parametrize("tokens", [[], ["--mhc-predictor", ","]])
def test_fresh_prediction_requires_a_nonempty_selection(factory, tokens):
    args = _optional_parser().parse_args(tokens)
    with pytest.raises(ValueError, match="MHC predictor required.*--mhc-predictor"):
        factory(args)


def test_optional_parser_still_validates_predictor_names():
    with pytest.raises(SystemExit) as exc:
        _optional_parser().parse_args(["--mhc-predictor", "not-a-predictor"])
    assert exc.value.code == 2


def test_optional_parser_still_requires_alleles_for_fresh_mhc_prediction():
    args = _optional_parser().parse_args(["--mhc-predictor", "random"])
    with pytest.raises(ValueError, match="MHC alleles required"):
        predictors_from_args(args)


def test_default_parser_requires_a_predictor(capsys):
    with pytest.raises(SystemExit) as exc:
        make_mhc_arg_parser().parse_args([])
    assert exc.value.code == 2
    assert "required: --mhc-predictor" in capsys.readouterr().err


def test_standalone_prediction_cli_requires_a_predictor(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--sequence", "SIINFEKL", "--mhc-alleles", "HLA-A*02:01"])
    assert exc.value.code == 2
    assert "required: --mhc-predictor" in capsys.readouterr().err
