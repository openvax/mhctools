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

import tempfile
from argparse import ArgumentParser
from os import remove

import pandas as pd
import pytest

from mhctools.cli.script import (
    add_output_args,
    format_predictions,
    main,
    parse_args,
    run_predictor,
)
from .common import eq_


def test_main_prints_prediction_table(capsys):
    """Without --output-csv the predictions must reach stdout (not a logger)."""
    main([
        "--mhc-predictor", "random",
        "--sequence", "SIINFEKL",
        "--mhc-alleles", "HLA-A*02:01"])
    stdout = capsys.readouterr().out
    assert "SIINFEKL" in stdout
    assert "peptide" in stdout
    assert "allele" in stdout


def test_main_with_output_csv_reports_the_file(capsys, tmp_path):
    output_csv = tmp_path / "predictions.csv"
    main([
        "--mhc-predictor", "random",
        "--sequence", "SIINFEKL",
        "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", str(output_csv)])
    assert "Wrote: %s" % output_csv in capsys.readouterr().out
    assert list(pd.read_csv(output_csv).peptide) == ["SIINFEKL"]


def test_format_predictions_has_no_index_column_and_short_floats():
    df = pd.DataFrame({
        "peptide": ["SIINFEKL"],
        "affinity": [11927.161441410432],
    })
    lines = format_predictions(df).splitlines()
    assert lines[0].split() == ["peptide", "affinity"]
    assert lines[1].split() == ["SIINFEKL", "11927.2"]


def test_format_predictions_when_empty():
    assert format_predictions(pd.DataFrame({"peptide": []})) == "No predictions."


def test_repeated_sequence_options_accumulate():
    args = parse_args([
        "--mhc-predictor", "random",
        "--sequence", "SIINFEKL",
        "--sequence", "GILGFVFTL", "AAAAAAAAA",
        "--mhc-alleles", "HLA-A*02:01"])
    assert args.sequence == ["SIINFEKL", "GILGFVFTL", "AAAAAAAAA"]
    assert [prediction.peptide for prediction in run_predictor(args)] == [
        "SIINFEKL", "GILGFVFTL", "AAAAAAAAA"]


@pytest.mark.parametrize("second_source", [
    ["--input-peptides-file", "peptides.txt"],
    ["--input-fasta-file", "proteins.fasta"],
])
def test_sequence_rejects_competing_input_source(second_source):
    with pytest.raises(SystemExit) as error:
        parse_args([
            "--mhc-predictor", "random",
            "--sequence", "SIINFEKL",
            *second_source,
            "--mhc-alleles", "HLA-A*02:01"])
    assert error.value.code == 2


def test_sequence_requires_at_least_one_value():
    with pytest.raises(SystemExit) as error:
        parse_args([
            "--mhc-predictor", "random",
            "--sequence",
            "--mhc-alleles", "HLA-A*02:01"])
    assert error.value.code == 2


def test_add_output_args_uses_supplied_parser():
    parser = ArgumentParser()
    add_output_args(parser)
    assert parser.parse_args(["--output-csv", "out.csv"]).output_csv == "out.csv"

def test_peptides_without_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)

def test_peptides_with_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])

def test_peptides_file_without_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)
    remove(f.name)


def test_peptides_file_ignores_blank_and_whitespace_only_lines(tmp_path):
    path = tmp_path / "peptides.txt"
    path.write_text("SIINFEKL\n\n   \nGILGFVFTL\n")
    args = parse_args([
        "--mhc-predictor", "random",
        "--input-peptides-file", str(path),
        "--mhc-alleles", "HLA-A*02:01"])
    binding_predictions = run_predictor(args)
    assert [prediction.peptide for prediction in binding_predictions] == [
        "SIINFEKL", "GILGFVFTL"]


def test_peptides_file_rejects_no_sequences(tmp_path):
    path = tmp_path / "peptides.txt"
    path.write_text("\n   \n")
    args = parse_args([
        "--mhc-predictor", "random",
        "--input-peptides-file", str(path),
        "--mhc-alleles", "HLA-A*02:01"])
    with pytest.raises(ValueError, match="No peptide sequences found"):
        run_predictor(args)

def test_peptides_file_with_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])
    remove(f.name)
