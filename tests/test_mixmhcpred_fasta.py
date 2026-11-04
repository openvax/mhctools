"""FASTA handling through the public allele-sequence prediction interface."""

from pathlib import Path
import sys

import pytest

from mhctools import MixMHCpred


def test_sequence_fasta_preserves_headers_order_and_alignment_gaps(tmp_path, monkeypatch):
    source = tmp_path / "alleles.fasta"
    source.write_text(
        "; input comment\n\n>novel-one full description\nACDE\nFGHI\n"
        "; between records\n>novel-two another description\nKLMN\nPQRS\n")
    predictor = MixMHCpred(program_name=sys.executable)
    predictor._version = "3.0"

    def fake_run(args, stdout_path):
        native_input = Path(args[args.index("-s") + 1]).read_text()
        assert native_input == (
            ">novel-one full description\nACDEFGHI\n"
            ">novel-two another description\nKLMNPQRS\n")
        output = Path(args[args.index("-o") + 1])
        output.mkdir()
        (output / "final_alignment.fasta").write_text(
            ">novel-one full description\nACD-E\nFGHI\n"
            "; alignment comment\n>novel-two another description\nKLMNPQRS-\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    result = predictor.predict_allele_sequences(source)
    assert result.aligned_sequences == (
        ("novel-one full description", "ACD-EFGHI"),
        ("novel-two another description", "KLMNPQRS-"),
    )
    assert [info.native_allele for info in result.allele_info] == [
        "novel-one full description", "novel-two another description"]


@pytest.mark.parametrize("content,error", [
    ("ACDE\n>allele\nFGHI\n", "before the first FASTA header"),
    ("\n; comment\nACDE\n>allele\nFGHI\n", "before the first FASTA header"),
    ("; only comment\n\n", "At least one"),
    (">same\nACDE\n>same\nFGHI\n", "non-empty and unique"),
    (">\nACDE\n", "non-empty and unique"),
    (">allele\n", "Invalid MHC-I sequence"),
])
def test_sequence_fasta_rejects_malformed_input_before_running(tmp_path, monkeypatch, content, error):
    source = tmp_path / "alleles.fasta"
    source.write_text(content)
    predictor = MixMHCpred(program_name=sys.executable)
    monkeypatch.setattr(predictor, "_run_command",
                        lambda *args: pytest.fail("invalid FASTA reached the predictor"))
    with pytest.raises(ValueError, match=error):
        predictor.predict_allele_sequences(source)


def test_native_alignment_without_records_is_an_error(monkeypatch):
    predictor = MixMHCpred(program_name=sys.executable)
    predictor._version = "3.0"

    def fake_run(args, stdout_path):
        output = Path(args[args.index("-o") + 1])
        output.mkdir()
        (output / "final_alignment.fasta").write_text("not a FASTA alignment\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    with pytest.raises(RuntimeError, match="no FASTA records"):
        predictor.predict_allele_sequences({"allele": "ACDEFGHI"})
