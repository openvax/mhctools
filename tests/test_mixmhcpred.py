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

import os
import shutil
import sys
from pathlib import Path

import pytest

from mhctools import MixMHCpred
from mhctools.mixmhcpred import (
    mixmhcpred_version,
    parse_mixmhcpred_output,
    parse_mixmhcpred_results,
    resolve_mixmhcpred_path,
)
from mhctools.pred import COLUMNS, Kind

V3_OUTPUT = """####################
# Output from MixMHCpred (v3.0)
# Alleles: A0201, A0102 - predicted motif for A0102
# Closest Allele (distance score): A0201 (0.0) -- A0101 (0.037)
# Input file: peptides.txt
####################
Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_A0201\t%Rank_A0201\tScore_A0102\t%Rank_A0102
SIINFEKL\t-1.119737\tA0201\t2.698623\t-1.119737\t2.698623\t-3.902923\t41.899867
MLDDFSAGA\t0.182093\tA0201\t0.3\t0.182093\t0.3\t-2.0\t20.0
"""

LEGACY_OUTPUT = """Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_A0201\t%Rank_A0201
MLDDFSAGA\t0.182093\tA0201\t0.3\t0.182093\t0.3
SPEGEETII\t-0.655341\tA0201\t51.0\t-0.655341\t51.0
ILDRIITNA\t0.203906\tA0201\t0.3\t0.203906\t0.3
"""

SINGLE_OUTPUT = """# Output from MixMHCpred (v3.0)
Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_A0201\t%Rank_A0201
SIINFEKL\t-1.119737\tA0201\t2.698623\t-1.119737\t2.698623
"""

SEQUENCE_OUTPUT = """####################
# Output from MixMHCpred (v3.0)
# Predictions for Alleles : novel-one, novel-two
####################
Peptide\tScore_bestAllele\tBestAllele\t%Rank_bestAllele\tScore_novel-one\t%Rank_novel-one\tScore_novel-two\t%Rank_novel-two
SIINFEKL\t-0.25\tnovel-two\t0.23\t-3.5\t32.6\t-0.25\t0.23
"""


def _write(path, text):
    Path(path).write_text(text)
    return path


def test_parse_v3_preserves_all_outputs(tmp_path):
    path = _write(tmp_path / "v3.txt", V3_OUTPUT)
    result = parse_mixmhcpred_output(
        path, alleles=["HLA-A*02:01", "HLA-A*01:02"])

    assert result.version == "3.0"
    assert list(result.table.columns) == [
        "Peptide", "Score_bestAllele", "BestAllele", "%Rank_bestAllele",
        "Score_A0201", "%Rank_A0201", "Score_A0102", "%Rank_A0102",
    ]
    assert result.table.loc[0, "Score_bestAllele"] == pytest.approx(-1.119737)
    assert result.table.loc[0, "BestAllele"] == "A0201"
    assert result.table.loc[0, "%Rank_bestAllele"] == pytest.approx(2.698623)

    a0201, a0102 = result.allele_info
    assert a0201.allele == "HLA-A*02:01"
    assert a0201.native_allele == "A0201"
    assert a0201.closest_training_allele == "HLA-A*02:01"
    assert a0201.distance == 0.0
    assert not a0201.pan_allele
    assert a0102.closest_training_allele == "HLA-A*01:01"
    assert a0102.distance == pytest.approx(0.037)
    assert a0102.pan_allele

    predictions = result.predictions
    assert len(predictions) == 2
    assert len(predictions[0].preds) == 2
    assert predictions[0].presentation.allele == "HLA-A*02:01"
    assert all(p.kind == Kind.pMHC_presentation for p in predictions[0].preds)
    assert all(p.predictor_version == "3.0" for p in predictions[0].preds)
    assert predictions[0].filter(allele="HLA-A*01:02")[0].score == pytest.approx(
        -3.902923)


def test_legacy_parser_remains_compatible(tmp_path):
    path = _write(tmp_path / "legacy.txt", LEGACY_OUTPUT)
    predictions = parse_mixmhcpred_results(path)
    assert len(predictions) == 3
    assert [p.peptide for p in predictions] == [
        "MLDDFSAGA", "SPEGEETII", "ILDRIITNA"]
    assert all(p.allele == "HLA-A*02:01" for p in predictions)

    v3_path = _write(tmp_path / "v3.txt", V3_OUTPUT)
    v3_predictions = parse_mixmhcpred_results(v3_path)
    assert len(v3_predictions) == 2
    assert all(p.allele == "HLA-A*02:01" for p in v3_predictions)


def test_parser_rejects_missing_per_allele_rank(tmp_path):
    bad = V3_OUTPUT.replace("\t%Rank_A0102", "")
    path = _write(tmp_path / "bad.txt", bad)
    with pytest.raises(ValueError, match="%Rank_A0102"):
        parse_mixmhcpred_output(path)


def test_resolve_path_from_environment_or_checkout(tmp_path, monkeypatch):
    executable = _write(tmp_path / "MixMHCpred", "#!/bin/sh\n")
    monkeypatch.setenv("MIXMHCPRED_PATH", str(tmp_path))
    assert resolve_mixmhcpred_path() == str(Path(executable).resolve())
    assert resolve_mixmhcpred_path(tmp_path) == str(Path(executable).resolve())


def test_resolve_path_error_has_install_guidance(monkeypatch):
    monkeypatch.delenv("MIXMHCPRED_PATH", raising=False)
    monkeypatch.setattr(shutil, "which", lambda _: None)
    with pytest.raises(FileNotFoundError, match="official v3.0 release"):
        resolve_mixmhcpred_path()


def test_version_detection(tmp_path):
    executable = _write(
        tmp_path / "MixMHCpred", "#!/bin/sh\necho MixMHCpred3.0\n")
    Path(executable).chmod(0o755)
    assert mixmhcpred_version(executable) == "3.0"


def test_predict_batches_alleles_and_filters_cysteine(tmp_path, monkeypatch):
    predictor = MixMHCpred(
        alleles=["HLA-A*02:01", "HLA-A*01:02"],
        program_name=sys.executable,
        exclude_peptides_with_cysteine=True,
    )
    calls = []

    def fake_run(args, stdout_path):
        calls.append(args)
        input_path = Path(args[args.index("-i") + 1])
        assert input_path.read_text() == "SIINFEKL\n"
        output_path = Path(args[args.index("-o") + 1])
        _write(output_path, V3_OUTPUT.splitlines()[0] + "\n" + "\n".join(
            line for line in V3_OUTPUT.splitlines()[1:]
            if not line.startswith("MLDDFSAGA")) + "\n")
        _write(stdout_path, "DONE\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    results = predictor.predict(["SIINFEKL", "CILGFVFTL"])

    assert len(calls) == 1
    assert calls[0][calls[0].index("-a") + 1] == (
        "HLA-A*02:01,HLA-A*01:02")
    assert "-c" not in calls[0]
    assert len(results) == 2
    assert len(results[0].preds) == 2
    assert results[1].preds == ()
    assert predictor.supported_kinds == (Kind.pMHC_presentation,)
    assert list(results[0].to_dataframe().columns) == list(COLUMNS)
    assert list(predictor.predict_dataframe([]).columns) == list(COLUMNS)


def test_predict_rejects_silent_or_stale_output(monkeypatch):
    predictor = MixMHCpred(
        alleles=["HLA-A*02:01"], program_name=sys.executable)

    def fake_run(args, stdout_path):
        _write(Path(args[args.index("-o") + 1]), LEGACY_OUTPUT)
        _write(stdout_path, "DONE\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    with pytest.raises(RuntimeError, match="do not match its inputs"):
        predictor.predict(["SIINFEKL"])


def test_motif_artifacts_are_returned_and_existing_path_is_safe(
        tmp_path, monkeypatch):
    predictor = MixMHCpred(
        alleles=["HLA-A*02:01"], program_name=sys.executable)
    predictor._version = "3.0"
    output_dir = tmp_path / "motifs"

    def fake_run(args, stdout_path):
        destination = Path(args[args.index("-o") + 1])
        destination.mkdir()
        (destination / "Motifs").mkdir()
        _write(destination / "Binding_predictions.txt", SINGLE_OUTPUT)
        _write(destination / "Data_overview.html", "<html></html>")
        _write(destination / "Motifs" / "A0201_PWM9.png", "fake")
        _write(stdout_path, "DONE\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    result = predictor.predict_detailed(
        ["SIINFEKL"], output_dir=output_dir, output_motifs=True)
    assert result.artifacts.output_dir == str(output_dir.resolve())
    assert result.artifacts.binding_predictions.endswith(
        "Binding_predictions.txt")
    assert result.artifacts.overview_html.endswith("Data_overview.html")
    assert len(result.artifacts.files) == 3
    with pytest.raises(FileExistsError, match="Refusing"):
        predictor.predict_detailed(
            ["SIINFEKL"], output_dir=output_dir, output_motifs=True)


def test_sequence_prediction_captures_alignment_quality_and_artifacts(
        tmp_path, monkeypatch):
    predictor = MixMHCpred(program_name=sys.executable)
    predictor._version = "3.0"
    output_dir = tmp_path / "sequence-output"
    sequences = {
        "novel-one": "ACDEFGHIKLMNPQRSTVWY",
        "novel-two": "YWVTSRQPNMLKIHGFEDCA",
    }

    def fake_run(args, stdout_path):
        destination = Path(args[args.index("-o") + 1])
        destination.mkdir()
        (destination / "Motifs").mkdir()
        sequence_path = Path(args[args.index("-s") + 1])
        shutil.copyfile(sequence_path, destination / "final_alignment.fasta")
        _write(destination / "Binding_predictions.txt", SEQUENCE_OUTPUT)
        _write(
            destination / "Data_overview.html",
            "<h2>novel-one</h2>"
            "<h2>closest Allele = A0201, distance = 0.1</h2>"
            "<h2>closest Allele from database = A0202, distance = 0.01</h2>"
            "<h2>novel-two</h2>"
            "<h2>closest Allele = B0702, distance = 0.2</h2>"
            "<h2>closest Allele from database = B0703, distance = 0.02</h2>",
        )
        _write(destination / "Motifs" / "PLD_sequences.txt", "fake")
        _write(stdout_path, "DONE\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    result = predictor.predict_allele_sequences(
        sequences,
        peptides=["SIINFEKL"],
        output_dir=output_dir,
        output_motifs=True,
    )

    assert result.aligned_sequences == tuple(sequences.items())
    assert result.artifacts.alignment.endswith("final_alignment.fasta")
    assert result.artifacts.overview_html.endswith("Data_overview.html")
    assert result.predictions[0].filter(allele="novel-two")[0].score == -0.25
    first, second = result.allele_info
    assert first.pan_allele and second.pan_allele
    assert first.closest_training_allele == "HLA-A*02:01"
    assert first.closest_database_allele == "HLA-A*02:02"
    assert first.distance == 0.1
    assert first.database_distance == 0.01


def test_alignment_only_retains_sequence_provenance(monkeypatch):
    predictor = MixMHCpred(program_name=sys.executable)
    predictor._version = "3.0"

    def fake_run(args, stdout_path):
        destination = Path(args[args.index("-o") + 1])
        destination.mkdir()
        shutil.copyfile(
            Path(args[args.index("-s") + 1]),
            destination / "final_alignment.fasta",
        )
        _write(stdout_path, "DONE\n")
        return "DONE\n"

    monkeypatch.setattr(predictor, "_run_command", fake_run)
    result = predictor.predict_allele_sequences(
        {"novel-one": "ACDEFGHIKLMNPQRSTVWY"})

    assert result.table.empty
    assert len(result.allele_info) == 1
    assert result.allele_info[0].allele == "novel-one"
    assert result.allele_info[0].native_allele == "novel-one"
    assert result.allele_info[0].pan_allele


MIXMHCPRED_V3 = os.environ.get("MIXMHCPRED_V3_PATH")
requires_mixmhcpred_v3 = pytest.mark.skipif(
    not MIXMHCPRED_V3,
    reason="set MIXMHCPRED_V3_PATH to an official MixMHCpred v3.0 checkout",
)


@requires_mixmhcpred_v3
def test_mixmhcpred_v3_end_to_end_known_and_pan_alleles():
    predictor = MixMHCpred(
        alleles=["HLA-A*02:01", "HLA-A*01:02"],
        program_name=MIXMHCPRED_V3,
    )
    result = predictor.predict_detailed(["SIINFEKL"])
    assert result.version == "3.0"
    assert result.predictions[0].filter(allele="HLA-A*02:01")[0].score == (
        pytest.approx(-1.119737, abs=1e-6))
    assert result.predictions[0].filter(allele="HLA-A*01:02")[0].score == (
        pytest.approx(-3.902923, abs=1e-6))
    assert [info.pan_allele for info in result.allele_info] == [False, True]
    assert result.allele_info[1].closest_training_allele == "HLA-A*01:01"


@requires_mixmhcpred_v3
def test_mixmhcpred_v3_end_to_end_sequence_prediction():
    executable = Path(resolve_mixmhcpred_path(MIXMHCPRED_V3))
    upstream_fixture = executable.parent / "input" / "To_align_sequences.fasta"
    if not upstream_fixture.is_file() or shutil.which("mafft") is None:
        pytest.skip("official sequence fixture and MAFFT are required")

    result = MixMHCpred(
        program_name=MIXMHCPRED_V3).predict_allele_sequences(
            upstream_fixture, peptides=["SIINFEKL"])
    assert [name for name, _ in result.aligned_sequences] == [
        "HLA-A*01:07", "Mafa-B*008:02", "Mamu-B*030:02:01:08"]
    assert len(result.predictions) == 1
    assert len(result.predictions[0].preds) == 3
    assert result.table.loc[0, "BestAllele"] == "Mamu-B*030:02:01:08"
