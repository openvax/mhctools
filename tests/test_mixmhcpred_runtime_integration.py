"""Compare isolated inference with the unmodified official launchers."""

import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from tempfile import TemporaryFile
import time

import pandas as pd
import pytest

from mhctools import MixMHCpred, PRIME
from mhctools.cli.script import main
from mhctools.mixmhcpred import (
    _collect_artifacts, _read_fasta, parse_mixmhcpred_output,
    resolve_mixmhcpred_path,
)
from mhctools.prime import parse_prime_results
from mhctools.process_helpers import run_command
from .common import process_is_running


PEPTIDES = ["SIINFEKL", "MLDDFSAGA", "SIINFEKL"]
ALLELES = ["HLA-A*02:01", "HLA-A*01:02"]


@pytest.fixture
def native_runtime():
    selected = os.environ.get("MIXMHCPRED_PYTHON")
    executable = os.environ.get("MIXMHCPRED_V3_PATH")
    if not selected or not executable:
        pytest.skip("provision the gfeller group with MIXMHCPRED_PYTHON")
    executable = Path(resolve_mixmhcpred_path(executable))
    assert executable.name == "MixMHCpred", "use the official launcher, not a test wrapper"
    assert Path(selected).is_file()
    # Native reference uses the venv's real python3 on PATH, with no mhctools shim.
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    environment.pop("PYTHONHOME", None)
    environment.update(PYTHONNOUSERSITE="1", MPLBACKEND="Agg")
    environment["PATH"] = str(Path(selected).parent) + os.pathsep + environment.get("PATH", os.defpath)
    before = (sys.executable, pd.__version__)
    before_environment = {key: value for key, value in os.environ.items() if key != "PYTEST_CURRENT_TEST"}
    yield str(executable), selected, environment
    assert (sys.executable, pd.__version__) == before
    after_environment = {key: value for key, value in os.environ.items() if key != "PYTEST_CURRENT_TEST"}
    # Report changed names only, so assertion output cannot expose secrets.
    changed = sorted(key for key in before_environment.keys() | after_environment.keys()
                     if before_environment.get(key) != after_environment.get(key))
    assert changed == []


def _native(args, environment, timeout=180):
    started = time.monotonic()
    with TemporaryFile() as output:
        try:
            run_command(
                args, env=environment, timeout=timeout,
                terminate_process_group=True,
                redirect_stdout_file=output, redirect_stderr_to_stdout=True)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
            output.seek(0, os.SEEK_END)
            output.seek(max(0, output.tell() - 2000))
            tail = output.read().decode("utf-8", errors="replace").strip()
            reason = (f"timed out (limit {timeout:g}s)"
                      if isinstance(error, subprocess.TimeoutExpired)
                      else f"exited with status {error.returncode}")
            pytest.fail(
                f"{args[0]} {reason} after {time.monotonic() - started:.2f}s.\n"
                f"Output tail:\n{tail or '(empty output)'}", pytrace=False)


def test_native_reference_failure_reports_bounded_output():
    command = [sys.executable, "-c",
               "import sys; print('discarded-prefix' + 'x' * 5000, flush=True); "
               "print('stdout marker', flush=True); "
               "print('stderr marker', file=sys.stderr, flush=True); sys.exit(7)"]
    with pytest.raises(pytest.fail.Exception) as error:
        _native(command, dict(os.environ))
    message = str(error.value)
    assert "exited with status 7" in message
    assert "stdout marker" in message and "stderr marker" in message
    assert "discarded-prefix" not in message
    assert len(message) < 2400


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_native_reference_timeout_stops_children(tmp_path):
    child = tmp_path / "child.py"
    child.write_text(
        "import os, sys, time\nfrom pathlib import Path\n"
        "print('reference child active', flush=True)\n"
        "Path(sys.argv[1]).write_text(str(os.getpid()))\ntime.sleep(60)\n")
    parent = tmp_path / "parent.py"
    parent.write_text(
        "import subprocess, sys\n"
        "subprocess.run([sys.executable, sys.argv[1], sys.argv[2]])\n")
    ready = tmp_path / "child.pid"
    try:
        with pytest.raises(pytest.fail.Exception, match="timed out \\(limit 3s\\)") as error:
            _native([sys.executable, str(parent), str(child), str(ready)],
                    dict(os.environ), timeout=3)
        assert ready.is_file(), "Reference child did not start before the timeout"
        assert "reference child active" in str(error.value)
        assert not process_is_running(int(ready.read_text()))
    finally:
        if ready.exists() and process_is_running(int(ready.read_text())):
            os.kill(int(ready.read_text()), signal.SIGKILL)


def _input(tmp_path):
    path = tmp_path / "peptides.txt"
    path.write_text("\n".join(PEPTIDES) + "\n")
    return str(path)


@pytest.mark.requires_external_tool
def test_isolated_multi_allele_rows_motifs_and_cli_match_native(native_runtime, tmp_path, monkeypatch):
    executable, selected, environment = native_runtime
    peptide_file = _input(tmp_path)
    reference_dir = tmp_path / "native-motifs"
    _native([executable, "-i", peptide_file, "-o", str(reference_dir),
             "-a", ",".join(ALLELES), "-m", "1"], environment)
    expected = parse_mixmhcpred_output(reference_dir / "Binding_predictions.txt", alleles=ALLELES)
    # Explicit API selection wins over a conflicting environment selection.
    monkeypatch.setenv("MIXMHCPRED_PYTHON", "/missing/conflicting/python")
    model = MixMHCpred(alleles=ALLELES, program_name=executable, python_executable=selected)
    actual = model.predict_detailed(PEPTIDES, output_dir=tmp_path / "isolated-motifs", output_motifs=True)
    pd.testing.assert_frame_equal(actual.table, expected.table)
    assert actual.allele_info == expected.allele_info
    assert actual.table.Peptide.tolist() == PEPTIDES
    assert actual.runtime_info["backend"]["python_executable"] == selected
    assert int(actual.runtime_info["backend"]["pandas_version"].split(".")[0]) < 3
    assert actual.runtime_info["host_pandas_version"] == pd.__version__
    assert Path(actual.artifacts.overview_html).stat().st_size > 100
    native_files = {p.relative_to(reference_dir) for p in reference_dir.rglob("*") if p.is_file()}
    actual_dir = Path(actual.artifacts.output_dir)
    assert {p.relative_to(actual_dir) for p in actual_dir.rglob("*") if p.is_file()} == native_files
    assert any(p.suffix == ".png" for p in native_files)
    for path in native_files:
        if path.name == "Binding_predictions.txt":
            continue  # Input filenames differ; all data columns were compared above.
        if path.suffix in {".txt", ".png"}:
            assert (actual_dir / path).read_bytes() == (reference_dir / path).read_bytes()

    csv_path = tmp_path / "cli.csv"
    main(["--mhc-predictor", "mixmhcpred", "--mhc-predictor-path", executable,
          "--mixmhcpred-python", selected, "--mhc-alleles", ",".join(ALLELES),
          "--input-peptides-file", peptide_file, "--output-csv", str(csv_path)])
    csv = pd.read_csv(csv_path)
    assert csv.peptide.tolist() == [peptide for peptide in PEPTIDES for _ in ALLELES]
    assert csv.allele.tolist() == ALLELES * len(PEPTIDES)
    assert csv.score.tolist() == pytest.approx([
        pred.score for row in expected.predictions for pred in row.preds], abs=5e-6)
    assert csv.percentile_rank.tolist() == pytest.approx([
        pred.percentile_rank for row in expected.predictions for pred in row.preds], abs=5e-5)


@pytest.mark.requires_external_tool
def test_isolated_sequence_alignment_prediction_and_artifacts_match_native(native_runtime, tmp_path):
    executable, selected, environment = native_runtime
    assert shutil.which("mafft"), "provision MAFFT for the sequence regression"
    sequences = Path(executable).parent / "input" / "To_align_sequences.fasta"
    assert sequences.is_file()
    reference_dir = tmp_path / "native-sequences"
    _native([executable, "-s", str(sequences), "-i", _input(tmp_path),
             "-o", str(reference_dir), "-p", "1", "-m", "1"], environment)
    aligned = _read_fasta(reference_dir / "final_alignment.fasta")
    expected = parse_mixmhcpred_output(
        reference_dir / "Binding_predictions.txt", alleles=[name for name, _ in aligned],
        artifacts=_collect_artifacts(reference_dir), sequence_mode=True, aligned_sequences=aligned)
    model = MixMHCpred(program_name=executable, python_executable=selected)
    actual = model.predict_allele_sequences(
        sequences, peptides=PEPTIDES, output_dir=tmp_path / "isolated-sequences", output_motifs=True)
    pd.testing.assert_frame_equal(actual.table, expected.table)
    assert actual.aligned_sequences == aligned
    assert actual.allele_info == expected.allele_info
    assert [name for name, _ in aligned] == ["HLA-A*01:07", "Mafa-B*008:02", "Mamu-B*030:02:01:08"]
    assert Path(actual.artifacts.overview_html).stat().st_size > 100
    assert any(Path(path).suffix == ".png" for path in actual.artifacts.files)
    alignment_only = model.predict_allele_sequences(sequences)
    assert alignment_only.aligned_sequences == aligned
    assert alignment_only.runtime_info["backend"]["python_executable"] == selected


@pytest.mark.requires_external_tool
def test_prime_nested_isolated_runtime_matches_native(native_runtime, tmp_path, monkeypatch):
    executable, selected, environment = native_runtime
    prime = os.environ.get("PRIME_EXECUTABLE")
    assert prime and Path(prime).is_file(), "provision PRIME alongside MixMHCpred"
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    output = tmp_path / "prime-native.txt"
    _native([prime, "-i", _input(tmp_path), "-o", str(output),
             "-a", ",".join(alleles), "-mix", executable], environment)
    expected = parse_prime_results(output, alleles)
    monkeypatch.setenv("MIXMHCPRED_PYTHON", "/missing/conflicting/python")
    model = PRIME(alleles=alleles, program_name=prime, mixmhcpred_path=executable,
                  mixmhcpred_python=selected)
    actual = model.predict(PEPTIDES)
    assert [{pred.allele: pred for pred in row.preds} for row in actual] == [
        {pred.allele: pred for pred in expected[peptide]} for peptide in PEPTIDES]
    assert model.runtime_info["backend"]["python_executable"] == selected
    assert model.runtime_info["host_pandas_version"] == pd.__version__
