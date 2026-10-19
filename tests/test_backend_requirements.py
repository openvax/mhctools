"""Installed but broken optional executables must never disappear as skips."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from .common import require_executable


def test_absent_executable_skips():
    with pytest.raises(pytest.skip.Exception, match="not found on PATH"):
        require_executable("mhctools-certainly-missing-executable")


def test_installed_executable_runs():
    require_executable(sys.executable, args=("-c", "print('ready')"))


def test_installed_executable_failure_is_not_a_skip():
    with pytest.raises(pytest.fail.Exception, match="python2: Too many arguments"):
        require_executable(sys.executable, args=(
            "-c", "import sys; sys.exit('python2: Too many arguments')"))


@pytest.mark.parametrize("runtime", ["missing", "broken", "working"])
def test_legacy_parameter_setup_distinguishes_absent_and_broken(tmp_path, runtime):
    """Exercise the actual epitope-test fixtures without licensed executables."""
    shutil.copyfile(Path(__file__).with_name("conftest.py"), tmp_path / "conftest.py")
    (tmp_path / "test_example.py").write_text('''
import pytest
from mhctools import NetMHC3, NetMHCcons
from tests.test_known_class1_epitopes import mhc_class

@pytest.mark.parametrize("mhc_class", [NetMHC3, NetMHCcons], indirect=True)
def test_runtime(mhc_class):
    pass
''')
    if runtime != "missing":
        for program in ("netMHC-3.4", "netMHCcons"):
            executable = tmp_path / program
            executable.write_text("#!/bin/sh\necho 'runtime diagnostic' >&2\nexit %d\n" %
                                  (1 if runtime == "broken" else 0))
            executable.chmod(0o755)
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("COV_CORE", "COVERAGE", "PYTEST"))}
    env.update(PATH=str(tmp_path), PYTEST_DISABLE_PLUGIN_AUTOLOAD="1",
               PYTHONPATH=str(Path(__file__).resolve().parent.parent))
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-rs"], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode == (1 if runtime == "broken" else 0), result.stdout + result.stderr
    expected = {"missing": "2 skipped", "broken": "2 errors", "working": "2 passed"}
    assert expected[runtime] in result.stdout, result.stdout
    if runtime == "broken":
        assert "runtime diagnostic" in result.stdout
