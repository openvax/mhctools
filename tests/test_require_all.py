"""The release gate must reject incomplete runs, including under xdist."""

import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.parametrize("body, complete", [
    ("def test_ok(): pass", True),
    ("def test_skip(): pytest.skip('missing model')", False),
    ("@pytest.mark.xfail\ndef test_xfail(): assert False", False),
    ("@pytest.mark.xfail\ndef test_xpass(): pass", False),
    ("pytest.skip('missing runtime', allow_module_level=True)", False),
])
def test_require_all_exit_status(tmp_path, body, complete):
    shutil.copyfile(Path(__file__).with_name("conftest.py"), tmp_path / "conftest.py")
    # A passing test keeps collection skips from returning NO_TESTS_COLLECTED.
    (tmp_path / "test_ok.py").write_text("def test_ok(): pass\n")
    (tmp_path / "test_example.py").write_text("import pytest\n" + body + "\n")
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("COV_CORE", "COVERAGE", "PYTEST"))}
    env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    worker_options = [[]]
    if importlib.util.find_spec("xdist"):
        worker_options.append(["-p", "xdist.plugin", "-n", "2"])
    for options in worker_options:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "--require-all", "-q", *options],
            cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)
        assert result.returncode == (0 if complete else 1), result.stdout + result.stderr
        assert ("incomplete test coverage" in result.stdout) is not complete
