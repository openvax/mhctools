"""The quality job must discover new offline tests and exclude only integrations."""

import os
from pathlib import Path
import shlex
import subprocess
import sys


def test_quality_command_discovers_new_files_and_keeps_offline_tests(tmp_path):
    root = Path(__file__).resolve().parents[1]
    workflow = (root / '.github/workflows/tests.yml').read_text()
    quality = workflow.split('\n  integration-', 1)[0]
    command = next(line.strip().removeprefix('run: ') for line in quality.splitlines()
                   if line.strip().startswith('run: python -m pytest '))
    args = shlex.split(command)
    test_dir = tmp_path / 'tests'
    test_dir.mkdir()
    (test_dir / 'test_new_file.py').write_text('def test_new(): pass\n')
    (test_dir / 'test_mixed.py').write_text(
        'import pytest\n'
        'def test_offline(): pass\n'
        '@pytest.mark.requires_external_tool\n'
        'def test_external(): assert False, "external tool must not run"\n')
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(('COV_CORE', 'COVERAGE', 'PYTEST'))}
    env['PYTEST_DISABLE_PLUGIN_AUTOLOAD'] = '1'
    result = subprocess.run(
        [sys.executable, *args[1:], '-c', str(root / 'pyproject.toml')],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert '2 passed, 1 deselected' in result.stdout
