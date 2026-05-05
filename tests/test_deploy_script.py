import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPLOY_SH = ROOT / "deploy.sh"


def test_deploy_script_is_valid_bash():
    subprocess.run(
        ["bash", "-n", str(DEPLOY_SH)],
        check=True,
        cwd=ROOT)


def test_deploy_script_matches_documented_release_guards():
    text = DEPLOY_SH.read_text()

    assert "git rev-parse --abbrev-ref HEAD" in text
    assert '"main" && "$branch" != "master"' in text
    assert "git diff --quiet" in text
    assert "git diff --cached --quiet" in text
    assert "python3 -m pip install" not in text
    assert 'PYTHON="${PYTHON:-python}"' in text
    assert "git tag -a" in text
    assert "git push origin \"$tag\"" in text


def test_deploy_script_is_executable():
    assert os.access(DEPLOY_SH, os.X_OK)
