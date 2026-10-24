"""Exercise deployment coordination with real Git and a local bare origin."""

from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import pytest


ROOT = Path(__file__).resolve().parents[1]
LOCK = "refs/heads/release-lock"

# Only replace expensive/external boundaries. Git, version edits, conditional
# pushes and exit traps execute unchanged.
PYTHON_STUB = r'''
import os
from pathlib import Path
import sys
import time

args = sys.argv[1:]
if args[:2] == ["-m", "build"]:
    if "--version" in args:
        sys.exit(0)
    stage = "build"
elif args[:2] == ["-m", "twine"]:
    if "--version" in args:
        sys.exit(0)
    stage = "upload"
elif args and args[0] == "scripts/check_pypi_release.py":
    stage = "check"
elif args and args[0] == "gate":
    stage = args[1]
else:
    os.execv(sys.executable, [sys.executable, *args])

state = Path(os.environ["RELEASE_TEST_STATE"])
with (state / "events").open("a") as output:
    output.write(stage + "\n")
(state / ("entered-" + stage)).touch()
if stage == "test" and (state / "block-test").exists():
    deadline = time.monotonic() + 20
    while not (state / "continue").exists():
        if time.monotonic() > deadline:
            sys.exit("test gate was not released")
        time.sleep(0.01)
if (state / ("fail-" + stage)).exists():
    sys.exit("failed " + stage)
if stage == "check" and (state / "published").exists():
    sys.exit("release already exists on PyPI")
if stage == "build":
    Path("dist").mkdir()
    Path("dist/test.whl").touch()
'''


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True,
        text=True, timeout=15).stdout.strip()


def remote_ref(repo, ref):
    output = git(repo, "ls-remote", "origin", ref)
    return output.split()[0] if output else None


def events(state):
    path = state / "events"
    return path.read_text().splitlines() if path.exists() else []


@pytest.fixture
def checkout(tmp_path, monkeypatch):
    # Disposable repositories must not inherit user hooks/signing/config.
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    origin = tmp_path / "origin.git"
    origin.mkdir()
    git(origin, "init", "--bare", "--initial-branch=main")
    seed = tmp_path / "seed"
    seed.mkdir()
    git(seed, "init", "--initial-branch=main")

    def identity(repo):
        git(repo, "config", "user.name", "Release tests")
        git(repo, "config", "user.email", "release-tests@example.invalid")

    identity(seed)
    shutil.copy2(ROOT / "deploy.sh", seed / "deploy.sh")
    (seed / "scripts").mkdir()
    shutil.copy2(ROOT / "scripts/check_pypi_release.py", seed / "scripts")
    (seed / "mhctools").mkdir()
    (seed / "mhctools/__init__.py").write_text('__version__ = "0.0.1"\n')
    (seed / ".gitignore").write_text("dist/\n")
    for stage in ("lint", "test"):
        script = seed / (stage + ".sh")
        script.write_text('#!/bin/sh\nexec "$PYTHON" gate ' + stage + '\n')
        script.chmod(0o755)
    git(seed, "add", ".")
    git(seed, "commit", "-m", "Release fixture")
    git(seed, "remote", "add", "origin", str(origin))
    git(seed, "push", "origin", "main")

    wrapper = tmp_path / "python-wrapper"
    wrapper.write_text(f"#!{sys.executable}\n" + PYTHON_STUB)
    wrapper.chmod(0o755)

    def make(name):
        repo = tmp_path / name
        git(tmp_path, "clone", str(origin), str(repo))
        identity(repo)
        state = tmp_path / (name + "-state")
        state.mkdir()
        env = dict(os.environ, PYTHON=str(wrapper), RELEASE_TEST_STATE=str(state))
        return repo, state, env

    return make


def deploy(repo, env, *args):
    return subprocess.run(
        ["./deploy.sh", *args], cwd=repo, env=env, text=True,
        capture_output=True, timeout=20)


@contextmanager
def paused_deploy(repo, state, env):
    (state / "block-test").touch()
    with (state / "output").open("w") as output:
        process = subprocess.Popen(
            ["./deploy.sh"], cwd=repo, env=env, stdout=output,
            stderr=subprocess.STDOUT, start_new_session=True)
        try:
            deadline = time.monotonic() + 15
            while not (state / "entered-test").exists():
                assert process.poll() is None, (state / "output").read_text()
                assert time.monotonic() < deadline, (state / "output").read_text()
                time.sleep(0.01)
            yield process
        finally:
            (state / "continue").touch()
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=5)
                raise


def test_competing_clones_have_one_owner_and_one_upload(checkout):
    first, first_state, first_env = checkout("first")
    second, second_state, second_env = checkout("second")
    with paused_deploy(first, first_state, first_env) as owner:
        owner_commit = remote_ref(first, LOCK)
        assert owner_commit and owner_commit != git(first, "rev-parse", "HEAD")
        result = deploy(second, second_env)
        assert result.returncode != 0
        assert "cannot acquire release lock" in result.stderr
        assert events(second_state) == ["check"]
        assert remote_ref(first, LOCK) == owner_commit
    assert owner.returncode == 0, (first_state / "output").read_text()
    assert events(first_state) == ["check", "lint", "test", "build", "check", "upload"]
    assert remote_ref(first, LOCK) is None
    assert remote_ref(first, "refs/tags/v0.0.1")
    result = deploy(second, second_env)
    assert result.returncode != 0 and "tag v0.0.1 already exists" in result.stderr
    assert events(second_state) == ["check"]


@pytest.mark.parametrize("stage", ["lint", "test", "build", "upload"])
def test_failure_releases_ownership(checkout, stage):
    repo, state, env = checkout("release")
    (state / ("fail-" + stage)).touch()
    result = deploy(repo, env)
    assert result.returncode != 0 and "failed " + stage in result.stderr
    assert remote_ref(repo, LOCK) is None
    assert remote_ref(repo, "refs/tags/v0.0.1") is None


@pytest.mark.parametrize("change", ["pypi", "tag", "lookup-failure"])
def test_publication_rechecked_after_long_gates(checkout, change):
    repo, state, env = checkout("release")
    with paused_deploy(repo, state, env) as process:
        if change == "tag":
            publisher, _, _ = checkout("publisher")
            git(publisher, "tag", "v0.0.1")
            git(publisher, "push", "origin", "v0.0.1")
        else:
            (state / ("published" if change == "pypi" else "fail-check")).touch()
    assert process.returncode != 0, (state / "output").read_text()
    assert "build" in events(state)
    assert "upload" not in events(state)
    assert remote_ref(repo, LOCK) is None


def test_changed_owner_is_not_deleted_or_allowed_to_publish(checkout):
    repo, state, env = checkout("release")
    other, _, _ = checkout("other")
    with paused_deploy(repo, state, env) as process:
        replacement = git(other, "rev-parse", "HEAD")
        git(other, "push", "--force", "origin", f"{replacement}:{LOCK}")
    assert process.returncode != 0
    assert "ownership changed" in (state / "output").read_text()
    assert "upload" not in events(state)
    assert remote_ref(repo, LOCK) == replacement


def test_interrupt_releases_ownership(checkout):
    repo, state, env = checkout("release")
    with paused_deploy(repo, state, env) as process:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=15)
    assert process.returncode != 0
    assert remote_ref(repo, LOCK) is None
    assert "upload" not in events(state)


@pytest.mark.parametrize("condition", ["branch", "dirty", "staged"])
def test_checkout_guards_run_before_remote_changes(checkout, condition):
    repo, state, env = checkout("release")
    if condition == "branch":
        git(repo, "checkout", "-b", "feature")
    else:
        (repo / "mhctools/__init__.py").write_text('__version__ = "0.0.2"\n')
        if condition == "staged":
            git(repo, "add", "mhctools/__init__.py")
    assert deploy(repo, env).returncode != 0
    assert events(state) == []
    assert remote_ref(repo, LOCK) is None


def test_requested_version_is_committed_and_tagged(checkout):
    repo, state, env = checkout("release")
    result = deploy(repo, env, "0.0.2")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (repo / "mhctools/__init__.py").read_text() == '__version__ = "0.0.2"\n'
    commit = git(repo, "rev-parse", "HEAD")
    assert git(repo, "rev-parse", "v0.0.2^{}") == commit
    assert remote_ref(repo, "refs/heads/main") == commit
    assert remote_ref(repo, LOCK) is None
