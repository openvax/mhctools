"""Exercise deployment coordination with real Git and a local bare origin."""

from contextlib import contextmanager
import json
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
STARTUP_IDLE_TIMEOUT = 20.0
STARTUP_MAX_TIMEOUT = 60.0
COMMAND_TIMEOUT = 60.0
GIT_TIMEOUT = 30.0
WAIT_TIMEOUT = 60.0

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
    deadline = time.monotonic() + 60
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
        text=True, timeout=GIT_TIMEOUT).stdout.strip()


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
        capture_output=True, timeout=COMMAND_TIMEOUT)


def progress_snapshot(paths):
    """Capture observable file progress (existence, size, mtime_ns) for paths."""
    snapshot = []
    for path in paths:
        try:
            stat = path.stat()
        except FileNotFoundError:
            snapshot.append((False, 0, 0))
        else:
            snapshot.append((True, stat.st_size, stat.st_mtime_ns))
    return tuple(snapshot)


def wait_for_stage(
        target,
        process,
        output_path,
        *,
        progress_paths=(),
        events_path=None,
        idle_timeout=STARTUP_IDLE_TIMEOUT,
        max_timeout=STARTUP_MAX_TIMEOUT,
        poll_interval=0.01,
        clock=time.monotonic,
        sleep=time.sleep):
    """Wait until *target* appears while tolerating slow but advancing startup.

    Extends the stall window whenever any path in *progress_paths* changes
    (for example captured logs or recorded deployment stages), while always
    enforcing the hard upper bound *max_timeout*.
    """
    watched = tuple(dict.fromkeys((output_path, *progress_paths)))
    started = clock()
    last_progress = started
    last_snapshot = progress_snapshot(watched)
    while not target.exists():
        returncode = process.poll()
        if returncode is not None:
            output = output_path.read_text() if output_path.exists() else ""
            raise AssertionError(
                f"deploy process exited early with status {returncode} before "
                f"creating {target.name}:\n{output}"
            )
        now = clock()
        current_snapshot = progress_snapshot(watched)
        if current_snapshot != last_snapshot:
            last_snapshot = current_snapshot
            last_progress = now
        elapsed = now - started
        idle = now - last_progress
        if idle >= idle_timeout or elapsed >= max_timeout:
            output = output_path.read_text() if output_path.exists() else ""
            observed = (
                events_path.read_text().splitlines()
                if events_path is not None and events_path.exists() else []
            )
            reason = (
                f"no progress for {idle:.2f}s (idle limit {idle_timeout:.2f}s)"
                if idle >= idle_timeout
                else f"exceeded max startup bound {max_timeout:.2f}s"
            )
            raise AssertionError(
                f"Timed out waiting for {target.name} after {elapsed:.2f}s "
                f"({reason}); observed events={observed!r}:\n{output}"
            )
        sleep(poll_interval)


@contextmanager
def paused_deploy(repo, state, env):
    (state / "block-test").touch()
    output_path = state / "output"
    with output_path.open("w") as output:
        process = subprocess.Popen(
            ["./deploy.sh"], cwd=repo, env=env, stdout=output,
            stderr=subprocess.STDOUT, start_new_session=True)
        try:
            wait_for_stage(
                state / "entered-test",
                process,
                output_path,
                progress_paths=(
                    output_path,
                    state / "events",
                    state / "entered-check",
                    state / "entered-lint",
                ),
                events_path=state / "events",
            )
            yield process
        finally:
            (state / "continue").touch()
            try:
                process.wait(timeout=WAIT_TIMEOUT)
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
        process.wait(timeout=WAIT_TIMEOUT)
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


def test_wait_for_stage_extends_idle_window_on_progress_within_max_bound(tmp_path):
    target = tmp_path / "entered-test"
    output = tmp_path / "output"
    events_file = tmp_path / "events"
    output.write_text("starting\n")
    now = [0.0]

    def fake_clock():
        return now[0]

    def fake_sleep(seconds):
        now[0] += seconds

    class AdvancingProcess:
        def __init__(self):
            self.ticks = 0

        def poll(self):
            self.ticks += 1
            if self.ticks == 4:
                events_file.write_text("check\n")
            elif self.ticks == 8:
                output.write_text("starting\nacquired release lock\n")
            elif self.ticks == 12:
                target.touch()
            return None

    wait_for_stage(
        target,
        AdvancingProcess(),
        output,
        progress_paths=(output, events_file),
        events_path=events_file,
        idle_timeout=5.0,
        max_timeout=20.0,
        poll_interval=1.0,
        clock=fake_clock,
        sleep=fake_sleep,
    )
    assert target.exists()
    assert now[0] > 5.0


def test_wait_for_stage_fails_on_stall_and_enforces_max_bound(tmp_path):
    target = tmp_path / "entered-test"
    output = tmp_path / "output"
    events_file = tmp_path / "events"
    output.write_text("acquired release lock\n")
    events_file.write_text("check\n")

    def make_clock():
        now = [0.0]

        def clock():
            return now[0]

        def sleep(seconds):
            now[0] += seconds

        return clock, sleep

    class RunningProcess:
        def poll(self):
            return None

    clock, sleep = make_clock()
    with pytest.raises(AssertionError, match="no progress") as exc_info:
        wait_for_stage(
            target,
            RunningProcess(),
            output,
            progress_paths=(output, events_file),
            events_path=events_file,
            idle_timeout=5.0,
            max_timeout=20.0,
            poll_interval=1.0,
            clock=clock,
            sleep=sleep,
        )
    assert "acquired release lock" in str(exc_info.value)
    assert "check" in str(exc_info.value)

    class NoisyProcess:
        def __init__(self):
            self.counter = 0

        def poll(self):
            self.counter += 1
            output.write_text(f"tick {self.counter}\n")
            return None

    clock, sleep = make_clock()
    with pytest.raises(AssertionError, match="max startup bound"):
        wait_for_stage(
            target,
            NoisyProcess(),
            output,
            progress_paths=(output,),
            idle_timeout=10.0,
            max_timeout=5.0,
            poll_interval=1.0,
            clock=clock,
            sleep=sleep,
        )

    class ExitedProcess:
        def poll(self):
            return 2

    clock, sleep = make_clock()
    with pytest.raises(AssertionError, match="exited early with status 2"):
        wait_for_stage(
            target,
            ExitedProcess(),
            output,
            idle_timeout=10.0,
            max_timeout=20.0,
            poll_interval=1.0,
            clock=clock,
            sleep=sleep,
        )


SSH_TRANSPORT_STUB = r'''
import json
import os
from pathlib import Path
import signal
import sys
import time

argv = sys.argv[1:]
remote_command = argv[-1]
state = Path(os.environ["RELEASE_TEST_STATE"])
log_path = state / "ssh-calls.jsonl"
existing = []
if log_path.exists():
    for line in log_path.read_text().splitlines():
        if line.strip():
            try:
                existing.append(json.loads(line))
            except ValueError:
                pass
if "git-upload-pack" in remote_command:
    kind = "upload-pack"
elif "git-receive-pack" in remote_command:
    kind = "receive-pack"
else:
    kind = "other"
kind_count = sum(1 for item in existing if item.get("kind") == kind) + 1
entry = {
    "index": len(existing) + 1,
    "kind": kind,
    "kind_count": kind_count,
    "argv": argv,
    "remote_command": remote_command,
}
payload = (json.dumps(entry) + "\n").encode()
fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
try:
    os.write(fd, payload)
finally:
    os.close(fd)

stall_spec = os.environ.get("SSH_STALL_MATCH", "")
should_stall = False
for token in [item.strip() for item in stall_spec.split(",") if item.strip()]:
    if token == f"{kind}:{kind_count}" or token == f"any:{entry['index']}":
        should_stall = True
        break

if should_stall:
    (state / "stalled-pids").open("a").write(f"{os.getpid()}\n")
    if os.environ.get("SSH_IGNORE_SIGTERM") == "1":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    time.sleep(300)
    sys.exit(0)

os.execvp("/bin/sh", ["/bin/sh", "-c", remote_command])
'''


def configure_ssh_transport(
        tmp_path,
        repo,
        state,
        env,
        *,
        stall_match="",
        ignore_sigterm=False,
        git_ssh_command=None,
        extra_env=None):
    """Route ``origin`` through an offline fake SSH transport script."""
    ssh_bin = tmp_path / (repo.name + "-ssh-bin")
    ssh_bin.mkdir(exist_ok=True)
    ssh_script = ssh_bin / "ssh"
    ssh_script.write_text(f"#!{sys.executable}\n" + SSH_TRANSPORT_STUB)
    ssh_script.chmod(0o755)

    origin_path = git(repo, "remote", "get-url", "origin")
    if not origin_path.startswith("ssh://"):
        git(
            repo,
            "remote",
            "set-url",
            "origin",
            f"ssh://release-test.invalid{origin_path}",
        )
    git(repo, "config", "ssh.variant", "ssh")

    updated = dict(env)
    updated["PATH"] = f"{ssh_bin}{os.pathsep}{updated.get('PATH', '')}"
    if git_ssh_command is None:
        updated.pop("GIT_SSH_COMMAND", None)
    else:
        updated["GIT_SSH_COMMAND"] = git_ssh_command.format(ssh=ssh_script)
    if stall_match:
        updated["SSH_STALL_MATCH"] = stall_match
    else:
        updated.pop("SSH_STALL_MATCH", None)
    if ignore_sigterm:
        updated["SSH_IGNORE_SIGTERM"] = "1"
    else:
        updated.pop("SSH_IGNORE_SIGTERM", None)
    if extra_env:
        updated.update(extra_env)
    return updated


def read_ssh_calls(state):
    """Return recorded invocations of the fake SSH transport."""
    log_path = state / "ssh-calls.jsonl"
    if not log_path.exists():
        return []
    return [json.loads(line) for line in log_path.read_text().splitlines()]


def process_is_alive(pid):
    """Return whether *pid* still exists in the OS process table."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def test_ssh_transport_applies_bounded_options_and_respects_user_command(
        checkout, tmp_path):
    repo, state, env = checkout("ssh-default")
    ssh_env = configure_ssh_transport(
        tmp_path,
        repo,
        state,
        env,
        extra_env={
            "DEPLOY_SSH_CONNECT_TIMEOUT": "7",
            "DEPLOY_SSH_SERVER_ALIVE_INTERVAL": "9",
            "DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX": "2",
        },
    )
    result = deploy(repo, ssh_env)
    assert result.returncode == 0, result.stdout + result.stderr
    calls = read_ssh_calls(state)
    assert calls, "expected git remote operations to invoke SSH transport"
    for call in calls:
        argv = call["argv"]
        assert "BatchMode=yes" in argv
        assert "ConnectTimeout=7" in argv
        assert "ServerAliveInterval=9" in argv
        assert "ServerAliveCountMax=2" in argv
        assert "-F" not in argv

    custom_repo, custom_state, custom_env = checkout("ssh-custom")
    custom_ssh_env = configure_ssh_transport(
        tmp_path,
        custom_repo,
        custom_state,
        custom_env,
        git_ssh_command="{ssh} -i /custom/deploy-key",
    )
    custom_result = deploy(custom_repo, custom_ssh_env, "0.0.2")
    assert custom_result.returncode == 0, (
        custom_result.stdout + custom_result.stderr
    )
    custom_calls = read_ssh_calls(custom_state)
    assert custom_calls
    for call in custom_calls:
        argv = call["argv"]
        assert argv[:2] == ["-i", "/custom/deploy-key"]
        assert "BatchMode=yes" in argv
        assert "ConnectTimeout=15" in argv
        assert "ServerAliveInterval=15" in argv
        assert "ServerAliveCountMax=3" in argv


def test_stalled_preupload_fetch_times_out_kills_child_and_releases_lock(
        checkout, tmp_path):
    # Reproduces #657: initial preflight, lock acquisition, gates, build, and
    # ownership recheck succeed, then the second `git fetch --tags origin`
    # (upload-pack #3) stalls on a hung SSH transport that ignores SIGTERM.
    repo, state, env = checkout("stalled-fetch")
    inspector, _, _ = checkout("inspector")
    ssh_env = configure_ssh_transport(
        tmp_path,
        repo,
        state,
        env,
        stall_match="upload-pack:3",
        ignore_sigterm=True,
        extra_env={"DEPLOY_GIT_REMOTE_TIMEOUT": "3"},
    )
    started = time.monotonic()
    result = deploy(repo, ssh_env)
    elapsed = time.monotonic() - started

    assert result.returncode == 124
    assert elapsed < 25.0
    assert "timed out after 3s" in result.stderr
    assert "fetch tags from origin (pre-upload recheck)" in result.stderr
    assert "build" in events(state)
    assert "upload" not in events(state)
    # Lock must be released cleanly by owner-checked cleanup once the stalled
    # fetch is terminated.
    assert remote_ref(inspector, LOCK) is None
    assert remote_ref(inspector, "refs/tags/v0.0.1") is None

    stalled_pids = [
        int(line)
        for line in (state / "stalled-pids").read_text().splitlines()
        if line.strip()
    ]
    assert len(stalled_pids) == 1
    assert not process_is_alive(stalled_pids[0])


def test_stalled_lock_release_times_out_bounded_with_recovery_diagnostic(
        checkout, tmp_path):
    # Both the pre-upload fetch (upload-pack #3) and the subsequent lock
    # release push (receive-pack #2) stall on an unresponsive transport.
    repo, state, env = checkout("stalled-cleanup")
    inspector, _, _ = checkout("inspector-cleanup")
    ssh_env = configure_ssh_transport(
        tmp_path,
        repo,
        state,
        env,
        stall_match="upload-pack:3,receive-pack:2",
        extra_env={"DEPLOY_GIT_REMOTE_TIMEOUT": "3"},
    )
    started = time.monotonic()
    result = deploy(repo, ssh_env)
    elapsed = time.monotonic() - started

    assert result.returncode == 124, result.stderr
    assert elapsed < 30.0
    assert "fetch tags from origin (pre-upload recheck)" in result.stderr, result.stderr
    assert "release lock (refs/heads/release-lock)" in result.stderr, result.stderr
    assert "could not release lock" in result.stderr
    assert "docs/testing.md#release-ownership" in result.stderr
    assert "upload" not in events(state)
    # Lock remains intact on origin for operator inspection; never stolen.
    assert remote_ref(inspector, LOCK) is not None


@pytest.mark.parametrize(
    "env_var,bad_value",
    [
        ("DEPLOY_SSH_CONNECT_TIMEOUT", "0"),
        ("DEPLOY_SSH_SERVER_ALIVE_INTERVAL", "abc"),
        ("DEPLOY_SSH_SERVER_ALIVE_COUNT_MAX", "-1"),
        ("DEPLOY_GIT_REMOTE_TIMEOUT", "0"),
        ("DEPLOY_GIT_REMOTE_TIMEOUT", "not-a-number"),
    ],
)
def test_invalid_timeout_env_vars_fail_before_remote_operations(
        checkout, env_var, bad_value):
    repo, state, env = checkout("bad-timeout")
    bad_env = dict(env, **{env_var: bad_value})
    result = deploy(repo, bad_env)
    assert result.returncode != 0
    assert env_var in result.stderr
    assert events(state) == []
    assert remote_ref(repo, LOCK) is None
