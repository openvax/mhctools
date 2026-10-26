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

"""Process lifecycle, concurrency limits and runtime probes."""

import errno
import os
import signal
import sys
from subprocess import Popen, TimeoutExpired
import tempfile
from threading import Event, Thread
import time
from unittest.mock import patch, MagicMock

import pytest

from .common import process_is_running
from mhctools.process_helpers import (
    AsyncProcess,
    legacy_keras_environment,
    legacy_keras_runtime_available,
    python_imports_available,
    run_multiple_commands_redirect_stdout,
)


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
@pytest.mark.parametrize("parent_ignores_term,child_ignores_term", [
    (False, False), (False, True), (True, True),
])
def test_timeout_stops_the_whole_process_group(
        tmp_path, parent_ignores_term, child_ignores_term):
    child = tmp_path / "child.py"
    child.write_text(
        "import os, signal, sys, time\nfrom pathlib import Path\n"
        "signal.signal(signal.SIGTERM, signal.SIG_IGN if sys.argv[2] == '1' "
        "else signal.SIG_DFL)\n"
        "Path(sys.argv[1]).write_text(str(os.getpid()))\ntime.sleep(60)\n")
    parent = tmp_path / "parent.py"
    parent.write_text(
        "import signal, subprocess, sys, time\n"
        "signal.signal(signal.SIGTERM, signal.SIG_IGN if sys.argv[4] == '1' "
        "else signal.SIG_DFL)\n"
        "subprocess.Popen([sys.executable, sys.argv[1], sys.argv[2], sys.argv[3]])\n"
        "time.sleep(60)\n")
    ready = tmp_path / "child.pid"
    proc = AsyncProcess(
        [sys.executable, str(parent), str(child), str(ready),
         str(int(child_ignores_term)), str(int(parent_ignores_term))],
        terminate_process_group=True)
    proc.start()
    try:
        deadline = time.monotonic() + 5
        while not ready.exists():
            assert proc.poll() is None
            assert time.monotonic() < deadline, "Child did not become ready"
            time.sleep(0.01)
        child_pid = int(ready.read_text())
        with pytest.raises(TimeoutExpired):
            proc.wait(timeout=0.05, terminate_grace_seconds=0.1)
        assert proc.process.returncode is not None
        deadline = time.monotonic() + 2
        while process_is_running(child_pid) and time.monotonic() < deadline:
            time.sleep(0.01)
        assert not process_is_running(child_pid), "Timeout left a running child"
    finally:
        # A failing regression must not leak its own sleeping child.
        try:
            os.killpg(proc.process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.process.wait()


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_group_cleanup_reaps_an_already_exited_leader():
    proc = AsyncProcess([sys.executable, "-c", "pass"], terminate_process_group=True)
    proc.start()
    try:
        deadline = time.monotonic() + 5
        # Deliberately leave the direct child unreaped, unlike Popen.poll().
        while process_is_running(proc.process.pid):
            assert time.monotonic() < deadline, "Child did not exit"
            time.sleep(0.01)
        proc._terminate(grace_seconds=0.1)
        assert proc.process.returncode == 0
    finally:
        proc.process.wait()


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_group_cleanup_waits_for_other_members_to_be_reaped(monkeypatch):
    proc = AsyncProcess(
        [sys.executable, "-c", "import time; time.sleep(30)"], terminate_process_group=True)
    # Own a group in this session so the test can parent both members and
    # control when its second, already-exited member gets reaped.
    proc.process = Popen(proc.args, process_group=0)
    member = Popen([sys.executable, "-c", "pass"], process_group=proc.process.pid)
    checked = Event()
    original = os.killpg

    def signal_group(group, signum):
        try:
            return original(group, signum)
        finally:
            if signum == 0 and proc.process.returncode is not None:
                checked.set()

    def reap_member():
        checked.wait(timeout=5)
        time.sleep(0.05)
        member.wait()

    reaper = Thread(target=reap_member)
    reaper.start()
    try:
        with monkeypatch.context() as patcher:
            patcher.setattr(os, "killpg", signal_group)
            proc._terminate(grace_seconds=1)
        assert proc.process.returncode is not None
    finally:
        checked.set()
        proc.process.kill()
        member.kill()
        proc.process.wait()
        reaper.join(timeout=5)
        assert not reaper.is_alive()


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
def test_group_cleanup_preserves_permission_errors_for_live_processes(monkeypatch):
    proc = AsyncProcess(
        [sys.executable, "-c", "import time; time.sleep(30)"], terminate_process_group=True)
    proc.start()
    denied = PermissionError("cannot signal live process")

    def deny(*args):
        raise denied

    try:
        with monkeypatch.context() as patcher:
            patcher.setattr(os, "killpg", deny)
            with pytest.raises(PermissionError) as error:
                proc._terminate(grace_seconds=0.1)
            assert error.value is denied
            assert proc.poll() is None
    finally:
        os.killpg(proc.process.pid, signal.SIGKILL)
        proc.process.wait()


def test_start_retries_on_blocking_io_error():
    """Popen should be retried on BlockingIOError (EAGAIN)."""
    proc = AsyncProcess(["echo", "hello"])
    mock_popen = MagicMock()
    # Fail twice with EAGAIN, then succeed
    mock_popen.side_effect = [
        BlockingIOError(errno.EAGAIN, "Resource temporarily unavailable"),
        BlockingIOError(errno.EAGAIN, "Resource temporarily unavailable"),
        MagicMock(),  # success
    ]
    with patch("mhctools.process_helpers.Popen", mock_popen), \
         patch("mhctools.process_helpers._POPEN_BASE_DELAY", 0.01):
        proc.start()
    assert mock_popen.call_count == 3
    assert proc.process is not None


def test_start_raises_after_max_retries():
    """After exhausting retries, the error should propagate."""
    proc = AsyncProcess(["echo", "hello"])
    mock_popen = MagicMock()
    mock_popen.side_effect = BlockingIOError(
        errno.EAGAIN, "Resource temporarily unavailable")
    with patch("mhctools.process_helpers.Popen", mock_popen), \
         patch("mhctools.process_helpers._POPEN_BASE_DELAY", 0.01), \
         patch("mhctools.process_helpers._POPEN_MAX_RETRIES", 3):
        with pytest.raises(BlockingIOError):
            proc.start()
    assert mock_popen.call_count == 3


def test_start_does_not_retry_other_os_errors():
    """Non-EAGAIN OSErrors should raise immediately."""
    proc = AsyncProcess(["nonexistent_command"])
    mock_popen = MagicMock()
    mock_popen.side_effect = FileNotFoundError("No such file")
    with patch("mhctools.process_helpers.Popen", mock_popen):
        with pytest.raises(FileNotFoundError):
            proc.start()
    assert mock_popen.call_count == 1


def test_start_succeeds_first_try():
    """Normal case — no retry needed."""
    proc = AsyncProcess(["echo", "hello"])
    mock_popen = MagicMock(return_value=MagicMock())
    with patch("mhctools.process_helpers.Popen", mock_popen):
        proc.start()
    assert mock_popen.call_count == 1


def test_run_multiple_with_low_process_limit():
    """Run several commands with process_limit=1 to exercise queuing."""
    n_commands = 5
    tmpdir = tempfile.mkdtemp()
    commands = {}
    for i in range(n_commands):
        path = os.path.join(tmpdir, "out_%d.txt" % i)
        f = open(path, "w")
        commands[f] = ["echo", "hello_%d" % i]
    run_multiple_commands_redirect_stdout(
        commands, print_commands=False, process_limit=1)
    for f in commands:
        f.close()
        with open(f.name) as r:
            content = r.read().strip()
        assert content.startswith("hello_")
        os.remove(f.name)
    os.rmdir(tmpdir)


def test_concurrency_never_exceeds_limit():
    """At no point should more than process_limit processes be running."""
    process_limit = 2
    n_commands = 6
    tmpdir = tempfile.mkdtemp()
    commands = {}
    for i in range(n_commands):
        path = os.path.join(tmpdir, "out_%d.txt" % i)
        f = open(path, "w")
        # sleep 0.1s so processes overlap
        commands[f] = ["sleep", "0.1"]

    max_concurrent = {"value": 0}
    current_concurrent = {"value": 0}
    original_popen = __import__("subprocess").Popen

    def tracking_popen(*args, **kwargs):
        current_concurrent["value"] += 1
        if current_concurrent["value"] > max_concurrent["value"]:
            max_concurrent["value"] = current_concurrent["value"]
        proc = original_popen(*args, **kwargs)
        original_wait = proc.wait

        def patched_wait(*a, **kw):
            result = original_wait(*a, **kw)
            current_concurrent["value"] -= 1
            return result
        proc.wait = patched_wait
        return proc

    with patch("mhctools.process_helpers.Popen", side_effect=tracking_popen):
        run_multiple_commands_redirect_stdout(
            commands, print_commands=False, process_limit=process_limit,
            polling_freq=0.02)

    assert max_concurrent["value"] <= process_limit, (
        "Had %d concurrent processes, limit was %d"
        % (max_concurrent["value"], process_limit))
    for f in commands:
        f.close()
        os.remove(f.name)
    os.rmdir(tmpdir)


def test_run_multiple_with_eagain_and_low_limit():
    """Simulate EAGAIN during a multi-command run with process_limit=2."""
    n_commands = 4
    tmpdir = tempfile.mkdtemp()
    commands = {}
    for i in range(n_commands):
        path = os.path.join(tmpdir, "out_%d.txt" % i)
        f = open(path, "w")
        commands[f] = ["echo", "hello_%d" % i]

    original_popen = __import__("subprocess").Popen
    call_count = {"n": 0}

    def flaky_popen(*args, **kwargs):
        call_count["n"] += 1
        # Fail on the 3rd call, then succeed on retry
        if call_count["n"] == 3:
            raise BlockingIOError(
                errno.EAGAIN, "Resource temporarily unavailable")
        return original_popen(*args, **kwargs)

    with patch("mhctools.process_helpers.Popen", side_effect=flaky_popen), \
         patch("mhctools.process_helpers._POPEN_BASE_DELAY", 0.01):
        run_multiple_commands_redirect_stdout(
            commands, print_commands=False, process_limit=2)
    # The 3rd call failed, so we expect n_commands + 1 total Popen calls
    assert call_count["n"] == n_commands + 1
    for f in commands:
        f.close()
        os.remove(f.name)
    os.rmdir(tmpdir)


def _probe_script():
    """The exact source legacy_keras_runtime_available sends to the interpreter."""
    from mhctools.process_helpers import _LEGACY_KERAS_PROBE
    return _LEGACY_KERAS_PROBE


def test_legacy_keras_runtime_missing_interpreter():
    assert not legacy_keras_runtime_available("/nonexistent/python")


def test_legacy_keras_runtime_interpreter_is_a_directory():
    assert not legacy_keras_runtime_available(tempfile.gettempdir())


def test_legacy_keras_runtime_timeout():
    # sys.executable is a real interpreter, so only the timeout can end this.
    assert not legacy_keras_runtime_available(sys.executable, timeout=0.001)


def test_legacy_keras_runtime_reports_probe_exit_status():
    # Stand in for the interpreter so the outcome is decided by the exit
    # status alone, with no TensorFlow installed anywhere in the test env.
    for exit_status, expected in ((0, True), (1, False)):
        with patch("mhctools.process_helpers.subprocess.run") as run:
            run.return_value = MagicMock(returncode=exit_status)
            assert legacy_keras_runtime_available("python") is expected


def test_legacy_keras_probe_rejects_keras_3_fallback():
    # TensorFlow >= 2.16 resolves `import tensorflow.keras` to Keras 3 when
    # tf-keras is absent, so the probe must import tf_keras itself instead.
    script = _probe_script()
    assert "import tf_keras" in script
    assert "(2, 16)" in script


def test_legacy_keras_runtime_sets_legacy_keras_env(monkeypatch):
    monkeypatch.delenv("TF_USE_LEGACY_KERAS", raising=False)
    monkeypatch.delenv("TF_CPP_MIN_LOG_LEVEL", raising=False)
    with patch("mhctools.process_helpers.subprocess.run") as run:
        run.return_value = MagicMock(returncode=0)
        legacy_keras_runtime_available("python")
    env = run.call_args.kwargs["env"]
    assert env["TF_USE_LEGACY_KERAS"] == "1"
    assert env["TF_CPP_MIN_LOG_LEVEL"] == "3"


def test_legacy_keras_runtime_preserves_explicit_environment(monkeypatch):
    monkeypatch.setenv("TF_USE_LEGACY_KERAS", "0")
    monkeypatch.setenv("TF_CPP_MIN_LOG_LEVEL", "1")
    with patch("mhctools.process_helpers.subprocess.run") as run:
        run.return_value = MagicMock(returncode=0)
        legacy_keras_runtime_available("python")
    env = run.call_args.kwargs["env"]
    assert env["TF_USE_LEGACY_KERAS"] == "0"
    assert env["TF_CPP_MIN_LOG_LEVEL"] == "1"


def test_legacy_keras_environment_does_not_modify_host_or_other_calls(monkeypatch):
    monkeypatch.delenv("TF_USE_LEGACY_KERAS", raising=False)
    monkeypatch.delenv("TF_CPP_MIN_LOG_LEVEL", raising=False)
    monkeypatch.setenv("MHCTOOLS_TEST_ENVIRONMENT", "preserved")
    environment = legacy_keras_environment()
    assert environment["MHCTOOLS_TEST_ENVIRONMENT"] == "preserved"
    assert "TF_USE_LEGACY_KERAS" not in os.environ
    assert "TF_CPP_MIN_LOG_LEVEL" not in os.environ
    environment["TF_USE_LEGACY_KERAS"] = "0"
    assert legacy_keras_environment()["TF_USE_LEGACY_KERAS"] == "1"


def test_python_imports_available_stdlib():
    # sys and os import everywhere, so a True here exercises the happy path
    # without depending on any optional scientific package.
    assert python_imports_available(sys.executable, ["sys", "os"])


def test_python_imports_available_missing_module():
    assert not python_imports_available(
        sys.executable, ["mhctools_definitely_not_a_real_module"])


def test_python_imports_available_requires_every_module():
    # One missing module fails the whole probe, so a partially provisioned
    # interpreter is not reported as runnable.
    assert not python_imports_available(
        sys.executable, ["sys", "mhctools_definitely_not_a_real_module"])


def test_python_imports_available_missing_interpreter():
    assert not python_imports_available("/nonexistent/python", ["sys"])


def test_python_imports_available_timeout():
    # Interpreter startup can beat a tiny timeout; exercise the timeout path
    # deterministically rather than depending on machine scheduling (#556).
    with patch("mhctools.process_helpers.subprocess.run",
               side_effect=TimeoutExpired("python", 7)) as run:
        assert not python_imports_available(sys.executable, ["sys"], timeout=7)
    assert run.call_args.kwargs["timeout"] == 7


def test_python_imports_available_empty_is_trivially_true():
    assert python_imports_available(sys.executable, [])


def test_python_imports_available_rejects_a_bare_string():
    # "torch" would otherwise become `import t / import o / ...` and quietly
    # report False, turning a caller's mistake into a permanent skip.
    with pytest.raises(TypeError):
        python_imports_available(sys.executable, "sys")


def test_python_imports_available_uses_the_supplied_environment():
    # The probe must see what the sidecar sees; PYTHONNOUSERSITE is the case
    # that matters, so check the env is passed through rather than dropped.
    with patch("mhctools.process_helpers.subprocess.run") as run:
        run.return_value = MagicMock(returncode=0)
        python_imports_available(
            "python", ["sys"], env={"PYTHONNOUSERSITE": "1"})
    assert run.call_args.kwargs["env"] == {"PYTHONNOUSERSITE": "1"}
