"""A failed command batch must not outlive its caller-owned output files."""

from contextlib import ExitStack
import os
import signal
import subprocess
import sys
import time

import pytest

from mhctools import process_helpers
from .common import process_is_running


@pytest.fixture
def batch_processes(monkeypatch, tmp_path):
    """Track real children and clean up even when testing the broken behavior."""
    processes = []
    original = process_helpers.Popen

    def start(*args, **kwargs):
        process = original(*args, **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr(process_helpers, "Popen", start)
    yield processes
    for process in processes:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except (ProcessLookupError, PermissionError):
            process.kill()
        process.wait()
    ready = tmp_path / "child.pid"
    if ready.exists() and ready.read_text():
        try:
            os.kill(int(ready.read_text()), signal.SIGKILL)
        except ProcessLookupError:
            pass


def _wait_for_file(path):
    deadline = time.monotonic() + 10
    while not path.exists() or not path.read_text():
        assert time.monotonic() < deadline, "Child did not become ready"
        time.sleep(0.01)


def _sleeper(tmp_path):
    ready = tmp_path / "child.pid"
    child = tmp_path / "child.py"
    child.write_text(
        "import os, signal, sys, time\nfrom pathlib import Path\n"
        "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "Path(sys.argv[1]).write_text(str(os.getpid()))\ntime.sleep(30)\n")
    parent = tmp_path / "parent.py"
    parent.write_text(
        "import subprocess, sys\n"
        "subprocess.run([sys.executable, sys.argv[1], sys.argv[2]])\n")
    return [sys.executable, str(parent), str(child), str(ready)], ready


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
@pytest.mark.parametrize("failure_first", [True, False])
@pytest.mark.parametrize("queued_command", [False, True])
def test_failed_batch_stops_siblings_and_descendants(
        tmp_path, batch_processes, failure_first, queued_command):
    sleeper, ready = _sleeper(tmp_path)
    failure = tmp_path / "fail.py"
    failure.write_text(
        "import sys, time\nfrom pathlib import Path\n"
        "deadline = time.monotonic() + 10\n"
        "while not Path(sys.argv[1]).exists():\n"
        "    if time.monotonic() > deadline: sys.exit(8)\n"
        "    time.sleep(.01)\n"
        "print('original failure', flush=True)\nsys.exit(7)\n")
    failing = [sys.executable, str(failure), str(ready)]
    arguments = [failing, sleeper] if failure_first else [sleeper, failing]
    not_started = tmp_path / "not-started"
    if queued_command:
        arguments.append([sys.executable, "-c",
                          "from pathlib import Path; Path(%r).touch()" % str(not_started)])
    with ExitStack() as stack:
        commands = {stack.enter_context((tmp_path / f"out-{i}").open("w+")): args
                    for i, args in enumerate(arguments)}
        started = time.monotonic()
        with pytest.raises(subprocess.CalledProcessError) as error:
            process_helpers.run_multiple_commands_redirect_stdout(
                commands, process_limit=2, polling_freq=0.01)
        assert error.value.returncode == 7
        assert time.monotonic() - started < 15, "Waited for a sleeping sibling to finish"
        assert not not_started.exists(), "Started queued work after observing failure"
        assert all(process.returncode is not None for process in batch_processes)
        assert not process_is_running(int(ready.read_text()))
        for output in commands:
            assert not output.closed
        failed_output = next(output for output, args in commands.items() if args == failing)
        failed_output.seek(0)
        assert "original failure" in failed_output.read()


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group behavior")
@pytest.mark.parametrize("failure", [FileNotFoundError("cannot start"), KeyboardInterrupt()])
def test_startup_failure_or_interrupt_stops_started_children(
        tmp_path, monkeypatch, batch_processes, failure):
    sleeper, ready = _sleeper(tmp_path)
    original = process_helpers.AsyncProcess.start

    def start(process):
        if process.args == ["second-command"]:
            _wait_for_file(ready)
            raise failure
        return original(process)

    monkeypatch.setattr(process_helpers.AsyncProcess, "start", start)
    with (tmp_path / "first").open("w") as first, (tmp_path / "second").open("w") as second:
        with pytest.raises(type(failure)) as error:
            process_helpers.run_multiple_commands_redirect_stdout(
                {first: sleeper, second: ["second-command"]},
                process_limit=2)
        assert error.value is failure
        assert all(process.returncode is not None for process in batch_processes)
        assert not process_is_running(int(ready.read_text()))


def test_logging_failure_stops_the_just_started_command(tmp_path, monkeypatch, batch_processes):
    def fail(*args):
        raise ValueError("command logging failed")

    monkeypatch.setattr(process_helpers.logger, "debug", fail)
    with (tmp_path / "output").open("w") as output:
        with pytest.raises(ValueError, match="command logging failed"):
            process_helpers.run_multiple_commands_redirect_stdout(
                {output: [sys.executable, "-c", "import time; time.sleep(30)"]},
                process_limit=1)
        assert all(process.returncode is not None for process in batch_processes)
        assert not output.closed


def test_cleanup_error_preserves_original_and_continues(tmp_path, monkeypatch, batch_processes):
    failure = ValueError("command logging failed")

    def fail_after_second_start(*args):
        if len(batch_processes) == 2:
            raise failure

    original = process_helpers.AsyncProcess._terminate
    cleaned = []

    def terminate(process, grace_seconds):
        original(process, grace_seconds)
        cleaned.append(process)
        if len(cleaned) == 1:
            raise OSError("cleanup diagnostic")

    monkeypatch.setattr(process_helpers.logger, "debug", fail_after_second_start)
    monkeypatch.setattr(process_helpers.AsyncProcess, "_terminate", terminate)
    with (tmp_path / "first").open("w") as first, (tmp_path / "second").open("w") as second:
        command = [sys.executable, "-c", "import time; time.sleep(30)"]
        with pytest.raises(ValueError) as error:
            process_helpers.run_multiple_commands_redirect_stdout(
                {first: command, second: command}, process_limit=2)
        assert error.value is failure
        assert any("cleanup diagnostic" in note for note in error.value.__notes__)
        assert len(cleaned) == 2
        assert all(process.returncode is not None for process in batch_processes)
