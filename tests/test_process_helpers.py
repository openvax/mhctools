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

"""Tests for process_helpers retry logic."""

import errno
import os
import tempfile
from unittest.mock import patch, MagicMock

import pytest

from mhctools.process_helpers import (
    AsyncProcess,
    run_multiple_commands_redirect_stdout,
)


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
