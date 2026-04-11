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
from unittest.mock import patch, MagicMock

import pytest

from mhctools.process_helpers import AsyncProcess


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
