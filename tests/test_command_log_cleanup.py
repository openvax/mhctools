"""Temporary command-log streams have an owner independent of output streams."""

import logging
import sys

import pytest

from mhctools import process_helpers


@pytest.mark.parametrize("logging_fails", [False, True])
def test_command_log_handler_is_closed_and_removed(tmp_path, monkeypatch, logging_fails):
    handlers = []
    original = logging.FileHandler

    def track(*args, **kwargs):
        handler = original(*args, **kwargs)
        handlers.append(handler)
        return handler

    monkeypatch.setattr(process_helpers.logging, "FileHandler", track)
    monkeypatch.setattr(process_helpers.logger, "debug", lambda *args: None)
    if logging_fails:
        def fail(*args):
            raise ValueError("logging failed")
        monkeypatch.setattr(process_helpers.logger, "debug", fail)
        # Keep this test about the handler's ownership; no child is left
        # running when the deliberate logging exception interrupts the helper.
        monkeypatch.setattr(process_helpers.AsyncProcess, "start", lambda self: None)
    with (tmp_path / "output.txt").open("w+") as output:
        commands = {output: [sys.executable, "-c", "print('prediction output')"]}
        if logging_fails:
            with pytest.raises(ValueError, match="logging failed"):
                process_helpers.run_multiple_commands_redirect_stdout(commands, process_limit=1)
        else:
            process_helpers.run_multiple_commands_redirect_stdout(commands, process_limit=1)
            output.seek(0)
            assert "prediction output" in output.read()
        assert not output.closed  # the caller still owns the command stream
    assert len(handlers) == 1
    assert handlers[0].stream is None
    assert handlers[0] not in process_helpers.logger.handlers
