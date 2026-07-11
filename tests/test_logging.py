import logging
import logging.config

import pytest

import mhctools
from mhctools import logging as mhctools_logging


@pytest.fixture(autouse=True)
def restore_mhctools_logger():
    logger = logging.getLogger("mhctools")
    old_level = logger.level
    old_handlers = list(logger.handlers)
    old_propagate = logger.propagate
    yield
    logger.setLevel(old_level)
    logger.handlers[:] = old_handlers
    logger.propagate = old_propagate


def test_get_logger_does_not_call_fileconfig(monkeypatch):
    def fail_fileconfig(*args, **kwargs):
        raise AssertionError("fileConfig must not be called by mhctools")

    monkeypatch.setattr(logging.config, "fileConfig", fail_fileconfig)

    logger = mhctools_logging.get_logger("mhctools.test")

    assert logger.name == "mhctools.test"


def test_get_logger_does_not_reconfigure_root_logger():
    root = logging.getLogger()
    old_level = root.level
    old_handlers = list(root.handlers)

    mhctools_logging.get_logger("mhctools.test")

    assert root.level == old_level
    assert list(root.handlers) == old_handlers


def test_package_logger_has_null_handler_by_default():
    handlers = logging.getLogger("mhctools").handlers

    assert any(isinstance(handler, logging.NullHandler)
               for handler in handlers)


def test_info_logs_are_quiet_by_default(capsys):
    logger = mhctools_logging.get_logger("mhctools.test.quiet")

    logger.info("should not print")

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_set_log_level_accepts_string_levels():
    logger = mhctools.set_log_level("WARNING")

    assert logger is logging.getLogger("mhctools")
    assert logger.level == logging.WARNING


def test_set_log_level_accepts_numeric_levels():
    logger = mhctools_logging.set_log_level(logging.ERROR)

    assert logger.level == logging.ERROR


def test_set_log_level_rejects_unknown_level():
    with pytest.raises(ValueError, match="NOTALEVEL"):
        mhctools.set_log_level("NOTALEVEL")
