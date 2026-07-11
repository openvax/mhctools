# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging


_PACKAGE_LOGGER_NAME = "mhctools"
_PACKAGE_LOGGER = logging.getLogger(_PACKAGE_LOGGER_NAME)
if not any(isinstance(handler, logging.NullHandler)
           for handler in _PACKAGE_LOGGER.handlers):
    _PACKAGE_LOGGER.addHandler(logging.NullHandler())

_LOG_LEVELS = {
    "CRITICAL": logging.CRITICAL,
    "ERROR": logging.ERROR,
    "WARNING": logging.WARNING,
    "WARN": logging.WARNING,
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "NOTSET": logging.NOTSET,
}


def get_logger(name):
    return logging.getLogger(name)


def set_log_level(level):
    """Set the log level for mhctools loggers.

    Parameters
    ----------
    level : str or int
        Logging level such as ``"WARNING"``, ``"INFO"``, or ``logging.DEBUG``.

    Returns
    -------
    logging.Logger
        The package logger whose level was set.
    """
    if isinstance(level, str):
        level_name = level.upper()
        if level_name not in _LOG_LEVELS:
            raise ValueError("Unknown logging level %r" % level)
        level = _LOG_LEVELS[level_name]
    _PACKAGE_LOGGER.setLevel(level)
    return _PACKAGE_LOGGER
