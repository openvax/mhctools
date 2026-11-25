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

"""Readable argparse help with optional terminal styling via ``datacache.cli_help``."""

import os
import sys

from datacache.cli_help import (
    ANSI_BOLD,
    ANSI_BOLD_CYAN,
    ANSI_CYAN,
    ANSI_DIM,
    ANSI_RESET,
    ANSI_YELLOW,
    DIM_STATUS_VALUES,
    HelpArgumentParser as _DatacacheHelpArgumentParser,
    HelpFormatter,
    OPTION_PATTERN,
    WARNING_STATUS_VALUES,
    classify_status,
    color_enabled as _datacache_color_enabled,
    enable_argcomplete,
    style_help as _datacache_style_help,
    style_status,
    style_table_header,
)

_OPTION = OPTION_PATTERN

__all__ = [
    "ANSI_BOLD",
    "ANSI_BOLD_CYAN",
    "ANSI_CYAN",
    "ANSI_DIM",
    "ANSI_RESET",
    "ANSI_YELLOW",
    "DIM_STATUS_VALUES",
    "HelpArgumentParser",
    "HelpFormatter",
    "OPTION_PATTERN",
    "WARNING_STATUS_VALUES",
    "classify_status",
    "color_enabled",
    "enable_argcomplete",
    "style_help",
    "style_status",
    "style_table_header",
]


def color_enabled(stream):
    """Return whether the destination supports optional terminal styling."""
    return "NO_COLOR" not in os.environ and _datacache_color_enabled(stream)


def style_help(text, stream):
    """Style headings and flags only when writing to a color-capable terminal."""
    if not color_enabled(stream):
        return text
    return _datacache_style_help(text, stream)


class HelpArgumentParser(_DatacacheHelpArgumentParser):
    """Keep format_help() plain; add color only when printing to a terminal."""

    def print_help(self, file=None):
        stream = sys.stdout if file is None else file
        self._print_message(style_help(self.format_help(), stream), stream)
