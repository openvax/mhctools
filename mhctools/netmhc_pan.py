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

import logging
import os
import re
import subprocess

from .netmhc_pan28 import NetMHCpan28
from .netmhc_pan3 import NetMHCpan3
from .netmhc_pan4 import NetMHCpan4
from .netmhc_pan41 import NetMHCpan41
from .netmhc_pan42 import NetMHCpan42

logger = logging.getLogger(__name__)

# Maps (major, minor) version tuple to the predictor class.
# Checked in order; the last entry is the fallback for any 4.x not
# explicitly listed (future-proofing).
_VERSION_MAP = [
    ((2, 8), NetMHCpan28),
    ((3, 0), NetMHCpan3),
    ((4, 0), NetMHCpan4),
    ((4, 1), NetMHCpan41),
    ((4, 2), NetMHCpan42),
]


def _parse_version(version_str):
    """Parse '4.2c' into (4, 2). Returns None on failure."""
    # strip trailing letter suffixes like 'b', 'c'
    cleaned = re.sub(r'[a-zA-Z]+$', '', version_str)
    parts = cleaned.split('.')
    try:
        return (int(parts[0]), int(parts[1]))
    except (IndexError, ValueError):
        return None


def NetMHCpan(
        alleles,
        program_name="netMHCpan",
        process_limit=-1,
        default_peptide_lengths=[9],
        extra_flags=[]):
    """
    Auto-detecting wrapper for any installed version of NetMHCpan.

    Runs ``netMHCpan --version`` to detect the installed version and returns
    the appropriate predictor class. For unrecognized versions >= 4.1,
    falls back to the latest known class (which uses the header-driven
    auto-detecting parser).
    """
    # Pass /dev/null as the input file so netMHCpan doesn't try to
    # open a nonexistent file (it hangs without any argument).
    # Use run() instead of check_output() because netMHCpan may exit
    # non-zero even though the version string is still in stdout.
    result = subprocess.run(
        [program_name, "--version", os.devnull],
        capture_output=True)
    output_str = result.stdout.decode("ascii", "ignore")

    match = re.search(r'# NetMHCpan version (\S+)', output_str)
    version_str = match.group(1) if match else ""
    version_tuple = _parse_version(version_str) if version_str else None

    common_kwargs = {
        "alleles": alleles,
        "default_peptide_lengths": default_peptide_lengths,
        "program_name": program_name,
        "process_limit": process_limit,
        "extra_flags": extra_flags,
    }

    # Exact match
    if version_tuple:
        for (major_minor, cls) in _VERSION_MAP:
            if version_tuple == major_minor:
                logger.info("Detected NetMHCpan %s, using %s", version_str, cls.__name__)
                return cls(**common_kwargs)

    # Fallback: use the latest known class (header-driven parser handles
    # any output format with Score_EL / Score_BA columns)
    fallback_cls = _VERSION_MAP[-1][1]
    logger.warning(
        "NetMHCpan version %s not explicitly supported, falling back to %s "
        "(header-driven auto-detecting parser)",
        version_str or "unknown", fallback_cls.__name__)
    return fallback_cls(**common_kwargs)
