# Copyright (c) 2016-2017. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import
import logging
from subprocess import check_output
import os

from .netmhc_pan28 import NetMHCpan28
from .netmhc_pan3 import NetMHCpan3
from .netmhc_pan4 import NetMHCpan4


logger = logging.getLogger(__name__)


def NetMHCpan(
        alleles,
        program_name="netMHCpan",
        process_limit=-1,
        default_peptide_lengths=[9],
        extra_flags=[]):
    """
    This function wraps NetMHCpan28 and NetMHCpan3 to automatically detect which class
    to use, with the help of the miraculous and strange '--version' netmhcpan argument.
    """
    # convert to str since Python3 returns a `bytes` object.
    # The '_MHCTOOLS_VERSION_SNIFFING' here is meaningless, but it is necessary
    # to call `netmhcpan --version` with some argument, otherwise it hangs.
    with open(os.devnull, 'w') as devnull:
        output = check_output([
            program_name, "--version", "_MHCTOOLS_VERSION_SNIFFING"],
            stderr=devnull)
    output_str = output.decode("ascii", "ignore")
    common_kwargs = {
        "alleles": alleles,
        "default_peptide_lengths": default_peptide_lengths,
        "program_name": program_name,
        "process_limit": process_limit,
        "extra_flags": extra_flags,
    }
    if "NetMHCpan version 2.8" in output_str:
        return NetMHCpan28(**common_kwargs)

    elif "NetMHCpan version 3.0" in output_str:
        return NetMHCpan3(**common_kwargs)

    elif "NetMHCpan version 4.0" in output_str:
        return NetMHCpan4(**common_kwargs)

    else:
        raise RuntimeError(
            "This software expects NetMHCpan version 2.8, 3.0, or 4.0")
