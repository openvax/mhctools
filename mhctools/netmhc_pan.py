# Copyright (c) 2016. Mount Sinai School of Medicine
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

from .netmhc_pan28 import NetMHCpan28
from .netmhc_pan3 import NetMHCpan3


logger = logging.getLogger(__name__)


def NetMHCpan(
        alleles,
        epitope_lengths=[9],
        program_name="netMHCpan",
        max_file_records=None,
        process_limit=0,
        extra_flags=[]):
    """
    This function wraps NetMHCpan28 and NetMHCpan3 to automatically detect which class
    to use, with the help of the miraculous and strange '--version' netmhcpan argument.
    """
    # convert to str since Python3 returns a `bytes` object. The '3' here is meaningless,
    # but it is necessary to call `netmhcpan --version` with some argument, otherwise
    # it hangs.
    output = check_output([
        program_name, "--version", "3",
    ])
    output_str = output.decode("ascii", "ignore")
    if "NetMHCpan version 2.8" in output_str:
        return NetMHCpan28(
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            program_name=program_name,
            max_file_records=max_file_records,
            process_limit=process_limit,
            extra_flags=extra_flags)

    elif "NetMHCpan version 3.0" in output_str:
        return NetMHCpan3(
            alleles=alleles,
            epitope_lengths=epitope_lengths,
            program_name=program_name,
            max_file_records=max_file_records,
            process_limit=process_limit,
            extra_flags=extra_flags)

    else:
        raise SystemError("This software expects NetMHCpan version 2.8 or 3.0")
