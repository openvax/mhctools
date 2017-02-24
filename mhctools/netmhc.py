# Copyright (c) 2015-2017. Mount Sinai School of Medicine
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

from subprocess import check_output
import os

from .netmhc3 import NetMHC3
from .netmhc4 import NetMHC4

def NetMHC(alleles,
           default_peptide_lengths=[9],
           program_name="netMHC"):
    """
    This function wraps NetMHC3 and NetMHC4 to automatically detect which class
    to use. Currently based on running the '-h' command and looking for
    discriminating substrings between the versions.
    """
    # run NetMHC's help command and parse discriminating substrings out of
    # the resulting str output
    with open(os.devnull, 'w') as devnull:
        help_output = check_output([program_name, "-h"], stderr=devnull)
    help_output_str = help_output.decode("ascii", "ignore")

    substring_to_netmhc_class = {
        "-listMHC": NetMHC4,
        "--Alleles": NetMHC3,
    }

    successes = []

    for substring, netmhc_class in substring_to_netmhc_class.items():
        if substring in help_output_str:
            successes.append(netmhc_class)

    if len(successes) > 1:
        raise SystemError("Command %s is valid for multiple NetMHC versions. "
                          "This is likely an mhctools bug." % program_name)
    if len(successes) == 0:
        raise SystemError("Command %s is not a valid way of calling any NetMHC software."
                          % program_name)

    netmhc_class = successes[0]
    return netmhc_class(
        alleles=alleles,
        default_peptide_lengths=default_peptide_lengths,
        program_name=program_name)
