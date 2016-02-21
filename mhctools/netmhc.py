# Copyright (c) 2015. Mount Sinai School of Medicine
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
from .netmhc3 import NetMHC3
from .netmhc4 import NetMHC4

def NetMHC(alleles,
           epitope_lengths=[9],
           program_name="netMHC",
           max_file_records=None):
    """
    This function wraps NetMHC3 and NetMHC4 to automatically detect which class
    to use. This is currently based on the results of supported_alleles_flag,
    and assumes that the flag for 3.x causes an error in 4.x and vice versa.

    TODO: Make this more robust.
    """
    netmhc_classes = [NetMHC3, NetMHC4]
    successes = []
    for netmhc_class in netmhc_classes:
        try:
            successes.append(netmhc_class(
                alleles=alleles,
                epitope_lengths=epitope_lengths,
                program_name=program_name,
                max_file_records=max_file_records))
        except SystemError:
            pass

    if len(successes) > 1:
        raise SystemError("Command %s is valid for multiple NetMHC versions. "
                          "This is likely an mhctools bug." % program_name)
    if len(successes) == 0:
        raise SystemError("Command %s is not a valid way of calling any NetMHC software."
                          % program_name)

    return successes[0]
