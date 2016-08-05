# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

def parse_int_list(string):
    """
    Parses a string of numbers and ranges into a list of integers. Ranges
    are separated by dashes and inclusive of both the start and end number.

    Example:
        parse_int_list("8 9 10,11-13") == [8,9,10,11,12,13]
    """
    integers = []
    for comma_part in string.split(","):
        for substring in comma_part.split(" "):
            if len(substring) == 0:
                continue
            if "-" in substring:
                left, right = substring.split("-")
                left_val = int(left.strip())
                right_val = int(right.strip())
                integers.extend(range(left_val, right_val + 1))
            else:
                integers.append(int(substring.strip()))
    return integers
