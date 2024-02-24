
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

from mhctools.cli.parsing_helpers import parse_int_list
from .common import eq_

def test_parse_int_list():
    # int by itself
    eq_(parse_int_list("1"), [1])
    # range of integers
    eq_(parse_int_list("1-3"), [1, 2, 3])
    # comma separated
    eq_(parse_int_list("9,10"), [9, 10])
    # comma separated with spaces
    eq_(parse_int_list("9, 10"), [9, 10])
    # comma separated and range together
    eq_(parse_int_list("9,10,12-13"), [9, 10, 12, 13])
