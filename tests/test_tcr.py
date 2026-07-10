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

"""Unit tests for the TCR input type (no NetTCR install required)."""

import json

import pandas as pd

from mhctools import TCR


def _tcr():
    return TCR(
        cdr1a="NSAFQY", cdr2a="TYSSGN", cdr3a="AMSGDGGSQGNLI",
        cdr1b="LNHDA", cdr2b="SQIVND", cdr3b="ASSIRAAYEQY",
        name="clone1")


def test_short_aliases_map_to_cdrs():
    t = _tcr()
    assert t.a1 == t.cdr1a
    assert t.a2 == t.cdr2a
    assert t.a3 == t.cdr3a
    assert t.b1 == t.cdr1b
    assert t.b2 == t.cdr2b
    assert t.b3 == t.cdr3b


def test_cdr_dict_keys_and_values():
    t = _tcr()
    assert t.cdr_dict() == {
        "a1": "NSAFQY", "a2": "TYSSGN", "a3": "AMSGDGGSQGNLI",
        "b1": "LNHDA", "b2": "SQIVND", "b3": "ASSIRAAYEQY"}


def test_identifier_uses_name_when_present():
    assert _tcr().identifier == "clone1"


def test_identifier_falls_back_to_cdr3_pair():
    t = TCR(cdr3a="AAA", cdr3b="BBB")
    assert t.identifier == "AAA/BBB"


def test_frozen():
    t = _tcr()
    try:
        t.cdr3b = "X"
        assert False, "TCR should be frozen"
    except AttributeError:
        pass


def test_to_dict_round_trip():
    t = _tcr()
    t2 = TCR.from_dict(t.to_dict())
    assert t == t2


def test_to_dict_json_serializable():
    t = _tcr()
    t2 = TCR.from_dict(json.loads(json.dumps(t.to_dict())))
    assert t == t2


def test_from_dict_ignores_unknown_keys():
    t = TCR.from_dict({"cdr3b": "ASSF", "bogus": 1})
    assert t.cdr3b == "ASSF"


def test_from_dict_accepts_short_aliases_case_insensitive():
    t = TCR.from_dict({
        "A1": "NSAFQY",
        "a2": "TYSSGN",
        "A3": "AMSGDGGSQGNLI",
        "b1": "LNHDA",
        "B2": "SQIVND",
        "b3": "ASSIRAAYEQY",
    })
    assert t.cdr1a == "NSAFQY"
    assert t.cdr2a == "TYSSGN"
    assert t.cdr3a == "AMSGDGGSQGNLI"
    assert t.cdr1b == "LNHDA"
    assert t.cdr2b == "SQIVND"
    assert t.cdr3b == "ASSIRAAYEQY"


def test_from_dict_canonical_keys_win_over_aliases():
    t = TCR.from_dict({"A1": "ALIAS", "cdr1a": "CANONICAL"})
    assert t.cdr1a == "CANONICAL"


def test_from_series_accepts_aliases():
    row = pd.Series({"A3": "CAVR", "B3": "CASS", "name": "clone1"})
    t = TCR.from_series(row)
    assert t.cdr3a == "CAVR"
    assert t.cdr3b == "CASS"
    assert t.identifier == "clone1"
