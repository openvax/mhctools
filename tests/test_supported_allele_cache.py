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

"""
Unit tests for the supported-allele determination in
BaseCommandlinePredictor: it should run the `-listMHC` subprocess and parse
its alleles once per (command, flag) rather than on every predictor
construction, and it should not spam the user at INFO level about alleles it
cannot parse. No predictor binary is required.
"""

import logging

import mhctools.base_commandline_predictor as bcp
from mhctools.base_commandline_predictor import BaseCommandlinePredictor

# A tiny stand-in for `netMHCpan -listMHC` output: two parseable HLA alleles,
# a comment line, and a name the legacy normalizer can't parse.
_CANNED_LISTMHC = (
    b"# comment line\n"
    b"HLA-A02:01\n"
    b"HLA-B07:02\n"
    b"BoLA-amani.1\n"
)

_KEY = ("faketool_cache_test", "-listMHC")


def _clear():
    bcp._supported_alleles_cache.pop(_KEY, None)


def test_supported_alleles_determined_once_and_cached(monkeypatch):
    calls = {"n": 0}

    def fake_check_output(args):
        calls["n"] += 1
        return _CANNED_LISTMHC

    monkeypatch.setattr(bcp, "check_output", fake_check_output)
    _clear()
    try:
        first = BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
        second = BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
        # The subprocess ran once; the second call hit the cache.
        assert calls["n"] == 1
        assert first is second
        assert "HLA-A*02:01" in first
        assert "HLA-B*07:02" in first
        # The unparseable name is skipped, not included.
        assert not any("BoLA" in a for a in first)
    finally:
        _clear()


def test_unparseable_alleles_not_logged_at_info(monkeypatch, caplog):
    monkeypatch.setattr(bcp, "check_output", lambda args: _CANNED_LISTMHC)
    _clear()
    try:
        with caplog.at_level(logging.INFO, logger=bcp.logger.name):
            BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
        # No per-allele "Skipping allele" spam at INFO — it is emitted at DEBUG.
        assert "Skipping allele" not in caplog.text
    finally:
        _clear()


def test_unparseable_alleles_available_at_debug(monkeypatch, caplog):
    monkeypatch.setattr(bcp, "check_output", lambda args: _CANNED_LISTMHC)
    _clear()
    try:
        with caplog.at_level(logging.DEBUG, logger=bcp.logger.name):
            BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
        # The detail is still there for anyone who opts into DEBUG.
        assert "BoLA-amani.1" in caplog.text
    finally:
        _clear()
