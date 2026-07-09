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
Unit tests for supported-allele handling in BaseCommandlinePredictor:

  * `-listMHC` output is read once per (command, flag) and returned verbatim
    (no mhcgnomes normalization of the ~13k names up front);
  * user alleles are validated with a fast verbatim path that only needs the
    normalized map (the expensive mhcgnomes pass) when the fast path misses;
  * when an allele is supported under a spelling our normalizer rewrites, the
    predictor's own -listMHC spelling is used on the command line so it
    round-trips (e.g. BoLA-1:00901 instead of the rejected BoLA-109:01).

No predictor binary is required.
"""

import logging

import pytest

import mhctools.base_commandline_predictor as bcp
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.unsupported_allele import UnsupportedAllele


# Raw `-listMHC`-style output: parseable HLA, a comment, non-human names, and
# two spellings of one BoLA allele (netMHCpan really does list both).
_CANNED_LISTMHC = (
    b"# comment line\n"
    b"HLA-A02:01\n"
    b"HLA-B07:02\n"
    b"BoLA-amani.1\n"
    b"BoLA-1:00901\n"
    b"BoLA-100901\n"
)

_KEY = ("faketool_cache_test", "-listMHC")


def _clear():
    bcp._supported_alleles_cache.pop(_KEY, None)
    bcp._normalized_supported_cache.pop(_KEY, None)


@pytest.fixture(autouse=True)
def _clean_caches():
    _clear()
    yield
    _clear()


class _StubResolve(BaseCommandlinePredictor):
    """Exercises _resolve_supported_alleles / _build_command without a real
    binary: sets the attributes they read and bypasses __init__."""
    def __init__(self, alleles, supported_names, prepare=None):
        self.program_name, self.supported_alleles_flag = _KEY
        self.alleles = list(alleles)              # already normalized
        self._supported_allele_names = set(supported_names)
        self._prepare = prepare
        # attributes _build_command needs
        self.peptide_mode_flags = ["-p"]
        self.allele_flag = "-a"
        self.length_flag = "-l"
        self.input_file_flag = "-f"
        self.tempdir_flag = None
        self.extra_flags = []

    def prepare_allele_name(self, allele_name):
        if self._prepare is not None:
            return self._prepare(allele_name)
        return allele_name.replace("*", "")


# ---------------------------------------------------------------------------
# _determine_supported_alleles: raw, cached, no normalization
# ---------------------------------------------------------------------------

def test_supported_alleles_are_raw_and_cached(monkeypatch):
    calls = {"n": 0}

    def fake_check_output(args):
        calls["n"] += 1
        return _CANNED_LISTMHC

    monkeypatch.setattr(bcp, "check_output", fake_check_output)
    first = BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
    second = BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
    assert calls["n"] == 1                 # subprocess ran once
    assert first is second                 # cached
    # names are verbatim, not normalized, and non-human names are kept
    assert first == {
        "HLA-A02:01", "HLA-B07:02",
        "BoLA-amani.1", "BoLA-1:00901", "BoLA-100901"}


def test_determine_only_probes_for_a_parseable_allele(monkeypatch):
    # It must not normalize the whole list up front — only probe until the
    # first parseable name (the wrong-executable sanity check).
    calls = {"n": 0}
    real = bcp.normalize_allele_name

    def counting(name, *a, **k):
        calls["n"] += 1
        return real(name, *a, **k)

    monkeypatch.setattr(bcp, "check_output", lambda args: _CANNED_LISTMHC)
    monkeypatch.setattr(bcp, "normalize_allele_name", counting)
    BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
    # 5 names; probe stops at the first that parses, so far fewer than 5.
    assert calls["n"] < 5


def test_garbage_output_raises_system_error(monkeypatch):
    # A mismatched executable's usage/error text has no parseable alleles and
    # must raise SystemError (the "incorrect executable version?" signal),
    # not silently pass through as a set of raw names.
    monkeypatch.setattr(
        bcp, "check_output",
        lambda args: b"usage: netMHC [-options]\nnot_an_allele_here\n")
    with pytest.raises(SystemError):
        BaseCommandlinePredictor._determine_supported_alleles(*_KEY)


# ---------------------------------------------------------------------------
# Fast path: verbatim match, no normalized map built
# ---------------------------------------------------------------------------

def test_fast_path_resolves_without_building_normalized_map():
    p = _StubResolve(
        ["HLA-A*02:01", "HLA-B*07:02"],
        {"HLA-A02:01", "HLA-B07:02", "BoLA-1:00901"})
    p._resolve_supported_alleles()
    assert p._allele_cli_names == {
        "HLA-A*02:01": "HLA-A02:01", "HLA-B*07:02": "HLA-B07:02"}
    # The expensive normalization of the whole list never ran.
    assert _KEY not in bcp._normalized_supported_cache


# ---------------------------------------------------------------------------
# Slow path: round-trip to the predictor's own spelling
# ---------------------------------------------------------------------------

def test_slow_path_uses_predictor_spelling_for_non_roundtripping_allele():
    # prepare("BoLA-1*09:01") -> "BoLA-109:01", which is NOT what netMHCpan
    # prints; the normalized map must map it back to "BoLA-1:00901".
    p = _StubResolve(["BoLA-1*09:01"], {"BoLA-1:00901", "BoLA-100901"})
    p._resolve_supported_alleles()
    assert p._allele_cli_names == {"BoLA-1*09:01": "BoLA-1:00901"}
    assert _KEY in bcp._normalized_supported_cache  # map was needed


def test_collision_prefers_colon_spelling():
    mapping = BaseCommandlinePredictor._normalized_supported_alleles(
        *_KEY, raw_names={"BoLA-100901", "BoLA-1:00901"})
    assert mapping["BoLA-1*09:01"] == "BoLA-1:00901"


def test_collision_tie_break_is_deterministic():
    # Two colon spellings for one Mamu allele: pick the one whose colon sits
    # right after the gene, deterministically regardless of set order.
    raw = {"Mamu-A1001:01", "Mamu-A1:00101"}
    a = BaseCommandlinePredictor._normalized_supported_alleles(*_KEY, raw_names=raw)
    _clear()
    b = BaseCommandlinePredictor._normalized_supported_alleles(
        *_KEY, raw_names=set(reversed(list(raw))))
    assert a == b == {"Mamu-A1*01:01": "Mamu-A1:00101"}


def test_resolved_spelling_used_in_command():
    # The predictor's own -listMHC spelling (not the mangled prepared form)
    # must reach the -a argument.
    p = _StubResolve(["BoLA-1*09:01"], {"BoLA-1:00901", "BoLA-100901"})
    p._resolve_supported_alleles()
    cmd = p._build_command(
        input_filename="pep.txt", alleles=["BoLA-1*09:01"], peptide_mode=True)
    i = cmd.index("-a")
    assert cmd[i + 1] == "BoLA-1:00901"
    assert "BoLA-109:01" not in cmd            # not the rewritten form


def test_prepare_failure_falls_back_to_slow_path():
    # prepare_allele_name that raises (as netMHCIIpan can on unexpected genes)
    # must not crash resolution — it falls through to the -listMHC spelling.
    def boom(allele):
        raise ValueError("prepare rejects %s" % allele)

    p = _StubResolve(["BoLA-1*09:01"], {"BoLA-1:00901"}, prepare=boom)
    p._resolve_supported_alleles()
    assert p._allele_cli_names == {"BoLA-1*09:01": "BoLA-1:00901"}


def test_unsupported_allele_raises():
    p = _StubResolve(["HLA-A*99:99"], {"HLA-A02:01", "HLA-B07:02"})
    with pytest.raises(UnsupportedAllele, match="HLA-A"):
        p._resolve_supported_alleles()


def test_no_supported_list_skips_validation():
    # Predictors that don't enumerate alleles keep working (no CLI-name map).
    p = _StubResolve(["HLA-A*02:01"], set())
    p._supported_allele_names = None
    p._resolve_supported_alleles()
    assert p._allele_cli_names == {}
    # falls back to prepare_allele_name for the command line
    assert p._cli_allele_name("HLA-A*02:01") == "HLA-A02:01"


# ---------------------------------------------------------------------------
# Logging: no INFO spam about unparseable names
# ---------------------------------------------------------------------------

def test_no_info_spam(monkeypatch, caplog):
    monkeypatch.setattr(bcp, "check_output", lambda args: _CANNED_LISTMHC)
    with caplog.at_level(logging.INFO, logger=bcp.logger.name):
        BaseCommandlinePredictor._determine_supported_alleles(*_KEY)
        # even building the normalized map (which skips BoLA-amani.1) is quiet
        BaseCommandlinePredictor._normalized_supported_alleles(
            *_KEY,
            raw_names={"HLA-A02:01", "BoLA-amani.1"})
    assert "Skipping allele" not in caplog.text
