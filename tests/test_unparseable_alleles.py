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
Unit tests for supporting alleles that netMHCpan lists but mhcgnomes can't
parse (issue #220): exotic non-human alleles (H-2-Qa1, BoLA-amani.1) and HLA
low/null-expression variants (HLA-A30:14L).

The key invariant is that a requested allele and the same allele echoed back
in a predictor's output resolve to a single identity, so `_check_results`
matches them. Both go through `normalize_allele_name_or_raw`; when parsing
fails it strips whitespace and the '*' gene/allele separator (which netMHCpan
adds or omits) so e.g. requested "HLA-A30:14L" and echoed "HLA-A*30:14L" agree.

No predictor binary is required.
"""

import pytest

from mhctools.allele_normalization import (
    normalize_allele_name_or_raw,
    AlleleParseError,
)
from mhctools.base_predictor import BasePredictor
from mhctools.parsing import parse_netmhc4_stdout


# ---------------------------------------------------------------------------
# normalize_allele_name_or_raw: parse when possible, canonical raw fallback
# ---------------------------------------------------------------------------

def test_parseable_allele_is_normalized():
    # Case-insensitive, and the canonical form keeps the '*'.
    assert normalize_allele_name_or_raw("hla-a*02:01") == "HLA-A*02:01"
    assert normalize_allele_name_or_raw("HLA-A02:01") == "HLA-A*02:01"


def test_unparseable_non_human_kept_verbatim():
    # These netMHCpan names don't parse at all; keep the tool's own spelling
    # (including original case, which the output echoes).
    for name in ["H-2-Qa1", "H-2-Qa2", "BoLA-amani.1", "BoLA-JSP.1",
                 "BoLA-T2a", "BoLA-gb1.7"]:
        assert normalize_allele_name_or_raw(name) == name


def test_whitespace_is_stripped():
    assert normalize_allele_name_or_raw("  H-2-Qa1 \n") == "H-2-Qa1"


def test_star_stripped_from_unparseable_fallback():
    # netMHCpan prints "HLA-A*30:14L" for a requested "HLA-A30:14L"; both
    # spellings must collapse to one identity so results match.
    bare = normalize_allele_name_or_raw("HLA-A30:14L")
    starred = normalize_allele_name_or_raw("HLA-A*30:14L")
    assert bare == starred == "HLA-A30:14L"


def test_request_and_echo_share_identity():
    # The invariant _check_results relies on: whatever the user typed and
    # whatever netMHCpan echoes back reduce to the same string.
    for requested, echoed in [
            ("HLA-A30:14L", "HLA-A*30:14L"),   # star added by netMHCpan
            ("H-2-Qa1", "H-2-Qa1"),            # verbatim
            ("BoLA-amani.1", "BoLA-amani.1")]:
        assert (normalize_allele_name_or_raw(requested) ==
                normalize_allele_name_or_raw(echoed))


# ---------------------------------------------------------------------------
# _check_hla_alleles: keep_unparseable gate
# ---------------------------------------------------------------------------

def test_check_hla_alleles_raises_by_default():
    # Default behavior (used by IEDB / in-process predictors) must still
    # reject names it can't normalize.
    with pytest.raises(AlleleParseError):
        BasePredictor._check_hla_alleles(["H-2-Qa1"])


def test_check_hla_alleles_keeps_unparseable_when_requested():
    result = set(BasePredictor._check_hla_alleles(
        ["H-2-Qa1", "BoLA-amani.1"], keep_unparseable=True))
    assert result == {"H-2-Qa1", "BoLA-amani.1"}


def test_check_hla_alleles_still_normalizes_parseable_when_keeping():
    # Mixing parseable and unparseable: parseable ones are still canonicalized.
    result = set(BasePredictor._check_hla_alleles(
        ["hla-a*02:01", "H-2-Qa1"], keep_unparseable=True))
    assert result == {"HLA-A*02:01", "H-2-Qa1"}


def test_check_hla_alleles_dedupes_star_variants_of_unparseable():
    # Two spellings of the same un-normalizable allele collapse to one entry.
    result = BasePredictor._check_hla_alleles(
        ["HLA-A30:14L", "HLA-A*30:14L"], keep_unparseable=True)
    assert result == ["HLA-A30:14L"]


# ---------------------------------------------------------------------------
# CLI --mhc-alleles must keep exotic alleles too, not crash on them (#220).
# ---------------------------------------------------------------------------

def test_cli_mhc_alleles_keep_exotic():
    from mhctools.cli.args import make_mhc_arg_parser, mhc_alleles_from_args
    args = make_mhc_arg_parser().parse_args([
        "--mhc-predictor", "netmhcpan",
        "--mhc-alleles", "BoLA-amani.1,HLA-A*02:01 H-2-Qa1"])
    # Parseable allele is canonicalized; exotic ones are kept verbatim rather
    # than raising an AlleleParseError at the CLI layer.
    assert mhc_alleles_from_args(args) == [
        "BoLA-amani.1", "HLA-A*02:01", "H-2-Qa1"]


# ---------------------------------------------------------------------------
# Output parser falls back to the raw name instead of raising (binary-free).
# ---------------------------------------------------------------------------

def test_parser_keeps_unparseable_output_allele():
    # A netMHCpan/netMHC-style table whose allele column is an exotic name
    # mhcgnomes can't parse. Before the fallback this raised in the parser;
    # now the raw name is carried through verbatim.
    stdout = """
# NetMHC version 4.0
-----------------------------------------------------------------------------------
  pos          HLA      peptide         Core Offset  I_pos  I_len  D_pos  D_len        iCore        Identity 1-log50k(aff) Affinity(nM)    %Rank  BindLevel
-----------------------------------------------------------------------------------
    0      H-2-Qa1    SIINFEKLL    SIINFEKLL      0      0      0      0      0    SIINFEKLL         SEQ_A           0.349      1147.39     4.50
-----------------------------------------------------------------------------------

Protein PEPLIST. Allele H-2-Qa1. Number of high binders 0. Number of weak binders 0. Number of peptides 1
-----------------------------------------------------------------------------------
"""
    preds = parse_netmhc4_stdout(stdout)
    assert len(preds) == 1
    assert preds[0].allele == "H-2-Qa1"       # raw name, not raised on
    assert preds[0].peptide == "SIINFEKLL"
