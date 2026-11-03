"""
Per-tool allele formatters (mhctools.allele_formats).

Each formatter takes an mhcgnomes parse result. The expected spellings below
are the ones the tools print in their own ``-listMHC`` output; the fixtures in
tests/data/allele_lists are those outputs with headers removed, captured from
netMHC 4.0, netMHCpan 4.2 and netMHCIIpan 4.3.
"""

from pathlib import Path

import pytest
from mhcgnomes import Allele, ParseError

from mhctools.allele_formats import (
    caphla_allele_key,
    deepimmuno_allele_name,
    netmhc4_allele_name,
    netmhc_allele_name,
    netmhciipan_allele_name,
    parse_allele,
    prefer_allele_spelling,
)

DATA = Path(__file__).parent / "data" / "allele_lists"


def _listed(filename):
    return {
        line.strip() for line in (DATA / filename).read_text().splitlines()
        if line.strip()
    }


def _format_listed(filename, formatter):
    """Format every listed name mhcgnomes parses; return (ok, missed, unparsed)."""
    listed = _listed(filename)
    ok, missed, unparsed = [], [], []
    for name in sorted(listed):
        try:
            parsed = parse_allele(name)
        except ParseError:
            unparsed.append(name)
            continue
        try:
            formatted = formatter(parsed)
        except ParseError:
            missed.append((name, None))
            continue
        (ok if formatted in listed else missed).append((name, formatted))
    return ok, missed, unparsed


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A02:01"),
    ("HLA-B*57:01", "HLA-B57:01"),
    ("H-2-Kb", "H-2-Kb"),
    ("H2-Kb", "H-2-Kb"),
])
def test_default_drops_star(name, expected):
    assert netmhc_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A0201"),
    ("HLA-B*57:01", "HLA-B5701"),
    ("H-2-Kb", "H-2-Kb"),
    ("Mamu-B*10:01", "Mamu-B1001"),
])
def test_netmhc4(name, expected):
    assert netmhc4_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-DRB1*01:01", "DRB1_0101"),
    ("HLA-DRA1*01:01/DRB1*01:01", "DRB1_0101"),
    ("HLA-DRB3*01:01", "DRB3_0101"),
    ("HLA-DQA1*05:01/DQB1*06:36", "HLA-DQA10501-DQB10636"),
    ("HLA-DPA1*01:05-DPB1*100:01", "HLA-DPA10105-DPB110001"),
    ("H-2-IAb", "H-2-IAb"),
    ("H2-IAd", "H-2-IAd"),
    ("H2-IEk", "H-2-IEk"),
])
def test_netmhciipan(name, expected):
    assert netmhciipan_allele_name(parse_allele(name)) == expected


def test_netmhciipan_rejects_class_i():
    with pytest.raises(ParseError):
        netmhciipan_allele_name(parse_allele("HLA-A*02:01"))


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A*0201"),
    ("HLA-B*57:01", "HLA-B*5701"),
])
def test_deepimmuno(name, expected):
    assert deepimmuno_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A*02:01"),
    ("HLA-DRB1*01:01", "HLA-DRB1*01:01"),
    ("HLA-DQA1*05:01/DQB1*06:36", "HLA-DQA1*05:01/DQB1*06:36"),
    ("H-2-Kb", "H2-Kb"),
    ("H-2-IAb", "H2-IAb"),
])
def test_caphla(name, expected):
    assert caphla_allele_key(parse_allele(name)) == expected


def test_caphla_keeps_expression_suffix():
    parsed = parse_allele("HLA-A*30:14L")
    assert isinstance(parsed, Allele)
    assert caphla_allele_key(parsed) == "HLA-A*30:14L"


@pytest.mark.parametrize("name", ["A2", "HLA-A", "not an allele"])
def test_parse_allele_rejects_non_alleles(name):
    with pytest.raises(ParseError):
        parse_allele(name)


# --- against each tool's own bundled list -------------------------------

def test_netmhciipan_formatter_matches_own_list():
    ok, missed, unparsed = _format_listed(
        "netmhciipan43_list.txt", netmhciipan_allele_name)
    assert len(ok) >= 10600
    # Everything mhcgnomes parses is either accepted or a documented non-DR/
    # non-pair name (non-human class II the formatter has no spelling for).
    assert len(missed) <= 10, missed
    # Unparseable names keep the verbatim path; make sure we notice growth.
    assert len(unparsed) <= 400


def test_netmhc4_formatter_matches_own_list():
    ok, missed, unparsed = _format_listed(
        "netmhc4_listmhc.txt", netmhc4_allele_name)
    assert len(ok) >= 120
    # Remaining misses are old-style Mamu-B names resolved by the -listMHC
    # mapping rather than by the formatter.
    assert len(missed) <= 4, missed
    assert len(unparsed) <= 8


def test_netmhc4_formatter_accepts_all_human_and_mouse_names():
    _, missed, _ = _format_listed("netmhc4_listmhc.txt", netmhc4_allele_name)
    assert [m for m in missed if m[0].startswith(("HLA-", "H-2-"))] == []


def test_default_formatter_accepts_all_human_and_mouse_names():
    listed = _listed("netmhcpan42_listmhc.txt")
    _, missed, _ = _format_listed(
        "netmhcpan42_listmhc.txt", netmhc_allele_name)
    # Some human alleles are listed only colon-less (HLA-A0111, no
    # HLA-A01:11). The -listMHC mapping supplies that spelling; the formatter
    # must at least name an allele the tool knows.
    unknown = [
        (name, out) for name, out in missed
        if name.startswith(("HLA-", "H-2-")) and out is not None
        and out.replace(":", "") not in listed]
    assert unknown == []
    # Expression variants (HLA-A30:14L) are declined and keep the verbatim
    # path rather than being silently reformatted.
    declined = [name for name, out in missed if out is None]
    assert all(name[-1] in "NLSCAQ" for name in declined), declined


def test_default_formatter_coverage_floor_on_netmhcpan():
    ok, missed, unparsed = _format_listed(
        "netmhcpan42_listmhc.txt", netmhc_allele_name)
    # The rest (non-human species whose families mhcgnomes pads differently)
    # resolve through the -listMHC mapping, not this formatter.
    assert len(ok) >= 12000
    assert len(unparsed) <= 25


# --- spelling preference ---------------------------------------------------

def test_prefer_allele_spelling_picks_colon_form():
    assert prefer_allele_spelling("BoLA-1:00901", "BoLA-100901") == \
        "BoLA-1:00901"
    assert prefer_allele_spelling("BoLA-100901", "BoLA-1:00901") == \
        "BoLA-1:00901"


def test_prefer_allele_spelling_earliest_colon_then_lexicographic():
    assert prefer_allele_spelling("Mamu-A1001:01", "Mamu-A1:00101") == \
        "Mamu-A1:00101"
    assert prefer_allele_spelling("b", "a") == "a"
