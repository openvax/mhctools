"""
Per-tool allele formatters (mhctools.allele_formats).

The spelling tables below were checked against the tools' own ``-listMHC``
output and against the behavior of the compatibility layer they replaced. The
fixtures in tests/data/allele_lists are those lists (headers removed, sampled
by scripts/capture_allele_lists.py): every formatter is checked against the
tool it serves.
"""

from pathlib import Path

import pytest
from mhcgnomes import Allele, AlleleWithoutGene, ParseError, parse

from mhctools.allele_normalization import AlleleIdentityError

from mhctools.allele_formats import (
    caphla_allele_key,
    deepimmuno_allele_name,
    format_allele,
    netmhc4_allele_name,
    netmhc_allele_name,
    netmhciipan_allele_name,
    parse_allele,
    prefer_allele_spelling,
)

DATA = Path(__file__).parent / "data" / "allele_lists"
HUMAN_OR_MOUSE = ("HLA-", "H-2-", "DRB")
# Gene-level names mhcgnomes parses as genes, not alleles (pirl-unc/mhcgnomes#190).
KNOWN_UNPARSED = {"H-2-Qa1", "H-2-Qa2", "H2-Qa1", "H2-Qa2"}


# --- spelling tables -------------------------------------------------------

@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A02:01"),
    ("HLA-B*57:01", "HLA-B57:01"),
    ("HLA-A*30:14L", "HLA-A30:14L"),
    ("H-2-Kb", "H-2-Kb"),
    ("H2-Kb", "H-2-Kb"),
    ("Mamu-A1*001:01", "Mamu-A1001:01"),
])
def test_default_drops_star(name, expected):
    assert netmhc_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A0201"),
    ("HLA-B*57:01", "HLA-B5701"),
    ("H-2-Kb", "H-2-Kb"),
    ("Mamu-B*10:01", "Mamu-B1001"),
    ("Mamu-B*001:01", "Mamu-B0101"),
])
def test_netmhc4(name, expected):
    assert netmhc4_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-DRB1*01:01", "DRB1_0101"),
    ("HLA-DRA1*01:01/DRB1*01:01", "DRB1_0101"),
    ("HLA-DRB3*01:01", "DRB3_0101"),
    ("HLA-DRB5*01:08N", "DRB5_0108N"),
    ("HLA-DQA1*05:01/DQB1*06:36", "HLA-DQA10501-DQB10636"),
    ("HLA-DQA1*01:02;DQB1*06:02", "HLA-DQA10102-DQB10602"),
    ("HLA-DPA1*01:05-DPB1*100:01", "HLA-DPA10105-DPB110001"),
    ("H-2-IAb", "H-2-IAb"),
    ("H2-IAd", "H-2-IAd"),
    ("H2-IEk", "H-2-IEk"),
])
def test_netmhciipan(name, expected):
    assert netmhciipan_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A*0201"),
    ("HLA-B*57:01", "HLA-B*5701"),
    ("HLA-C*07:02", "HLA-C*0702"),
])
def test_deepimmuno(name, expected):
    assert deepimmuno_allele_name(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("HLA-A*02:01", "HLA-A*02:01"),
    ("HLA-A*30:14L", "HLA-A*30:14L"),
    ("HLA-DRB1*01:01", "HLA-DRB1*01:01"),
    ("HLA-DRB5*01:08N", "HLA-DRB5*01:08N"),
    ("HLA-DQA1*05:01/DQB1*06:36", "HLA-DQA1*05:01/DQB1*06:36"),
    ("H-2-Kb", "H2-Kb"),
    ("H-2-IAb", "H2-IAb"),
    ("H-2-IEk", "H2-IEk"),
])
def test_caphla(name, expected):
    assert caphla_allele_key(parse_allele(name)) == expected


@pytest.mark.parametrize("name,expected", [
    ("H2-Qa1", "H2-Qa1"),
    ("H2-Qa2", "H2-Qa2"),
])
def test_caphla_non_classical_mouse_genes(name, expected):
    assert caphla_allele_key(parse(name, raise_on_error=True)) == expected


def test_caphla_rejects_other_genes():
    with pytest.raises(ParseError, match="gene-level"):
        caphla_allele_key(parse("H2-K", raise_on_error=True))


# --- declining: the tool has no spelling for these -------------------------

@pytest.mark.parametrize("formatter,name", [
    (netmhciipan_allele_name, "HLA-A*02:01"),            # class I
    (netmhciipan_allele_name, "H-2-Kb"),                 # mouse class I
    (netmhc4_allele_name, "HLA-DRB1*01:01"),             # class II
    (deepimmuno_allele_name, "HLA-DRB1*01:01"),          # class II
    (netmhc4_allele_name, "HLA-A*30:14L"),               # expression suffix
    (deepimmuno_allele_name, "HLA-A*30:14L"),            # expression suffix
])
def test_formatter_declines(formatter, name):
    with pytest.raises(ParseError):
        formatter(parse_allele(name))


@pytest.mark.parametrize("name", ["A2", "HLA-A", "not an allele", ""])
def test_parse_allele_rejects_non_alleles(name):
    with pytest.raises((ParseError, AlleleIdentityError)):
        parse_allele(name)


def test_parse_allele_requires_a_string():
    with pytest.raises(TypeError):
        parse_allele(None)


def test_parse_allele_returns_alleles_and_pairs():
    assert isinstance(parse_allele("HLA-A*02:01"), Allele)
    assert isinstance(parse_allele("HLA-DRB1*01:01"), Allele)
    assert not isinstance(parse_allele("HLA-DRA*01:01/DRB1*01:01"), Allele)


# --- format_allele ---------------------------------------------------------

def test_format_allele_parses_and_formats():
    assert format_allele("HLA-A*02:01", netmhc4_allele_name) == "HLA-A0201"


def test_format_allele_propagates_errors_without_fallback():
    with pytest.raises(ParseError):
        format_allele("native-Xfoo", netmhc4_allele_name)
    with pytest.raises(ParseError):
        format_allele("HLA-A*02:01", netmhciipan_allele_name)


def test_format_allele_fallback_gets_the_original_name():
    assert format_allele(
        "native-Xfoo", netmhc4_allele_name, fallback=lambda name: name + "!") == \
        "native-Xfoo!"


# --- against each tool's own list -------------------------------------------

def _listed(filename):
    return [
        line.strip() for line in (DATA / filename).read_text().splitlines()
        if line.strip()
    ]


def _sweep(filename, formatter):
    """
    Format every listed name that mhcgnomes parses.

    Returns ``(formatted, unparsed)``: ``formatted`` maps listed name to its
    formatted spelling (``None`` where the formatter declined).
    """
    formatted, unparsed = {}, []
    for name in _listed(filename):
        try:
            parsed = parse_allele(name)
        except (ParseError, AlleleIdentityError):
            unparsed.append(name)
            continue
        try:
            formatted[name] = formatter(parsed)
        except ParseError:
            formatted[name] = None
    return formatted, unparsed


def _colonless(names):
    return {name.replace(":", "") for name in names}


@pytest.mark.parametrize("filename,formatter", [
    ("netmhc4_listmhc.txt", netmhc4_allele_name),
    ("netmhciipan43_list.txt", netmhciipan_allele_name),
])
def test_formatter_reproduces_listed_names_exactly(filename, formatter):
    """For these tools every parseable listed name round-trips to a listed one."""
    listed = set(_listed(filename))
    formatted, _ = _sweep(filename, formatter)
    assert formatted
    assert {n: o for n, o in formatted.items()
            if o is not None and o not in listed
            and not isinstance(parse_allele(n), AlleleWithoutGene)} == {}
    # A listed allele without a generic formatter still round-trips through
    # the native catalog; no naming rule is guessed for it.
    from mhctools.base_commandline_predictor import BaseCommandlinePredictor
    mapping = BaseCommandlinePredictor._normalized_supported_alleles(
        filename, "-list", listed)
    for name, output in formatted.items():
        if output is None or output not in listed:
            assert parse_allele(name).to_string() in mapping


@pytest.mark.parametrize("filename,formatter", [
    ("netmhc4_listmhc.txt", netmhc4_allele_name),
    ("netmhciipan43_list.txt", netmhciipan_allele_name),
    ("netmhcpan42_listmhc.txt", netmhc_allele_name),
])
def test_unparseable_listed_names_are_never_human_or_mouse(filename, formatter):
    """Names mhcgnomes cannot parse keep the tool-list path; none of them may
    be a human or mouse allele, where the formatters must be complete."""
    _, unparsed = _sweep(filename, formatter)
    assert [n for n in unparsed
            if n.startswith(HUMAN_OR_MOUSE) and n not in KNOWN_UNPARSED] == []


def test_default_formatter_names_every_human_and_mouse_allele_the_tool_lists():
    """
    netMHCpan lists some human alleles only without a colon (``HLA-A0111`` with
    no ``HLA-A01:11``). The ``-listMHC`` mapping supplies that spelling, so
    the formatter must only name an allele the tool knows.
    """
    listed = _listed("netmhcpan42_listmhc.txt")
    formatted, _ = _sweep("netmhcpan42_listmhc.txt", netmhc_allele_name)
    known = _colonless(listed)
    human = {n: o for n, o in formatted.items()
             if n.startswith(HUMAN_OR_MOUSE)}
    assert human
    assert {n: o for n, o in human.items()
            if o is None or o.replace(":", "") not in known} == {}


def test_default_formatter_is_only_a_fallback_for_other_species():
    """
    mhcgnomes renames species and pads families (``Chi-B0401`` is
    ``Gaga-BF2*04:01``), so for non-human species the tool's own list, not
    the formatter, decides the spelling. Pin which species this affects so a
    new one is noticed.
    """
    listed = _listed("netmhcpan42_listmhc.txt")
    formatted, _ = _sweep("netmhcpan42_listmhc.txt", netmhc_allele_name)
    known = _colonless(listed)
    differing_species = {
        n.split("-")[0] for n, o in formatted.items()
        if o is not None and o.replace(":", "") not in known}
    assert differing_species <= {"Mamu", "SLA", "Chi"}


# --- spelling preference ----------------------------------------------------

def test_prefer_allele_spelling_picks_colon_form():
    assert prefer_allele_spelling("BoLA-1:00901", "BoLA-100901") == \
        "BoLA-1:00901"
    assert prefer_allele_spelling("BoLA-100901", "BoLA-1:00901") == \
        "BoLA-1:00901"


def test_prefer_allele_spelling_earliest_colon_then_lexicographic():
    assert prefer_allele_spelling("Mamu-A1001:01", "Mamu-A1:00101") == \
        "Mamu-A1:00101"
    assert prefer_allele_spelling("b", "a") == "a"
    assert prefer_allele_spelling("a", "b") == "a"
