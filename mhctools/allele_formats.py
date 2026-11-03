"""
Per-tool allele spellings.

mhcgnomes owns parsing and the generic spelling of an MHC name, but each
external predictor accepts its own idiosyncratic form on the command line
(``HLA-A0201`` for NetMHC 4, ``DRB1_0101`` for NetMHCIIpan, ...). That is
knowledge about the tool, so it lives here: one public formatter per tool,
each taking an mhcgnomes parse result (see :func:`parse_allele`) and returning
the spelling that tool expects.

Names mhcgnomes cannot parse yet (e.g. ``H-2-Qa1``, ``BoLA-amani.1``) have no
parse result to format; callers keep passing those through verbatim.
"""

from functools import lru_cache

from mhcgnomes import Allele, Pair, ParseError, parse

__all__ = [
    "parse_allele",
    "netmhc_allele_name",
    "netmhc4_allele_name",
    "netmhciipan_allele_name",
    "deepimmuno_allele_name",
    "caphla_allele_key",
    "prefer_allele_spelling",
]


@lru_cache(maxsize=None)
def parse_allele(name, infer_class2_pairing=True):
    """
    Parse an allele name with mhcgnomes, tolerating the ``;`` pair separator
    used in some publications and ``_`` standing in for ``*``.

    Returns an mhcgnomes ``Allele`` or ``Pair``; raises ``ParseError`` for
    anything else (serotypes, genes, unparseable names).
    """
    if not isinstance(name, str):
        raise TypeError("Expected allele name string but got %s" % (type(name),))
    cleaned = name.strip().replace(";", "/").replace("_", "*")
    parsed = parse(
        cleaned,
        infer_class2_pairing=infer_class2_pairing,
        raise_on_error=True)
    if not isinstance(parsed, (Allele, Pair)):
        raise ParseError(
            "%r is a %s, not a specific allele or allele pair" % (
                name, type(parsed).__name__))
    if isinstance(parsed, Pair) and not (
            isinstance(parsed.alpha, Allele) and
            isinstance(parsed.beta, Allele)):
        raise ParseError("Malformed MHC pair %r" % (name,))
    return parsed


def _species(allele):
    prefix = allele.species_prefix
    return "H-2" if prefix in ("H2", "H-2") else prefix


def _gene(allele):
    return str(allele.gene_name).upper()


def _fields(allele, legacy_family=False):
    fields = [str(field) for field in allele.allele_fields]
    if len(fields) == 0:
        raise ParseError("Allele %s has no allele fields" % (allele,))
    # mhcgnomes pads non-human families to three digits (Mamu-B*010:01);
    # NetMHC 4 still lists the older two-digit spelling (Mamu-B1001).
    if (legacy_family and len(fields) > 1 and len(fields[0]) == 3
            and fields[0].isdigit() and fields[0][0] == "0"):
        fields[0] = fields[0][1:]
    return fields


def _require_plain(allele):
    if allele.annotations or allele.mutations:
        raise ParseError(
            "Annotated or mutant allele %s has no tool spelling" % (allele,))


def _is_mouse(allele):
    return _species(allele) == "H-2"


def _gene_and_fields(allele, sep=":", star="*", legacy_family=False):
    """``A*02:01`` (default), or e.g. ``A0201`` with ``sep="", star=""``.

    Single-field alleles (mouse ``Kb``) have no separator at all.
    """
    _require_plain(allele)
    fields = _fields(allele, legacy_family)
    if len(fields) == 1:
        code = fields[0].lower() if _is_mouse(allele) else fields[0]
        return "%s%s" % (_gene(allele), code)
    return "%s%s%s" % (_gene(allele), star, sep.join(fields))


def _single_name(allele, sep=":", star="*", legacy_family=False):
    return "%s-%s" % (
        _species(allele), _gene_and_fields(allele, sep, star, legacy_family))


def _mouse_class2_locus(pair):
    """``IA``/``IE`` plus the shared allele code for H-2 class II pairs."""
    alpha, beta = pair.alpha, pair.beta
    if not (_is_mouse(alpha) and _is_mouse(beta)):
        return None
    loci = {("AA", "AB"): "IA", ("EA", "EB"): "IE"}
    locus = loci.get((_gene(alpha), _gene(beta)))
    alpha_fields, beta_fields = _fields(alpha), _fields(beta)
    if (locus is None or len(alpha_fields) != 1 or len(beta_fields) != 1
            or alpha_fields[0].lower() != beta_fields[0].lower()):
        return None
    return locus, alpha_fields[0].lower()


def _human_drb_gene(allele):
    return _species(allele) == "HLA" and "DRB" in _gene(allele)


def netmhc_allele_name(parsed):
    """
    The default command-line spelling: the mhcgnomes name with ``*`` dropped
    (``HLA-A02:01``). Used by the NetMHCpan family and most other wrappers.
    Pairs are joined with ``-`` as these tools have no pair syntax of their own.
    """
    if isinstance(parsed, Pair):
        return "-".join([
            netmhc_allele_name(parsed.alpha),
            _gene_and_fields(parsed.beta, star="")])
    return _single_name(parsed, star="")


def netmhc4_allele_name(parsed):
    """NetMHC 4.0 spelling: no ``*`` and no ``:`` (``HLA-A0201``)."""
    return _single_name(parsed, sep="", star="", legacy_family=True)


def netmhciipan_allele_name(parsed):
    """
    NetMHCIIpan spelling:

    - ``DRB1_0101`` for human DR (DRA is invariant, so it is dropped),
    - ``HLA-DQA10501-DQB10636`` for other alpha/beta pairs,
    - ``H-2-IAb`` / ``H-2-IEk`` for the mouse class II shorthand.

    Raises ``ParseError`` for anything else, as NetMHCIIpan accepts nothing
    else.
    """
    if isinstance(parsed, Allele):
        if _is_mouse(parsed):
            return "H-2-%s%s" % (_gene(parsed), _fields(parsed)[-1].lower())
        return _drb_name(parsed)
    mouse = _mouse_class2_locus(parsed)
    if mouse is not None:
        locus, code = mouse
        return "H-2-%s%s" % (locus, code)
    alpha, beta = parsed.alpha, parsed.beta
    if _species(alpha) == "HLA" and _gene(alpha) in ("DRA", "DRA1"):
        return _drb_name(beta)
    return "HLA-%s-%s" % (
        _gene_and_fields(alpha, sep="", star=""),
        _gene_and_fields(beta, sep="", star=""))


def _drb_name(allele):
    _require_plain(allele)
    if "DRB" not in _gene(allele):
        raise ParseError("Unexpected allele %s" % (allele,))
    return "%s_%s" % (_gene(allele), "".join(_fields(allele)))


def deepimmuno_allele_name(parsed):
    """DeepImmuno paratope-table spelling: ``HLA-A*0201``."""
    if isinstance(parsed, Pair):
        raise ParseError("DeepImmuno takes class I alleles, got %s" % (parsed,))
    return _single_name(parsed, sep="")


def caphla_allele_key(parsed):
    """
    CapHLA key: DRA is implicit for DR, DP/DQ chains are joined with ``/``
    (beta without species prefix), mouse alleles use the ``H2-`` prefix with
    the ``IA``/``IE`` shorthand for class II, and expression suffixes
    (``N``, ``L``) are kept because CapHLA has pseudosequences for them.
    """
    if isinstance(parsed, Allele):
        return _caphla_component(parsed, include_species=True)
    mouse = _mouse_class2_locus(parsed)
    if mouse is not None:
        return "H2-%s%s" % mouse
    alpha, beta = parsed.alpha, parsed.beta
    if _gene(alpha) in ("DRA", "DRA1"):
        return _caphla_component(beta, include_species=True)
    return "%s/%s" % (
        _caphla_component(alpha, include_species=True),
        _caphla_component(beta, include_species=False))


def _caphla_component(allele, include_species):
    if allele.mutations:
        raise ParseError("Mutant allele %s has no CapHLA key" % (allele,))
    species = _species(allele)
    if species == "H-2":
        species = "H2"
    fields = _fields(allele)
    body = _gene(allele) + (
        "*" + ":".join(fields) if len(fields) > 1 else fields[0])
    body += "".join(allele.annotations)
    return "%s-%s" % (species, body) if include_species else body


def prefer_allele_spelling(new, existing):
    """
    When several raw ``-listMHC`` spellings normalize to the same allele, pick
    the one most likely accepted on ``-a``, deterministically (independent of
    set iteration order):

      1. a colon-containing spelling (e.g. ``BoLA-1:00901``, which netMHCpan
         accepts, over ``BoLA-100901``, which it rejects);
      2. among those, the spelling whose colon sits earliest, i.e. right after
         the gene (``Mamu-A1:00101`` over ``Mamu-A1001:01``);
      3. the lexicographically smaller string as a final stable tie-break.
    """
    def rank(name):
        colon = name.find(":")
        # higher is preferred: has-colon first, then earliest colon
        return (colon >= 0, -colon if colon >= 0 else 1)

    rank_new, rank_existing = rank(new), rank(existing)
    if rank_new != rank_existing:
        return new if rank_new > rank_existing else existing
    return min(new, existing)
