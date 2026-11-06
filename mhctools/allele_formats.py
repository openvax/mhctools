"""
Per-tool allele spellings.

mhcgnomes owns parsing and the generic spelling of an MHC name. Each external
predictor, however, accepts its own idiosyncratic form on the command line
(``HLA-A0201`` for NetMHC 4, ``DRB1_0101`` for NetMHCIIpan, ...). That is
knowledge about the tool, not about MHC nomenclature, so it lives here.

API
---
:func:`parse_allele`
    Parse a validated canonical ``Allele``, ``AlleleWithoutGene`` or ``Pair``.
:func:`format_allele`
    Parse a name and format it for a tool, with an optional fallback for names
    mhcgnomes cannot parse.
Formatters
    One per tool, each taking an mhcgnomes parse result and returning that
    tool's spelling: :func:`netmhc_allele_name`, :func:`netmhc4_allele_name`,
    :func:`netmhciipan_allele_name`, :func:`deepimmuno_allele_name`,
    :func:`caphla_allele_key`.
:func:`prefer_allele_spelling`
    Choose between several spellings of one allele in a tool's own allele list.

Contract
--------
A formatter either returns a spelling or raises ``mhcgnomes.ParseError`` when
the tool has no spelling for that name (wrong MHC class for the tool, a mutant
allele, ...). Formatters never guess.

Command-line predictors prefer the spelling found in the tool's own
``-listMHC`` output (see ``BaseCommandlinePredictor``); a formatter is only
used when the name is not found there. Names mhcgnomes cannot parse
have no parse result to format and require a verified native-list match.
Parsed genes, groups and invalid pairs are rejected even with a raw fallback.
"""

from mhcgnomes import Allele, AlleleWithoutGene, Gene, Pair, ParseError

from .allele_normalization import parse_allele, validate_allele_identity

__all__ = [
    "parse_allele",
    "format_allele",
    "netmhc_allele_name",
    "netmhc4_allele_name",
    "netmhciipan_allele_name",
    "deepimmuno_allele_name",
    "caphla_allele_key",
    "prefer_allele_spelling",
]


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def format_allele(name, formatter, fallback=None):
    """
    Format an allele name for one tool.

    Parameters
    ----------
    name : str
    formatter : callable
        One of the formatters in this module.
    fallback : callable, optional
        Called with ``name`` when the name cannot be parsed or the formatter
        has no spelling for it. Without it, ``ParseError`` propagates.
    """
    try:
        return formatter(parse_allele(name))
    except ParseError:
        if fallback is None:
            raise
        return fallback(name)


# --------------------------------------------------------------------------
# Shared pieces
# --------------------------------------------------------------------------

def _species(allele):
    prefix = allele.species_prefix
    return "H-2" if prefix in ("H2", "H-2") else prefix


def _is_mouse(allele):
    return _species(allele) == "H-2"


def _gene(allele):
    return str(allele.gene_name).upper()


def _fields(allele):
    fields = [str(field) for field in allele.allele_fields]
    if not fields or allele.mutations:
        raise ParseError("%s has no tool spelling" % (allele,))
    return fields


def _body(allele, fields, sep=":", star="*", annotations=True):
    """
    ``A*02:01`` by default; ``A0201`` with ``sep="", star=""``.

    A single field (mouse ``Kb``, ``Mamu-B*001``) is appended to the gene with
    no separator. Expression suffixes (``N``, ``L``, ``Q``) are kept unless
    ``annotations`` is false, in which case they are an error.
    """
    if len(fields) == 1:
        code = fields[0].lower() if _is_mouse(allele) else fields[0]
        body = _gene(allele) + code
    else:
        body = _gene(allele) + star + sep.join(fields)
    if allele.annotations:
        if not annotations:
            raise ParseError("%s has an unsupported suffix" % (allele,))
        body += "".join(allele.annotations)
    return body


def _with_species(allele, body):
    return "%s-%s" % (_species(allele), body)


def _mouse_class2(pair):
    """
    ``("IA" | "IE", code)`` for a mouse class II pair whose chains share an
    allele (mhcgnomes ``H2-AA*b/AB*b`` is ``H-2-IAb``), else ``None``.
    """
    alpha, beta = pair.alpha, pair.beta
    if not (_is_mouse(alpha) and _is_mouse(beta)):
        return None
    locus = {("AA", "AB"): "IA", ("EA", "EB"): "IE"}.get(
        (_gene(alpha), _gene(beta)))
    alpha_fields, beta_fields = _fields(alpha), _fields(beta)
    if (locus is None or len(alpha_fields) != 1 or len(beta_fields) != 1
            or alpha_fields[0].lower() != beta_fields[0].lower()):
        return None
    return locus, alpha_fields[0].lower()


def _is_human_dr_alpha(allele):
    return _species(allele) == "HLA" and _gene(allele) in ("DRA", "DRA1")


# --------------------------------------------------------------------------
# Formatters
# --------------------------------------------------------------------------

def netmhc_allele_name(parsed):
    """
    Default command-line spelling: the allele with ``*`` dropped
    (``HLA-A02:01``, ``H-2-Kb``, ``HLA-A30:14L``).

    Used by the NetMHCpan family and the other wrappers that do not override
    ``prepare_allele_name``. Class II pairs have no syntax in these tools and
    are joined with ``-``.
    """
    validate_allele_identity(parsed)
    if isinstance(parsed, AlleleWithoutGene):
        return parsed.to_string()
    if isinstance(parsed, Pair):
        return "%s-%s" % (
            netmhc_allele_name(parsed.alpha),
            _body(parsed.beta, _fields(parsed.beta), star=""))
    return _with_species(
        parsed, _body(parsed, _fields(parsed), star=""))


def netmhc4_allele_name(parsed):
    """
    NetMHC 4.0 spelling: no ``*`` and no ``:`` (``HLA-A0201``, ``H-2-Kb``).

    mhcgnomes pads non-human allele families to three digits
    (``Mamu-B*010:01``, ``Mamu-B*001``); NetMHC 4 lists the older two-digit
    spelling (``Mamu-B1001``, ``Mamu-B01``), so a leading zero on a
    three-digit family is dropped.
    Class II and expression suffixes are not supported (its list has none).
    """
    validate_allele_identity(parsed)
    if isinstance(parsed, AlleleWithoutGene):
        return parsed.to_string()
    if isinstance(parsed, Pair) or not parsed.is_class1:
        raise ParseError("NetMHC 4 is class I only, got %s" % (parsed,))
    fields = _fields(parsed)
    if len(fields[0]) == 3 and fields[0].startswith("0"):
        fields[0] = fields[0][1:]
    return _with_species(
        parsed, _body(parsed, fields, sep="", star="", annotations=False))


def netmhciipan_allele_name(parsed):
    """
    NetMHCIIpan spelling. The tool accepts only:

    - ``DRB1_0101`` for human DR beta chains (the DR alpha chain
      is implicit in this tool and is omitted), including expression variants (``DRB5_0108N``);
    - ``HLA-DQA10501-DQB10636`` for other human alpha/beta pairs;
    - ``H-2-IAb`` / ``H-2-IEk`` for mouse class II.

    Anything else, including class I, raises ``ParseError``.
    """
    validate_allele_identity(parsed)
    if not isinstance(parsed, (Allele, Pair)):
        raise ParseError("NetMHCIIpan requires a class-II allele or explicit pair")
    if isinstance(parsed, Allele):
        # A bare beta chain only arrives here when mhcgnomes did not pair it.
        return _drb_name(parsed)
    mouse = _mouse_class2(parsed)
    if mouse is not None:
        return "H-2-%s%s" % mouse
    alpha, beta = parsed.alpha, parsed.beta
    if _species(alpha) != "HLA" or _species(beta) != "HLA":
        raise ParseError("NetMHCIIpan has no spelling for %s" % (parsed,))
    if _is_human_dr_alpha(alpha):
        return _drb_name(beta)
    return "HLA-%s-%s" % (
        _body(alpha, _fields(alpha), sep="", star=""),
        _body(beta, _fields(beta), sep="", star=""))


def _drb_name(allele):
    if _species(allele) != "HLA" or not _gene(allele).startswith("DRB"):
        raise ParseError("NetMHCIIpan has no spelling for %s" % (allele,))
    fields = _fields(allele)
    return "%s_%s%s" % (
        _gene(allele), "".join(fields), "".join(allele.annotations))


def deepimmuno_allele_name(parsed):
    """
    DeepImmuno paratope-table spelling: ``HLA-A*0201``.

    Class I only, no expression suffixes (the table has none).
    """
    validate_allele_identity(parsed)
    if not isinstance(parsed, Allele) or not parsed.is_class1:
        raise ParseError("DeepImmuno is class I only, got %s" % (parsed,))
    return _with_species(
        parsed,
        _body(parsed, _fields(parsed), sep="", annotations=False))


#: CapHLA keeps the historical names for two non-classical mouse genes that
#: mhcgnomes resolves to their official names.
_CAPHLA_GENE_KEYS = {"H2-T23": "H2-Qa1", "H2-Q7": "H2-Qa2"}


def caphla_allele_key(parsed):
    """
    CapHLA library key.

    - DR alpha is implicit: ``HLA-DRB1*15:01``.
    - Other class II chains are joined by ``/`` with the beta chain's species
      prefix dropped: ``HLA-DPA1*01:03/DPB1*04:01``.
    - Mouse uses the ``H2-`` prefix (``H2-Db``) and the class II shorthand
      (``H2-IAb``).
    - Expression suffixes are kept (``HLA-A*30:14L``); CapHLA has
      pseudosequences for some.
    - The two non-classical mouse genes accept ``H2-Qa1`` / ``H2-Qa2``
      (mhcgnomes ``Gene`` results for ``H2-T23`` / ``H2-Q7``).
    """
    if isinstance(parsed, Gene):
        name = parsed.to_string()
        try:
            return _CAPHLA_GENE_KEYS[name]
        except KeyError:
            raise ParseError(
                "CapHLA does not support gene-level MHC name %s" % name)
    validate_allele_identity(parsed)
    if isinstance(parsed, AlleleWithoutGene):
        raise ParseError("CapHLA has no spelling for %s" % parsed.to_string())
    if isinstance(parsed, Allele):
        return _caphla_chain(parsed, include_species=True)
    mouse = _mouse_class2(parsed)
    if mouse is not None:
        return "H2-%s%s" % mouse
    alpha, beta = parsed.alpha, parsed.beta
    if _is_human_dr_alpha(alpha):
        return _caphla_chain(beta, include_species=True)
    return "%s/%s" % (
        _caphla_chain(alpha, include_species=True),
        _caphla_chain(beta, include_species=False))


def _caphla_chain(allele, include_species):
    body = _body(allele, _fields(allele))
    if not include_species:
        return body
    species = "H2" if _is_mouse(allele) else _species(allele)
    return "%s-%s" % (species, body)


# --------------------------------------------------------------------------
# Choosing among a tool's own spellings
# --------------------------------------------------------------------------

def prefer_allele_spelling(new, existing):
    """
    Choose between two spellings of one allele from a tool's allele list.

    When several raw ``-listMHC`` names normalize to the same allele, pick the
    one most likely accepted on ``-a``, deterministically (independent of set
    iteration order):

    1. a colon-containing spelling (``BoLA-1:00901``, which netMHCpan accepts,
       over ``BoLA-100901``, which it rejects);
    2. among those, the spelling whose colon sits earliest, i.e. right after
       the gene (``Mamu-A1:00101`` over ``Mamu-A1001:01``);
    3. the lexicographically smaller string as a final stable tie-break.

    Returns
    -------
    str
        Either ``new`` or ``existing``.
    """
    def rank(name):
        colon = name.find(":")
        # higher is preferred: has-colon first, then earliest colon
        return (colon >= 0, -colon if colon >= 0 else 1)

    rank_new, rank_existing = rank(new), rank(existing)
    if rank_new != rank_existing:
        return new if rank_new > rank_existing else existing
    return min(new, existing)
