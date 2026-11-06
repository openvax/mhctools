"""Lossless canonical MHC identities, with parsing owned by mhcgnomes.

Canonical names describe requests. Predictor-specific spellings and model
resolution belong at the adapter boundary, never in these identity keys.
"""

from functools import lru_cache
from importlib.metadata import version

from mhcgnomes import Allele, AlleleWithoutGene, Pair, ParseError, parse


AlleleParseError = ParseError
ALLELE_IDENTITY_POLICY = "mhctools-canonical-v1"
PARSER_OPTIONS = dict(
    use_allele_aliases=False,
    infer_class2_pairing=False,
    collapse_singleton_haplotypes=True,
    collapse_singleton_serotypes=False,
    max_allele_fields=None,
    raise_on_error=True,
)


class AlleleIdentityError(ValueError):
    """A parsed MHC name does not identify one specific allele or valid pair.

    Deliberately distinct from ParseError: raw-name fallback must never rescue
    an ambiguous group, gene, or invalid pair.
    """


def validate_allele_identity(parsed):
    """Return a specific mhcgnomes identity, or raise AlleleIdentityError.

    Pairs must contain a concrete alpha chain followed by a concrete beta
    chain of the same species. Supported locus combinations are a model
    capability, not a nomenclature rule. G/P typings and allele groups cannot
    stand in for one member; expression suffixes and full resolution survive.
    """
    if isinstance(parsed, Pair):
        alpha, beta = parsed.alpha, parsed.beta
        if not (
                isinstance(alpha, Allele) and isinstance(beta, Allele)
                and alpha.species == beta.species
                and alpha.is_class2_alpha and beta.is_class2_beta):
            raise AlleleIdentityError(
                "%s is not a valid class-II alpha/beta pair; provide explicit "
                "alleles of the same species in alpha/beta order" % parsed.to_string())
        validate_allele_identity(alpha)
        validate_allele_identity(beta)
    elif isinstance(parsed, Allele):
        if parsed.is_allele_group or set(parsed.annotations).intersection("GP"):
            raise AlleleIdentityError(
                "%s is an ambiguous allele group; provide a specific allele. "
                "mhctools does not choose a group member." % parsed.to_string())
    elif not isinstance(parsed, AlleleWithoutGene):
        raise AlleleIdentityError(
            "%s is a %s, not a specific allele or explicit class-II pair; "
            "provide a specific allele. mhctools does not choose a group member."
            % (parsed.to_string(), type(parsed).__name__))
    return parsed


@lru_cache(maxsize=16384)
def parse_allele(name):
    """Parse and validate a lossless canonical allele identity.

    Accept mhcgnomes Allele, AlleleWithoutGene, and ordered class-II Pair
    identities. Preserve fields, zeros, and expression annotations; never infer
    a missing chain or resolve historical allele aliases. The documented ``;``
    pair separator and ``_`` spelling of ``*`` remain accepted.

    Raises
    ------
    mhcgnomes.ParseError
        The spelling cannot be parsed.
    AlleleIdentityError
        The parsed name is ambiguous or is not a valid allele identity.
    """
    if not isinstance(name, str):
        raise TypeError("Expected allele name string but got %s" % (type(name),))
    parsed = parse(name.strip().replace(";", "/"),
                   **PARSER_OPTIONS)
    return validate_allele_identity(parsed)


def normalize_allele_name(raw_allele, omit_dra1=False, infer_class2_pair=False):
    """Return mhcgnomes' canonical, lossless spelling of a specific identity.

    The legacy projection/inference keywords are accepted only when false
    during the 4.9 deprecation cycle. True values raise a migration error:
    explicit pairs and native formatters replace those identity-changing modes.
    """
    if omit_dra1 or infer_class2_pair:
        raise AlleleIdentityError(
            "Canonical identities no longer infer class-II partners or omit "
            "DRA. Supply an explicit alpha/beta pair, or use a native allele "
            "formatter for tool-specific spellings; leave infer_class2_pair "
            "and omit_dra1 false.")
    return parse_allele(raw_allele).to_string()


def normalize_allele_name_or_raw(name):
    """Canonicalize a name, retaining raw spelling only on a parsing failure.

    Raw names preserve case and strip whitespace and ``*``. Callers must
    validate these against the native tool's catalog before prediction.
    Parsed groups, genes and invalid pairs always raise AlleleIdentityError.
    """
    stripped = str(name).strip()
    try:
        parsed = parse_allele(stripped)
    except ParseError:
        return stripped.replace("*", "")
    return parsed.to_string()


def native_allele_identity(name):
    """Resolve one *verified native catalog entry* to a canonical identity.

    Only an allele-group spelling may consult mhcgnomes' curated historical
    aliases. Specific identities are never renamed. This function is not a
    user-input normalizer: an ambiguous request must fail even if a catalog
    uses that spelling for a historical allele.
    """
    normalized = name.strip().replace(";", "/")
    try:
        parsed = parse(normalized, **PARSER_OPTIONS)
    except ParseError:
        return name.strip().replace("*", "")
    if isinstance(parsed, Allele) and parsed.is_allele_group:
        parsed = parse(normalized, **dict(PARSER_OPTIONS, use_allele_aliases=True))
    return validate_allele_identity(parsed).to_string()


def allele_identity_provenance():
    """Return the canonicalization policy and parser version for cache keys."""
    return dict(policy=ALLELE_IDENTITY_POLICY, mhcgnomes=version("mhcgnomes"),
                parser_options=dict(PARSER_OPTIONS))
