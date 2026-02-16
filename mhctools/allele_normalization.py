"""
Compatibility helpers for allele normalization/parsing backed by mhcgnomes.

This module preserves the subset of mhcnames behavior relied on by mhctools:
- normalize_allele_name
- parse_classi_or_classii_allele_name
- AlleleParseError
"""

from __future__ import print_function, division, absolute_import

from collections import namedtuple

from mhcgnomes import Allele, Pair, ParseError, parse


AlleleParseError = ParseError

AlleleName = namedtuple("AlleleName", [
    "species",
    "gene",
    "allele_family",
    "allele_code",
])

_DRA1_0101 = AlleleName(
    species="HLA",
    gene="DRA1",
    allele_family="01",
    allele_code="01")

_parse_cache = {}
_normalized_allele_cache = {}


def _normalize_species(species_prefix):
    if species_prefix in ("H2", "H-2"):
        return "H-2"
    return species_prefix


def _normalize_gene(species, gene_name):
    gene_name = str(gene_name).upper()
    # Preserve historical mhcnames behavior of expanding human DRA to DRA1.
    if species == "HLA" and gene_name == "DRA":
        return "DRA1"
    return gene_name


def _normalize_family(family):
    family = str(family)
    if family.isdigit():
        if len(family) == 1:
            family = "0" + family
        elif len(family) == 3 and family[0] == "0":
            family = family[1:]
    return family


def _normalize_code(code):
    code = str(code)
    if len(code) == 0:
        return "01"
    if code.isdigit() and len(code) == 3 and code[0] == "0":
        return code[1:]
    return code


def _annotation_parse_error(parsed_allele):
    original = parsed_allele.raw_string or parsed_allele.to_string()
    suffix = "".join(parsed_allele.annotations)
    if len(suffix) == 0:
        suffix = "".join(str(m) for m in parsed_allele.mutations)
    return AlleleParseError(
        "The suffix '%s' of '%s' was not parsed" % (suffix, original))


def _allele_name_from_parsed_allele(parsed_allele):
    if not isinstance(parsed_allele, Allele):
        raise AlleleParseError("Malformed MHC type %s" % (parsed_allele,))

    if parsed_allele.annotations or parsed_allele.mutations:
        raise _annotation_parse_error(parsed_allele)

    species = _normalize_species(parsed_allele.species_prefix)
    gene = _normalize_gene(species, parsed_allele.gene_name)
    allele_fields = [str(field) for field in parsed_allele.allele_fields]

    if len(allele_fields) == 0:
        raise AlleleParseError("Malformed MHC type %s" % (parsed_allele,))

    if len(allele_fields) == 1:
        allele_family = ""
        allele_code = allele_fields[0]
    else:
        allele_family = _normalize_family(allele_fields[0])
        allele_code = _normalize_code(allele_fields[1])

    if species == "H-2" and len(allele_family) == 0:
        allele_code = allele_code.lower()

    return AlleleName(
        species=species,
        gene=gene,
        allele_family=allele_family,
        allele_code=allele_code)


def _mouse_shorthand_from_pair(parsed_pair):
    """
    Convert mhcgnomes' H2-AA/AB or H2-EA/EB class-II pair forms to
    historical mhcnames shorthand (H-2-IAx / H-2-IEx).
    """
    alpha = parsed_pair.alpha
    beta = parsed_pair.beta

    if not isinstance(alpha, Allele) or not isinstance(beta, Allele):
        return None

    alpha_species = _normalize_species(alpha.species_prefix)
    beta_species = _normalize_species(beta.species_prefix)
    if alpha_species != "H-2" or beta_species != "H-2":
        return None

    alpha_gene = str(alpha.gene_name).upper()
    beta_gene = str(beta.gene_name).upper()
    if (alpha_gene, beta_gene) == ("AA", "AB"):
        locus = "IA"
    elif (alpha_gene, beta_gene) == ("EA", "EB"):
        locus = "IE"
    else:
        return None

    alpha_fields = [str(x) for x in alpha.allele_fields]
    beta_fields = [str(x) for x in beta.allele_fields]
    if len(alpha_fields) != 1 or len(beta_fields) != 1:
        return None

    alpha_code = alpha_fields[0].lower()
    beta_code = beta_fields[0].lower()
    if alpha_code != beta_code:
        return None

    return AlleleName(
        species="H-2",
        gene=locus,
        allele_family="",
        allele_code=alpha_code)


def _choose_representative_allele(parsed_result):
    if isinstance(parsed_result, Allele):
        return parsed_result

    representative = getattr(parsed_result, "representative", None)
    if isinstance(representative, Allele):
        return representative

    alleles = getattr(parsed_result, "alleles", None)
    if alleles:
        first = alleles[0]
        if isinstance(first, Allele):
            return first

    raise AlleleParseError("Could not parse '%s'" % (parsed_result.raw_string,))


def _serotype_or_supertype_to_allele_name(parsed_result):
    """
    mhcnames maps serotypes/supertypes (e.g. A2, B7, A02) to a two-field
    allele ending in :01 (e.g. HLA-A*02:01, HLA-B*07:01).
    """
    species = _normalize_species(getattr(parsed_result, "species_prefix", ""))
    gene = str(getattr(parsed_result, "gene_name", "")).upper()
    name = str(getattr(parsed_result, "name", ""))

    if not species or not gene or not name:
        return None

    upper_name = name.upper()
    if not upper_name.startswith(gene):
        return None

    family = upper_name[len(gene):]
    if not family.isdigit():
        return None

    return AlleleName(
        species=species,
        gene=gene,
        allele_family=_normalize_family(family),
        allele_code="01")


def parse_classi_or_classii_allele_name(name, infer_pair=True):
    cache_key = (name, infer_pair)
    if cache_key in _parse_cache:
        return _parse_cache[cache_key]

    if not isinstance(name, str):
        raise TypeError("Expected allele name string but got %s" % (type(name),))

    # Match legacy mhcnames behavior for separators.
    normalized_name = name.strip().replace("/", "-").replace("_", "*")

    try:
        parsed = parse(
            normalized_name,
            infer_class2_pairing=infer_pair,
            raise_on_error=True)
    except ParseError as e:
        raise AlleleParseError(str(e))

    if isinstance(parsed, Pair):
        mouse_shorthand = _mouse_shorthand_from_pair(parsed)
        if mouse_shorthand is not None:
            result = (mouse_shorthand,)
        else:
            if not isinstance(parsed.alpha, Allele) or not isinstance(parsed.beta, Allele):
                raise AlleleParseError("Malformed MHC type %s" % (name,))
            result = (
                _allele_name_from_parsed_allele(parsed.alpha),
                _allele_name_from_parsed_allele(parsed.beta),
            )
    else:
        serotype_or_supertype = _serotype_or_supertype_to_allele_name(parsed)
        if serotype_or_supertype is not None:
            result = (serotype_or_supertype,)
        else:
            representative_allele = _choose_representative_allele(parsed)
            result = (_allele_name_from_parsed_allele(representative_allele),)

    _parse_cache[cache_key] = result
    return result


def normalize_allele_name(raw_allele, omit_dra1=False, infer_class2_pair=True):
    cache_key = (raw_allele, omit_dra1, infer_class2_pair)
    if cache_key in _normalized_allele_cache:
        return _normalized_allele_cache[cache_key]

    parsed_alleles = parse_classi_or_classii_allele_name(
        raw_allele,
        infer_pair=infer_class2_pair)

    species = parsed_alleles[0].species

    # Optionally omit the alpha allele, e.g. for IEDB predictors.
    if omit_dra1 and len(parsed_alleles) == 2:
        alpha, beta = parsed_alleles
        if alpha == _DRA1_0101:
            parsed_alleles = (beta,)

    normalized_parts = [species]
    for parsed_allele in parsed_alleles:
        if len(parsed_allele.allele_family) > 0:
            normalized_parts.append("%s*%s:%s" % (
                parsed_allele.gene,
                parsed_allele.allele_family,
                parsed_allele.allele_code))
        else:
            normalized_parts.append("%s%s" % (
                parsed_allele.gene,
                parsed_allele.allele_code))

    normalized = "-".join(normalized_parts)
    _normalized_allele_cache[cache_key] = normalized
    return normalized
