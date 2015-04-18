
SPECIES_PREFIXES = {
    "BoLA",
    "DLA"    # dogs
    "HLA",   # human
    "OVA",   # ovine (sheep)
    "SLA",   # swine (pigs)
    "H-2",   # mice
    "Mamu",  # Rhesus macaques
}

def _parse_hla_allele_name(hla):
    """Takes an allele name and splits it into four parts:
        1) species prefix
        2) gene name
        3) allele family
        4) allele code

    For example, in all of the following inputs:
        "HLA-A*02:01"
        "A0201"
        "A00201"
    The result should be a four-tuple:
        (
            "HLA",  # species prefix
            "A",    # gene name
            "02",   # allele family
            "01",   # allele code
        )
    """
    original = hla
    hla = hla.strip().upper()

    for
    if "H-2" in hla:
        # mouse allele
        species_species = "H-2"
        parts = hla.split("H-2")
    else:
        parts = hla.split("-")

    if len(parts) == 1:
        # no prefix, assume it's a human allele
        hla = parts[0]
    elif len(parts) == 2:
        prefix, hla = parts
    else:
        raise ValueError("Invalid MHC allele: %s" % (hla,))

    if prefix in {"SLA", "BoLA"}:
        # parse gene names like SLA-1, BoLA-2, BoLA-N
        gene, hla = _parse_until(hla, ":")
    elif prefix == "H-2":
        print(hla)
        assert False
    else:
        # gene name is sequence of letters at start of HLA string
        gene, hla = _parse_letters(hla)

    if len(gene) == 0:
        raise ValueError("No HLA gene name given in %s" % original)
    if len(hla) == 0:
        raise ValueError("Malformed HLA type %s" % original)

    gene = gene.upper()

    # skip initial separator
    sep, hla = _parse_not_numbers(hla)
    if sep not in ("", ":", "*", "-"):
        raise ValueError(
            "Malformed separator %s in HLA type %s" % (sep, original))

    family, hla = _parse_numbers(hla, max_len=2)

    sep, hla = _parse_not_numbers(hla)

    if sep not in ("", ":"):
        raise ValueError(
            "Malformed separator %s in HLA type %s" % (sep, original))

    allele, hla = _parse_numbers(hla)

    if len(family) == 1:
        family = "0" + family
    if len(allele) == 0:
        allele = "01"
    elif len(allele) == 1:
        allele = "0" + allele
    return species, gene, family, allele_code

_normalized_allele_cache
def normalize_allele_name(raw_allele):
    """MHC alleles are named with a frustatingly loose system. It's not uncommon
    to see dozens of different forms for the same allele.

    For example, these all refer to the same MHC sequence:
        - HLA-A*02:01
        - HLA-A02:01
        - HLA-A:02:01
        - HLA-A0201
        - HLA-A00201

    Additionally, for human alleles, the species prefix is often omitted:
        - A*02:01
        - A*00201
        - A*0201
        - A02:01
        - A:02:01
        - A:002:01
        - A0201
        - A00201

    We might also encounter "6 digit" and "8 digit" MHC types (which specify
    variants that don't affect amino acid sequence), for our purposes these
    should be truncated to their "4-digit" forms:
        - A*02:01:01
        - A*02:01:01:01
    There are also suffixes which we're going to ignore:
        - HLA-A*02:01:01G

    And lastly, for human alleles, there are serotypes which we'll treat
    as approximately equal to a 4-digit type.
        - HLA-A2
        - A2

    These should all be normalized to:
        HLA-A*02:01
    """
    if raw_allele in _normalized_allele_cache:
        return _normalized_allele_cache[raw_allele]
    (species, gene, family, allele_code) = _parse_allele_name(raw_allele)
    normalized = "%s-%s*%s:%s" % (species, gene, family, allele)
    _normalized_allele_cache[raw_allele] = normalized
    return normalized

MHC_1_GENE_SET = set(["A", "B", "C", "E", "F", "G", "K", "L"])
MHC_2_GENE_SET = set(["DM", "DO", "DP", "DQ", "DR"])

def compact_hla_allele_name(hla):
    long_name = normalize_allele_name(hla)
    # turn HLA-A*02:01 into A0201
    return long_name[4:].replace("*", "").replace(":", "")

def mhc_class_from_normalized_allele_name(normalized_hla):
    """
    Given a normalized HLA allele name, returns 1 or 2 (corresponding to the
    MHC class).

    Returns 1 for: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, HLA-K, HLA-L
    Returns 2 for: HLA-DM, HLA-DO, HLA-DP, HLA-DQ, HLA-DR
    """
    assert normalized_hla.startswith("HLA-") and all(
        delim in normalized_hla for delim in set(["*", ":", "-"])), \
        "Expected normalized HLA allele name, but received %s" % normalized_hla

    gene_end_pos = normalized_hla.index("*")
    gene = normalized_hla[4:gene_end_pos]

    if gene in MHC_1_GENE_SET:
        return 1
    elif gene in MHC_2_GENE_SET:
        return 2
    raise ValueError(
        "HLA gene %s is not a part of the MHC 1 or MHC 2 gene sets." % gene)
    return 1


def _parse_substring(hla, pred, max_len=None):
    """
    Extract substring of letters for which predicate is True
    """
    result = ""
    pos = 0
    if max_len is None:
        max_len = len(hla)
    else:
        max_len = min(max_len, len(hla))
    while pos < max_len and pred(hla[pos]):
        result += hla[pos]
        pos += 1
    return result, hla[pos:]

def _parse_letters(hla, max_len=None):
    return _parse_substring(hla, lambda c: c.isalpha(), max_len=max_len)

def _parse_numbers(hla, max_len=None):
    return _parse_substring(hla, lambda c: c.isdigit(), max_len=max_len)

def _parse_not_numbers(hla, max_len=None):
    return _parse_substring(hla, lambda c: not c.isdigit(), max_len=max_len)

def _parse_until(hla, sep):
    return _parse_substring(hla, lambda c: c != sep)
