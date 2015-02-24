# Copyright (c) 2014. Mount Sinai School of Medicine
#
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
import logging
import re

MHC_1_GENE_SET = set(["A", "B", "C", "E", "F", "G", "K", "L"])
MHC_2_GENE_SET = set(["DM", "DO", "DP", "DQ", "DR"])

def seq_to_str(obj):
    """
    Given a sequence convert it to a comma separated string.
    If, however, the argument is a single object, return its string
    representation.
    """
    if isinstance(obj, (unicode, str)):
        return obj
    elif isinstance(obj, (list, tuple)):
        return  ",".join([str(x) for x in obj])
    else:
        return str(obj)

def convert_str(obj):
    """
    Given a string, convert it to an int or float if possible.
    """
    if obj is None:
        return obj
    try:
        try:
            return int(obj)
        except:
            return float(obj)
    except:
        return str(obj)

def _parse_substring(hla, pred, max_len = None):
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
        pos +=1
    return result, hla[pos:]

def _parse_letters(hla, max_len = None):
    return _parse_substring(hla, lambda c: c.isalpha(), max_len  = max_len)

def _parse_numbers(hla, max_len = None):
    return _parse_substring(hla, lambda c: c.isdigit(), max_len = max_len)

def _parse_not_numbers(hla, max_len = None):
    return _parse_substring(hla, lambda c: not c.isdigit(), max_len = max_len)

def normalize_hla_allele_name(hla):
    """
    HLA allele names can look like:
        - HLA-A*03:02
        - HLA-A02:03
        - HLA-A:02:03
        - HLA-A2
        - A2
        - A*03:02
        - A02:02
        - A:02:03
        - A*02:01:01
        - HLA-A*02:01:01G
    ...should all be normalized to:
        HLA-A*03:02
    """
    original = hla
    hla = hla.strip()

    if hla.startswith("HLA-"):
        hla = hla[4:]

    # gene name is sequence of letters at start of HLA string
    gene, hla = _parse_letters(hla)

    assert len(gene) > 0, "No HLA gene name given in %s" % original
    assert len(hla) > 0, "Malformed HLA type %s" % original

    gene = gene.upper()

    # skip initial separator
    sep, hla = _parse_not_numbers(hla)
    assert sep in ("", ":", "*"), \
        "Malformed separator %s in HLA type %s" % (sep, original)

    family, hla = _parse_numbers(hla, max_len = 2)

    sep, hla = _parse_not_numbers(hla)

    assert sep in ("", ":"), \
        "Malformed separator %s in HLA type %s" % (sep, original)

    allele, hla = _parse_numbers(hla)

    #assert len(hla) == 0, \
    #    "Unexpected suffix %s in HLA type %s" % (hla, original)

    if len(family) == 1:
        family = "0" + family
    if len(allele) == 0:
        allele = "01"
    elif len(allele) == 1:
        allele = "0" + allele
    return "HLA-%s*%s:%s" % (gene, family, allele )

def compact_hla_allele_name(hla):
    long_name = normalize_hla_allele_name(hla)
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
