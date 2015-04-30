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

from __future__ import print_function, division, absolute_import
from collections import namedtuple

# copied from https://www.ebi.ac.uk/ipd/mhc/species.html
SPECIES_PREFIXES = dict(
    human="HLA",
    cattle="BoLA",
    bison="Bibi",
    dog="DLA",
    sheep=["OVA", "Ovar", "Ovca"],
    swine="SLA",
    mouse="H2",
    rainbow_trout="Onmy",
    rat=["Rano", "Rara", "RT1"],
    salmon="Sasa",
    cat="FLA",
    horse=["ELA", "Eqca"],
    chimp=["Patr", "ChLA"],
    bonobo="Papa",
    white_handed_gibbon="Hyla",
    gorilla="Gogo",
    orangutan=["Popy", "OrLA"],
    blue_monkey="Cemi",
    de_brazzas_monkey="Cene",
    vervet_monkey="Chae",
    mantled_colobus="Cogu",
    black_mangabey="Loat",
    stump_tailed_macaque="Maar",
    crab_eating_macaque="Mafa",
    japanese_macaque="Mafu",
    rhesus_macaque=["Mamu", "RhLA"],
    pig_tailed_macaque="Mane",
    lion_tailed_macaque="Masi",
    drill="Male",
    mandrill="Masp",
    olive_baboon="Paan",
    yellow_baboon="Pacy",
    hamadryas_baboon="Paha",
    guinea_baboon="Papp",
    chacma_baboon="Paur",
    entelus_langur="Pren",
    gelada_baboon="Thge",
    owl_monkey=["Aoaz", "Aovo"],
    northern_night_owl_monkey=["Aona", "Aoni", "OmLA"],
    long_haired_spider_monkey="Atbe",
    brown_headed_spider_monkey="Atfu",
    marmoset=["Caja", "MaLA"],
    pygmy_marmoset="Cepy",
    dusk_titi_monkey="Camo",
    tufted_capuchin="Ceap",
    golden_lion_tamarin="Lero",
    white_faced_saki="Pipi",
    saddle_backed_tamarin="Safu",
    red_crested_tamarin="Sage",
    moustached_tamarin="Samy",
    cotton_top_tamarin="Saoe",
    squirrel_monkey="Sasc",
    lemur="Leca")

AlleleName = namedtuple("AlleleName", [
    "species",
    "gene",
    "allele_family",
    "allele_code"])

def _parse_mouse_allele_name(name):
    """Parses mouse MHc alleles such as H2-Kd, H-2-Db, H2-IAb.
    Returns pair of (gene, allele_code).
    """
    original = name

    if name.startswith("H2"):
        name = name[2:]
    elif name.startswith("H-2"):
        name = name[3:]

    _, name = _parse_separator(name)

    # special logic for mouse alleles
    if name.startswith("I"):
        # class II mouse allele
        if len(name) < 3:
            raise ValueError("Incomplete mouse MHC allele: %s" % original)
        elif len(name) > 3:
            raise ValueError("Malformed mouse MHC allele: %s" % original)
        # mice don't seem to have allele families, only a small list of
        # alleles per gene code as single lowercase letters
        return name[:2].upper(), name[2].lower()

    else:
        # class I mouse allele
        if len(name) < 2:
            raise ValueError("Incomplete mouse MHC allele: %s" % original)
        elif len(name) > 2:
            raise ValueError("Malformed mouse MHC allele: %s" % original)
        return name[0].upper(), name[1].lower()

def parse_allele_name(name):
    """Takes an allele name and splits it into four parts:
        1) species prefix
        2) gene name
        3) allele family
        4) allele code

    For example, in all of the following inputs:
        "HLA-A*02:01"
        "A0201"
        "A00201"
    The result is a AlleleName object. Example:
        AlleleName(
            species="HLA",  # species prefix
            gene="A",    # gene name
            allele_family="02",   # allele family
            allele_code="01",   # allele code
        )

    The logic for other species mostly resembles the naming system for humans,
    except for mice, rats, and swine, which have archaic nomenclature.
    """
    original = name
    name = name.strip()

    if len(name) == 0:
        raise ValueError("Can't normalize empty MHC allele name")

    if name.startswith("H2") or name.startswith("H-2"):
        gene, allele_code = _parse_mouse_allele_name(name)
        # mice don't have allele families
        return AlleleName("H2", gene, "", allele_code)
    species = None
    for species_list in SPECIES_PREFIXES.values():
        if isinstance(species_list, str):
            species_list = [species_list]
        for curr_species in species_list:
            prefix = name[:len(curr_species) + 1].upper()
            if prefix == (curr_species.upper() + "-"):
                species = curr_species
                name = name[len(curr_species) + 1:]
                break

    if len(name) == 0:
        raise ValueError("Incomplete MHC allele name: %s" % (original,))

    elif not species:
        # assume that a missing species name means we're dealing with a
        # human HLA allele
        if "-" in name:
            raise ValueError("Can't parse allele name: %s" % original)
        species = "HLA"

    if name[0] == "D":
        # MHC class II genes like "DQA1" need to be parsed with both
        # letters and numbers
        gene, name = _parse_alphanum(name)
    elif name[0].isalpha():
        # if there are more separators to come, then assume the gene names
        # can have the form "DQA1"
        gene, name = _parse_letters(name)
    elif name[0].isdigit():
        gene, name = _parse_numbers(name)
    else:
        raise ValueError("Can't parse gene name from allele: %s" % original)

    if len(gene) == 0:
        raise ValueError("No MHC gene name given in %s" % original)
    if len(name) == 0:
        raise ValueError("Malformed MHC type %s" % original)

    gene = gene.upper()

    # skip initial separator
    sep, name = _parse_separator(name)
    # if all that's left is e.g. "0201" then only parse the
    # first two digits as the family code
    if len(name) == 4:
        family, name = _parse_numbers(name, max_len=2)
    else:
        family, name = _parse_numbers(name, max_len=3)

    sep, name = _parse_separator(name)

    allele_code, name = _parse_numbers(name)

    if len(family) == 1:
        family = "0" + family
    elif len(family) == 3 and family[0] == "0":
        family = family[1:]

    if len(allele_code) == 0:
        allele_code = "01"
    elif len(allele_code) == 1:
        # change HLA-A*2:01 into HLA-A*02:01
        allele_code = "0" + allele_code
    elif len(allele_code) == 3 and allele_code[0] == "0":
        # normalize HLA-A*002:01 into HLA-A*02:01
        allele_code = allele_code[1:]
    return AlleleName(species, gene, family, allele_code)

_normalized_allele_cache = {}
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
    parsed_allele = parse_allele_name(raw_allele)
    if len(parsed_allele.allele_family) > 0:
        normalized = "%s-%s*%s:%s" % (
            parsed_allele.species,
            parsed_allele.gene,
            parsed_allele.allele_family,
            parsed_allele.allele_code)
    else:
        # mice don't have allele families
        # e.g. H2-Kd
        # species = H2
        # gene = K
        # allele = d
        normalized = "%s-%s%s" % (
            parsed_allele.species,
            parsed_allele.gene,
            parsed_allele.allele_code)
    _normalized_allele_cache[raw_allele] = normalized
    return normalized

def compact_allele_name(allele):
    long_name = normalize_allele_name(allele)
    # turn HLA-A*02:01 into A0201
    return long_name.split("-")[1].replace("*", "").replace(":", "")

def _parse_substring(allele, pred, max_len=None):
    """
    Extract substring of letters for which predicate is True
    """
    result = ""
    pos = 0
    if max_len is None:
        max_len = len(allele)
    else:
        max_len = min(max_len, len(allele))
    while pos < max_len and pred(allele[pos]):
        result += allele[pos]
        pos += 1
    return result, allele[pos:]

def _parse_alphanum(allele, max_len=None):
    return _parse_substring(allele, lambda c: c.isalnum(), max_len=max_len)

def _parse_letters(allele, max_len=None):
    return _parse_substring(allele, lambda c: c.isalpha(), max_len=max_len)

def _parse_numbers(allele, max_len=None):
    return _parse_substring(allele, lambda c: c.isdigit(), max_len=max_len)

def _parse_not_numbers(allele, max_len=None):
    return _parse_substring(allele, lambda c: not c.isdigit(), max_len=max_len)

def _parse_until(allele, sep):
    return _parse_substring(allele, lambda c: c != sep)

SEPARATORS = {":", "*", "-"}
def _parse_separator(allele, max_len=None):
    return _parse_substring(allele, lambda c: c in SEPARATORS, max_len=max_len)
