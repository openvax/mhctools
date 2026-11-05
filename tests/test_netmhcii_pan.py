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

from .common import eq_

from mhctools import NetMHCIIpan
from mhctools.allele_normalization import normalize_allele_name

import pytest


pytestmark = pytest.mark.requires_external_tool

def test_netmhcii_pan_DRB():
    alleles = [normalize_allele_name("HLA-DRB1*01:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DRB1*01:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_alpha_beta():
    alleles = [normalize_allele_name("HLA-DPA1*01:05/DPB1*100:01")]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }

    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])
    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"HLA-DPA1*01:05/DPB1*100:01"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_multiple_alleles():
    alleles = [
        normalize_allele_name("HLA-DPA1*01:05/DPB1*100:01"),
        normalize_allele_name("HLA-DQA1*05:11/DQB1*03:02"),
        normalize_allele_name("HLA-DRA*01:01/DRB1*01:01")
    ]
    ii_pan_predictor = NetMHCIIpan(
        alleles=alleles)
    fasta_dictionary = {
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {
        "HLA-DPA1*01:05/DPB1*100:01",
        "HLA-DQA1*05:11/DQB1*03:02",
        "HLA-DRA*01:01/DRB1*01:01"
    })

    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect 3 * ((21-15+1) + (21-16+1)) = 39 entries
    assert len(binding_predictions) == 39, \
        "Expected 39 epitopes from %s" % (binding_predictions,)

def test_netmhcii_pan_mouse():
    alleles = [normalize_allele_name("H2-IAb")]
    ii_pan_predictor = NetMHCIIpan(alleles=alleles)
    fasta_dictionary = {
        "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
        "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
    }
    binding_predictions = ii_pan_predictor.predict_subsequences(
        sequence_dict=fasta_dictionary,
        peptide_lengths=[15, 16])

    unique_lengths = {x.length for x in binding_predictions}
    eq_(unique_lengths, {15, 16})

    unique_alleles = {x.allele for x in binding_predictions}
    eq_(unique_alleles, {"H2-AA*b/AB*b"})

    # length of "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT" is 34
    # length of "SQAMDDLMLSPDDIEQWFTED" is 21
    # Expect (34-15+1) + (34-16+1) + (21-15+1) + (21-16+1) = 52 entries
    assert len(binding_predictions) == 52, \
        "Expected 52 epitopes from %s" % (binding_predictions,)

# Tests will also be run using the following program names if they are installed.
# In any case, a program called "netMHCIIpan" MUST be installed and working for
# this test suite to succeeed.
OPTIONAL_NETMHCIIPAN_PROGRAM_NAMES = [
    "netMHCIIpan-4.0",
]

DEFAULT_ALLELE = "HLA-DPA1*01:05/DPB1*100:01"

protein_sequence_dict = {
    "SMAD4-001": "PAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGT",
    "TP53-001": "SQAMDDLMLSPDDIEQWFTED"
}

# follows similar model to test_netmhc_pan.py
def test_netmhcii_pan():
    check_netmhcii_pan("netMHCIIpan", True)  # required

    for program_name in OPTIONAL_NETMHCIIPAN_PROGRAM_NAMES:
        check_netmhcii_pan(program_name, False)  # optional

# TODO: add more actual tests to this
def check_netmhcii_pan(program_name, fail_if_no_such_program=True):
    try:
        NetMHCIIpan(
            alleles=[DEFAULT_ALLELE], program_name=program_name)
    except FileNotFoundError:
        if fail_if_no_such_program:
            raise
        print("Skipping because no such program: %s" % program_name)
        return


def test_beta_and_explicit_pair_keep_distinct_identities_in_one_native_batch():
    """One DR model input must produce both explicit request identities."""
    from mhctools import NetMHCIIpan
    beta = 'HLA-DRB1*01:01'
    pair = 'HLA-DRA*01:01/DRB1*01:01'
    predictor = NetMHCIIpan(alleles=[beta, pair], program_name='netMHCIIpan-4.3')
    peptides = ['YARFASHNGALIQRR']
    results = predictor.predict(peptides)
    assert {p.allele for r in results for p in r.preds} == {beta, pair}
    binding = list(predictor.predict_peptides(peptides))
    assert {p.allele for p in binding} == {beta, pair}
    assert len(binding) == 2
    assert binding[0].score == binding[1].score
    provenance = predictor.allele_provenance()
    assert provenance['native_alleles'] == {beta: 'DRB1_0101', pair: 'DRB1_0101'}
