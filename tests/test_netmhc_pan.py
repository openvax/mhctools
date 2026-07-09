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

from mhctools import NetMHCpan


DEFAULT_ALLELE = 'HLA-A*02:01'

protein_sequence_dict = {
    "SMAD4-001": "ASIINFKELA",
    "TP53-001": "ASILLLVFYW"
}

# Tests will also be run using the following program names if they are installed.
# In any case, a program called "netMHCpan" MUST be installed and working for
# this test suite to succeeed.
OPTIONAL_NETMHCPAN_PROGRAM_NAMES = [
    "netMHCpan-2.8",
    "netMHCpan-3.0",
    "netMHCpan-4.0",
    "netMHCpan-4.1",
]


def test_netmhc_pan():
    check_netmhc_pan("netMHCpan", True)  # required

    for program_name in OPTIONAL_NETMHCPAN_PROGRAM_NAMES:
        check_netmhc_pan(program_name, False)  # optional


def check_netmhc_pan(program_name, fail_if_no_such_program=True):
    try:
        predictor = NetMHCpan(
            alleles=[DEFAULT_ALLELE], program_name=program_name)
    except FileNotFoundError:
        if fail_if_no_such_program:
            raise
        print("Skipping because no such program: %s" % program_name)
        return

    binding_predictions = predictor.predict_subsequences(
        protein_sequence_dict,
        peptide_lengths=[9])
    assert len(binding_predictions) == 4, \
        "Expected 4 binding predictions from %s" % (binding_predictions,)
    for x in binding_predictions:
        # recompute the peptide from the offset and starting sequence, and make sure it matches.
        # this is currently wrong in netMHCpan-3.0 and we want to test our wrapper fix to that
        offset = x.offset
        length = x.length
        seq_name = x.source_sequence_name
        expected_peptide = protein_sequence_dict[seq_name][offset:offset + length]
        eq_(expected_peptide, x.peptide,
            "Peptide mismatch: expected %s but got %s in binding prediction '%s'" % (
                expected_peptide, x.peptide, x,))


def test_netmhc_pan_multiple_lengths():
    predictor = NetMHCpan(alleles=["A6801"])
    binding_predictions = predictor.predict_peptides(
        ["A" * 8, "A" * 9, "A" * 10, "A" * 11])
    assert len(binding_predictions) == 4, \
        "Expected 4 epitopes from %s" % (binding_predictions,)

def test_netmhc_pan_multiple_alleles():
    alleles = 'A*02:01,B*35:02'
    predictor = NetMHCpan(
        alleles=alleles,
        default_peptide_lengths=[9])
    sequence_dict = {
        "SMAD4-001": "ASIINFKELA",
        "TP53-001": "ASILLLVFYW"
    }
    binding_predictions = predictor.predict_subsequences(
        sequence_dict=sequence_dict)
    assert len(binding_predictions) == 8, \
        "Expected 8 binding predictions from %s" % (binding_predictions,)
    # With allele batching both alleles are produced by a single netMHCpan
    # invocation (-a A02:01,B35:02); make sure both are attributed correctly.
    observed_alleles = {bp.allele for bp in binding_predictions}
    assert observed_alleles == {"HLA-A*02:01", "HLA-B*35:02"}, \
        "Expected both alleles, got %s" % (observed_alleles,)


def test_netmhc_pan_exotic_unnormalizable_alleles():
    """netMHCpan lists non-human alleles (H-2-Qa1, BoLA-amani.1) and HLA
    low/null-expression variants (HLA-A30:14L) that mhcgnomes can't normalize
    (issue #220). They must round-trip: request them, run, and get predictions
    back keyed by the requested identity. netMHCpan echoes HLA-A30:14L back
    with a '*' (HLA-A*30:14L), so both spellings must resolve to one identity.
    """
    cases = [
        ("H-2-Qa1", "H-2-Qa1"),
        ("BoLA-amani.1", "BoLA-amani.1"),
        ("HLA-A30:14L", "HLA-A30:14L"),   # requested without '*'
        ("HLA-A*30:14L", "HLA-A30:14L"),  # requested with '*', same identity
    ]
    peptides = ["SIINFEKLL", "GILGFVFTL"]
    for requested, expected_identity in cases:
        predictor = NetMHCpan(alleles=[requested])
        binding_predictions = predictor.predict_peptides(peptides)
        assert len(binding_predictions) == len(peptides), \
            "Expected %d predictions for %s, got %s" % (
                len(peptides), requested, binding_predictions)
        observed = {bp.allele for bp in binding_predictions}
        assert observed == {expected_identity}, \
            "Expected identity %r for requested %r, got %s" % (
                expected_identity, requested, observed)


def test_netmhc_pan_batched_matches_per_allele():
    """Batching alleles into one `-a A,B,C` invocation must produce exactly
    the same scores as running one process per allele. Uses the real binary's
    own output on both paths as non-circular ground truth."""
    alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02", "HLA-A*01:01"]
    peptides = [
        "SIINFEKLL", "GILGFVFTL", "NLVPMVATV", "LLWNGPMAV", "AAAWYLWEV"]

    def scores(max_alleles_per_command):
        predictor = NetMHCpan(
            alleles=alleles,
            max_alleles_per_command=max_alleles_per_command)
        return {
            (bp.allele, bp.peptide): bp.value
            for bp in predictor.predict_peptides(peptides)
        }

    batched = scores(None)     # all alleles in a single process
    per_allele = scores(1)     # one process per allele
    assert set(batched) == set(per_allele), (
        "Batched and per-allele runs produced different (allele, peptide) "
        "pairs")
    assert len(batched) == len(alleles) * len(peptides)
    for key, value in per_allele.items():
        eq_(value, batched[key],
            "Score mismatch for %s: per-allele=%s batched=%s" % (
                key, value, batched[key]))
