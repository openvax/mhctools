# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the PRIME immunogenicity wrapper.

The parser tests are binary-free (they parse a captured PRIME output file).
The end-to-end tests near the bottom run only when a PRIME executable is
available (``PRIME_EXECUTABLE`` env var or ``PRIME`` on ``PATH``); PRIME also
needs MixMHCpred, whose path may be given via ``MIXMHCPRED_PATH``.
"""

import os
import shutil
from tempfile import NamedTemporaryFile

import pytest

from mhctools import PRIME
from mhctools.prime import parse_prime_results
from mhctools.pred import Kind


# A captured PRIME (v2.1) output for alleles A0201, B0702 on three peptides.
_PRIME_OUTPUT = """\
####################
# Output from PRIME (v2.1)
# Alleles: A0201, B0702
####################
Peptide\t%Rank_bestAllele\tScore_bestAllele\t%RankBinding_bestAllele\tBestAllele\t%Rank_A0201\tScore_A0201\t%RankBinding_A0201\t%Rank_B0702\tScore_B0702\t%RankBinding_B0702
SIINFEKL\t40.396\t0.001067\t10.000\tA0201\t40.396\t0.001067\t10.000\t80.982\t0.000332\t52.000
GILGFVFTL\t0.011\t0.236432\t0.100\tA0201\t0.011\t0.236432\t0.100\t2.556\t0.013242\t16.000
NLVPMVATV\t0.087\t0.128786\t0.200\tA0201\t0.087\t0.128786\t0.200\t12.663\t0.003802\t32.000
"""


def _write_output():
    f = NamedTemporaryFile("w", suffix="_prime.txt", delete=False)
    f.write(_PRIME_OUTPUT)
    f.close()
    return f.name


# --- parser (binary-free) ---------------------------------------------------

def test_parse_prime_results_shape():
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    path = _write_output()
    try:
        by_peptide = parse_prime_results(path, alleles)
    finally:
        os.remove(path)

    assert set(by_peptide) == {"SIINFEKL", "GILGFVFTL", "NLVPMVATV"}
    for peptide, preds in by_peptide.items():
        assert len(preds) == 2, "one prediction per allele"
        assert {p.allele for p in preds} == set(alleles)
        for p in preds:
            assert p.kind == Kind.immunogenicity
            assert p.peptide == peptide
            assert p.predictor_name == "prime"


def test_parse_prime_results_values_and_position_mapping():
    # Position mapping: first triple -> alleles[0], second -> alleles[1].
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    path = _write_output()
    try:
        by_peptide = parse_prime_results(path, alleles)
    finally:
        os.remove(path)

    by_allele = {p.allele: p for p in by_peptide["GILGFVFTL"]}
    a = by_allele["HLA-A*02:01"]
    b = by_allele["HLA-B*07:02"]
    assert a.score == 0.236432 and a.percentile_rank == 0.011
    assert b.score == 0.013242 and b.percentile_rank == 2.556


def test_parse_prime_results_allele_identity_is_caller_supplied():
    # The returned allele is the normalized string we passed, not PRIME's
    # short column form (A0201).
    alleles = ["HLA-A*02:01", "HLA-B*07:02"]
    path = _write_output()
    try:
        by_peptide = parse_prime_results(path, alleles)
    finally:
        os.remove(path)
    alleles_seen = {p.allele for preds in by_peptide.values() for p in preds}
    assert alleles_seen == {"HLA-A*02:01", "HLA-B*07:02"}
    assert "A0201" not in alleles_seen


def test_parse_prime_results_allele_count_mismatch_raises():
    # Two alleles' worth of columns but three alleles claimed.
    path = _write_output()
    try:
        with pytest.raises(ValueError, match="per-allele columns"):
            parse_prime_results(path, ["A", "B", "C"])
    finally:
        os.remove(path)


def test_prime_default_kind_and_support():
    predictor = PRIME(alleles=["HLA-A*02:01"], program_name="PRIME")
    assert predictor._default_pred_kind() == Kind.immunogenicity
    support = predictor.kind_support()[Kind.immunogenicity]
    assert support["mhc_dependence"] == "single_allele"
    assert support["mhc_class"] == "I"


# --- end-to-end (requires a PRIME install) ----------------------------------

PRIME_EXECUTABLE = os.environ.get("PRIME_EXECUTABLE") or shutil.which("PRIME")
MIXMHCPRED_PATH = os.environ.get("MIXMHCPRED_PATH") or shutil.which("MixMHCpred")

# PRIME shells out to MixMHCpred, so both must be present; otherwise skip
# (rather than fail) so CI without MixMHCpred stays green.
requires_prime = pytest.mark.skipif(
    not (PRIME_EXECUTABLE and MIXMHCPRED_PATH),
    reason="PRIME + MixMHCpred not both installed (set PRIME_EXECUTABLE and "
           "MIXMHCPRED_PATH, or put PRIME and MixMHCpred on PATH)")


@requires_prime
def test_prime_end_to_end():
    predictor = PRIME(
        alleles=["HLA-A*02:01", "HLA-B*07:02"],
        program_name=PRIME_EXECUTABLE,
        mixmhcpred_path=MIXMHCPRED_PATH)
    peptides = ["GILGFVFTL", "NLVPMVATV", "SIINFEKL"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 2, "one immunogenicity pred per allele"
        for pred in result.preds:
            assert pred.kind == Kind.immunogenicity
            assert pred.allele in {"HLA-A*02:01", "HLA-B*07:02"}
            assert pred.score is not None
            assert pred.percentile_rank is not None
        # a best-immunogenicity accessor resolves
        assert result.immunogenicity is not None


@requires_prime
def test_prime_dataframe_end_to_end():
    predictor = PRIME(
        alleles=["HLA-A*02:01"],
        program_name=PRIME_EXECUTABLE,
        mixmhcpred_path=MIXMHCPRED_PATH)
    df = predictor.predict_dataframe(["GILGFVFTL", "SIINFEKL"])
    assert len(df) == 2
    assert set(df["kind"]) == {Kind.immunogenicity}
    assert (df["allele"] == "HLA-A*02:01").all()


def test_predict_dataframe_inherits_flank_signature():
    # Regression: PRIME used to override predict_dataframe with a narrower
    # signature that dropped n_flanks/c_flanks (a TypeError trap). It should
    # inherit the base method, which forwards flanks to predict().
    import inspect
    params = inspect.signature(PRIME.predict_dataframe).parameters
    assert "n_flanks" in params and "c_flanks" in params
