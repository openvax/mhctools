# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for BasePredictor._check_results, including the large-input path
that previously silently accepted duplicate/missing predictions."""

import pytest

from mhctools import RandomBindingPredictor
from mhctools.binding_prediction import BindingPrediction


ALLELES = ["HLA-A*02:01", "HLA-B*07:02"]
PEPTIDES = ["SIINFEKLA", "ASILLLVFY", "KKLLLFFAA"]


def _bp(allele, peptide):
    return BindingPrediction(peptide=peptide, allele=allele, affinity=100.0)


def _predictor():
    # RandomBindingPredictor subclasses BasePredictor and doesn't need
    # a network or external binary to instantiate.
    return RandomBindingPredictor(alleles=ALLELES)


def test_check_results_accepts_complete():
    p = _predictor()
    preds = [_bp(a, pep) for a in ALLELES for pep in PEPTIDES]
    p._check_results(preds, PEPTIDES, ALLELES)


def test_check_results_detects_missing():
    p = _predictor()
    preds = [_bp(a, pep) for a in ALLELES for pep in PEPTIDES]
    preds.pop()  # drop one
    with pytest.raises(ValueError, match="Missing"):
        p._check_results(preds, PEPTIDES, ALLELES)


def test_check_results_detects_unexpected_peptide():
    p = _predictor()
    preds = [_bp(a, pep) for a in ALLELES for pep in PEPTIDES]
    preds.append(_bp(ALLELES[0], "ZZZZZZZZZ"))
    with pytest.raises(ValueError, match="Unexpected"):
        p._check_results(preds, PEPTIDES, ALLELES)


def test_check_results_detects_unexpected_allele():
    p = _predictor()
    preds = [_bp(a, pep) for a in ALLELES for pep in PEPTIDES]
    preds.append(_bp("HLA-C*07:01", PEPTIDES[0]))
    with pytest.raises(ValueError, match="Unexpected"):
        p._check_results(preds, PEPTIDES, ALLELES)


def test_check_results_detects_duplicate_and_missing_at_large_scale():
    """At sizes >100k the previous fast path only warned on count mismatch
    and silently accepted duplicates that masked missing pairs. This test
    ensures the correct behavior: duplicate-plus-missing raises."""
    # 500 alleles * 500 peptides = 250,000 expected pairs (>100k threshold)
    alleles = ["A%04d" % i for i in range(500)]
    peptides = ["PEPT%05d" % i for i in range(500)]
    preds = [_bp(a, pep) for a in alleles for pep in peptides]
    # Replace one prediction with a duplicate of another — count still
    # matches n_expected but one (allele, peptide) is missing.
    preds[-1] = _bp(alleles[0], peptides[0])
    p = RandomBindingPredictor(alleles=alleles[:1])  # alleles arg doesn't matter here
    with pytest.raises(ValueError, match="Missing"):
        p._check_results(preds, peptides, alleles)


def test_check_results_accepts_complete_at_large_scale():
    alleles = ["A%04d" % i for i in range(400)]
    peptides = ["PEPT%05d" % i for i in range(300)]  # 120k pairs
    preds = [_bp(a, pep) for a in alleles for pep in peptides]
    p = RandomBindingPredictor(alleles=alleles[:1])
    p._check_results(preds, peptides, alleles)
