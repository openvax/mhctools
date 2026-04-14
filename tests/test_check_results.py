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


def test_check_results_error_example_is_deterministic():
    """The 'example' missing pair in the error message should follow the
    caller's input order, not set iteration order."""
    alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"]
    peptides = ["PEP1AAAAA", "PEP2BBBBB", "PEP3CCCCC"]
    # Drop exactly one pair: (alleles[1], peptides[2])
    preds = [
        _bp(a, pep)
        for a in alleles for pep in peptides
        if not (a == alleles[1] and pep == peptides[2])
    ]
    p = RandomBindingPredictor(alleles=alleles)
    with pytest.raises(ValueError) as exc_info:
        p._check_results(preds, peptides, alleles)
    # Must name the exact missing pair, not an arbitrary one
    assert "peptide='PEP3CCCCC'" in str(exc_info.value)
    assert "allele='HLA-B*07:02'" in str(exc_info.value)


def test_check_results_detects_duplicate_masking_missing_pair():
    """Regression: the previous >100k fast path only warned on count
    mismatch, silently accepting duplicates that masked missing pairs.
    A duplicate replacing a missing pair now raises."""
    alleles = ["A%03d" % i for i in range(50)]
    peptides = ["PEPT%04d" % i for i in range(50)]
    preds = [_bp(a, pep) for a in alleles for pep in peptides]
    # Replace the last prediction with a duplicate of the first: count
    # still matches n_expected but one (allele, peptide) pair is missing.
    preds[-1] = _bp(alleles[0], peptides[0])
    with pytest.raises(ValueError, match="Missing"):
        _predictor()._check_results(preds, peptides, alleles)


def test_check_results_accepts_empty_inputs():
    """n_expected = 0 with no predictions is trivially valid."""
    _predictor()._check_results([], [], [])
