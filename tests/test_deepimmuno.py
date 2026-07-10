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

"""Tests for the DeepImmuno class-I immunogenicity wrapper.

Parser and validation tests are binary-free. The end-to-end tests run only when
a DeepImmuno checkout is available (``DEEPIMMUNO_HOME`` pointing at a clone,
with ``DEEPIMMUNO_PYTHON`` optionally naming an interpreter that has TensorFlow
with Keras 2, or newer TensorFlow plus the ``tf-keras`` shim).
"""

import os
from os.path import isfile, join
from tempfile import NamedTemporaryFile

import pytest

from mhctools import DeepImmuno, Kind
from mhctools.deepimmuno import _deepimmuno_allele, parse_deepimmuno_results


# Captured DeepImmuno (multiple-mode) output — tab separated.
_OUTPUT = (
    "peptide\tHLA\timmunogenicity\n"
    "NLVPMVATV\tHLA-A*0201\t0.95676666\n"
    "GILGFVFTL\tHLA-A*0201\t0.8871707\n")


def _write(text):
    f = NamedTemporaryFile("w", suffix="_deepimmuno.txt", delete=False)
    f.write(text)
    f.close()
    return f.name


# --- parser / helpers (binary-free) -----------------------------------------

def test_parse_deepimmuno_results_ordered_scores():
    path = _write(_OUTPUT)
    try:
        scores = parse_deepimmuno_results(path)
    finally:
        os.remove(path)
    assert scores == pytest.approx([0.95676666, 0.8871707])


def test_parse_deepimmuno_rejects_missing_column():
    path = _write("peptide\tHLA\nNLVPMVATV\tHLA-A*0201\n")
    try:
        with pytest.raises(ValueError, match="immunogenicity"):
            parse_deepimmuno_results(path)
    finally:
        os.remove(path)


def test_deepimmuno_allele_format():
    # mhctools canonical HLA-A*02:01 -> DeepImmuno's colon-free HLA-A*0201.
    assert _deepimmuno_allele("HLA-A*02:01") == "HLA-A*0201"
    assert _deepimmuno_allele("A0201") == "HLA-A*0201"
    assert _deepimmuno_allele("HLA-B*35:01") == "HLA-B*3501"


# --- construction / validation (no binary needed) ---------------------------

def _stub_home(tmp_path):
    """A minimal fake DeepImmuno checkout (just a deepimmuno-cnn.py marker)."""
    (tmp_path / "deepimmuno-cnn.py").write_text("# stub\n")
    return str(tmp_path)


def test_default_kind_and_support(tmp_path):
    predictor = DeepImmuno(
        alleles=["HLA-A*02:01"], deepimmuno_home=_stub_home(tmp_path))
    assert predictor._default_pred_kind() == Kind.immunogenicity
    support = predictor.kind_support()[Kind.immunogenicity]
    assert support["mhc_dependence"] == "single_allele"
    assert support["mhc_class"] == "I"
    assert predictor.supported_kinds == (Kind.immunogenicity,)


def test_rejects_missing_home(tmp_path):
    # An empty directory is not a DeepImmuno checkout.
    with pytest.raises(FileNotFoundError, match="deepimmuno-cnn.py"):
        DeepImmuno(alleles=["HLA-A*02:01"], deepimmuno_home=str(tmp_path))


def test_rejects_no_alleles(tmp_path):
    with pytest.raises(ValueError, match="at least one allele"):
        DeepImmuno(alleles=[], deepimmuno_home=_stub_home(tmp_path))


def test_peptide_length_validation(tmp_path):
    predictor = DeepImmuno(
        alleles=["HLA-A*02:01"], deepimmuno_home=_stub_home(tmp_path))
    # DeepImmuno only scores 9- and 10-mers; an 8mer is rejected before any
    # subprocess runs.
    with pytest.raises(ValueError, match="9- and 10-mers"):
        predictor.predict(["SIINFEK"])
    with pytest.raises(ValueError, match="9- and 10-mers"):
        predictor.predict(["A" * 11])
    # Empty input short-circuits to an empty result list.
    assert predictor.predict([]) == []


# --- end-to-end (requires a DeepImmuno checkout) ----------------------------

DEEPIMMUNO_HOME = os.environ.get("DEEPIMMUNO_HOME")
_has_deepimmuno = bool(DEEPIMMUNO_HOME) and isfile(
    join(DEEPIMMUNO_HOME, "deepimmuno-cnn.py"))

requires_deepimmuno = pytest.mark.skipif(
    not _has_deepimmuno,
    reason="DeepImmuno not installed (set DEEPIMMUNO_HOME to a clone; "
           "optionally DEEPIMMUNO_PYTHON to an interpreter with TensorFlow "
           "and Keras 2 / tf-keras)")


@requires_deepimmuno
def test_deepimmuno_end_to_end():
    predictor = DeepImmuno(
        alleles=["HLA-A*02:01"], deepimmuno_home=DEEPIMMUNO_HOME)
    peptides = ["NLVPMVATV", "GILGFVFTL"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.immunogenicity
        assert pred.allele == "HLA-A*02:01"
        assert pred.predictor_name == "deepimmuno"
        assert 0.0 <= pred.score <= 1.0
        assert result.immunogenicity is pred

    # DeepImmuno is deterministic; these match a local reference run.
    by_peptide = {r.peptide: r.preds[0].score for r in results}
    assert by_peptide["NLVPMVATV"] == pytest.approx(0.9568, abs=1e-3)
    assert by_peptide["GILGFVFTL"] == pytest.approx(0.8872, abs=1e-3)


@requires_deepimmuno
def test_deepimmuno_multiple_alleles_per_peptide():
    predictor = DeepImmuno(
        alleles=["HLA-A*02:01", "HLA-B*07:02"], deepimmuno_home=DEEPIMMUNO_HOME)
    results = predictor.predict(["NLVPMVATV"])
    assert len(results) == 1
    preds = results[0].preds
    assert len(preds) == 2
    assert {p.allele for p in preds} == {"HLA-A*02:01", "HLA-B*07:02"}
