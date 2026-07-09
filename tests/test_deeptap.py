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

"""Tests for the DeepTAP TAP-transport wrapper.

Parser and validation tests are binary-free. The end-to-end tests run only when
a DeepTAP checkout is available (``DEEPTAP_HOME`` pointing at a clone, with
``DEEPTAP_PYTHON`` optionally naming an interpreter that has torch +
pytorch-lightning).
"""

import os
from os.path import isfile, join
from tempfile import NamedTemporaryFile

import pytest

from mhctools import DeepTAP, Kind
from mhctools.deeptap import parse_deeptap_results
from mhctools.pred import VALUE_BEST_DIRECTIONS, best_direction


# Captured DeepTAP (classification) output, matching the committed demo.
_CLA_OUTPUT = (
    "peptide,pred_score,pred_label\n"
    "KADDDKPGA,0.2964,0\n"
    "AEASAAAAY,0.9795,1\n"
    "ALAKAGAAV,0.3405,0\n")

# Captured DeepTAP (regression) output — note the extra pred_affinity column.
_REG_OUTPUT = (
    "peptide,pred_score,pred_affinity,pred_label\n"
    "KADDDKPGA,0.2527,101416.8532,1\n"
    "AEASAAAAY,0.5676,788.6656,1\n"
    "ALAKAGAAV,0.3499,22637.9888,1\n")


def _write(text):
    f = NamedTemporaryFile("w", suffix="_deeptap.csv", delete=False)
    f.write(text)
    f.close()
    return f.name


# --- parser (binary-free) ---------------------------------------------------

def test_parse_deeptap_cla():
    path = _write(_CLA_OUTPUT)
    try:
        by_peptide = parse_deeptap_results(path, "cla")
    finally:
        os.remove(path)

    assert set(by_peptide) == {"KADDDKPGA", "AEASAAAAY", "ALAKAGAAV"}
    pred = by_peptide["AEASAAAAY"][0]
    assert pred.kind == Kind.tap_transport
    assert pred.predictor_name == "deeptap_cla"
    assert pred.peptide == "AEASAAAAY"
    assert pred.allele == ""
    assert pred.score == 0.9795
    assert pred.value is None


def test_parse_deeptap_reg_fills_value():
    path = _write(_REG_OUTPUT)
    try:
        by_peptide = parse_deeptap_results(path, "reg")
    finally:
        os.remove(path)

    pred = by_peptide["AEASAAAAY"][0]
    assert pred.kind == Kind.tap_transport
    assert pred.predictor_name == "deeptap_reg"
    assert pred.score == 0.5676
    assert pred.value == 788.6656  # predicted affinity, nM


def test_parse_deeptap_reg_requires_affinity_column():
    # A cla-shaped file has no pred_affinity, so reg parsing must reject it.
    path = _write(_CLA_OUTPUT)
    try:
        with pytest.raises(ValueError, match="pred_affinity"):
            parse_deeptap_results(path, "reg")
    finally:
        os.remove(path)


def test_tap_transport_value_direction():
    # Regression-mode affinity is nM: lower is stronger binding.
    assert VALUE_BEST_DIRECTIONS[Kind.tap_transport] == "min"
    assert best_direction(Kind.tap_transport, "value") == "min"
    assert best_direction(Kind.tap_transport, "score") == "max"


# --- construction / validation (no binary needed) ---------------------------

def _stub_home(tmp_path):
    """A minimal fake DeepTAP checkout (just a deeptap.py marker file)."""
    (tmp_path / "deeptap.py").write_text("# stub\n")
    return str(tmp_path)


def test_default_kind_and_support(tmp_path):
    predictor = DeepTAP(deeptap_home=_stub_home(tmp_path))
    assert predictor._default_pred_kind() == Kind.tap_transport
    support = predictor.kind_support()[Kind.tap_transport]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "none"
    assert predictor.supported_kinds == (Kind.tap_transport,)


def test_rejects_bad_task_type(tmp_path):
    with pytest.raises(ValueError, match="task_type"):
        DeepTAP(task_type="bogus", deeptap_home=_stub_home(tmp_path))


def test_rejects_oversized_max_length(tmp_path):
    with pytest.raises(ValueError, match="pads peptides"):
        DeepTAP(deeptap_home=_stub_home(tmp_path), max_peptide_length=25)


def test_rejects_missing_home(tmp_path):
    # An empty directory is not a DeepTAP checkout.
    with pytest.raises(FileNotFoundError, match="deeptap.py"):
        DeepTAP(deeptap_home=str(tmp_path))


def test_peptide_validation(tmp_path):
    predictor = DeepTAP(deeptap_home=_stub_home(tmp_path))
    with pytest.raises(ValueError, match="up to 17"):
        predictor.predict(["A" * 18])
    with pytest.raises(ValueError, match="cannot encode"):
        predictor.predict(["SIINFEKL1"])
    with pytest.raises(ValueError, match="Empty peptide"):
        predictor.predict([""])


# --- end-to-end (requires a DeepTAP checkout) -------------------------------

DEEPTAP_HOME = os.environ.get("DEEPTAP_HOME")
_has_deeptap = bool(DEEPTAP_HOME) and isfile(join(DEEPTAP_HOME, "deeptap.py"))

requires_deeptap = pytest.mark.skipif(
    not _has_deeptap,
    reason="DeepTAP not installed (set DEEPTAP_HOME to a clone; optionally "
           "DEEPTAP_PYTHON to an interpreter with torch + pytorch-lightning)")


@requires_deeptap
def test_deeptap_end_to_end_cla():
    predictor = DeepTAP(task_type="cla", deeptap_home=DEEPTAP_HOME)
    peptides = ["KADDDKPGA", "AEASAAAAY", "ALAKAGAAV"]
    results = predictor.predict(peptides)

    assert len(results) == len(peptides)
    for result, peptide in zip(results, peptides):
        assert result.peptide == peptide
        assert len(result.preds) == 1
        pred = result.preds[0]
        assert pred.kind == Kind.tap_transport
        assert pred.allele == ""
        assert 0.0 <= pred.score <= 1.0
        assert result.tap_transport is pred

    # DeepTAP is deterministic; these match the committed demo scores.
    by_peptide = {r.peptide: r.preds[0].score for r in results}
    assert by_peptide["AEASAAAAY"] == pytest.approx(0.9795, abs=1e-3)
    assert by_peptide["KADDDKPGA"] == pytest.approx(0.2964, abs=1e-3)


@requires_deeptap
def test_deeptap_end_to_end_reg_sets_value():
    predictor = DeepTAP(task_type="reg", deeptap_home=DEEPTAP_HOME)
    results = predictor.predict(["AEASAAAAY"])
    pred = results[0].preds[0]
    assert pred.kind == Kind.tap_transport
    assert pred.value is not None and pred.value > 0  # affinity in nM
    assert pred.predictor_name == "deeptap_reg"
