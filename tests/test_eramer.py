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

"""Tests for the ERAMER ERAP1-trimming predictor.

The trimming-cascade math is tested binary-free with a synthetic PWM. The
end-to-end tests run only when an ERAMER checkout with ``PWM.xlsx`` is available
(``ERAMER_HOME``) and openpyxl is installed.
"""

import os
from os.path import isfile, join

import pytest

from mhctools import ERAMER, Kind
from mhctools.eramer import _specificity, eramer_score


# --- cascade math (synthetic PWM, no xlsx) ----------------------------------

def test_specificity_is_position_mean():
    aa_weights = {"A": [0.2, 0.4, 0.6], "C": [0.0, 0.0, 0.0]}
    assert _specificity("AAA", aa_weights) == 0.4      # (0.2+0.4+0.6)/3
    assert _specificity("CCC", aa_weights) == 0.0
    assert _specificity("AAC", aa_weights) == 0.2      # (0.2+0.4+0.0)/3


def test_eramer_score_cascade_average():
    # Two synthetic precursor lengths: score('AAA' via len-3 sheet)=0.3,
    # score('AA' via len-2 sheet)=0.5.
    wbl = {3: {"A": [0.3, 0.3, 0.3]}, 2: {"A": [0.5, 0.5]}}
    # epitope_length 1 -> trim 'AAA'(l=3) then 'AA'(l=2); mean(0.3, 0.5) = 0.4
    assert eramer_score("AAA", 1, wbl) == pytest.approx(0.4)
    # epitope_length 2 -> only the full 3-mer (l=3): 0.3
    assert eramer_score("AAA", 2, wbl) == pytest.approx(0.3)
    # nothing to trim
    assert eramer_score("AAA", 3, wbl) is None


# --- construction / validation (needs a PWM path but not a real PWM) --------

@pytest.fixture
def stub_predictor(tmp_path):
    # ERAMER validates peptides before it lazily loads the PWM, so a stub file
    # is enough to exercise construction + validation without a real workbook.
    pwm = tmp_path / "PWM.xlsx"
    pwm.write_bytes(b"not a real workbook")
    return ERAMER(pwm_path=str(pwm))


def test_default_kind_and_support(stub_predictor):
    assert stub_predictor._default_pred_kind() == Kind.erap_trimming
    support = stub_predictor.kind_support()[Kind.erap_trimming]
    assert support["mhc_dependence"] == "none"
    assert support["mhc_class"] == "I"
    assert stub_predictor.supported_kinds == (Kind.erap_trimming,)


def test_rejects_bad_epitope_length(tmp_path):
    pwm = tmp_path / "PWM.xlsx"
    pwm.write_bytes(b"x")
    with pytest.raises(ValueError, match="epitope_length"):
        ERAMER(epitope_length=7, pwm_path=str(pwm))
    with pytest.raises(ValueError, match="epitope_length"):
        ERAMER(epitope_length=16, pwm_path=str(pwm))


def test_rejects_missing_pwm(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        ERAMER(pwm_path=str(tmp_path / "nope.xlsx"))


def test_peptide_validation(stub_predictor):
    with pytest.raises(ValueError, match="9-16 residues"):
        stub_predictor.predict(["AAAAAAAA"])          # 8-mer, too short
    with pytest.raises(ValueError, match="9-16 residues"):
        stub_predictor.predict(["A" * 17])            # 17-mer, too long
    with pytest.raises(ValueError, match="non-standard amino acid"):
        stub_predictor.predict(["GGGGGVVVVBX"])       # 11-mer, bad residues


# --- end-to-end (requires an ERAMER checkout + openpyxl) --------------------

ERAMER_HOME = os.environ.get("ERAMER_HOME")
_has_eramer = bool(ERAMER_HOME) and isfile(join(ERAMER_HOME, "PWM.xlsx"))
try:
    import openpyxl  # noqa: F401
    _has_openpyxl = True
except ImportError:
    _has_openpyxl = False

requires_eramer = pytest.mark.skipif(
    not (_has_eramer and _has_openpyxl),
    reason="ERAMER not installed (set ERAMER_HOME to a clone with PWM.xlsx; "
           "needs openpyxl)")


@requires_eramer
def test_eramer_end_to_end_matches_reference():
    # Reference values are ERAMER v1.0's own sample output (README/test.fasta).
    predictor = ERAMER(epitope_length=8, eramer_home=ERAMER_HOME)
    results = predictor.predict(["GGGGGVVVVVVAAAEE", "LLLLLLLLLLLAAAAA"])

    by_peptide = {r.peptide: r.preds[0] for r in results}
    assert by_peptide["GGGGGVVVVVVAAAEE"].kind == Kind.erap_trimming
    assert by_peptide["GGGGGVVVVVVAAAEE"].allele == ""
    assert by_peptide["GGGGGVVVVVVAAAEE"].score == pytest.approx(
        0.085322875, abs=1e-6)
    assert by_peptide["LLLLLLLLLLLAAAAA"].score == pytest.approx(
        0.273310375, abs=1e-6)
    assert results[0].erap_trimming is by_peptide["GGGGGVVVVVVAAAEE"]


@requires_eramer
def test_eramer_epitope_length_changes_cascade():
    # epitope_length 9 drops the shortest (9-mer) trimming step vs default 8.
    e8 = ERAMER(epitope_length=8, eramer_home=ERAMER_HOME)
    e9 = ERAMER(epitope_length=9, eramer_home=ERAMER_HOME)
    (r8,) = e8.predict(["GGGGGVVVVVVAAAEE"])
    (r9,) = e9.predict(["GGGGGVVVVVVAAAEE"])
    assert r8.preds[0].score == pytest.approx(0.085322875, abs=1e-6)
    assert r9.preds[0].score == pytest.approx(0.078594142857, abs=1e-6)
