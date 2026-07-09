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

import os

import pytest

from mhctools.netcleave import NetCleave, NetCleave_I, NetCleave_II
from mhctools.pred import COLUMNS, Kind


# ---------------------------------------------------------------------------
# Constructor / validation — no NetCleave install required.
# ---------------------------------------------------------------------------

def test_init_invalid_class():
    with pytest.raises(ValueError, match="mhc_class"):
        NetCleave(mhc_class="bad")


def test_init_missing_path_raises():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        NetCleave(netcleave_path="/nonexistent/netcleave")


# ---------------------------------------------------------------------------
# Model tests — require a cloned NetCleave (weights) and a Python with
# NetCleave's own deps (tensorflow/keras, scikit-learn, biopython).
# ---------------------------------------------------------------------------

NETCLEAVE_DIR = None
for candidate in [
    os.environ.get("NETCLEAVE_DIR", ""),
    os.path.join(os.path.expanduser("~"), "NetCleave"),
    os.path.join(os.path.expanduser("~"), "code", "NetCleave"),
]:
    if candidate and os.path.isfile(os.path.join(candidate, "NetCleave.py")) \
            and os.path.isfile(os.path.join(
                candidate, "data", "models", "I_mass-spectrometry_HLA",
                "I_mass-spectrometry_HLA_model.h5")):
        NETCLEAVE_DIR = candidate
        break

requires_netcleave = pytest.mark.skipif(
    NETCLEAVE_DIR is None,
    reason="NetCleave not installed (set NETCLEAVE_DIR or clone to ~/NetCleave)")


# Reference C-terminal cleavage scores from NetCleave's OWN CLI
# (`NetCleave.py --predict ... --pred_input 3`) on the pan-allele models,
# with protein_seq = peptide + c_flank. Upstream code, not this wrapper.
# Keys: (peptide, c_flank).
CLASS_I_REF = {
    ("SIINFEKL", "DGH"): 0.842449,
    ("GILGFVFTL", "QRS"): 0.908024,
    ("NLVPMVATV", "KKK"): 0.582904,
    ("LLWNGPMAV", "QRS"): 0.696255,
}
CLASS_II_REF = {
    ("SIINFEKL", "DGH"): 0.360996,
    ("GILGFVFTL", "QRS"): 0.260713,
    ("NLVPMVATV", "KKK"): 0.245104,
    ("LLWNGPMAV", "QRS"): 0.227258,
}


@pytest.fixture(scope="module")
def predictor_I():
    return NetCleave_I(netcleave_path=NETCLEAVE_DIR)


@pytest.fixture(scope="module")
def predictor_II():
    return NetCleave_II(netcleave_path=NETCLEAVE_DIR)


def _predict_ref(predictor, ref):
    pairs = list(ref)
    peptides = [p for p, _ in pairs]
    c_flanks = [c for _, c in pairs]
    return pairs, predictor.predict(peptides, c_flanks=c_flanks)


@requires_netcleave
def test_reproduces_cli_class_I(predictor_I):
    pairs, results = _predict_ref(predictor_I, CLASS_I_REF)
    for (pep, cflank), pp in zip(pairs, results):
        got = pp.cleavage.score
        assert got == pytest.approx(CLASS_I_REF[(pep, cflank)], abs=1e-3)


@requires_netcleave
def test_reproduces_cli_class_II(predictor_II):
    pairs, results = _predict_ref(predictor_II, CLASS_II_REF)
    for (pep, cflank), pp in zip(pairs, results):
        got = pp.endolysosomal_cleavage.score
        assert got == pytest.approx(CLASS_II_REF[(pep, cflank)], abs=1e-3)


@requires_netcleave
def test_class_I_emits_proteasome_kind(predictor_I):
    pp = predictor_I.predict(["SIINFEKL"], c_flanks=["DGH"])[0]
    pred = pp.preds[0]
    assert pred.kind == Kind.proteasome_cleavage
    assert pred.predictor_name == "netcleave"
    assert pred.c_flank == "DGH"
    assert pred.allele == ""
    assert 0.0 <= pred.score <= 1.0


@requires_netcleave
def test_class_II_emits_endolysosomal_kind(predictor_II):
    pp = predictor_II.predict(["SIINFEKL"], c_flanks=["DGH"])[0]
    assert pp.preds[0].kind == Kind.endolysosomal_cleavage


@requires_netcleave
def test_kind_support_metadata(predictor_I, predictor_II):
    assert predictor_I.kind_support() == {
        Kind.proteasome_cleavage: {"mhc_dependence": "none", "mhc_class": "I"}}
    assert predictor_II.kind_support() == {
        Kind.endolysosomal_cleavage: {"mhc_dependence": "none", "mhc_class": "II"}}


@requires_netcleave
def test_class_I_stronger_than_class_II(predictor_I, predictor_II):
    """NetCleave's class-I C-terminal signal is much stronger than class-II
    (paper AUC 0.91 vs 0.66); the pan models reflect this on these peptides."""
    pairs = list(CLASS_I_REF)
    peptides = [p for p, _ in pairs]
    c_flanks = [c for _, c in pairs]
    s1 = [pp.cleavage.score
          for pp in predictor_I.predict(peptides, c_flanks=c_flanks)]
    s2 = [pp.endolysosomal_cleavage.score
          for pp in predictor_II.predict(peptides, c_flanks=c_flanks)]
    for a, b in zip(s1, s2):
        assert a > b


@requires_netcleave
def test_predict_requires_c_flanks(predictor_I):
    with pytest.raises(ValueError, match="c_flanks"):
        predictor_I.predict(["SIINFEKL"])


@requires_netcleave
def test_short_c_flank_yields_empty_result(predictor_I):
    # 2-residue c_flank can't build the 4+3 site -> empty PeptideResult,
    # but the result list stays aligned 1:1 with the input.
    results = predictor_I.predict(["SIINFEKL"], c_flanks=["AB"])
    assert len(results) == 1
    assert results[0].preds == ()


@requires_netcleave
def test_mixed_valid_and_short_flanks_stay_aligned(predictor_I):
    results = predictor_I.predict(
        ["SIINFEKL", "GILGFVFTL"], c_flanks=["AB", "QRS"])
    assert len(results) == 2
    assert results[0].preds == ()             # short flank
    assert results[1].cleavage is not None    # valid


@requires_netcleave
def test_predict_proteins(predictor_I):
    protein = "MASIINFEKLDGHKQRLLWNGPMAVQRSTTT"  # SIINFEKL@2, LLWNGPMAV@16
    results = predictor_I.predict_proteins({"p": protein}, peptide_lengths=[8, 9])
    assert "p" in results
    scored = {pp.cleavage.peptide: pp.cleavage
              for pp in results["p"] if pp.cleavage}
    # peptide scored in real protein context matches the CLI reference
    assert scored["SIINFEKL"].score == pytest.approx(0.842449, abs=1e-3)
    assert scored["SIINFEKL"].offset == 2
    assert scored["SIINFEKL"].source_sequence_name == "p"


@requires_netcleave
def test_predict_proteins_repeated_peptide(predictor_I):
    """Regression: a peptide occurring more than once in the protein must not
    crash (NetCleave emits one row per occurrence). Each occurrence is scored
    in its own downstream context."""
    # SIINFEKL appears at offset 2 (downstream DGH) and offset 16 (downstream QRS)
    protein = "MASIINFEKLDGHKQRSIINFEKLQRSTTT"
    results = predictor_I.predict_proteins({"p": protein}, peptide_lengths=[8])
    hits = [pp.cleavage for pp in results["p"]
            if pp.cleavage and pp.cleavage.peptide == "SIINFEKL"]
    assert len(hits) == 2
    by_offset = {h.offset: h for h in hits}
    assert by_offset[2].c_flank == "DGH"
    assert by_offset[16].c_flank == "QRS"
    # Each occurrence matches predict() with the same downstream flank, and
    # the differing downstream context gives differing scores.
    assert by_offset[2].score == pytest.approx(
        predictor_I.predict(["SIINFEKL"], c_flanks=["DGH"])[0].cleavage.score,
        abs=1e-6)
    assert by_offset[16].score == pytest.approx(
        predictor_I.predict(["SIINFEKL"], c_flanks=["QRS"])[0].cleavage.score,
        abs=1e-6)
    assert by_offset[2].score != by_offset[16].score


@requires_netcleave
def test_predict_recurring_peptide_in_flank(predictor_I):
    """Regression: a peptide that recurs within peptide+flank must not crash."""
    pp = predictor_I.predict(["AAAAAAAA"], c_flanks=["AAA"])[0]
    assert pp.cleavage is not None
    assert 0.0 <= pp.cleavage.score <= 1.0


@requires_netcleave
def test_two_instances_do_not_interfere(predictor_I, predictor_II):
    """Regression: distinct instances must not collide on the shared output
    dir — each keeps its own class-correct score and Kind."""
    a = predictor_I.predict(["SIINFEKL"], c_flanks=["DGH"])[0].cleavage
    b = predictor_II.predict(["SIINFEKL"], c_flanks=["DGH"])[0].endolysosomal_cleavage
    assert a.kind == Kind.proteasome_cleavage
    assert b.kind == Kind.endolysosomal_cleavage
    assert a.score == pytest.approx(CLASS_I_REF[("SIINFEKL", "DGH")], abs=1e-3)
    assert b.score == pytest.approx(CLASS_II_REF[("SIINFEKL", "DGH")], abs=1e-3)


@requires_netcleave
def test_predict_dataframe_schema(predictor_II):
    df = predictor_II.predict_dataframe(["SIINFEKL"], c_flanks=["DGH"])
    assert list(df.columns) == list(COLUMNS)
    assert df["kind"].iloc[0] == Kind.endolysosomal_cleavage
    assert df["c_flank"].iloc[0] == "DGH"


@requires_netcleave
def test_repeated_calls_consistent(predictor_I):
    a = predictor_I.predict(["SIINFEKL"], c_flanks=["DGH"])[0].cleavage.score
    b = predictor_I.predict(["SIINFEKL"], c_flanks=["DGH"])[0].cleavage.score
    assert a == b


@requires_netcleave
def test_predict_empty_returns_empty(predictor_I):
    assert predictor_I.predict([], c_flanks=[]) == []
    assert predictor_I.predict([]) == []   # no c_flanks needed for empty input


@requires_netcleave
def test_subclasses_set_class(predictor_I, predictor_II):
    assert predictor_I.mhc_class == "I"
    assert predictor_II.mhc_class == "II"
    assert isinstance(predictor_I, NetCleave)
    assert isinstance(predictor_II, NetCleave)
