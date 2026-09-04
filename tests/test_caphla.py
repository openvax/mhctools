# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

import csv
import os
from pathlib import Path

import pandas as pd
import pytest

from mhctools.caphla import (
    CAPHLA_REVISION,
    CapHLA,
    CapHLA_BA,
    CapHLA_EL,
    _affinity_nm,
    normalize_caphla_allele,
)
from mhctools.pred import COLUMNS, Kind


def test_artifact_revision_matches_wrapper():
    from mhctools.artifacts import _SNAPSHOTS
    assert _SNAPSHOTS["caphla"].revision == CAPHLA_REVISION


@pytest.fixture
def caphla_home(tmp_path):
    for filename in ("BA_model.py", "EL_model.py"):
        (tmp_path / filename).touch()
    with (tmp_path / "HLA_library.csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(("Allele Name", "MHC pseudo-seq"))
        writer.writerow(("HLA-A*02:01", "A" * 34))
        writer.writerow(("HLA-DRB1*15:01", "D" * 34))
        writer.writerow(("HLA-DPA1*01:03/DPB1*04:01", "P" * 34))
        writer.writerow(("HLA-A*30:14L", "L" * 34))
        writer.writerow(("HLA-DRB5*01:08N", "N" * 34))
        writer.writerow(("H2-Db", "M" * 34))
    return tmp_path


def test_allele_normalization_class_i():
    assert normalize_caphla_allele("A0201") == (
        "HLA-A*02:01", "HLA-A*02:01", "I")


def test_allele_normalization_class_ii_dr():
    assert normalize_caphla_allele("HLA-DRB1*15:01") == (
        "HLA-DRA1*01:01-DRB1*15:01", "HLA-DRB1*15:01", "II")


def test_allele_normalization_class_ii_heterodimer():
    assert normalize_caphla_allele(
        "HLA-DPA1*01:03-DPB1*04:01") == (
            "HLA-DPA1*01:03-DPB1*04:01",
            "HLA-DPA1*01:03/DPB1*04:01",
            "II")


def test_allele_normalization_mouse():
    assert normalize_caphla_allele("H-2-Db") == (
        "H-2-Db", "H2-Db", "I")


def test_allele_normalization_expression_variants():
    assert normalize_caphla_allele("HLA-A*30:14L") == (
        "HLA-A*30:14L", "HLA-A*30:14L", "I")
    assert normalize_caphla_allele("HLA-DRB5*01:08N") == (
        "HLA-DRA1*01:01-DRB5*01:08N", "HLA-DRB5*01:08N", "II")


def test_allele_normalization_nonclassical_mouse_aliases():
    assert normalize_caphla_allele("H2-Qa1") == (
        "H-2-T23", "H2-Qa1", "I")
    assert normalize_caphla_allele("H2-Qa2") == (
        "H-2-Q7", "H2-Qa2", "I")


def test_init_and_kind_support(caphla_home):
    predictor = CapHLA(
        ["HLA-A*02:01", "HLA-DRB1*15:01"],
        caphla_path=caphla_home)
    assert predictor.alleles == [
        "HLA-A*02:01", "HLA-DRA1*01:01-DRB1*15:01"]
    assert predictor._models == {}
    assert predictor.kind_support() == {
        Kind.pMHC_presentation: {
            "mhc_dependence": "single_allele", "mhc_class": "both"},
        Kind.pMHC_affinity: {
            "mhc_dependence": "single_allele", "mhc_class": "both"},
    }


def test_mode_subclasses(caphla_home):
    el = CapHLA_EL("HLA-A*02:01", caphla_path=caphla_home)
    ba = CapHLA_BA("HLA-A*02:01", caphla_path=caphla_home)
    assert tuple(el.kind_support()) == (Kind.pMHC_presentation,)
    assert tuple(ba.kind_support()) == (Kind.pMHC_affinity,)


def test_invalid_init(caphla_home):
    with pytest.raises(ValueError, match="at least one allele"):
        CapHLA([], caphla_path=caphla_home)
    with pytest.raises(ValueError, match="mode"):
        CapHLA("HLA-A*02:01", mode="bad", caphla_path=caphla_home)
    with pytest.raises(ValueError, match="batch_size"):
        CapHLA("HLA-A*02:01", batch_size=0, caphla_path=caphla_home)
    with pytest.raises(ValueError, match="unsupported"):
        CapHLA("HLA-A*99:99", caphla_path=caphla_home)


class _FakeCapHLA(CapHLA):
    def _predict_raw(self, peptides, allele_keys):
        self.calls = (list(peptides), list(allele_keys))
        count = len(peptides)
        outputs = {}
        if self.mode in ("el", "both"):
            outputs["el"] = [0.1 + 0.1 * i for i in range(count)]
        if self.mode in ("ba", "both"):
            outputs["ba"] = [0.2 + 0.1 * i for i in range(count)]
        return outputs


def test_predict_pairs_preserves_order_duplicates_and_outputs(caphla_home):
    predictor = _FakeCapHLA("HLA-A*02:01", caphla_path=caphla_home)
    results = predictor.predict_pairs([
        ("siinfekl", "A0201"),
        ("SIINFEKL", "HLA-A*02:01"),
    ])
    assert predictor.calls == (
        ["SIINFEKL", "SIINFEKL"],
        ["HLA-A*02:01", "HLA-A*02:01"],
    )
    assert len(results) == 2
    assert [prediction.kind for prediction in results[0].preds] == [
        Kind.pMHC_presentation, Kind.pMHC_affinity]
    affinity = results[0].affinity
    assert affinity.score == pytest.approx(0.2)
    assert affinity.value == pytest.approx(_affinity_nm(0.2))
    assert affinity.percentile_rank is None
    assert affinity.predictor_version == CAPHLA_REVISION


def test_predict_cross_product(caphla_home):
    predictor = _FakeCapHLA(
        ["HLA-A*02:01", "HLA-DRB1*15:01"],
        mode="el",
        caphla_path=caphla_home)
    results = predictor.predict(["SIINFEKL", "GILGFVFTL"])
    assert len(results) == 2
    assert [prediction.allele for prediction in results[0].preds] == [
        "HLA-A*02:01", "HLA-DRA1*01:01-DRB1*15:01"]
    assert predictor.calls[0] == [
        "SIINFEKL", "SIINFEKL", "GILGFVFTL", "GILGFVFTL"]


def test_empty_and_invalid_peptides(caphla_home):
    predictor = _FakeCapHLA("HLA-A*02:01", caphla_path=caphla_home)
    assert predictor.predict_pairs([]) == []
    with pytest.raises(ValueError, match="lengths 7-25"):
        predictor.predict_pairs([("SHORT", "HLA-A*02:01")])
    with pytest.raises(ValueError, match="invalid amino acid"):
        predictor.predict_pairs([("SIINFEKZ", "HLA-A*02:01")])


def test_predict_pairs_dataframe(caphla_home):
    predictor = _FakeCapHLA("HLA-A*02:01", caphla_path=caphla_home)
    dataframe = predictor.predict_pairs_dataframe(
        [("SIINFEKL", "HLA-A*02:01")], sample_name="sample")
    assert list(dataframe.columns) == list(COLUMNS)
    assert len(dataframe) == 2
    assert dataframe["sample_name"].tolist() == ["sample", "sample"]


def _installed_caphla_home():
    candidates = [
        os.environ.get("CAPHLA_HOME", ""),
        str(Path("~/CapHLA").expanduser()),
    ]
    for candidate in candidates:
        if candidate and (Path(candidate) / "test_out.csv").is_file():
            return candidate
    from mhctools.artifacts import artifact_status
    status = artifact_status("caphla")
    if status.status == "ready":
        return status.path
    return ""


CAPHLA_HOME = _installed_caphla_home()
requires_caphla = pytest.mark.skipif(
    not CAPHLA_HOME,
    reason="CapHLA not installed (run `mhctools fetch caphla`)")


@requires_caphla
def test_all_upstream_alleles_normalize_to_library_keys():
    library = pd.read_csv(Path(CAPHLA_HOME) / "HLA_library.csv")
    keys = set(library["Allele Name"])
    classes = []
    for allele in library["Allele Name"]:
        _, key, mhc_class = normalize_caphla_allele(allele)
        assert key in keys
        classes.append(mhc_class)
    assert set(classes) == {"I", "II"}


@requires_caphla
def test_reproduces_official_el_and_ba_fixture():
    pytest.importorskip("torch")
    inputs = pd.read_csv(Path(CAPHLA_HOME) / "test.csv", header=None)
    expected = pd.read_csv(Path(CAPHLA_HOME) / "test_out.csv")
    predictor = CapHLA(
        expected["Allele Name"].tolist(), caphla_path=CAPHLA_HOME)
    results = predictor.predict_pairs(inputs[0], inputs[1])
    presentation = [result.presentation.score for result in results]
    affinity = [result.affinity.score for result in results]
    assert presentation == pytest.approx(
        expected["presentation_score"].tolist(), abs=3e-7)
    assert affinity == pytest.approx(
        expected["affinity_score"].tolist(), abs=3e-7)
