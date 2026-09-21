"""Real local inference against the official pinned SMM/SMM-PMBEC models."""

import csv
import os
from pathlib import Path
import shutil

import pytest

from mhctools import SMM, SMMPMBEC


ROOT = Path(__file__).parent / "data/osteosarc"
PROGRAM = os.environ.get("IEDB_MHCI_EXECUTABLE", "iedb-mhci")
pytestmark = pytest.mark.skipif(
    not shutil.which(PROGRAM), reason="Install local SMM: see docs/testing.md")


@pytest.mark.parametrize("cls,method", [(SMM, "smm"), (SMMPMBEC, "smmpmbec")])
def test_local_vaccine_predictions_match_recording(cls, method):
    peptides = (ROOT / "inputs/class_i.txt").read_text().splitlines()
    predictor = cls(alleles=["HLA-A*01:01", "HLA-B*08:01"])
    predictions = predictor.predict_peptides(peptides)
    expected = {(r["peptide"], r["allele"]): (float(r["ic50"]), float(r["rank"]))
                for length in (9, 10) for r in csv.DictReader(
                    (ROOT / "smm" / ("%s-%d.tsv" % (method, length))).read_text().splitlines(),
                    delimiter="\t")}
    assert len(predictions) == len(expected) == 38
    for prediction in predictions:
        assert (prediction.affinity, prediction.percentile_rank) == pytest.approx(
            expected[prediction.peptide, prediction.allele])


def test_unsupported_allele_length_is_an_error():
    with pytest.raises(ValueError):
        SMM(alleles=["HLA-A*24:01"]).predict_peptides(["SIINFEKLA"])
