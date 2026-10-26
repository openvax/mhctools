"""Public window metadata and native scanning must describe the same inputs."""

import pytest

from mhctools import Kind, PeptideResult, Prediction
from mhctools.base_predictor import BasePredictor


class ContextPredictor(BasePredictor):
    uses_flanking_sequences = True

    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        return [PeptideResult(preds=(Prediction(
            kind=Kind.pMHC_affinity, peptide=peptide, allele=self.alleles[0],
            score=0.5, value=100., predictor_name="context-test",
            n_flank=n, c_flank=c),))
            for peptide, n, c in zip(peptides, n_flanks, c_flanks)]


@pytest.mark.parametrize("n_width,c_width", [(2, 3), (0, 3), (2, 0), (5, 5)])
def test_public_window_limits_drive_native_scan_at_and_between_boundaries(n_width, c_width):
    model = ContextPredictor(alleles=["HLA-A*02:01"], default_peptide_lengths=[9])
    model.n_flank_length, model.c_flank_length = n_width, c_width
    assert model.protein_flank_lengths() == (n_width, c_width)
    proteins = {"short": "MASIINFEKLAGG", "long": "WWMASIINFEKLAGGWW"}
    result = model.predict_proteins_dataframe(proteins)
    for row in result.to_dict("records"):
        sequence, offset = proteins[row["source_sequence_name"]], row["offset"]
        assert row["n_flank"] == sequence[max(0, offset - n_width):offset]
        assert row["c_flank"] == sequence[offset + 9:offset + 9 + c_width]
    assert model.protein_flank_lengths() == (n_width, c_width)


def test_peptide_only_predictors_declare_zero_windows():
    model = BasePredictor(alleles=["HLA-A*02:01"])
    model.flank_length = 7
    assert model.protein_flank_lengths() == (0, 0)


def test_public_override_is_the_native_scanners_contract():
    class LoadedModel(ContextPredictor):
        def protein_flank_lengths(self):
            return 1, 4
    model = LoadedModel(alleles=["HLA-A*02:01"], default_peptide_lengths=[9])
    row = model.predict_proteins_dataframe({"p": "MASIINFEKLAGGWW"}).iloc[2]
    assert row.n_flank == "A"
    assert row.c_flank == "GGWW"
