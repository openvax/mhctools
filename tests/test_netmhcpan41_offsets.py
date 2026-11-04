"""Native 4.1b output must map back onto zero-based protein coordinates."""

from pathlib import Path

import pytest

from mhctools import NetMHCpan41, NetMHCpan42
from mhctools.base_commandline_predictor import BaseCommandlinePredictor


SEQUENCE = "ACDEFGHIKLMN"


@pytest.mark.parametrize("predictor_class", [NetMHCpan41, NetMHCpan42])
@pytest.mark.parametrize("mode,suffix", [("binding_affinity", "ba"), ("elution_score", "el")])
def test_native_fixture_offsets_and_scores(monkeypatch, predictor_class, mode, suffix):
    # Parser selection belongs to the wrapper; avoid requiring a binary just
    # to construct it in the offline regression.
    monkeypatch.setattr(BaseCommandlinePredictor, "__init__",
                        lambda self, **kwargs: self.__dict__.update(kwargs))
    predictor = predictor_class(mode=mode)
    output = (Path(__file__).parent / "data" / f"netmhcpan41_offsets_{suffix}.txt").read_text()
    predictions = predictor.parse_output_fn(
        output, prediction_method_name="netMHCpan-4.1",
        sequence_key_mapping={"offset_probe": "original protein"})
    assert [p.offset for p in predictions] == [0, 1, 2, 3]
    assert [p.peptide for p in predictions] == [SEQUENCE[i:i + 9] for i in range(4)]
    assert all(p.source_sequence_name == "original protein" for p in predictions)
    assert all(p.allele == "HLA-A*02:01" for p in predictions)
    first = predictions[0]
    assert first.affinity == (pytest.approx(35192.55) if suffix == "ba" else None)
    assert first.score == pytest.approx(0.032458 if suffix == "ba" else 0.000417)
    assert first.percentile_rank == pytest.approx(56.398 if suffix == "ba" else 24.979)


@pytest.mark.requires_external_tool
@pytest.mark.parametrize("mode", ["binding_affinity", "elution_score"])
def test_real_netmhcpan41_protein_coordinates(mode):
    predictor = NetMHCpan41(alleles=["HLA-A*02:01"], mode=mode)
    predictions = predictor.predict_subsequences({"protein": SEQUENCE}, peptide_lengths=[9])
    assert [p.offset for p in predictions] == [0, 1, 2, 3]
    assert all(p.peptide == SEQUENCE[p.offset:p.offset + p.length] for p in predictions)
