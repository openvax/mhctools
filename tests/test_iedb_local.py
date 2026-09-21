"""Historical CLI names select local executables and affinity endpoints."""

from pathlib import Path
import socket
import subprocess

import pytest

from mhctools import base_commandline_predictor as commandline
from mhctools.allele_normalization import normalize_allele_name
from mhctools.cli.args import make_mhc_arg_parser, predictors_from_args


@pytest.mark.parametrize("name,executable,allele,raw,peptide,affinity,rank", [
    ("netmhcpan-iedb", "netMHCpan-4.1", "HLA-A*01:01", "netmhcpan.stdout",
     "AAPVATPAL", 32539.90, 35.431),
    ("netmhciipan-iedb", "netMHCIIpan-4.3", "HLA-DRB1*03:01", "netmhciipan.stdout",
     "AAKAVKPKVVKPKKA", 6998.92, 51.30),
    ("netmhccons-iedb", "netMHCcons", "HLA-A*01:01", "netmhccons-a0101.stdout",
     "AAPVATPAL", 26123.53, 50.00),
])
def test_historical_names_use_local_affinity_parser(
        name, executable, allele, raw, peptide, affinity, rank, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Recorded migration regression must stay offline")

    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket.socket, "connect_ex", forbidden)
    monkeypatch.setattr(commandline.BaseCommandlinePredictor, "_determine_supported_alleles",
                        lambda *args: [allele])
    monkeypatch.setattr(commandline, "run_command", lambda *args: "local help")
    args = make_mhc_arg_parser().parse_args([
        "--mhc-predictor", name, "--mhc-alleles", allele])
    predictor, = predictors_from_args(args)
    assert predictor.program_name == executable
    if name != "netmhccons-iedb":
        assert "-BA" in predictor.extra_flags
    text = (Path(__file__).parent / "data/osteosarc/outputs" / raw).read_text()
    predictions = predictor.parse_output_fn(text)
    anchor = next(p for p in predictions if p.peptide == peptide and p.allele == normalize_allele_name(allele))
    assert anchor.affinity == affinity
    assert anchor.percentile_rank == rank


@pytest.mark.parametrize("name", [
    "netmhcpan-iedb", "netmhccons-iedb", "smm-iedb", "smm-pmbec-iedb",
])
@pytest.mark.parametrize("api", ["predict_subsequences", "predict_proteins"])
def test_legacy_protein_scans_preserve_all_default_lengths(name, api, monkeypatch):
    import sys
    from mhctools import BindingPrediction, BindingPredictionCollection

    allele = "HLA-A*02:01"
    monkeypatch.setattr(commandline.BaseCommandlinePredictor, "_determine_supported_alleles",
                        lambda *args: [allele])
    monkeypatch.setattr(commandline, "run_command", lambda *args: "local help")
    args = make_mhc_arg_parser().parse_args([
        "--mhc-predictor", name, "--mhc-alleles", allele,
        "--mhc-predictor-path", sys.executable])
    predictor, = predictors_from_args(args)

    def score(peptides, **kwargs):
        # Synthetic scores isolate window enumeration from model inference.
        return BindingPredictionCollection([
            BindingPrediction(peptide=p, allele=allele, affinity=100) for p in peptides])

    monkeypatch.setattr(predictor, "predict_peptides", score)
    if isinstance(predictor, commandline.BaseCommandlinePredictor):
        monkeypatch.setattr(predictor, "_predict_binding_predictions_for_alleles", score)

    protein = "ACDEFGHIKLMN"
    scan = getattr(predictor, api)
    for lengths in (None, [10]):
        results = scan({"protein": protein}, peptide_lengths=lengths)
        predictions = results if api == "predict_subsequences" else [
            p for result in results["protein"] for p in result.preds]
        expected = {(protein[i:i + n], i, "protein")
                    for n in (lengths if lengths is not None else [8, 9, 10, 11])
                    for i in range(len(protein) - n + 1)}
        assert len(predictions) == len(expected)
        assert {(p.peptide, p.offset, p.source_sequence_name) for p in predictions} == expected


def test_legacy_netmhcpan_predict_and_pairs_emit_only_native_affinity(monkeypatch):
    from mhctools import IedbNetMHCpan
    from mhctools.binding_prediction_collection import BindingPredictionCollection
    from mhctools.pred import Kind

    allele = "HLA-A*01:01"
    monkeypatch.setattr(commandline.BaseCommandlinePredictor, "_determine_supported_alleles",
                        lambda *args: [allele])
    predictor = IedbNetMHCpan(alleles=[allele])
    text = (Path(__file__).parent / "data/osteosarc/outputs/netmhcpan.stdout").read_text()

    def replay(commands, input_filenames, temp_dir_list, sequence_key_mapping=None):
        return BindingPredictionCollection([
            p for p in predictor.parse_output_fn(text)
            if p.peptide == "AAPVATPAL" and p.allele == allele])

    monkeypatch.setattr(predictor, "_build_peptide_commands", lambda **kwargs: ({}, [], []))
    monkeypatch.setattr(predictor, "_run_commands_and_collect_predictions", replay)
    # A regression to the native multi-endpoint path must fail before any process starts.
    def forbidden(*args, **kwargs):
        raise AssertionError("Legacy interface selected the multi-endpoint runner")
    monkeypatch.setattr(predictor, "_run_commands_and_collect_preds", forbidden)
    for results in (predictor.predict(["AAPVATPAL"]),
                    predictor.predict_pairs([("AAPVATPAL", allele)])):
        predictions = [p for result in results for p in result.preds]
        assert len(predictions) == 1
        assert predictions[0].kind == Kind.pMHC_affinity
        assert predictions[0].value == 32539.90
        assert predictions[0].percentile_rank == 35.431
    assert set(predictor.kind_support()) == {Kind.pMHC_affinity}
