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
