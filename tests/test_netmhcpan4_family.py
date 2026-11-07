import pytest

from mhctools.base_commandline_predictor import BaseCommandlinePredictor

from mhctools import NetMHCpan4, NetMHCpan41, NetMHCpan42
from mhctools.pred import Kind


@pytest.fixture(autouse=True)
def fake_supported_alleles(monkeypatch):
    monkeypatch.setattr(
        BaseCommandlinePredictor,
        "_determine_supported_alleles",
        staticmethod(lambda command, flag: ["HLA-A*02:01"]))


@pytest.mark.parametrize("cls, kinds", [
    (NetMHCpan4, {Kind.pMHC_affinity}),
    (NetMHCpan41, {Kind.pMHC_affinity, Kind.pMHC_presentation}),
    (NetMHCpan42, {Kind.pMHC_affinity, Kind.pMHC_presentation}),
])
def test_kind_support_by_mode(cls, kinds):
    predictor = cls(alleles=["HLA-A*02:01"], program_name="true")
    assert set(predictor.kind_support()) == kinds
    el = cls(alleles=["HLA-A*02:01"], program_name="true", mode="elution_score")
    assert set(el.kind_support()) == {Kind.pMHC_presentation}


@pytest.mark.parametrize("cls", [NetMHCpan4, NetMHCpan41, NetMHCpan42])
def test_mode_sets_ba_flag_and_rejects_unknown(cls):
    assert "-BA" in cls(alleles=["HLA-A*02:01"], program_name="true").extra_flags
    with pytest.raises(ValueError):
        cls(alleles=["HLA-A*02:01"], program_name="true", mode="bogus")
