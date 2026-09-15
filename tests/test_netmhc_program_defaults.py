"""Regression tests for version-specific NetMHC executable defaults."""

import inspect

from mhctools import (
    NetMHC,
    NetMHC3,
    NetMHC4,
    NetMHCIIpan,
    NetMHCIIpan3,
    NetMHCIIpan4,
    NetMHCIIpan4_BA,
    NetMHCIIpan4_EL,
    NetMHCIIpan43,
    NetMHCIIpan43_BA,
    NetMHCIIpan43_EL,
    NetMHCpan,
    NetMHCpan3,
    NetMHCpan4,
    NetMHCpan4_BA,
    NetMHCpan4_EL,
    NetMHCpan28,
    NetMHCpan41,
    NetMHCpan41_BA,
    NetMHCpan41_EL,
    NetMHCpan42,
    NetMHCpan42_BA,
    NetMHCpan42_EL,
)


def _program_default(predictor):
    return inspect.signature(predictor).parameters["program_name"].default


def test_version_specific_wrappers_default_to_versioned_executables():
    expected = {
        NetMHC3: "netMHC-3.4",
        NetMHC4: "netMHC-4.0",
        NetMHCpan28: "netMHCpan-2.8",
        NetMHCpan3: "netMHCpan-3.0",
        NetMHCpan4: "netMHCpan-4.0",
        NetMHCpan4_BA: "netMHCpan-4.0",
        NetMHCpan4_EL: "netMHCpan-4.0",
        NetMHCpan41: "netMHCpan-4.1",
        NetMHCpan41_BA: "netMHCpan-4.1",
        NetMHCpan41_EL: "netMHCpan-4.1",
        NetMHCpan42: "netMHCpan-4.2",
        NetMHCpan42_BA: "netMHCpan-4.2",
        NetMHCpan42_EL: "netMHCpan-4.2",
        NetMHCIIpan3: "netMHCIIpan-3.0",
        NetMHCIIpan4: "netMHCIIpan-4.0",
        NetMHCIIpan4_BA: "netMHCIIpan-4.0",
        NetMHCIIpan4_EL: "netMHCIIpan-4.0",
        NetMHCIIpan43: "netMHCIIpan-4.3",
        NetMHCIIpan43_BA: "netMHCIIpan-4.3",
        NetMHCIIpan43_EL: "netMHCIIpan-4.3",
    }
    assert {_cls: _program_default(_cls) for _cls in expected} == expected


def test_auto_detecting_wrappers_keep_unversioned_executables():
    assert _program_default(NetMHC) == "netMHC"
    assert _program_default(NetMHCpan) == "netMHCpan"
    assert _program_default(NetMHCIIpan) == "netMHCIIpan"

