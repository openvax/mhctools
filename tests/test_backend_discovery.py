# Copyright (c) 2026 Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Production and integration tests must share optional-backend discovery."""

from pathlib import Path

import pytest

from mhctools.caphla import _find_caphla_home
from mhctools.deepimmuno import _find_deepimmuno_home
from mhctools.deeptap import _find_deeptap_home
from mhctools.eramer import _find_pwm_path
from mhctools.nettcr import _find_nettcr_dir
from mhctools.tlimmuno2 import _find_tlimmuno2_home


@pytest.mark.parametrize(
    "resolver,env_var,directory_name,marker,returns_marker",
    [
        (_find_caphla_home, "CAPHLA_HOME", "CapHLA", "EL_model.py", False),
        (
            _find_deepimmuno_home,
            "DEEPIMMUNO_HOME",
            "DeepImmuno",
            "deepimmuno-cnn.py",
            False,
        ),
        (_find_deeptap_home, "DEEPTAP_HOME", "DeepTAP", "deeptap.py", False),
        (_find_pwm_path, "ERAMER_HOME", "ERAMER", "PWM.xlsx", True),
        (
            _find_tlimmuno2_home,
            "TLIMMUNO2_HOME",
            "TLimmuno2",
            "Python/TLimmuno2.py",
            False,
        ),
    ],
)
def test_resolver_finds_common_code_checkout(
        monkeypatch, tmp_path, resolver, env_var, directory_name, marker,
        returns_marker):
    monkeypatch.delenv(env_var, raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    checkout = tmp_path / "code" / directory_name
    marker_path = checkout / marker
    marker_path.parent.mkdir(parents=True)
    marker_path.touch()
    expected = marker_path if returns_marker else checkout
    assert Path(resolver()) == expected


def test_nettcr_resolver_accepts_lowercase_common_checkout(monkeypatch, tmp_path):
    monkeypatch.delenv("NETTCR_DIR", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    checkout = tmp_path / "code" / "nettcr"
    checkout.mkdir(parents=True)
    assert Path(_find_nettcr_dir()) == checkout
