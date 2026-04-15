# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for mhctools.mhcflurry._normalize_models_path.

Two different strings that point at the same directory should produce
the same cache key so the models aren't loaded twice into memory.
"""

import os

from mhctools.mhcflurry import _normalize_models_path


def test_none_stays_none():
    assert _normalize_models_path(None) is None


def test_expands_tilde(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    assert _normalize_models_path("~/models") == str(tmp_path / "models")


def test_relative_and_absolute_dedupe(tmp_path, monkeypatch):
    target = tmp_path / "models"
    target.mkdir()
    monkeypatch.chdir(tmp_path)
    assert _normalize_models_path("models") == _normalize_models_path(str(target))


def test_symlink_resolves_to_target(tmp_path):
    target = tmp_path / "real_models"
    target.mkdir()
    link = tmp_path / "link_models"
    os.symlink(target, link)
    assert _normalize_models_path(str(link)) == _normalize_models_path(str(target))


def test_dotdot_and_redundant_separators_dedupe(tmp_path):
    target = tmp_path / "models"
    target.mkdir()
    messy = "%s/./sub/../models" % tmp_path
    assert _normalize_models_path(messy) == _normalize_models_path(str(target))
