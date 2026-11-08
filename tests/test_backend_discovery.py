# Copyright (c) 2026 Mount Sinai School of Medicine
# Licensed under the Apache License, Version 2.0.

"""Inventory and constructors must agree on the selected installation."""

from pathlib import Path
import shutil

import pytest

import mhctools
from mhctools import artifacts
from mhctools.artifacts import artifact_status, resolve_artifact_path
from tests.snapshot_fixtures import populate_snapshot


# Constructor, override argument, selected-path attribute, required inputs.
WRAPPERS = {
    "bigmhc": ("BigMHC", "bigmhc_path", "bigmhc_dir", {"alleles": ["A0201"]}),
    "caphla": ("CapHLA", "caphla_path", "caphla_home", {"alleles": ["A0201"]}),
    "cleavenet": ("CleaveNet", "cleavenet_home", "home", {}),
    "deepimmuno": ("DeepImmuno", "deepimmuno_home", "deepimmuno_home", {"alleles": ["A0201"]}),
    "deeptap": ("DeepTAP", "deeptap_home", "deeptap_home", {}),
    "eramer": ("ERAMER", "eramer_home", "pwm_path", {}),
    "mixtcrpred": ("MixTCRpred", "mixtcrpred_path", "mixtcrpred_home", {"model": "A0201_GILGFVFTL"}),
    "netcleave": ("NetCleave", "netcleave_path", "netcleave_dir", {}),
    "nettcr": ("NetTCR", "nettcr_path", "nettcr_dir", {}),
    "tlimmuno2": ("TLimmuno2", "tlimmuno2_home", "tlimmuno2_home", {"alleles": ["DRB1_0803"]}),
    "tulip": ("Tulip", "tulip_home", "tulip_home", {}),
}


@pytest.fixture
def isolated_snapshots(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setenv("MHCTOOLS_DATA_DIR", str(tmp_path / "managed"))
    for snapshot in artifacts._SNAPSHOTS.values():
        monkeypatch.delenv(snapshot.environment_variable, raising=False)
    monkeypatch.delenv("ERAMER_PWM", raising=False)
    return home


def complete_checkout(root, name):
    populate_snapshot(root, name)
    if name == "caphla":
        (root / "HLA_library.csv").write_text(
            "Allele Name,MHC pseudo-seq\nHLA-A*02:01," + "A" * 34 + "\n")
    elif name == "mixtcrpred":
        (root / "pretrained_models/info_models.csv").write_text(
            "MixTCRpred_model_name,Peptide,MHC,MHC_class,Host_species,Origin,Number_training_abTCR,AUC_5fold\n"
            "A0201_GILGFVFTL,GILGFVFTL,HLA-A*02:01,MHCI,HomoSapiens,Influenza,2079,0.94\n")
    elif name == "tulip":
        (root / "configs/shallow.config.json").write_text("{}")
        (root / "model_weights/pytorch_model.bin").touch()
    return root


@pytest.fixture(params=sorted(WRAPPERS))
def installation(request, isolated_snapshots, monkeypatch):
    name = request.param
    snapshot = artifacts._SNAPSHOTS[name]
    target = artifacts.managed_path(name)
    populate_snapshot(target, name, managed=True)
    complete_checkout(target, name)
    if name == "cleavenet":
        from mhctools import cleavenet
        from mhctools.optional_backend import sha256_file
        # Keep checksum verification while using a tiny fixture asset.
        monkeypatch.setattr(cleavenet, "_ARTIFACTS", {
            "LICENSE": ("fixture", sha256_file(target / "LICENSE"), "text")})
    return name, snapshot, target, isolated_snapshots


def construct(name, **overrides):
    constructor, argument, attribute, inputs = WRAPPERS[name]
    instance = getattr(mhctools, constructor)(**inputs, **overrides)
    return Path(getattr(instance, attribute))


def expected_path(name, root):
    return root / "PWM.xlsx" if name == "eramer" else root


def test_inventory_covers_every_snapshot_wrapper():
    assert set(WRAPPERS) == set(artifacts._SNAPSHOTS)


def test_constructors_use_managed_snapshot(installation):
    name, snapshot, target, home = installation
    assert construct(name) == expected_path(name, target)


@pytest.mark.parametrize("location", ["", "code"])
@pytest.mark.parametrize("complete", [False, True])
def test_partial_home_never_shadows_complete_install(installation, location, complete):
    name, snapshot, target, home = installation
    for directory in snapshot.checkout_names or ("cleavenet",):
        checkout = home / location / directory
        checkout.mkdir(parents=True, exist_ok=True)
        if complete:
            complete_checkout(checkout, name)
        else:
            marker = checkout / snapshot.required_paths[0]
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.touch()
        expected = (checkout if snapshot.checkout_names and (complete or name == "eramer")
                    else target)
        status = artifact_status(name)
        assert status.status == "ready"
        assert Path(status.path) == expected
        assert construct(name) == expected_path(name, expected)
        shutil.rmtree(checkout)


@pytest.mark.parametrize("override", ["argument", "environment"])
def test_explicit_selection_and_invalid_override_do_not_fall_back(installation, monkeypatch, override):
    name, snapshot, target, home = installation
    explicit = home / "custom"
    argument = WRAPPERS[name][1]
    kwargs = {argument: explicit} if override == "argument" else {}
    monkeypatch.setenv(snapshot.environment_variable, str(explicit))
    with pytest.raises(FileNotFoundError, match="incomplete"):
        construct(name, **kwargs)
    status = artifact_status(name)
    assert status.status == "missing"
    assert status.path == str(explicit)
    complete_checkout(explicit, name)
    if override == "argument":
        monkeypatch.setenv(snapshot.environment_variable, str(home / "absent"))
    assert construct(name, **kwargs) == expected_path(name, explicit)


@pytest.mark.parametrize("name,pattern", [("nettcr", "*.tflite"), ("caphla", "ba_fold*.params")])
def test_explicit_install_requires_model_files(isolated_snapshots, name, pattern):
    root = complete_checkout(isolated_snapshots / "custom", name)
    for path in root.rglob(pattern):
        path.unlink()
    with pytest.raises(FileNotFoundError, match="missing"):
        resolve_artifact_path(name, root)


def test_custom_managed_root_and_provenance(isolated_snapshots, tmp_path):
    root = artifacts.managed_path("deeptap", data_dir=tmp_path / "other-data")
    populate_snapshot(root, "deeptap")
    with pytest.raises(FileNotFoundError, match="invalid provenance"):
        resolve_artifact_path("deeptap", data_dir=tmp_path / "other-data")
    populate_snapshot(root, "deeptap", managed=True)
    assert resolve_artifact_path("deeptap", data_dir=tmp_path / "other-data") == root


def test_eramer_file_and_directory_overrides(isolated_snapshots, monkeypatch):
    root = complete_checkout(isolated_snapshots / "ERAMER", "eramer")
    pwm = isolated_snapshots / "custom.xlsx"
    pwm.touch()
    monkeypatch.setenv("ERAMER_PWM", str(pwm))
    assert Path(mhctools.ERAMER().pwm_path) == pwm
    assert Path(mhctools.ERAMER(eramer_home=root).pwm_path) == root / "PWM.xlsx"
    assert Path(mhctools.ERAMER(eramer_home=root, pwm_path=pwm).pwm_path) == pwm
    pwm.unlink()
    assert artifact_status("eramer").status == "missing"
    with pytest.raises(FileNotFoundError, match="ERAMER_PWM"):
        mhctools.ERAMER()


def test_fetch_does_not_hide_invalid_override(isolated_snapshots, monkeypatch):
    monkeypatch.setenv("DEEPTAP_HOME", str(isolated_snapshots / "absent"))
    monkeypatch.setattr(artifacts, "_fetch_snapshot", lambda *a, **kw: pytest.fail("download attempted"))
    with pytest.raises(FileNotFoundError, match="DEEPTAP_HOME"):
        artifacts.fetch("deeptap")
