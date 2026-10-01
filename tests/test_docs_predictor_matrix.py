"""The docs predictor matrix must cover every exported predictor and CLI name."""

import importlib.util
import inspect
from pathlib import Path

import mhctools
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools import Kind
from mhctools.artifacts import list_artifacts
from mhctools.cli.args import mhc_predictors

_SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "predictor_matrix.py"
_spec = importlib.util.spec_from_file_location("predictor_matrix", _SCRIPT)
matrix = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(matrix)

_FACTORY_FUNCTIONS = {"NetMHC", "NetMHCIIpan", "NetMHCpan"}
_KIND_NAMES = {k for k in dir(Kind) if not k.startswith("_")}


def _exported_predictor_names():
    names = set(mhctools.__all__) | set(mhctools._LAZY_IMPORTS)
    found = set(_FACTORY_FUNCTIONS)
    for name in names:
        obj = getattr(mhctools, name, None)
        if (inspect.isclass(obj) and hasattr(obj, "predict")
                and hasattr(obj, "kind_support")
                and name not in matrix.NOT_PREDICTORS):
            found.add(name)
    return found


def test_every_exported_predictor_has_a_row():
    missing = _exported_predictor_names() - matrix.covered_classes()
    assert not missing, (
        "add these predictors to scripts/predictor_matrix.py: %s" % sorted(missing))


def test_matrix_names_only_real_classes():
    unknown = {c for c in matrix.covered_classes() if not hasattr(mhctools, c)}
    assert not unknown, "matrix lists classes mhctools does not export: %s" % sorted(unknown)


def test_every_cli_name_has_a_row():
    missing = set(mhc_predictors) - matrix.covered_cli_names()
    assert not missing, (
        "add these CLI names to scripts/predictor_matrix.py: %s" % sorted(missing))
    unknown = matrix.covered_cli_names() - set(mhc_predictors)
    assert not unknown, "matrix lists CLI names that do not exist: %s" % sorted(unknown)


def test_every_fetchable_artifact_has_a_row():
    artifact_names = {a.name for a in list_artifacts()}
    missing = artifact_names - matrix.covered_artifacts()
    assert not missing, "artifacts with no matrix row: %s" % sorted(missing)
    unknown = matrix.covered_artifacts() - artifact_names
    assert not unknown, "matrix names unknown artifacts: %s" % sorted(unknown)


def test_row_values_are_valid():
    tiers = set(matrix.LICENSE_TIERS)
    families = {key for key, _, _ in matrix.FAMILIES}
    for row in matrix.ROWS:
        assert row["family"] in families, row["name"]
        assert row["license"] in tiers, row["name"]
        for kind in row["kinds"]:
            assert kind in _KIND_NAMES, (row["name"], kind)


def _try_build(factory):
    for alleles in (["HLA-A*02:01"], ["HLA-DRB1*15:01"], None):
        kwargs = {"program_name": "/usr/bin/true"}
        if alleles:
            kwargs["alleles"] = alleles
        try:
            return factory(**kwargs)
        except Exception:
            continue
    return None


def test_rows_agree_with_instances(monkeypatch):
    """Where a predictor can be built without its tool, the kinds it reports
    must be among the kinds the row documents. Predictors that need their
    external tool to construct are checked by their own test modules."""
    monkeypatch.setattr(
        BaseCommandlinePredictor, "_determine_supported_alleles",
        staticmethod(lambda *args: {"HLA-A*02:01", "HLA-DRB1*15:01"}))
    checked = 0
    for cli_name, factory in sorted(mhc_predictors.items()):
        row = next(r for r in matrix.ROWS if cli_name in r["cli"])
        instance = _try_build(factory)
        if instance is None:
            continue
        checked += 1
        live_kinds = set(instance.kind_support())
        assert live_kinds <= set(row["kinds"]), (cli_name, live_kinds, row["kinds"])
    assert checked >= 20, "too few predictors could be constructed to cross-check"


def test_generated_page_is_current():
    assert matrix.PAGE.read_text() == matrix.render(), (
        "docs/predictor-matrix.md is stale; run python scripts/predictor_matrix.py")


def test_every_cleavage_model_is_documented():
    """docs/cleavage/models.md must mention every model `cleavage --list-models` reports."""
    import json
    from mhctools.cli.script import main
    import contextlib
    import io
    import sys

    buf = io.StringIO()
    old_argv = sys.argv
    sys.argv = ["mhctools", "cleavage", "--list-models", "--json"]
    try:
        with contextlib.redirect_stdout(buf):
            try:
                main()
            except SystemExit:
                pass
    finally:
        sys.argv = old_argv
    catalog = json.loads(buf.getvalue())
    models = catalog["models"] if isinstance(catalog, dict) else catalog
    page = (matrix.PAGE.parent / "cleavage" / "models.md").read_text()
    missing = [m["name"] for m in models if "`%s`" % m["name"] not in page]
    assert not missing, "docs/cleavage/models.md does not mention: %s" % missing
