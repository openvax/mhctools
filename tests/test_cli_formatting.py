"""Catalog formatting stays consistent across CLI commands."""

import pytest

from mhctools.cli.formatting import format_table
from mhctools.cli.script import main


def test_table_aligns_headers_and_rows_without_trailing_spaces():
    rows = iter([("a", "long value"), ("longer name", "-")])
    assert format_table(("NAME", "VALUE"), rows) == (
        "NAME         VALUE\n"
        "a            long value\n"
        "longer name  -"
    )


@pytest.mark.parametrize("command,provider,headers", [
    (["ls"], "mhctools.cli.artifacts.list_artifacts",
     "NAME  STATUS  MANAGER  VERSION  FETCHABLE"),
    (["predictors"], "mhctools.cli.integrations.list_integrations",
     "NAME  LOCATED  RUNNABLE  REPRODUCED  CAPABILITY  PATH"),
    (["integrations"], "mhctools.cli.integrations.list_integrations",
     "NAME  LOCATED  RUNNABLE  REPRODUCED  CAPABILITY  PATH"),
    (["cleavage", "--list-models"], "mhctools.cli.cleavage.cleavage_models",
     "NAME  ENZYME  COMPARTMENTS  EVIDENCE  SCORE UNITS"),
])
def test_empty_catalog_prints_headers(command, provider, headers, monkeypatch, capsys):
    monkeypatch.setattr(provider, lambda *args, **kwargs: [])
    assert main(command) in (None, 0)
    captured = capsys.readouterr()
    assert captured.out == headers + "\n"
    assert captured.err == ""
