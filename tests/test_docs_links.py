"""Relative links and anchors in README.md and docs/ must resolve."""

import importlib.util
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "check_docs_links.py"
_spec = importlib.util.spec_from_file_location("check_docs_links", _SCRIPT)
check_docs_links = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_docs_links)


def test_docs_links_resolve():
    errors = check_docs_links.check()
    assert not errors, "\n".join(errors)


def test_slug_matches_github_anchors():
    slug = check_docs_links.slug
    assert slug("score, value, and percentile_rank") == "score-value-and-percentile_rank"
    assert slug("Annotate a table (`predict-table`)") == "annotate-a-table-predict-table"
    assert slug("SMM and SMM-PMBEC") == "smm-and-smm-pmbec"
