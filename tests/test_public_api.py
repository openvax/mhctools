"""The public import surface of ``mhctools`` and the lists that define it.

``mhctools/__init__.py`` names every public object in ``_LAZY_IMPORTS``, the
one list that decides what ``from mhctools import X`` can load. The
``TYPE_CHECKING`` block repeats those names for editors and mkdocstrings, and
``__all__`` is derived from ``_LAZY_IMPORTS``. These tests keep the repeated
block in step with the list (#656) and compare the surface with the names
recorded in ``tests/data/public_api_names.json``.

Adding, renaming or removing a public name is a deliberate API change: edit
``_LAZY_IMPORTS`` and the ``TYPE_CHECKING`` block, then update the JSON file.
"""

import ast
import json
from pathlib import Path

import mhctools

PUBLIC_API_NAMES_PATH = Path(__file__).parent / "data" / "public_api_names.json"


def type_checking_imports(source):
    """Collect the imports inside a module's ``if TYPE_CHECKING:`` block.

    Parameters
    ----------
    source : str
        Source code of the module.

    Returns
    -------
    dict
        Maps each imported alias to ``(module, attribute)``, with the module
        written as a relative path such as ``".pred"``, which is the form
        ``_LAZY_IMPORTS`` uses.
    """
    imports = {}
    for node in ast.parse(source).body:
        if not (isinstance(node, ast.If)
                and isinstance(node.test, ast.Name)
                and node.test.id == "TYPE_CHECKING"):
            continue
        for statement in node.body:
            assert isinstance(statement, ast.ImportFrom), ast.dump(statement)
            module = "." * statement.level + (statement.module or "")
            for alias in statement.names:
                exported_name = alias.asname or alias.name
                assert exported_name not in imports, exported_name
                imports[exported_name] = (module, alias.name)
    return imports


def recorded_public_names():
    """Load the recorded public names.

    Returns
    -------
    dict
        ``"lazy"`` holds every name importable from ``mhctools`` and ``"all"``
        holds the names in ``mhctools.__all__``, both sorted.
    """
    return json.loads(PUBLIC_API_NAMES_PATH.read_text())


def test_type_checking_imports_match_lazy_imports():
    source = Path(mhctools.__file__).read_text()
    static_imports = type_checking_imports(source)
    assert static_imports, "no TYPE_CHECKING imports found"
    assert static_imports == mhctools._LAZY_IMPORTS


def test_type_checking_imports_use_explicit_reexport_aliases():
    # ``from .module import Name as Name`` is the form type checkers treat as
    # a re-export, so the alias must repeat the attribute name.
    source = Path(mhctools.__file__).read_text()
    for name, (_module, attribute) in type_checking_imports(source).items():
        assert name == attribute


def test_star_export_exclusions_are_lazy_imports():
    excluded = mhctools._NOT_STAR_EXPORTED
    assert len(excluded) == len(set(excluded))
    assert set(excluded) <= set(mhctools._LAZY_IMPORTS)


def test_all_is_lazy_imports_without_exclusions():
    assert len(mhctools.__all__) == len(set(mhctools.__all__))
    assert set(mhctools.__all__) == (
        set(mhctools._LAZY_IMPORTS) - set(mhctools._NOT_STAR_EXPORTED))


def test_importable_names_match_record():
    assert sorted(mhctools._LAZY_IMPORTS) == recorded_public_names()["lazy"]


def test_all_matches_record():
    assert sorted(mhctools.__all__) == recorded_public_names()["all"]


def test_star_import_binds_exactly_the_recorded_names():
    namespace = {}
    exec("from mhctools import *", namespace)
    namespace.pop("__builtins__")
    assert sorted(namespace) == recorded_public_names()["all"]
