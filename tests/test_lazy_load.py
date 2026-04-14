# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Regression tests for lazy loading of heavy predictor dependencies.

Importing the CLI argument helpers or the top-level mhctools package must not
pull in torch or tensorflow. These libraries are only needed when an instance
of MHCflurry or BigMHC is actually constructed.
"""

import subprocess
import sys


def _import_and_check(import_stmt):
    """Run `import_stmt` in a fresh subprocess and return which of
    {torch, tensorflow} were imported as a side effect."""
    code = (
        "import sys; " + import_stmt + "; "
        "loaded = [m for m in ('torch', 'tensorflow') if m in sys.modules]; "
        "print(','.join(loaded))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
    )
    out = result.stdout.strip()
    return set(out.split(",")) - {""}


def test_importing_mhctools_cli_args_does_not_load_torch_or_tensorflow():
    loaded = _import_and_check("import mhctools.cli.args")
    assert not loaded, (
        "Importing mhctools.cli.args pulled in: %s. "
        "Heavy predictors (BigMHC/MHCflurry) must stay behind _LazyPredictor."
        % sorted(loaded)
    )


def test_importing_top_level_mhctools_does_not_load_torch_or_tensorflow():
    loaded = _import_and_check("import mhctools")
    assert not loaded, (
        "Importing mhctools pulled in: %s. "
        "BigMHC/MHCflurry must stay behind the PEP 562 __getattr__ in __init__."
        % sorted(loaded)
    )


def test_accessing_mhcflurry_attribute_does_load_wrapper_module():
    """Sanity: touching `mhctools.MHCflurry` loads the wrapper module.
    (The heavy mhcflurry package stays deferred to MHCflurry.__init__.)"""
    code = (
        "import sys; import mhctools; _ = mhctools.MHCflurry; "
        "print('mhctools.mhcflurry' in sys.modules)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == "True"
