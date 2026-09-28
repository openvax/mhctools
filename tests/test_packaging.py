# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Packaging declares what the package actually ships and uses (#351).

varcode was required by every install for years without being imported once,
and logging.conf was shipped as package data in two locations, one of which
did not exist. Neither is visible from the code, so pin both here.
"""

from pathlib import Path
import re
import subprocess

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    tomllib = pytest.importorskip("tomli", reason="no TOML parser available")

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"


@pytest.fixture(scope="module")
def pyproject():
    with open(PYPROJECT, "rb") as handle:
        return tomllib.load(handle)


def _distribution_name(requirement):
    """'numpy>=1.26; python_version < "3.10"' -> 'numpy'."""
    return re.split(r"[<>=!;\[ ]", requirement.strip())[0]


def test_every_runtime_dependency_is_imported(pyproject):
    """A dependency nothing imports is weight on every install."""
    unused = []
    for requirement in pyproject["project"]["dependencies"]:
        name = _distribution_name(requirement)
        pattern = r"^\s*(import|from)\s+%s\b" % re.escape(name)
        found = subprocess.run(
            ["grep", "-rEl", pattern, "mhctools"],
            cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()
        if not found:
            unused.append(name)
    assert not unused, (
        "dependency declared but never imported: %s" % ", ".join(unused))


def test_declared_package_data_exists(pyproject):
    """Shipping a path that isn't in the tree silently ships nothing."""
    package_data = pyproject["tool"]["setuptools"]["package-data"]
    missing = []
    for package, patterns in package_data.items():
        directory = REPO_ROOT / Path(*package.split("."))
        for pattern in patterns:
            if not list(directory.glob(pattern)):
                missing.append("%s: %s" % (package, pattern))
    assert not missing, (
        "package-data pattern matches no file: %s" % "; ".join(missing))


def test_no_unreferenced_requirements_file():
    """requirements.txt was read by nothing and had drifted from pyproject."""
    assert not (REPO_ROOT / "requirements.txt").exists()
