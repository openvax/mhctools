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

"""The hosted SMM subset must stay traceable to one upstream IEDB release.

Three files pin IEDB 3.1.7 independently: the subset builder, the fixture
recorder, and the backend installer. If one is bumped to a new upstream
release without the others, the recorded fixtures silently stop describing
what CI installs.
"""

import importlib.util
from pathlib import Path
import re
import sys

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name):
    path = SCRIPTS / ("%s.py" % name)
    if not path.exists():
        pytest.skip("%s not present" % path)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    # These scripts import each other as top-level modules when run from
    # scripts/, so reproduce that import path here.
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


def test_subset_builder_and_fixture_recorder_pin_the_same_release():
    builder = _load("build_iedb_smm_subset")
    recorder = _load("record_smm_fixtures")

    assert builder.ARCHIVE_SHA256 == recorder.ARCHIVE_SHA256


def test_installer_pins_a_subset_checksum():
    installer = _load("setup_test_backends")

    assert re.fullmatch(r"[0-9a-f]{64}", installer.SMM_SUBSET_SHA256)


def test_installer_downloads_the_release_the_builder_produces():
    """The installed artifact must be the one build_iedb_smm_subset writes."""
    installer = _load("setup_test_backends")

    assert installer.SMM_SUBSET_URL.startswith(
        "https://github.com/openvax/mhctools/releases/download/")
    assert installer.SMM_SUBSET_URL.endswith(
        "IEDB_MHC_I-3.1.7-smm-subset.tar.gz")


def test_subset_allowlist_keeps_the_license_and_percentile_data():
    """Both are load-bearing and easy to drop while pruning for size."""
    builder = _load("build_iedb_smm_subset")

    assert "LIAI_license.txt" in builder.INCLUDED_PREFIXES
    # Percentile ranks read distribution_consensus_bin.cpickle from here.
    assert any(
        prefix.endswith("consensus/")
        for prefix in builder.INCLUDED_PREFIXES)
    # The bulk of the upstream archive must stay out.
    assert not any(
        name in prefix
        for prefix in builder.INCLUDED_PREFIXES
        for name in ("netmhcpan", "netmhccons", "netmhcstabpan", "pickpocket"))
