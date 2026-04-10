
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

from mhctools.cli import mhc_predictors
from mhctools.cli.args import _no_allele_predictors


def test_all_models_in_registry():
    """Verify all major model classes are available via CLI."""
    expected = [
        "netmhcpan42", "netmhcpan42-ba", "netmhcpan42-el",
        "netmhciipan43", "netmhciipan43-ba", "netmhciipan43-el",
        "bigmhc", "netmhcstabpan", "netchop", "pepsickle",
        # legacy entries that should still be present
        "netmhcpan", "netmhcpan4", "netmhcpan41",
        "netmhciipan", "netmhciipan4",
        "mhcflurry", "mixmhcpred", "random",
    ]
    for name in expected:
        assert name in mhc_predictors, f"{name} missing from mhc_predictors registry"


def test_no_allele_predictors_in_registry():
    """Verify processing predictors are marked as not needing alleles."""
    for name in _no_allele_predictors:
        assert name in mhc_predictors, (
            f"{name} in _no_allele_predictors but not in mhc_predictors"
        )


def test_all_registry_values_are_callable():
    """Every registry value should be callable (class or factory function)."""
    for name, cls in mhc_predictors.items():
        assert callable(cls), (
            f"mhc_predictors[{name!r}] = {cls!r} is not callable"
        )
