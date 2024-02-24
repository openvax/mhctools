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

import tempfile
from os import remove

from mhctools.cli.script import parse_args, run_predictor
from .common import eq_

def test_peptides_without_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)

def test_peptides_with_subsequences():
    peptide = "SIINFEKLQY"
    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--sequence", peptide,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])

def test_peptides_file_without_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = run_predictor(args)
    eq_(len(binding_predictions), 1, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide)
    remove(f.name)

def test_peptides_file_with_subsequences():
    peptide = "SIINFEKLQY"
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
        f.write("%s\n" % peptide)

    args = parse_args([
        "--mhc-predictor", "netmhc",
        "--mhc-peptide-lengths", "9",
        "--input-peptides-file", f.name,
        "--extract-subsequences",
        "--mhc-alleles", "H-2-Kb"])
    binding_predictions = sorted(run_predictor(args), key=lambda bp: bp.offset)
    eq_(len(binding_predictions), 2, binding_predictions)
    eq_(binding_predictions[0].peptide, peptide[:9])
    eq_(binding_predictions[1].peptide, peptide[1:10])
    remove(f.name)

