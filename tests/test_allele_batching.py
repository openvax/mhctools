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

"""
Unit tests for allele batching in BaseCommandlinePredictor.

These exercise the command-construction logic (which alleles end up in which
process) WITHOUT running any external predictor binary, so they run on every
platform. End-to-end validation against the real netMHCpan output lives in
tests/test_netmhc_pan.py (which needs the binary) and the multi-allele parser
tests in tests/test_mhc_formats.py.
"""

import inspect
import os

import pytest

from mhctools import (
    NetMHCpan,
    NetMHCpan28,
    NetMHCpan3,
    NetMHCpan4,
    NetMHCpan41,
    NetMHCpan42,
)
from mhctools.base_commandline_predictor import BaseCommandlinePredictor
from mhctools.binding_prediction_collection import BindingPredictionCollection


class _StubPredictor(BaseCommandlinePredictor):
    """
    Exercises the pure command-construction logic without a real binary.

    Bypasses BaseCommandlinePredictor.__init__ (which shells out to netMHCpan
    to list supported alleles) and sets only the attributes that
    _allele_groups(), _build_command() and predict_peptides() read.
    """
    def __init__(
            self,
            alleles,
            max_alleles_per_command=1,
            prepare=None,
            max_peptides_per_file=10 ** 4,
            process_limit=-1):
        self.program_name = "netMHCpan"
        self.peptide_mode_flags = ["-p"]
        self.allele_flag = "-a"
        self.length_flag = "-l"
        self.input_file_flag = "-f"
        self.tempdir_flag = None
        self.extra_flags = []
        self.alleles = list(alleles)
        self.max_alleles_per_command = max_alleles_per_command
        self.max_peptides_per_file = max_peptides_per_file
        self.process_limit = process_limit
        self.group_peptides_by_length = False
        # peptide-validation attributes touched by _check_peptide_inputs
        self.allow_X_in_peptides = False
        self.allow_lowercase_in_peptides = False
        self.min_peptide_length = None
        self.max_peptide_length = None
        self.parse_output_fn = None
        self.parse_to_preds_fn = None
        self._prepare = prepare
        self.captured_commands = None

    def prepare_allele_name(self, allele_name):
        if self._prepare is not None:
            return self._prepare(allele_name)
        return super().prepare_allele_name(allele_name)

    def _run_commands_and_collect_predictions(
            self, commands, input_filenames, temp_dir_list,
            sequence_key_mapping=None):
        # Capture the commands instead of running them, and clean up the temp
        # files predict_peptides() created.
        self.captured_commands = dict(commands)
        for output_file in commands:
            output_file.close()
            _silent_remove(output_file.name)
        for fname in input_filenames:
            _silent_remove(fname)
        return BindingPredictionCollection([])

    def _check_results(self, binding_predictions, peptides, alleles):
        pass  # result completeness is not what these tests check


def _silent_remove(path):
    try:
        os.remove(path)
    except OSError:
        pass


def _allele_arg(command):
    """The single value passed after the -a flag in a built command."""
    assert command.count("-a") == 1, \
        "expected exactly one -a flag, got: %s" % (command,)
    return command[command.index("-a") + 1]


# ---------------------------------------------------------------------------
# _allele_groups(): how self.alleles is partitioned into per-process groups
# ---------------------------------------------------------------------------

def test_allele_groups_one_per_command_is_default_behavior():
    p = _StubPredictor(["A", "B", "C"], max_alleles_per_command=1)
    assert p._allele_groups() == [["A"], ["B"], ["C"]]


def test_allele_groups_none_batches_all():
    p = _StubPredictor(["A", "B", "C"], max_alleles_per_command=None)
    assert p._allele_groups() == [["A", "B", "C"]]


def test_allele_groups_chunk_size():
    p = _StubPredictor(["A", "B", "C", "D", "E"], max_alleles_per_command=2)
    assert p._allele_groups() == [["A", "B"], ["C", "D"], ["E"]]


def test_allele_groups_non_positive_batches_all():
    p = _StubPredictor(["A", "B", "C"], max_alleles_per_command=0)
    assert p._allele_groups() == [["A", "B", "C"]]


def test_allele_groups_larger_than_alleles_is_single_group():
    p = _StubPredictor(["A", "B"], max_alleles_per_command=10)
    assert p._allele_groups() == [["A", "B"]]


def test_allele_groups_empty():
    assert _StubPredictor([], max_alleles_per_command=None)._allele_groups() == []
    assert _StubPredictor([], max_alleles_per_command=1)._allele_groups() == []


# ---------------------------------------------------------------------------
# "auto": batch alleles but keep enough processes to use the cores
# ---------------------------------------------------------------------------

def test_auto_batches_all_when_files_already_saturate_cores():
    # 8 input files already reach the 4-process target, so all alleles can
    # share one command per file (fewest network reloads).
    p = _StubPredictor(
        ["A", "B", "C", "D", "E"],
        max_alleles_per_command="auto",
        process_limit=4)
    groups = p._allele_groups(n_input_files=8)
    assert groups == [["A", "B", "C", "D", "E"]]


def test_auto_keeps_parallelism_with_a_single_input_file():
    # 1 input file, 4-process target -> need 4 allele groups to keep cores
    # busy, so a single file must not collapse to one process.
    p = _StubPredictor(
        ["A", "B", "C", "D"],
        max_alleles_per_command="auto",
        process_limit=4)
    groups = p._allele_groups(n_input_files=1)
    assert len(groups) == 4
    assert [a for g in groups for a in g] == ["A", "B", "C", "D"]


def test_auto_never_more_groups_than_alleles():
    # target far exceeds allele count: at most one allele per group.
    p = _StubPredictor(
        ["A", "B"], max_alleles_per_command="auto", process_limit=32)
    groups = p._allele_groups(n_input_files=1)
    assert groups == [["A"], ["B"]]


def test_auto_partial_batching_with_few_files():
    # 2 files, target 8 -> 4 allele groups per file; 8 alleles -> size-2 groups.
    p = _StubPredictor(
        ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8"],
        max_alleles_per_command="auto",
        process_limit=8)
    groups = p._allele_groups(n_input_files=2)
    assert groups == [
        ["A1", "A2"], ["A3", "A4"], ["A5", "A6"], ["A7", "A8"]]


def test_auto_caps_group_size():
    # Many files would let "auto" batch every allele into one call; the cap
    # keeps groups bounded so a single call never grows unbounded.
    from mhctools.base_commandline_predictor import (
        AUTO_MAX_ALLELES_PER_COMMAND as CAP)
    n_alleles = CAP * 3 + 5
    p = _StubPredictor(
        ["A%d" % i for i in range(n_alleles)],
        max_alleles_per_command="auto",
        process_limit=2)
    # 100 input files easily saturate the 2-process target, so absent a cap
    # auto would put all alleles in one group.
    groups = p._allele_groups(n_input_files=100)
    assert max(len(g) for g in groups) <= CAP
    assert sum(len(g) for g in groups) == n_alleles


def test_none_still_batches_all_without_cap():
    # The explicit "unbounded" escape hatch is not subject to the auto cap.
    from mhctools.base_commandline_predictor import (
        AUTO_MAX_ALLELES_PER_COMMAND as CAP)
    n_alleles = CAP * 3
    p = _StubPredictor(
        ["A%d" % i for i in range(n_alleles)], max_alleles_per_command=None)
    groups = p._allele_groups(n_input_files=100)
    assert len(groups) == 1
    assert len(groups[0]) == n_alleles


# ---------------------------------------------------------------------------
# _build_command(): a group becomes one comma-separated -a argument
# ---------------------------------------------------------------------------

def test_build_command_joins_alleles_into_single_flag():
    p = _StubPredictor(["HLA-A*02:01", "HLA-B*35:02"])
    cmd = p._build_command(
        input_filename="peptides.txt",
        alleles=["HLA-A*02:01", "HLA-B*35:02"],
        peptide_mode=True)
    # default prepare_allele_name strips the '*'
    assert _allele_arg(cmd) == "HLA-A02:01,HLA-B35:02"
    assert cmd[:2] == ["netMHCpan", "-p"]
    assert cmd[-2:] == ["-f", "peptides.txt"]


def test_build_command_single_allele_string_still_works():
    p = _StubPredictor(["HLA-A*02:01"])
    cmd = p._build_command(
        input_filename="peptides.txt",
        alleles="HLA-A*02:01",
        peptide_mode=True)
    assert _allele_arg(cmd) == "HLA-A02:01"


def test_build_command_applies_prepare_per_allele_then_joins():
    # netMHC4-style prepare that also strips ':'; must be applied to each
    # allele individually before joining, not to the joined string.
    def prepare(a):
        return a.replace("*", "").replace(":", "")
    p = _StubPredictor(["HLA-A*02:01", "HLA-B*35:02"], prepare=prepare)
    cmd = p._build_command(
        input_filename="peptides.txt",
        alleles=["HLA-A*02:01", "HLA-B*35:02"],
        peptide_mode=True)
    assert _allele_arg(cmd) == "HLA-A0201,HLA-B3502"


# ---------------------------------------------------------------------------
# predict_peptides() wiring: number of processes and their -a arguments
# ---------------------------------------------------------------------------

def test_predict_peptides_batches_all_alleles_into_one_process():
    alleles = ["HLA-A*02:01", "HLA-B*35:02", "HLA-C*07:02"]
    p = _StubPredictor(alleles, max_alleles_per_command=None)
    p.predict_peptides(["SIINFEKLL", "GILGFVFTL"])
    commands = list(p.captured_commands.values())
    # one input file (2 peptides) x one allele group = one process
    assert len(commands) == 1
    got = set(_allele_arg(commands[0]).split(","))
    assert got == {"HLA-A02:01", "HLA-B35:02", "HLA-C07:02"}


def test_predict_peptides_one_process_per_allele_when_not_batched():
    alleles = ["HLA-A*02:01", "HLA-B*35:02", "HLA-C*07:02"]
    p = _StubPredictor(alleles, max_alleles_per_command=1)
    p.predict_peptides(["SIINFEKLL", "GILGFVFTL"])
    commands = list(p.captured_commands.values())
    # one input file x three alleles = three processes, each single-allele
    assert len(commands) == 3
    per_command = [_allele_arg(c) for c in commands]
    assert all("," not in a for a in per_command)
    assert set(per_command) == {"HLA-A02:01", "HLA-B35:02", "HLA-C07:02"}


def test_predict_peptides_command_count_is_files_times_groups():
    alleles = ["HLA-A*02:01", "HLA-B*35:02"]
    peptides = ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE"]
    # max_peptides_per_file=2 -> 4 peptides split into 2 input files
    batched = _StubPredictor(
        alleles, max_alleles_per_command=None, max_peptides_per_file=2)
    batched.predict_peptides(peptides)
    assert len(batched.captured_commands) == 2  # 2 files x 1 group

    per_allele = _StubPredictor(
        alleles, max_alleles_per_command=1, max_peptides_per_file=2)
    per_allele.predict_peptides(peptides)
    assert len(per_allele.captured_commands) == 4  # 2 files x 2 alleles


# ---------------------------------------------------------------------------
# The netMHCpan family batches all alleles by default (perf win out of the box)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("factory", [
    NetMHCpan, NetMHCpan28, NetMHCpan3, NetMHCpan4, NetMHCpan41, NetMHCpan42,
])
def test_pan_wrappers_default_to_auto_batching(factory):
    default = inspect.signature(factory).parameters[
        "max_alleles_per_command"].default
    assert default == "auto"


@pytest.mark.parametrize("factory", [
    NetMHCpan, NetMHCpan28, NetMHCpan3, NetMHCpan4, NetMHCpan41, NetMHCpan42,
])
def test_pan_wrappers_expose_max_peptides_per_file(factory):
    assert "max_peptides_per_file" in inspect.signature(factory).parameters


def test_base_predictor_defaults_to_one_allele_per_command():
    # Non-pan command-line predictors keep the historical behavior.
    default = inspect.signature(
        BaseCommandlinePredictor).parameters["max_alleles_per_command"].default
    assert default == 1
