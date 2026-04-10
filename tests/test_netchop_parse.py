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
Unit tests for NetChop parsing and error handling that do NOT require the
netChop binary to be installed.
"""

import pytest
from unittest.mock import patch

from mhctools import NetChop


# ---------------------------------------------------------------------------
# Sample netChop output (trimmed from real output)
# ---------------------------------------------------------------------------

GOOD_OUTPUT = b"""\
some header junk
 pos  AA  C/N  score
------
   1   M    C   0.547588
   2   D    N   0.233000
   3   S    N   0.100000
---------
"""

TWO_SEQUENCE_OUTPUT = b"""\
 pos  AA  C/N  score
------
   1   M    C   0.500000
   2   D    N   0.200000
---------
 pos  AA  C/N  score
------
   1   A    C   0.600000
   2   R    N   0.300000
   3   G    N   0.100000
---------
"""

EMPTY_OUTPUT = b"Number of cleavage sites 0\n"

BAD_SCORE_OUTPUT = b"""\
 pos  AA  C/N  score
------
   1   M    C   notanumber
---------
"""

MISSING_DASHES_OUTPUT = b"""\
 pos  AA  C/N  score
no dashes here
   1   M    C   0.5
---------
"""

SHORT_LINE_OUTPUT = b"""\
 pos  AA  C/N  score
------
   1   M
---------
"""


class TestParseNetchop:

    def test_good_output(self):
        result = NetChop.parse_netchop(GOOD_OUTPUT)
        assert len(result) == 1
        assert len(result[0]) == 3
        assert abs(result[0][0] - 0.547588) < 1e-6
        assert abs(result[0][1] - 0.233000) < 1e-6

    def test_two_sequences(self):
        result = NetChop.parse_netchop(TWO_SEQUENCE_OUTPUT)
        assert len(result) == 2
        assert len(result[0]) == 2
        assert len(result[1]) == 3

    def test_empty_output_returns_no_sequences(self):
        result = NetChop.parse_netchop(EMPTY_OUTPUT)
        assert result == []

    def test_bad_score_raises_valueerror(self):
        with pytest.raises(ValueError, match="Could not parse score"):
            NetChop.parse_netchop(BAD_SCORE_OUTPUT)

    def test_missing_dashes_raises_valueerror(self):
        with pytest.raises(ValueError, match="Expected dashes"):
            NetChop.parse_netchop(MISSING_DASHES_OUTPUT)

    def test_short_line_raises_valueerror(self):
        with pytest.raises(ValueError, match="Unexpected netChop output line"):
            NetChop.parse_netchop(SHORT_LINE_OUTPUT)


class TestNetChopInit:

    def test_missing_executable_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError, match="Could not find"):
            NetChop(program_name="surely_not_installed_netchop_xyz")

    @patch("shutil.which", return_value="/usr/local/bin/netChop")
    def test_found_executable_succeeds(self, _mock_which):
        obj = NetChop(program_name="netChop")
        assert obj.program_name == "netChop"


class TestNetChopCleavageProbs:

    @patch("shutil.which", return_value="/usr/local/bin/netChop")
    def test_wrong_number_of_sequences_raises(self, _mock_which):
        """If netChop returns 0 sequences, cleavage_probs should raise."""
        obj = NetChop()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = EMPTY_OUTPUT
            mock_run.return_value.stderr = b""
            mock_run.return_value.returncode = 0
            with pytest.raises(ValueError, match="Expected 1 result"):
                obj.cleavage_probs("MDS")

    @patch("shutil.which", return_value="/usr/local/bin/netChop")
    def test_wrong_length_raises_with_pepsickle_hint(self, _mock_which):
        """Length mismatch should mention pepsickle as alternative."""
        obj = NetChop()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = GOOD_OUTPUT  # 3 scores
            mock_run.return_value.stderr = b""
            mock_run.return_value.returncode = 0
            with pytest.raises(ValueError, match="pepsickle"):
                obj.cleavage_probs("MDSH")  # 4 chars, but output has 3

    @patch("shutil.which", return_value="/usr/local/bin/netChop")
    def test_nonzero_returncode_raises(self, _mock_which):
        obj = NetChop()
        with patch("subprocess.run") as mock_run:
            mock_run.return_value.stdout = b""
            mock_run.return_value.stderr = b"some error"
            mock_run.return_value.returncode = 1
            with pytest.raises(RuntimeError, match="exited with code 1"):
                obj.cleavage_probs("MDS")

    @patch("shutil.which", return_value="/usr/local/bin/netChop")
    def test_timeout_raises(self, _mock_which):
        import subprocess
        obj = NetChop()
        with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("netChop", 120)):
            with pytest.raises(RuntimeError, match="timed out"):
                obj.cleavage_probs("MDS")
