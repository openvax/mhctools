# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import shutil
import subprocess
import tempfile

from .proteasome_predictor import ProteasomePredictor

logger = logging.getLogger(__name__)

# Timeout in seconds for a single netChop invocation.
NETCHOP_TIMEOUT_SECONDS = 120


class NetChop(ProteasomePredictor):
    """
    Wrapper around the netChop command-line tool.

    Assumes ``netChop`` (or a custom *program_name*) is on your PATH.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : callable, optional
        See :class:`ProcessingPredictor`.  Default:
        ``score_cterm_anti_max_internal``.

    program_name : str
        Name or path of the netChop executable (default ``"netChop"``).
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            program_name="netChop"):
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.program_name = program_name
        if not shutil.which(self.program_name):
            raise FileNotFoundError(
                "Could not find '%s' on PATH. Is netChop installed? "
                "You can pass a full path via program_name=."
                % self.program_name)

    def __str__(self):
        return "%s(program_name=%r, scoring=%s)" % (
            self.__class__.__name__,
            self.program_name,
            getattr(self.scoring, "__name__", repr(self.scoring)))

    def _predictor_name(self):
        return "netchop"

    def cleavage_probs(self, sequence):
        """
        Run netChop on a single sequence.

        Returns
        -------
        list of float
            Per-position cleavage probabilities.
        """
        with tempfile.NamedTemporaryFile(suffix=".fsa", mode="w") as fd:
            fd.write("> seq\n")
            fd.write(sequence)
            fd.write("\n")
            fd.flush()
            try:
                result = subprocess.run(
                    [self.program_name, fd.name],
                    capture_output=True,
                    timeout=NETCHOP_TIMEOUT_SECONDS,
                )
            except subprocess.TimeoutExpired:
                raise RuntimeError(
                    "%s timed out after %d seconds on a sequence of "
                    "length %d"
                    % (self.program_name, NETCHOP_TIMEOUT_SECONDS,
                       len(sequence)))
            except FileNotFoundError:
                raise FileNotFoundError(
                    "Could not execute '%s'. Is netChop installed and on "
                    "your PATH?" % self.program_name)
        stderr_text = result.stderr.decode("utf-8", errors="replace").strip()
        if stderr_text:
            logger.warning("%s stderr:\n%s", self.program_name, stderr_text)
        if result.returncode != 0:
            raise RuntimeError(
                "%s exited with code %d.\nstdout: %s\nstderr: %s"
                % (self.program_name, result.returncode,
                   result.stdout.decode("utf-8", errors="replace").strip(),
                   stderr_text))
        parsed = self.parse_netchop(result.stdout)
        if len(parsed) != 1:
            raise ValueError(
                "Expected 1 result sequence from %s, got %d. "
                "stdout: %s\nstderr: %s"
                % (self.program_name, len(parsed),
                   result.stdout.decode("utf-8", errors="replace").strip(),
                   stderr_text))
        if len(parsed[0]) != len(sequence):
            raise ValueError(
                "Expected %d per-position scores from %s, got %d. "
                "This usually means netChop's internal temp files are "
                "broken — try reinstalling netChop or use pepsickle as "
                "an alternative (--mhc-predictor pepsickle).\n"
                "stdout: %s\nstderr: %s"
                % (len(sequence), self.program_name, len(parsed[0]),
                   result.stdout.decode("utf-8", errors="replace").strip(),
                   stderr_text))
        return parsed[0]

    @staticmethod
    def parse_netchop(netchop_output):
        """
        Parse netChop stdout.

        Returns
        -------
        list of list of float
            One inner list per input sequence, with per-position
            cleavage scores.
        """
        text = netchop_output.decode("utf-8", errors="replace")
        lines = text.split("\n")
        line_iterator = iter(lines)
        scores = []
        for line in line_iterator:
            if "pos" in line and 'AA' in line and 'score' in line:
                scores.append([])
                dashes_line = next(line_iterator, "")
                if "----" not in dashes_line:
                    raise ValueError(
                        "Expected dashes after netChop header, got: %r"
                        % dashes_line)
                line = next(line_iterator, "-------")
                while '-------' not in line:
                    parts = line.split()
                    if len(parts) < 4:
                        raise ValueError(
                            "Unexpected netChop output line: %r" % line)
                    try:
                        score = float(parts[3])
                    except ValueError:
                        raise ValueError(
                            "Could not parse score from netChop "
                            "output line: %r" % line)
                    scores[-1].append(score)
                    line = next(line_iterator, "-------")
        return scores
