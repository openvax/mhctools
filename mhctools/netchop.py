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
import subprocess
import tempfile

from .processing_predictor import ProcessingPredictor

logger = logging.getLogger(__name__)


class NetChop(ProcessingPredictor):
    """
    Wrapper around the netChop command-line tool.

    Assumes ``netChop`` (or a custom *program_name*) is on your PATH.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : str
        Aggregation method (see :class:`ProcessingPredictor`).

    program_name : str
        Name or path of the netChop executable (default ``"netChop"``).
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring="nterm_cterm_max_internal",
            program_name="netChop"):
        ProcessingPredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.program_name = program_name

    def __str__(self):
        return "%s(program_name=%r, scoring=%r)" % (
            self.__class__.__name__, self.program_name, self.scoring)

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
                output = subprocess.check_output([self.program_name, fd.name])
            except subprocess.CalledProcessError as e:
                logger.error(
                    "Error calling %s: %s:\n%s",
                    self.program_name, e, e.output)
                raise
        parsed = self.parse_netchop(output)
        assert len(parsed) == 1, \
            "Expected 1 result from netChop, got %d" % len(parsed)
        assert len(parsed[0]) == len(sequence), \
            "Expected %d scores, got %d" % (len(sequence), len(parsed[0]))
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
        line_iterator = iter(netchop_output.decode().split("\n"))
        scores = []
        for line in line_iterator:
            if "pos" in line and 'AA' in line and 'score' in line:
                scores.append([])
                if "----" not in next(line_iterator):
                    raise ValueError("Dashes expected")
                line = next(line_iterator)
                while '-------' not in line:
                    score = float(line.split()[3])
                    scores[-1].append(score)
                    line = next(line_iterator)
        return scores
