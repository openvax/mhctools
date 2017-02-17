# Copyright (c) 2014-2017. Mount Sinai School of Medicine
#
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


import subprocess
import logging
import tempfile


class NetChop(object):
    """
    Wrapper around netChop tool. Assumes netChop is in your PATH.
    """

    def predict(self, sequences):
        """
        Return netChop predictions for each position in each sequence.

        Parameters
        -----------
        sequences : list of string
            Amino acid sequences to predict cleavage for

        Returns
        -----------
        list of list of float

        The i'th list corresponds to the i'th sequence. Each list gives
        the cleavage probability for each position in the sequence.
        """
        with tempfile.NamedTemporaryFile(suffix=".fsa", mode="w") as input_fd:
            for (i, sequence) in enumerate(sequences):
                input_fd.write("> %d\n" % i)
                input_fd.write(sequence)
                input_fd.write("\n")
            input_fd.flush()
            try:
                output = subprocess.check_output(["netChop", input_fd.name])
            except subprocess.CalledProcessError as e:
                logging.error("Error calling netChop: %s:\n%s" % (e, e.output))
                raise

        parsed = self.parse_netchop(output)
        assert len(parsed) == len(sequences), \
            "Expected %d results but got %d" % (
                len(sequences), len(parsed))
        assert [len(x) for x in parsed] == [len(x) for x in sequences]
        return parsed

    @staticmethod
    def parse_netchop(netchop_output):
        """
        Parse netChop stdout.
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
