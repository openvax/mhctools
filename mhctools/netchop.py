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
        assert len(parsed) == len(sequences)
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
