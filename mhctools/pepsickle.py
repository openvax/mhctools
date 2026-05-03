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

import json
import logging
import subprocess
import sys

from .proteasome_predictor import ProteasomePredictor

# Module-level cache for loaded pepsickle models. Keyed by human_only.
_model_cache = {}

logger = logging.getLogger(__name__)

PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS = 300

_PEPSICKLE_SUBPROCESS_SCRIPT = r"""
import json
import sys

from pepsickle.model_functions import (
    initialize_epitope_model,
    predict_protein_cleavage_locations,
)

request = json.loads(sys.stdin.read())
model = initialize_epitope_model(human_only=request["human_only"])
results = {}
for sequence in request["sequences"]:
    preds_raw = predict_protein_cleavage_locations(
        sequence,
        model,
        mod_type="epitope",
        proteasome_type="C",
        threshold=request["threshold"],
    )
    results[sequence] = [entry[2] for entry in preds_raw]
json.dump({"results": results}, sys.stdout)
"""


class Pepsickle(ProteasomePredictor):
    """
    Proteasomal cleavage predictor using pepsickle's epitope model.

    Uses the in-vivo epitope model from Weeder et al. (Bioinformatics
    2021), which the paper shows outperforms the in-vitro alternatives
    and NetChop for neoantigen identification.

    Parameters
    ----------
    default_peptide_lengths : list of int, optional
        Peptide lengths used when scanning proteins. Default ``[9]``.

    scoring : callable, optional
        See :class:`ProcessingPredictor`.  Default:
        ``score_cterm_anti_max_internal``.

    threshold : float
        Cleavage probability threshold used by pepsickle internally
        (default 0.5).

    human_only : bool
        If True, use human-only trained model instead of all-mammal.

    isolate_subprocess : bool
        If True, run pepsickle inference in a short-lived Python subprocess.
        This avoids macOS duplicate OpenMP runtime crashes when the parent
        process has already imported packages such as pandas, numpy, or
        pyarrow.

    subprocess_timeout : int
        Timeout in seconds for isolated pepsickle inference.
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            threshold=0.5,
            human_only=False,
            isolate_subprocess=False,
            subprocess_timeout=PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS):
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.threshold = threshold
        self.human_only = human_only
        self.isolate_subprocess = isolate_subprocess
        self.subprocess_timeout = subprocess_timeout
        self._model = None

    def __str__(self):
        return "%s(scoring=%s, isolate_subprocess=%s)" % (
            self.__class__.__name__,
            getattr(self.scoring, "__name__", repr(self.scoring)),
            self.isolate_subprocess)

    def _predictor_name(self):
        return "pepsickle"

    def _load_model(self):
        if self._model is None:
            cache_key = self.human_only
            if cache_key not in _model_cache:
                from pepsickle.model_functions import initialize_epitope_model
                _model_cache[cache_key] = initialize_epitope_model(
                    human_only=self.human_only)
            self._model = _model_cache[cache_key]
        return self._model

    def cleavage_probs(self, sequence):
        return self.cleavage_probs_many([sequence])[sequence]

    def cleavage_probs_many(self, sequences):
        unique_sequences = list(dict.fromkeys(sequences))
        if not unique_sequences:
            return {}
        if self.isolate_subprocess:
            return self._cleavage_probs_many_subprocess(unique_sequences)
        return {
            sequence: self._cleavage_probs_in_process(sequence)
            for sequence in unique_sequences
        }

    def _cleavage_probs_in_process(self, sequence):
        from pepsickle.model_functions import predict_protein_cleavage_locations
        model = self._load_model()
        preds_raw = predict_protein_cleavage_locations(
            sequence,
            model,
            mod_type="epitope",
            proteasome_type="C",
            threshold=self.threshold,
        )
        return [entry[2] for entry in preds_raw]

    def _cleavage_probs_many_subprocess(self, sequences):
        payload = json.dumps({
            "human_only": bool(self.human_only),
            "threshold": float(self.threshold),
            "sequences": sequences,
        })
        try:
            result = subprocess.run(
                [sys.executable, "-c", _PEPSICKLE_SUBPROCESS_SCRIPT],
                input=payload,
                text=True,
                capture_output=True,
                timeout=self.subprocess_timeout,
            )
        except subprocess.TimeoutExpired as e:
            msg = (
                "pepsickle subprocess timed out after %d seconds "
                "while scoring %d sequences"
                % (self.subprocess_timeout, len(sequences)))
            logger.warning(msg)
            raise RuntimeError(msg) from e
        except OSError as e:
            msg = "Could not start pepsickle subprocess: %s" % e
            logger.warning(msg)
            raise RuntimeError(msg) from e

        stderr_text = result.stderr.strip()
        if stderr_text:
            logger.warning("pepsickle subprocess stderr:\n%s", stderr_text)
        if result.returncode != 0:
            logger.warning(
                "pepsickle subprocess exited with code %d",
                result.returncode)
            raise RuntimeError(
                "pepsickle subprocess exited with code %d.\nstdout: %s\n"
                "stderr: %s"
                % (result.returncode, result.stdout.strip(), stderr_text))

        try:
            parsed = json.loads(result.stdout)
        except ValueError as e:
            logger.warning("Could not parse pepsickle subprocess JSON output")
            raise RuntimeError(
                "Could not parse pepsickle subprocess JSON output: %s"
                % result.stdout.strip()) from e

        results = parsed.get("results")
        if not isinstance(results, dict):
            logger.warning(
                "pepsickle subprocess output is missing a results object")
            raise RuntimeError(
                "pepsickle subprocess output is missing a results object")
        missing = [sequence for sequence in sequences if sequence not in results]
        if missing:
            logger.warning(
                "pepsickle subprocess omitted %d sequences",
                len(missing))
            raise RuntimeError(
                "pepsickle subprocess omitted %d sequences" % len(missing))

        output = {}
        for sequence in sequences:
            probs = results[sequence]
            if len(probs) != len(sequence):
                logger.warning(
                    "pepsickle subprocess returned %d scores for a "
                    "%d-residue sequence",
                    len(probs),
                    len(sequence))
                raise ValueError(
                    "Expected %d pepsickle scores for sequence, got %d"
                    % (len(sequence), len(probs)))
            output[sequence] = [float(p) for p in probs]
        return output
