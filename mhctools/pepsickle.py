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

import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import subprocess
import sys

from .cleavage import CleavageModel, CleavageResult, CleavageSite, coerce_peptide
from .proteasome_predictor import ProteasomePredictor

# Module-level cache for loaded pepsickle models. Keyed by human_only.
_model_cache = {}
_identity_cache = {}

logger = logging.getLogger(__name__)

PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS = 300


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _pepsickle_identity(human_only):
    """Return verified code, feature, and weight identity for pepsickle."""
    cache_key = bool(human_only)
    if cache_key in _identity_cache:
        return _identity_cache[cache_key]

    import pepsickle.model_functions as model_functions
    import pepsickle.sequence_featurization_tools as feature_functions

    package_dir = Path(model_functions.__file__).resolve().parent
    weights_path = package_dir / "trained_model_dict.pickle"
    if not weights_path.is_file():
        raise RuntimeError(
            "Installed pepsickle is missing trained_model_dict.pickle at %s"
            % weights_path)
    identity = {
        "package_version": importlib.metadata.version("pepsickle"),
        "weights_path": str(weights_path),
        "weights_sha256": _sha256(weights_path),
        "inference_path": str(Path(model_functions.__file__).resolve()),
        "inference_sha256": _sha256(model_functions.__file__),
        "features_path": str(Path(feature_functions.__file__).resolve()),
        "features_sha256": _sha256(feature_functions.__file__),
        "model_key": (
            "human_epitope_sequence_mod+human_epitope_motif_mod"
            if human_only else
            "all_mammal_epitope_sequence_mod+all_mammal_epitope_motif_mod"
        ),
    }
    _identity_cache[cache_key] = identity
    return identity

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

    @staticmethod
    def _cleavage_model_from_identity(human_only, identity):
        population = "human-only" if human_only else "all-mammal"
        if identity is None:
            version = "unresolved: pepsickle package/assets not located"
        else:
            version = (
                "package:%s;model:%s;weights-sha256:%s;inference-sha256:%s;"
                "features-sha256:%s"
                % (
                    identity["package_version"],
                    identity["model_key"],
                    identity["weights_sha256"],
                    identity["inference_sha256"],
                    identity["features_sha256"],
                )
            )
        return CleavageModel(
            name="pepsickle-in-vivo-%s" % population,
            version=version,
            enzyme="proteasome epitope proxy",
            uniprot="",
            species="Homo sapiens" if human_only else "Mammalia",
            compartments=("cytosol",),
            evidence="quantitative_model",
            references=(
                "https://doi.org/10.1093/bioinformatics/btab628",
                "https://github.com/pdxgx/pepsickle",
            ),
            assay=(
                "Neural ensemble trained from observed epitope C termini; %s "
                "training subset" % population
            ),
            limitations=(
                "Native dimensionless model output, not an empirical cleavage, "
                "degradation, presentation, or vaccine efficacy probability. "
                "Uses an eight-residue context on each side with terminal "
                "padding. The upstream package deserializes a pickle model; "
                "isolate_subprocess controls process isolation but does not make "
                "an untrusted pickle safe."
            ),
            score_name="pepsickle_epitope_model_output",
            score_units="dimensionless native model output",
            scored_endpoint="site_cleavage",
        )

    @classmethod
    def catalog_cleavage_model(cls, human_only):
        """Describe a known optional model without claiming absent assets."""
        try:
            identity = _pepsickle_identity(human_only)
        except (ImportError, OSError, RuntimeError,
                importlib.metadata.PackageNotFoundError):
            identity = None
        return cls._cleavage_model_from_identity(human_only, identity)

    def cleavage_model(self):
        """Return canonical per-bond metadata tied to installed assets."""
        return self._cleavage_model_from_identity(
            self.human_only, _pepsickle_identity(self.human_only))

    def predict_cleavage(self, peptide):
        """Return every internal bond through the canonical cleavage contract.

        Pepsickle emits one score after each residue and forces its last score
        to zero because a sequence endpoint is not an internal peptide bond.
        This adapter maps array index ``i`` to bond ``i + 1`` and deliberately
        excludes that endpoint sentinel.
        """
        peptide = coerce_peptide(peptide)
        model = self.cleavage_model()
        if peptide.n_term != "free" or peptide.c_term != "free":
            return CleavageResult(
                peptide,
                model,
                unsupported_reason=(
                    "Pepsickle does not model modified terminal chemistry"
                ),
            )
        scores = self.cleavage_probs(peptide.sequence)
        if len(scores) != len(peptide.sequence):
            raise ValueError(
                "Expected %d pepsickle scores for sequence, got %d"
                % (len(peptide.sequence), len(scores)))
        sites = tuple(
            CleavageSite(
                bond,
                "scored",
                "Pepsickle epitope-model output after residue %d" % bond,
                float(score),
            )
            for bond, score in enumerate(scores[:-1], start=1)
        )
        identity = _pepsickle_identity(self.human_only)
        conditions = (
            ("human_only", str(bool(self.human_only)).lower()),
            ("model_mode", "epitope"),
            ("model_key", identity["model_key"]),
            ("threshold", "%.17g" % self.threshold),
            ("isolate_subprocess", str(bool(self.isolate_subprocess)).lower()),
            ("subprocess_timeout_seconds", str(self.subprocess_timeout)),
        )
        return CleavageResult(peptide, model, sites, conditions=conditions)

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


class PepsickleCleavage:
    """Canonical facade that leaves legacy peptide-level ``predict`` intact."""

    def __init__(self, **kwargs):
        kwargs.setdefault("isolate_subprocess", True)
        self.predictor = Pepsickle(**kwargs)

    @property
    def model(self):
        return self.predictor.cleavage_model()

    def predict(self, peptide):
        return self.predictor.predict_cleavage(peptide)
