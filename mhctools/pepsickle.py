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
import os
from pathlib import Path
import shutil
import subprocess
import sys

from .cleavage import CleavageModel, CleavageResult, CleavageSite, coerce_peptide
from .proteasome_predictor import ProteasomePredictor
from .pepsickle_runtime import runtime_identity

# Models are cached by family/population; C/I is an inference input.
_model_cache = {}
_identity_cache = {}

logger = logging.getLogger(__name__)

PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS = 300

PEPSICKLE_MODELS = {
    "pepsickle-in-vivo-" + population: dict(human_only=human)
    for population, human in (("human-only", True), ("all-mammal", False))
}
for _family in ("in-vitro", "in-vitro-2"):
    for _population, _human in (("human-only", True), ("all-mammal", False)):
        if _family == "in-vitro" and _human:
            continue  # Upstream's gradient-boosted model ignores human_only.
        for _label, _type in (("constitutive", "C"), ("immunoproteasome", "I")):
            PEPSICKLE_MODELS["pepsickle-%s-%s-%s" % (_family, _population, _label)] = dict(
                model_type=_family, human_only=_human, proteasome_type=_type)


def _runtime_executable(model_type, python_executable=None):
    value = python_executable
    if value is None and model_type == "in-vitro":
        value = os.environ.get("PEPSICKLE_GB_PYTHON")
    value = value or os.environ.get("PEPSICKLE_PYTHON")
    if not value:
        return None
    resolved = shutil.which(os.path.expanduser(str(value)))
    if resolved is None:
        raise RuntimeError("Pepsickle Python executable is not runnable: %s" % value)
    return os.path.abspath(resolved)


def _pepsickle_identity(human_only, model_type="epitope", python_executable=None,
                       timeout=PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS):
    """Return provenance from the interpreter that will actually predict."""
    cache_key = (bool(human_only), model_type, python_executable)
    if cache_key not in _identity_cache:
        if python_executable is None:
            identity = runtime_identity(human_only, model_type)
        else:
            request = dict(operation="identity", human_only=bool(human_only), model_type=model_type)
            try:
                completed = subprocess.run(
                    [python_executable, "-c", _PEPSICKLE_SUBPROCESS_SCRIPT],
                    input=json.dumps(request), text=True, capture_output=True, timeout=timeout)
                if completed.returncode:
                    raise RuntimeError("Pepsickle runtime inspection failed: %s" % completed.stderr.strip())
                identity = json.loads(completed.stdout)["identity"]
            except (OSError, subprocess.TimeoutExpired, ValueError, KeyError) as error:
                raise RuntimeError("Could not inspect Pepsickle runtime: %s" % error) from error
        _identity_cache[cache_key] = identity
    return _identity_cache[cache_key]


_PEPSICKLE_SUBPROCESS_SCRIPT = Path(__file__).with_name("pepsickle_runtime.py").read_text(encoding="utf-8")


class Pepsickle(ProteasomePredictor):
    """
    Proteasomal cleavage predictor using an explicit upstream model family.

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

    model_type : {"epitope", "in-vitro", "in-vitro-2"}
        Epitope-trained neural ensemble (default), digestion-trained gradient
        boosting, or digestion-trained neural ensemble. The gradient-boosted
        artifact requires its compatible upstream scikit-learn runtime.

    python_executable : str or None
        Separate prediction interpreter or Python-compatible launcher. Implies
        subprocess isolation. Defaults to PEPSICKLE_GB_PYTHON for gradient
        boosting, then PEPSICKLE_PYTHON, then the current interpreter.

    proteasome_type : {"C", "I"} or None
        Constitutive or immunoproteasome; required for digestion models.
        Must be None for the proteasome-type-agnostic epitope model.
    """

    def __init__(
            self,
            default_peptide_lengths=None,
            scoring=None,
            threshold=0.5,
            human_only=False,
            isolate_subprocess=False,
            subprocess_timeout=PEPSICKLE_SUBPROCESS_TIMEOUT_SECONDS,
            model_type="epitope",
            proteasome_type=None,
            python_executable=None):
        if model_type not in ("epitope", "in-vitro", "in-vitro-2"):
            raise ValueError("Unknown Pepsickle model_type %r" % model_type)
        if model_type == "epitope" and proteasome_type is not None:
            raise ValueError("The epitope model is proteasome-type agnostic")
        if model_type != "epitope" and proteasome_type not in ("C", "I"):
            raise ValueError("Digestion models require proteasome_type C or I")
        if model_type == "in-vitro" and human_only:
            raise ValueError("The gradient-boosted model does not support human_only")
        ProteasomePredictor.__init__(
            self,
            default_peptide_lengths=default_peptide_lengths,
            scoring=scoring,
        )
        self.threshold = threshold
        self.human_only = human_only
        self.python_executable = _runtime_executable(model_type, python_executable)
        self.isolate_subprocess = bool(isolate_subprocess or self.python_executable)
        self.subprocess_timeout = subprocess_timeout
        self.model_type = model_type
        self.proteasome_type = proteasome_type
        self._model = None

    def __str__(self):
        return "%s(scoring=%s, isolate_subprocess=%s)" % (
            self.__class__.__name__,
            getattr(self.scoring, "__name__", repr(self.scoring)),
            self.isolate_subprocess)

    def _predictor_name(self):
        if self.model_type != "epitope":
            return self._cleavage_model_from_identity(
                self.human_only, None, self.model_type, self.proteasome_type).name
        return "pepsickle"

    @staticmethod
    def _cleavage_model_from_identity(human_only, identity, model_type="epitope",
                                      proteasome_type=None):
        population = "human-only" if human_only else "all-mammal"
        if identity is None:
            version = "unresolved: pepsickle runtime assets not inspected or not available"
        else:
            version = (
                "package:%s;model:%s;weights-sha256:%s;inference-sha256:%s;"
                "features-sha256:%s;runtime-sha256:%s"
                % (
                    identity["package_version"],
                    identity["model_key"],
                    identity["weights_sha256"],
                    identity["inference_sha256"],
                    identity["features_sha256"],
                    hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
                )
            )
        if model_type == "epitope":
            name = "pepsickle-in-vivo-%s" % population
            assay = "Neural ensemble trained from observed epitope C termini"
            context = 8
        else:
            name = "pepsickle-%s-%s-%s" % (
                model_type, population,
                "constitutive" if proteasome_type == "C" else "immunoproteasome")
            assay = ("20S in-vitro digestion-trained %s; proteasome input %s" % (
                "gradient boosting" if model_type == "in-vitro" else "neural ensemble",
                proteasome_type))
            context = 3
        return CleavageModel(
            name=name,
            version=version,
            enzyme="proteasome epitope proxy" if model_type == "epitope" else "20S proteasome",
            uniprot="",
            species="Homo sapiens" if human_only else "Mammalia",
            compartments=("cytosol",),
            evidence="quantitative_model",
            references=(
                "https://doi.org/10.1093/bioinformatics/btab628",
                "https://github.com/pdxgx/pepsickle",
            ),
            assay=assay + "; %s training subset" % population,
            limitations=(
                "Native dimensionless model output, not an empirical cleavage, "
                "degradation, presentation, or vaccine efficacy probability. "
                "Uses a %d-residue context before P1 and after P1 with terminal "
                "padding. The upstream package deserializes a pickle model; "
                "isolate_subprocess controls process isolation but does not make "
                "an untrusted pickle safe." % context
            ),
            score_name="pepsickle_%s_model_output" % model_type.replace("-", "_"),
            score_units="dimensionless native model output",
            scored_endpoint="site_cleavage",
        )

    @classmethod
    def catalog_cleavage_model(cls, human_only=False, model_type="epitope", proteasome_type=None):
        """Describe a known optional model without claiming absent assets."""
        try:
            # A catalog/motif-only batch must not boot configured containers.
            # Selected predictors resolve their real identity before inference.
            identity = (None if _runtime_executable(model_type) else
                        _pepsickle_identity(human_only, model_type))
        except (ImportError, OSError, RuntimeError,
                importlib.metadata.PackageNotFoundError):
            identity = None
        return cls._cleavage_model_from_identity(human_only, identity, model_type, proteasome_type)

    def cleavage_model(self):
        """Return canonical per-bond metadata tied to installed assets."""
        return self._cleavage_model_from_identity(
            self.human_only, self._identity(),
            self.model_type, self.proteasome_type)

    def predict_cleavage(self, peptide):
        """Return every internal bond through the canonical cleavage contract.

        Pepsickle emits one score after each residue and forces its last score
        to zero because a sequence endpoint is not an internal peptide bond.
        This adapter maps array index ``i`` to bond ``i + 1`` and deliberately
        excludes that endpoint sentinel.
        """
        return self.predict_cleavage_many([peptide])[0]

    def predict_cleavage_many(self, peptides):
        """Assess canonical inputs with one model load and deduplicated inference."""
        peptides = tuple(coerce_peptide(p) for p in peptides)
        model = self.cleavage_model()
        eligible = {p.sequence for p in peptides if p.n_term == p.c_term == "free"}
        scores = self.cleavage_probs_many(sorted(eligible))
        return tuple(self._canonical_result(peptide, model, scores.get(peptide.sequence))
                     for peptide in peptides)

    def _canonical_result(self, peptide, model, scores):
        if peptide.n_term != "free" or peptide.c_term != "free":
            return CleavageResult(
                peptide,
                model,
                unsupported_reason=(
                    "Pepsickle does not model modified terminal chemistry"
                ),
            )
        if len(scores) != len(peptide.sequence):
            raise ValueError(
                "Expected %d pepsickle scores for sequence, got %d"
                % (len(peptide.sequence), len(scores)))
        sites = tuple(
            CleavageSite(
                bond,
                "scored",
                "Pepsickle %s model output after residue %d" % (self.model_type, bond),
                float(score),
            )
            for bond, score in enumerate(scores[:-1], start=1)
        )
        identity = self._identity()
        conditions = (
            ("human_only", str(bool(self.human_only)).lower()),
            ("model_mode", self.model_type),
            ("proteasome_type", self.proteasome_type or "agnostic"),
            ("model_key", identity["model_key"]),
            ("threshold", "%.17g" % self.threshold),
            ("isolate_subprocess", str(bool(self.isolate_subprocess)).lower()),
            ("subprocess_timeout_seconds", str(self.subprocess_timeout)),
            ("runtime", json.dumps(identity, sort_keys=True)),
        )
        return CleavageResult(peptide, model, sites, conditions=conditions)

    def _identity(self):
        return _pepsickle_identity(
            self.human_only, self.model_type, self.python_executable, self.subprocess_timeout)

    def _load_model(self):
        if self._model is None:
            cache_key = (self.human_only, self.model_type)
            if cache_key not in _model_cache:
                from pepsickle import model_functions
                initialize = {
                    "epitope": model_functions.initialize_epitope_model,
                    "in-vitro": model_functions.initialize_digestion_gb_model,
                    "in-vitro-2": model_functions.initialize_digestion_model,
                }[self.model_type]
                _model_cache[cache_key] = initialize(
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
            mod_type=self.model_type,
            proteasome_type=self.proteasome_type or "C",
            threshold=self.threshold,
        )
        return [entry[2] for entry in preds_raw]

    def _cleavage_probs_many_subprocess(self, sequences):
        expected_identity = self._identity()
        payload = json.dumps({
            "model_type": self.model_type,
            "proteasome_type": self.proteasome_type,
            "human_only": bool(self.human_only),
            "threshold": float(self.threshold),
            "sequences": sequences,
        })
        try:
            result = subprocess.run(
                [self.python_executable or sys.executable, "-c", _PEPSICKLE_SUBPROCESS_SCRIPT],
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
            if self.model_type == "in-vitro" and "sklearn.ensemble._gb_losses" in stderr_text:
                raise RuntimeError(
                    "Pepsickle's gradient-boosted artifact requires a compatible legacy "
                    "scikit-learn runtime; configure PEPSICKLE_GB_PYTHON (mhctools #471). "
                    "No alternative model was substituted.")
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

        if parsed.get("identity") != expected_identity:
            raise RuntimeError("Pepsickle runtime/assets changed between inspection and inference")
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

    def predict_many(self, peptides):
        """Return canonical results with one isolated inference call per batch."""
        return self.predictor.predict_cleavage_many(peptides)
