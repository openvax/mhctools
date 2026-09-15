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

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
from functools import lru_cache
from typing import Optional

import pandas as pd

from .peptide_input import PeptideInput


MHC_DEPENDENCE_VALUES = frozenset((
    "none",
    "single_allele",
    "haplotype",
))
"""Allowed ``kind_support()[kind]["mhc_dependence"]`` values."""


MHC_CLASS_VALUES = frozenset((
    "none",
    "I",
    "II",
    "both",
))
"""Allowed ``kind_support()[kind]["mhc_class"]`` values."""


class Kind:
    """String constants for prediction kinds.

    You can use ``Kind.pMHC_affinity`` or just ``"pMHC_affinity"`` —
    they're the same string. These constants name what is measured, but
    predictor instances define the MHC context required for their supported
    kinds through ``kind_support()``.
    """
    pMHC_affinity = "pMHC_affinity"
    pMHC_presentation = "pMHC_presentation"
    pMHC_stability = "pMHC_stability"
    pMHC_TCR_binding = "pMHC_TCR_binding"
    immunogenicity = "immunogenicity"
    antigen_processing = "antigen_processing"
    proteasome_cleavage = "proteasome_cleavage"
    endolysosomal_cleavage = "endolysosomal_cleavage"
    tap_transport = "tap_transport"
    erap_trimming = "erap_trimming"
    # Half-life of the parent peptide. MeasurementContext distinguishes a
    # defined solution, serum/plasma/whole blood, a cellular compartment, and
    # systemic in-vivo PK. This is deliberately distinct from pMHC_stability,
    # whose analyte is the assembled peptide-MHC complex.
    peptide_half_life = "peptide_half_life"

    # Deprecated input spellings retained for old callers and serialized data.
    # Prediction canonicalizes them to peptide_half_life and recovers the
    # matrix or systemic scope that used to be embedded in the kind.
    serum_half_life = "serum_half_life"
    plasma_half_life = "plasma_half_life"
    blood_half_life = "blood_half_life"
    systemic_elimination_half_life = "systemic_elimination_half_life"
    # Pharmacokinetic quantities whose units and interpretation depend on the
    # study and therefore live in MeasurementContext rather than VALUE_UNITS.
    systemic_clearance = "systemic_clearance"
    distribution_volume = "distribution_volume"
    systemic_exposure = "systemic_exposure"
    # CPP class confidence and quantitative uptake are deliberately distinct.
    cpp_classification = "cpp_classification"
    cellular_uptake = "cellular_uptake"
    tissue_concentration = "tissue_concentration"


CONTEXT_DEPENDENT_KINDS = frozenset((
    Kind.peptide_half_life,
    Kind.systemic_clearance,
    Kind.distribution_volume,
    Kind.systemic_exposure,
    Kind.cpp_classification,
    Kind.cellular_uptake,
    Kind.tissue_concentration,
))
"""Kinds for which mhctools intentionally defines no universal ordering."""

PHYSICAL_VALUE_KINDS = CONTEXT_DEPENDENT_KINDS - {
    Kind.cpp_classification,
    # May expose only a native score when conversion to a duration is unknown.
    Kind.peptide_half_life,
}
"""Context-dependent endpoints that require a units-bearing value."""


ESTIMATE_TYPE_VALUES = frozenset((
    "observed",
    "fitted",
    "simulated",
    "ml_predicted",
))
"""How an endpoint value was obtained."""


RESULT_STATUS_VALUES = frozenset((
    "available",
    "unsupported",
    "missing",
    "out_of_domain",
    "failed",
))
"""Availability states for a result, independent of estimate type."""


CONCENTRATION_BASIS_VALUES = frozenset(("total", "unbound"))
PK_SCOPE_VALUES = frozenset(("systemic", "apparent"))


@dataclass(frozen=True)
class MeasurementContext:
    """Versioned semantics shared by every prediction.

    ``unit`` and ``transform`` describe :attr:`Prediction.value`; the stored
    physical value remains linear, so ``transform`` must currently be
    ``"linear"`` when a value is present. Predictor-native or confidence
    outputs belong in ``score`` and are identified by ``score_semantics``.
    Ordinary model outputs receive a small cached default; assay-specific
    wrappers fill only the fields they actually know. Unknown descriptive
    fields are ``None``, never a plausible biological default.
    """

    estimate_type: str
    status: str = "available"
    analyte: str | None = None
    compartment: str | None = None
    matrix: str | None = None
    unit: str | None = None
    transform: str | None = None
    score_semantics: str | None = None
    class_label: str | None = None
    concentration_basis: str | None = None
    pk_scope: str | None = None
    timepoint: float | None = None
    time_unit: str | None = None
    time_origin: str | None = None
    series_id: str | None = None
    detail: str | None = None
    schema_version: int = 1

    def __post_init__(self):
        if self.schema_version != 1:
            raise ValueError(
                f"Unsupported MeasurementContext schema_version "
                f"{self.schema_version!r}")
        if self.estimate_type not in ESTIMATE_TYPE_VALUES:
            raise ValueError(
                f"estimate_type must be one of "
                f"{sorted(ESTIMATE_TYPE_VALUES)}, got {self.estimate_type!r}")
        if self.status not in RESULT_STATUS_VALUES:
            raise ValueError(
                f"status must be one of {sorted(RESULT_STATUS_VALUES)}, "
                f"got {self.status!r}")
        if (self.concentration_basis is not None and
                self.concentration_basis not in CONCENTRATION_BASIS_VALUES):
            raise ValueError(
                "concentration_basis must be 'total' or 'unbound', got "
                f"{self.concentration_basis!r}")
        if self.pk_scope is not None and self.pk_scope not in PK_SCOPE_VALUES:
            raise ValueError(
                "pk_scope must be 'systemic' or 'apparent', got "
                f"{self.pk_scope!r}")
        time_fields = (self.timepoint, self.time_unit, self.time_origin,
                       self.series_id)
        if any(v is not None for v in time_fields) and not all(
                v is not None for v in time_fields):
            raise ValueError(
                "timepoint, time_unit, time_origin and series_id must be "
                "provided together")

    def to_dict(self):
        """Serialize to a JSON-friendly dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, value):
        """Deserialize a context while ignoring forward-compatible fields."""
        if isinstance(value, cls):
            return intern_measurement_context(value)
        valid = {f.name for f in fields(cls)}
        return intern_measurement_context(
            cls(**{k: v for k, v in value.items() if k in valid}))


@lru_cache(maxsize=4096)
def intern_measurement_context(context):
    """Return a shared instance for an immutable measurement context."""
    if not isinstance(context, MeasurementContext):
        raise TypeError("context must be a MeasurementContext")
    return context


# Canonical "best direction" for each prediction field. Used by
# downstream aggregators (e.g. "best across alleles" or "best across
# methods") and by the :class:`PeptideResult` ``.best_*`` accessors.
#
# - ``score``: every kind uses higher = better as its numerical ordering
#   convention (binding strength, presentation likelihood, immunogenicity,
#   ...). Its scale and units remain predictor-specific.
# - ``percentile_rank``: 0 means best, smaller is better, every kind.
# - ``value``: kind-dependent — see :data:`VALUE_BEST_DIRECTIONS`.
FIELD_BEST_DIRECTIONS = {
    "score": "max",
    "percentile_rank": "min",
}

# The canonical unit of ``value`` for each kind that has one.
#
# ``value`` is always a physical quantity on a LINEAR scale in the unit named
# here — never a log, never a rescaling, never "whatever upstream printed".
# Wrappers convert. That rule predates these constants: affinity predictors
# commonly work in 1-log50k space internally and every affinity wrapper here
# inverts it to nM before filling ``value`` (see ``parse_stdout``'s
# ``50000 ** (1 - score)`` salvage and ``caphla._affinity_nm``), so a consumer
# can compare a NetMHCpan IC50 with an MHCflurry one without asking which
# transform each applied. PeptiVerse reports the hours produced by its
# upstream sequence model. PlifePred2's target transform and units are not
# established, so its wrapper leaves ``value`` empty unless a caller opts in
# to the inferred log10-seconds conversion.
#
# A predictor's native output is not lost, it just does not belong in a
# units-bearing field: wrappers keep it in their ``last_qc`` frame
# (``PlifePred2.last_qc["log10_seconds"]``, for instance).
#
# Kind and unit are independent: every prediction has a kind because every
# prediction measures something, but only some kinds have a canonical unit for
# ``value``. ``score`` has predictor-specific semantics and may itself carry
# units. Kinds absent from this mapping are exactly those without a canonical
# ``value`` unit.
VALUE_UNITS = {
    Kind.pMHC_affinity: "nM",
    Kind.pMHC_stability: "hours",
    Kind.tap_transport: "nM",
    Kind.peptide_half_life: "hours",
}

# Per-kind direction for ``value``. Whether higher or lower is "better" depends
# on what the unit in :data:`VALUE_UNITS` means. Add an entry alongside a
# ``VALUE_UNITS`` entry when introducing a new ``value``-bearing kind.
VALUE_BEST_DIRECTIONS = {
    Kind.pMHC_affinity: "min",   # IC50: tighter binding is a smaller number
    Kind.pMHC_stability: "max",  # pMHC complex dissociation half-life
    Kind.tap_transport: "min",   # predicted TAP-binding affinity
    Kind.peptide_half_life: "max",
}


_LEGACY_HALF_LIFE_CONTEXT = {
    Kind.serum_half_life: {"matrix": "serum"},
    Kind.plasma_half_life: {"matrix": "plasma"},
    Kind.blood_half_life: {"matrix": "whole blood"},
    Kind.systemic_elimination_half_life: {
        "compartment": "systemic circulation",
        "pk_scope": "systemic",
    },
}


def canonical_kind(kind):
    """Return the canonical kind for a current or legacy spelling."""
    if kind in _LEGACY_HALF_LIFE_CONTEXT:
        return Kind.peptide_half_life
    return kind


def _matches_legacy_half_life_context(prediction, kind):
    expected = _LEGACY_HALF_LIFE_CONTEXT[kind]
    context = prediction.measurement_context
    matrix = expected.get("matrix")
    if matrix is not None:
        observed = context.matrix
        if observed is None:
            return False
        observed = observed.lower()
        if observed != matrix and not observed.endswith(" " + matrix):
            return False
    return all(
        key == "matrix" or getattr(context, key) == value
        for key, value in expected.items())


@lru_cache(maxsize=256)
def _default_measurement_context(kind, value_is_present, status="available"):
    """Return the small shared context used by ordinary model predictions."""
    unit = VALUE_UNITS.get(kind) if value_is_present else None
    return intern_measurement_context(MeasurementContext(
        estimate_type="ml_predicted",
        status=status,
        unit=unit,
        transform="linear" if unit is not None else None,
    ))


@lru_cache(maxsize=64)
def _legacy_half_life_measurement_context(
        kind, value_is_present, status="available"):
    """Recover context that an old half-life kind encoded in its name."""
    return intern_measurement_context(MeasurementContext(
        estimate_type="ml_predicted",
        status=status,
        analyte="parent peptide",
        unit="hours" if status == "available" and value_is_present else None,
        transform=(
            "linear" if status == "available" and value_is_present else None),
        **_LEGACY_HALF_LIFE_CONTEXT[kind],
    ))


def value_unit(kind) -> Optional[str]:
    """Canonical unit of ``value`` for *kind*, or ``None`` if it has no value.

    The unit is always linear — ``"nM"``, ``"hours"`` — never a log-transformed
    scale. See :data:`VALUE_UNITS`.

    Examples
    --------
    >>> value_unit(Kind.pMHC_affinity)
    'nM'
    >>> value_unit(Kind.peptide_half_life)
    'hours'
    >>> value_unit(Kind.immunogenicity) is None
    True
    """
    return VALUE_UNITS.get(canonical_kind(kind))


def best_direction(kind, field) -> str:
    """Canonical "best" direction (``"max"`` or ``"min"``) for ``(kind, field)``.

    ``score`` (max) and ``percentile_rank`` (min) are uniform across kinds;
    ``value`` is kind-dependent (see :data:`VALUE_BEST_DIRECTIONS`). Raises
    ``ValueError`` for an unknown ``field``, or for ``value`` on a kind with no
    registered direction.
    """
    kind = canonical_kind(kind)
    if kind in CONTEXT_DEPENDENT_KINDS:
        raise ValueError(
            f"best_direction is context-dependent for {kind!r}; mhctools "
            "does not define a biological ranking policy for this endpoint")
    direction = FIELD_BEST_DIRECTIONS.get(field)
    if direction is not None:
        return direction
    if field == "value":
        if kind not in VALUE_BEST_DIRECTIONS:
            raise ValueError(
                f"best_direction undefined for ({kind!r}, 'value') — "
                f"`value` semantics depend on the kind. Add an entry "
                f"to mhctools.pred.VALUE_BEST_DIRECTIONS."
            )
        return VALUE_BEST_DIRECTIONS[kind]
    known = sorted(FIELD_BEST_DIRECTIONS) + ["value"]
    raise ValueError(
        f"best_direction undefined for field {field!r}. Known: {known}."
    )


def reduce_op(kind, field):
    """``max`` or ``min`` — the reducer that picks the best ``(kind, field)``."""
    return max if best_direction(kind, field) == "max" else min


COLUMNS = (
    "sample_name",
    "peptide",
    "n_flank",
    "c_flank",
    "source_sequence_name",
    "offset",
    "predictor_name",
    "predictor_version",
    "allele",
    "tcr",
    "kind",
    "score",
    "value",
    "percentile_rank",
    "measurement_context",
    "peptide_input",
    "cache_key",
)


@dataclass(frozen=True, repr=False)
class Prediction:
    """Single prediction from one model on one peptide."""
    kind: str
    score: float | None
    peptide: str = ""
    allele: str = ""
    tcr: str = ""
    n_flank: str = ""
    c_flank: str = ""
    value: Optional[float] = None
    percentile_rank: Optional[float] = None
    source_sequence_name: Optional[str] = None
    offset: int = 0
    predictor_name: str = ""
    predictor_version: str = ""
    measurement_context: MeasurementContext | None = None
    peptide_input: PeptideInput | None = None
    cache_key: str | None = None

    def __post_init__(self):
        original_kind = self.kind
        object.__setattr__(self, "kind", canonical_kind(original_kind))

        context = self.measurement_context
        if isinstance(context, Mapping):
            context = MeasurementContext.from_dict(context)
        elif context is not None and not isinstance(context, MeasurementContext):
            raise TypeError(
                "measurement_context must be MeasurementContext, a mapping, "
                "or None")
        elif context is not None:
            context = intern_measurement_context(context)
        else:
            status = (
                "available"
                if self.score is not None or self.value is not None
                else "missing")
            if original_kind in _LEGACY_HALF_LIFE_CONTEXT:
                context = _legacy_half_life_measurement_context(
                    original_kind, self.value is not None, status)
            else:
                context = _default_measurement_context(
                    self.kind, self.value is not None, status)
        object.__setattr__(self, "measurement_context", context)

        peptide_input = self.peptide_input
        if isinstance(peptide_input, Mapping):
            peptide_input = PeptideInput.from_dict(peptide_input)
            object.__setattr__(self, "peptide_input", peptide_input)
        elif peptide_input is not None and not isinstance(
                peptide_input, PeptideInput):
            raise TypeError(
                "peptide_input must be PeptideInput, a mapping, or None")
        if peptide_input is not None:
            if self.peptide and self.peptide != peptide_input.sequence:
                raise ValueError("Prediction peptide differs from peptide_input")
            if not self.peptide:
                object.__setattr__(self, "peptide", peptide_input.sequence)
            if (self.source_sequence_name is not None and
                    peptide_input.source_sequence_name is not None and
                    self.source_sequence_name !=
                    peptide_input.source_sequence_name):
                raise ValueError(
                    "Prediction source_sequence_name differs from peptide_input")
            if self.source_sequence_name is None:
                object.__setattr__(
                    self, "source_sequence_name",
                    peptide_input.source_sequence_name)
            if self.offset not in (0, peptide_input.source_start):
                raise ValueError("Prediction offset differs from peptide_input")
            if self.offset == 0 and peptide_input.source_start:
                object.__setattr__(self, "offset", peptide_input.source_start)
        if self.cache_key is not None:
            if not isinstance(self.cache_key, str) or not self.cache_key:
                raise ValueError("cache_key must be a nonempty string or None")
            if peptide_input is None:
                raise ValueError("cache_key requires peptide_input")
        if context.status == "available":
            if self.score is None and self.value is None:
                raise ValueError(
                    "available predictions require score or value")
            if self.kind in PHYSICAL_VALUE_KINDS and self.value is None:
                raise ValueError(
                    f"{self.kind} requires a quantitative value")
        elif self.score is not None or self.value is not None:
            raise ValueError(
                f"{context.status} predictions cannot carry score or value")
        if self.value is not None:
            if not context.unit:
                raise ValueError(
                    "measurement_context.unit is required when value is set")
            if context.transform != "linear":
                raise ValueError(
                    "Prediction.value must use the linear transform")
        if self.kind in CONTEXT_DEPENDENT_KINDS and self.allele:
            raise ValueError(
                f"{self.kind} is MHC-independent and cannot carry an allele")
        if self.kind == Kind.cpp_classification and self.value is not None:
            raise ValueError(
                "cpp_classification confidence belongs in score, not value")
        if (self.kind == Kind.cpp_classification and
                context.status == "available" and
                (not context.class_label or not context.score_semantics)):
            raise ValueError(
                "cpp_classification requires class_label and score_semantics")

    def __repr__(self):
        parts = [self.peptide or "?", self.kind]
        if self.allele:
            parts.insert(1, self.allele)
        if self.tcr:
            parts.insert(1, self.tcr)
        if self.score is not None:
            parts.append("score=%.4g" % self.score)
        if self.value is not None:
            parts.append("value=%.4g" % self.value)
        if self.percentile_rank is not None:
            parts.append("rank=%.2f%%" % self.percentile_rank)
        if self.predictor_name:
            parts.append(self.predictor_name)
        return "Prediction(%s)" % " | ".join(parts)

    def __str__(self):
        return repr(self)

    def to_row(self, sample_name=""):
        return {
            "sample_name": sample_name,
            "peptide": self.peptide,
            "n_flank": self.n_flank,
            "c_flank": self.c_flank,
            "source_sequence_name": self.source_sequence_name,
            "offset": self.offset,
            "predictor_name": self.predictor_name,
            "predictor_version": self.predictor_version,
            "allele": self.allele,
            "tcr": self.tcr,
            "kind": self.kind,
            "score": self.score,
            "value": self.value,
            "percentile_rank": self.percentile_rank,
            "measurement_context": (
                self.measurement_context.to_dict()
                if self.measurement_context is not None else None),
            "peptide_input": (
                self.peptide_input.to_dict()
                if self.peptide_input is not None else None),
            "cache_key": self.cache_key,
        }

    def to_dict(self):
        """Serialize to a JSON-friendly dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a dict (as produced by :meth:`to_dict`)."""
        valid = {f.name for f in fields(cls)}
        values = {k: v for k, v in d.items() if k in valid}
        context = values.get("measurement_context")
        if context is not None:
            values["measurement_context"] = MeasurementContext.from_dict(
                context)
        peptide_input = values.get("peptide_input")
        if peptide_input is not None:
            values["peptide_input"] = PeptideInput.from_dict(peptide_input)
        return cls(**values)


@dataclass(repr=False)
class PeptideResult:
    """All predictions for one peptide (a tuple of ``Prediction`` objects)."""
    preds: tuple[Prediction, ...] = ()

    def __repr__(self):
        if not self.preds:
            return "PeptideResult(empty)"
        from collections import Counter
        kinds = Counter(p.kind for p in self.preds)
        kind_str = ", ".join(
            "%d\u00d7%s" % (n, k) for k, n in kinds.items())
        return "PeptideResult(%s, %d preds: %s)" % (
            self.peptide or "?", len(self.preds), kind_str)

    def __str__(self):
        return repr(self)

    # --- shared fields (same across all Preds in this result) ---

    @property
    def peptide(self) -> str:
        return self.preds[0].peptide if self.preds else ""

    @property
    def offset(self) -> int:
        return self.preds[0].offset if self.preds else 0

    @property
    def source_sequence_name(self) -> Optional[str]:
        return self.preds[0].source_sequence_name if self.preds else None

    @property
    def kinds(self) -> set:
        """Set of Kind values present in this result."""
        return {p.kind for p in self.preds}

    @property
    def alleles(self) -> set:
        """Set of allele strings present in this result."""
        return {p.allele for p in self.preds if p.allele}

    @property
    def tcrs(self) -> set:
        """Set of TCR identifiers present in this result."""
        return {p.tcr for p in self.preds if p.tcr}

    # --- kind accessors (best by score, wrapped for safe field access) ---

    @property
    def affinity(self) -> Optional[Prediction]:
        """Best affinity prediction, or None."""
        return self.best_by_score(Kind.pMHC_affinity)

    @property
    def presentation(self) -> Optional[Prediction]:
        """Best presentation prediction, or None."""
        return self.best_by_score(Kind.pMHC_presentation)

    @property
    def stability(self) -> Optional[Prediction]:
        """Best stability prediction, or None."""
        return self.best_by_score(Kind.pMHC_stability)

    @property
    def immunogenicity(self) -> Optional[Prediction]:
        """Best immunogenicity prediction, or None."""
        return self.best_by_score(Kind.immunogenicity)

    @property
    def processing(self) -> Optional[Prediction]:
        """Best antigen-processing prediction, or None."""
        return self.best_by_score(Kind.antigen_processing)

    @property
    def cleavage(self) -> Optional[Prediction]:
        """Best proteasomal cleavage prediction, or None."""
        return self.best_by_score(Kind.proteasome_cleavage)

    @property
    def endolysosomal_cleavage(self) -> Optional[Prediction]:
        """Best endolysosomal (MHC-II) cleavage prediction, or None."""
        return self.best_by_score(Kind.endolysosomal_cleavage)

    @property
    def tap_transport(self) -> Optional[Prediction]:
        """Best TAP transport prediction, or None."""
        return self.best_by_score(Kind.tap_transport)

    @property
    def erap_trimming(self) -> Optional[Prediction]:
        """Best ERAP1 trimming prediction, or None."""
        return self.best_by_score(Kind.erap_trimming)

    def _half_life_for_matrix(self, matrix) -> Optional[Prediction]:
        """Best peptide half-life restricted to one specimen matrix."""
        matrix = matrix.lower()
        candidates = [
            pred for pred in self.preds
            if pred.kind == Kind.peptide_half_life
            and pred.score is not None
            and pred.measurement_context.matrix is not None
            and (pred.measurement_context.matrix.lower() == matrix
                 or pred.measurement_context.matrix.lower().endswith(
                     " " + matrix))
        ]
        return max(candidates, key=lambda pred: pred.score, default=None)

    @property
    def peptide_half_life(self) -> Optional[Prediction]:
        """Peptide half-life when all available results share one context.

        Comparing different matrices or systemic scopes is intentionally
        rejected; callers can select those contexts explicitly instead.
        """
        candidates = [
            pred for pred in self.preds
            if pred.kind == Kind.peptide_half_life and pred.score is not None]
        if not candidates:
            return None
        contexts = {pred.measurement_context for pred in candidates}
        if len(contexts) != 1:
            raise ValueError(
                "peptide_half_life results span multiple measurement contexts")
        return max(candidates, key=lambda pred: pred.score)

    @property
    def serum_half_life(self) -> Optional[Prediction]:
        """Best peptide half-life measured in serum, or None."""
        return self._half_life_for_matrix("serum")

    @property
    def plasma_half_life(self) -> Optional[Prediction]:
        """Best peptide half-life measured in plasma, or None."""
        return self._half_life_for_matrix("plasma")

    @property
    def blood_half_life(self) -> Optional[Prediction]:
        """Best peptide half-life measured in whole blood, or None."""
        return self._half_life_for_matrix("whole blood")

    @property
    def tcr_binding(self) -> Optional[Prediction]:
        """Best pMHC:TCR binding prediction, or None."""
        return self.best_by_score(Kind.pMHC_TCR_binding)

    # backward compat aliases
    @property
    def best_affinity(self) -> Optional[Prediction]:
        return self.affinity

    @property
    def best_presentation(self) -> Optional[Prediction]:
        return self.presentation

    @property
    def best_stability(self) -> Optional[Prediction]:
        return self.stability

    # --- best by rank (return raw Prediction) ---

    @property
    def best_affinity_by_rank(self) -> Optional[Prediction]:
        return self.best_by_rank(Kind.pMHC_affinity)

    @property
    def best_presentation_by_rank(self) -> Optional[Prediction]:
        return self.best_by_rank(Kind.pMHC_presentation)

    @property
    def best_stability_by_rank(self) -> Optional[Prediction]:
        return self.best_by_rank(Kind.pMHC_stability)

    # --- filtering ---

    def filter(self, kind=None, allele=None):
        """Filter preds. None means don't filter on that field."""
        legacy_kind = (
            kind if kind in _LEGACY_HALF_LIFE_CONTEXT else None)
        if kind is not None:
            kind = canonical_kind(kind)
        return [p for p in self.preds
                if (kind is None or p.kind == kind)
                and (legacy_kind is None
                     or _matches_legacy_half_life_context(p, legacy_kind))
                and (allele is None or p.allele == allele)]

    # --- serialization ---

    def to_dict(self):
        """Serialize to a JSON-friendly dict."""
        return {"preds": [p.to_dict() for p in self.preds]}

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a dict (as produced by :meth:`to_dict`)."""
        return cls(preds=tuple(Prediction.from_dict(p) for p in d["preds"]))

    # --- dataframe ---

    def to_dataframe(self, sample_name=""):
        rows = [p.to_row(sample_name) for p in self.preds]
        if not rows:
            return pd.DataFrame(columns=COLUMNS)
        return pd.DataFrame(rows, columns=COLUMNS)

    # --- best_by (public, direction-aware) ---

    def best_by(self, kind, field) -> Optional[Prediction]:
        """Return the prediction of ``kind`` with the best ``field``, or None.

        "Best" direction comes from :func:`best_direction` —
        ``score`` is max-better, ``percentile_rank`` is min-better, and
        ``value`` is kind-dependent (IC50 lower-better, half-life higher-better).

        Predictions with ``None`` for ``field`` are skipped. Predictions
        with an allele are preferred; if none of the matching kind have an
        allele, falls back to allele-less predictions (e.g. processing
        predictors that emit allele-independent scores).
        """
        original_kind = kind
        kind = canonical_kind(kind)
        if original_kind in _LEGACY_HALF_LIFE_CONTEXT:
            op = min if field == "percentile_rank" else max
        else:
            op = reduce_op(kind, field)

        def has_value(p):
            return getattr(p, field) is not None

        matching = self.filter(kind=original_kind)
        with_allele = [p for p in matching if p.allele and has_value(p)]
        if with_allele:
            return op(with_allele, key=lambda p: getattr(p, field))
        without_allele = [p for p in matching if has_value(p)]
        if without_allele:
            return op(without_allele, key=lambda p: getattr(p, field))
        return None

    def best_by_score(self, kind) -> Optional[Prediction]:
        """Best prediction of ``kind`` by ``score`` (max-better)."""
        return self.best_by(kind, "score")

    def best_by_rank(self, kind) -> Optional[Prediction]:
        """Best prediction of ``kind`` by ``percentile_rank`` (min-better)."""
        return self.best_by(kind, "percentile_rank")

    def best_by_value(self, kind) -> Optional[Prediction]:
        """Best prediction of ``kind`` by ``value``. Direction is kind-specific
        (see :data:`VALUE_BEST_DIRECTIONS`). Raises ``ValueError`` for kinds
        without a registered ``value`` direction."""
        return self.best_by(kind, "value")


def preds_from_rows(rows, **shared):
    """Build a PeptideResult from dicts, with shared fields filled in.

    Example::

        preds_from_rows(
            [
                dict(kind=Kind.pMHC_affinity, allele="HLA-A*02:01",
                     score=0.85, value=120.5, percentile_rank=0.8),
                dict(kind=Kind.pMHC_presentation, allele="HLA-A*02:01",
                     score=0.92, percentile_rank=0.3),
            ],
            peptide="SIINFEKL",
            predictor_name="netMHCpan",
            predictor_version="4.1",
        )
    """
    return PeptideResult(preds=tuple(
        Prediction(**{**shared, **row}) for row in rows
    ))


# Backward compatibility
Pred = Prediction
PeptidePreds = PeptideResult
