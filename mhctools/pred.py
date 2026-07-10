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

from dataclasses import asdict, dataclass, fields
from typing import Optional

import pandas as pd


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


# Canonical "best direction" for each prediction field. Used by
# downstream aggregators (e.g. "best across alleles" or "best across
# methods") and by the :class:`PeptideResult` ``.best_*`` accessors.
#
# - ``score``: every kind that uses score normalizes higher = better
#   (binding strength, presentation likelihood, immunogenicity, ...).
# - ``percentile_rank``: 0 means best, smaller is better, every kind.
# - ``value``: kind-dependent — see :data:`VALUE_BEST_DIRECTIONS`.
FIELD_BEST_DIRECTIONS = {
    "score": "max",
    "percentile_rank": "min",
}

# Per-kind direction for ``value``. ``value`` carries the predictor's
# raw output unit (IC50 nM for affinity, half-life for stability, ...);
# whether higher or lower is "better" depends on what the unit means.
# Add an entry when introducing a new ``value``-bearing kind.
VALUE_BEST_DIRECTIONS = {
    Kind.pMHC_affinity: "min",   # IC50 nM
    Kind.pMHC_stability: "max",  # half-life
    Kind.tap_transport: "min",   # predicted TAP-binding affinity, nM
}


def best_direction(kind, field) -> str:
    """Canonical "best" direction (``"max"`` or ``"min"``) for ``(kind, field)``.

    ``score`` (max) and ``percentile_rank`` (min) are uniform across kinds;
    ``value`` is kind-dependent (see :data:`VALUE_BEST_DIRECTIONS`). Raises
    ``ValueError`` for an unknown ``field``, or for ``value`` on a kind with no
    registered direction.
    """
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
)


@dataclass(frozen=True, repr=False)
class Prediction:
    """Single prediction from one model on one peptide."""
    kind: str
    score: float
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

    def __repr__(self):
        parts = [self.peptide or "?", self.kind]
        if self.allele:
            parts.insert(1, self.allele)
        if self.tcr:
            parts.insert(1, self.tcr)
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
        }

    def to_dict(self):
        """Serialize to a JSON-friendly dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a dict (as produced by :meth:`to_dict`)."""
        valid = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in valid})


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
        return [p for p in self.preds
                if (kind is None or p.kind == kind)
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
        op = reduce_op(kind, field)

        def has_value(p):
            return getattr(p, field) is not None

        with_allele = [p for p in self.preds
                       if p.kind == kind and p.allele and has_value(p)]
        if with_allele:
            return op(with_allele, key=lambda p: getattr(p, field))
        without_allele = [p for p in self.preds
                          if p.kind == kind and has_value(p)]
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
