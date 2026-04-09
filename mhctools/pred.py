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
from enum import Enum
from typing import Optional

import pandas as pd


class Kind(Enum):
    """What biological quantity is being predicted."""
    # Peptide-MHC
    pMHC_affinity = "pMHC_affinity"
    pMHC_presentation = "pMHC_presentation"
    pMHC_stability = "pMHC_stability"
    immunogenicity = "immunogenicity"
    # Processing pathway
    antigen_processing = "antigen_processing"
    proteasome_cleavage = "proteasome_cleavage"
    tap_transport = "tap_transport"
    erap_trimming = "erap_trimming"


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
    "kind",
    "score",
    "value",
    "percentile_rank",
)


@dataclass(frozen=True, repr=False)
class Pred:
    """Single prediction from one model on one peptide. Self-contained."""
    kind: Kind
    score: float
    peptide: str = ""
    allele: str = ""
    n_flank: str = ""
    c_flank: str = ""
    value: Optional[float] = None
    percentile_rank: Optional[float] = None
    source_sequence_name: Optional[str] = None
    offset: int = 0
    predictor_name: str = ""
    predictor_version: str = ""

    def __repr__(self):
        parts = [self.peptide or "?", self.kind.value]
        if self.allele:
            parts.insert(1, self.allele)
        parts.append("score=%.4g" % self.score)
        if self.value is not None:
            parts.append("value=%.4g" % self.value)
        if self.percentile_rank is not None:
            parts.append("rank=%.2f%%" % self.percentile_rank)
        if self.predictor_name:
            parts.append(self.predictor_name)
        return "Pred(%s)" % " | ".join(parts)

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
            "kind": self.kind.value,
            "score": self.score,
            "value": self.value,
            "percentile_rank": self.percentile_rank,
        }

    def to_dict(self):
        """Serialize to a JSON-friendly dict."""
        d = asdict(self)
        d["kind"] = self.kind.value
        return d

    @classmethod
    def from_dict(cls, d):
        """Deserialize from a dict (as produced by :meth:`to_dict`)."""
        d = dict(d)
        d["kind"] = Kind(d["kind"])
        # Only pass fields that Pred knows about
        valid = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in valid})


@dataclass(repr=False)
class PeptideResult:
    """All predictions for one peptide. Contains a tuple of Pred objects."""
    preds: tuple[Pred, ...] = ()

    def __repr__(self):
        if not self.preds:
            return "PeptideResult(empty)"
        peptide = self.preds[0].peptide or "?"
        from collections import Counter
        kinds = Counter(p.kind.value for p in self.preds)
        kind_str = ", ".join(
            "%d\u00d7%s" % (n, k) for k, n in kinds.items())
        return "PeptideResult(%s, %d preds: %s)" % (
            peptide, len(self.preds), kind_str)

    def __str__(self):
        return repr(self)

    # --- best allele by score (higher = better) ---

    @property
    def best_affinity(self) -> Optional[Pred]:
        return self._best_by_score(Kind.pMHC_affinity)

    @property
    def best_presentation(self) -> Optional[Pred]:
        return self._best_by_score(Kind.pMHC_presentation)

    @property
    def best_stability(self) -> Optional[Pred]:
        return self._best_by_score(Kind.pMHC_stability)

    # --- best allele by rank (lower = better) ---

    @property
    def best_affinity_by_rank(self) -> Optional[Pred]:
        return self._best_by_rank(Kind.pMHC_affinity)

    @property
    def best_presentation_by_rank(self) -> Optional[Pred]:
        return self._best_by_rank(Kind.pMHC_presentation)

    @property
    def best_stability_by_rank(self) -> Optional[Pred]:
        return self._best_by_rank(Kind.pMHC_stability)

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
        return cls(preds=tuple(Pred.from_dict(p) for p in d["preds"]))

    # --- dataframe ---

    def to_dataframe(self, sample_name=""):
        rows = [p.to_row(sample_name) for p in self.preds]
        if not rows:
            return pd.DataFrame(columns=COLUMNS)
        return pd.DataFrame(rows, columns=COLUMNS)

    # --- internals ---

    def _best_by_score(self, kind) -> Optional[Pred]:
        candidates = [p for p in self.preds if p.kind == kind and p.allele]
        return max(candidates, key=lambda p: p.score) if candidates else None

    def _best_by_rank(self, kind) -> Optional[Pred]:
        candidates = [p for p in self.preds
                      if p.kind == kind and p.allele
                      and p.percentile_rank is not None]
        return min(candidates, key=lambda p: p.percentile_rank) if candidates else None


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
        Pred(**{**shared, **row}) for row in rows
    ))


# Backward compatibility
PeptidePreds = PeptideResult
