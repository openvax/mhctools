"""Per-bond peptidase evidence with explicit chemistry and native scores.

Bond ``b`` splits ``sequence[:b] | sequence[b:]``; only 1..length-1 are
peptide bonds. Motif matches are recognition rules, never probabilities.
"""

from dataclasses import asdict, dataclass
import math
from typing import Optional, Tuple


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


@dataclass(frozen=True)
class CleavageInput:
    """Linear, canonical L-peptide; defaults explicitly assume free termini.

    Other chemistry (D-residues, cyclization, conjugation, disulfides, etc.)
    is outside this input type. No sequence normalization is performed.
    ``source_start`` is a zero-based offset in an optional parent sequence.
    """

    sequence: str
    n_term: str = "free"
    c_term: str = "free"
    source_id: Optional[str] = None
    source_start: int = 0

    def __post_init__(self):
        if (not isinstance(self.sequence, str) or not self.sequence or
                not set(self.sequence) <= AMINO_ACIDS):
            raise ValueError("Expected a nonempty canonical uppercase L-peptide")
        if self.n_term not in ("free", "acetylated", "unknown"):
            raise ValueError("n_term must be free, acetylated or unknown")
        if self.c_term not in ("free", "amidated", "unknown"):
            raise ValueError("c_term must be free, amidated or unknown")
        if not _integer(self.source_start) or self.source_start < 0:
            raise ValueError("source_start must be a nonnegative integer")
        if self.source_id is not None and not isinstance(self.source_id, str):
            raise ValueError("source_id must be a string or None")

    def fragment(self, start, end, *, n_term, c_term):
        """Describe a conditional fragment, explicitly supplying its chemistry.

        This does not predict whether or when the fragment is produced.
        Retained parent termini must retain their original chemistry.
        """
        if (not _integer(start) or not _integer(end) or
                not 0 <= start < end <= len(self.sequence)):
            raise ValueError("Invalid fragment interval")
        if start == 0 and n_term != self.n_term:
            raise ValueError("Retained N terminus must preserve chemistry")
        if end == len(self.sequence) and c_term != self.c_term:
            raise ValueError("Retained C terminus must preserve chemistry")
        return CleavageInput(
            self.sequence[start:end], n_term, c_term, self.source_id,
            self.source_start + start)


@dataclass(frozen=True)
class CleavageModel:
    """Model provenance and applicability; compartments are enzyme locations.

    Compartment annotation is not evidence of calibration in that matrix.
    ``assay`` describes the source of the specificity evidence.
    """

    name: str
    version: str
    enzyme: str
    uniprot: str
    species: str
    compartments: Tuple[str, ...]
    evidence: str
    references: Tuple[str, ...]
    assay: str
    limitations: str
    score_name: Optional[str] = None
    score_units: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "compartments", tuple(self.compartments))
        object.__setattr__(self, "references", tuple(self.references))
        if self.evidence not in ("motif_rule", "quantitative_model"):
            raise ValueError("Unknown cleavage evidence type")
        if self.evidence == "quantitative_model":
            if not self.score_name or not self.score_units:
                raise ValueError("Quantitative models must name their native score and units")
        elif self.score_name is not None or self.score_units is not None:
            raise ValueError("Motif rules do not have numerical scores")


@dataclass(frozen=True)
class CleavageSite:
    """An assessed peptide bond, with a motif decision or a native score."""

    bond: int
    status: str
    reason: str
    score: Optional[float] = None

    def __post_init__(self):
        if not _integer(self.bond) or self.bond < 1:
            raise ValueError("bond must be a positive integer")
        if self.status not in ("matched", "not_matched", "scored"):
            raise ValueError("Unknown cleavage site status")
        if self.status == "scored":
            if (isinstance(self.score, bool) or
                    not isinstance(self.score, (float, int)) or
                    not math.isfinite(self.score)):
                raise ValueError("A scored site requires a finite native score")
        elif self.score is not None:
            raise ValueError("A motif decision cannot carry a numerical score")


@dataclass(frozen=True)
class CleavageResult:
    """Evidence for one input/model; unsupported inputs never receive zeros.

    ``not_matched`` means the limited recognition rule did not match. It
    does not establish resistance. Empty sites with ``unsupported_reason``
    mean no assessment could be made, rather than absence of cleavage.
    """

    peptide: CleavageInput
    model: CleavageModel
    sites: Tuple[CleavageSite, ...] = ()
    unsupported_reason: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "sites", tuple(self.sites))
        if self.unsupported_reason is not None and self.sites:
            raise ValueError("Unsupported results cannot contain scored sites")
        bonds = [site.bond for site in self.sites]
        if len(set(bonds)) != len(bonds) or any(
                b >= len(self.peptide.sequence) for b in bonds):
            raise ValueError("Sites must identify distinct internal peptide bonds")
        for site in self.sites:
            if (site.status == "scored") != (self.model.evidence == "quantitative_model"):
                raise ValueError("Site values must agree with model evidence semantics")

    def to_dict(self):
        """Return JSON-compatible data including original source bond offsets."""
        result = asdict(self)
        result["sites"] = [
            dict(asdict(site), source_bond=self.peptide.source_start + site.bond)
            for site in self.sites]
        return result
