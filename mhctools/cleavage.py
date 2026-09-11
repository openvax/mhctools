"""Per-bond peptidase evidence with explicit chemistry and native scores.

Bond ``b`` splits ``sequence[:b] | sequence[b:]``; only 1..length-1 are
peptide bonds. Motif matches are recognition rules, never probabilities.
"""

from dataclasses import asdict, dataclass
import math
from typing import Optional, Tuple


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# The single source of truth for which CleavageModel.evidence values exist
# and which CleavageSite.status values each one legally carries. Keying the
# evidence-membership check off these same keys (below, in CleavageModel)
# means a future evidence type only needs adding here once.
_ALLOWED_SITE_STATUSES = {"quantitative_model": {"scored"}, "motif_rule": {"matched", "not_matched"},
                          "substrate_reference": {"reported"}}


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def coerce_peptide(peptide):
    """Accept a canonical peptide string or pass a :class:`CleavageInput` through.

    Every predictor's ``predict`` entry point does exactly this coercion;
    sharing it keeps the accepted-input contract identical everywhere.
    """
    if isinstance(peptide, str):
        peptide = CleavageInput(peptide)
    if not isinstance(peptide, CleavageInput):
        raise TypeError("Expected CleavageInput or canonical peptide string")
    return peptide


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

    ``motif_strictness`` grades how much of the enzyme's real specificity a
    motif rule captures, so a caller can tell how far to trust a decision:

    ``required``
        The source establishes the pattern as necessary for this activity,
        so a non-match is meaningful evidence against cleavage by this
        enzyme through this route. It is still not proof of resistance.
    ``preferred``
        The source reports a favoured context. Non-matching bonds can be
        cleaved, usually more slowly, so a non-match is weak evidence.
    ``permissive``
        A broad flag with many exceptions. A match constrains little and a
        non-match almost nothing.

    ``strictness_basis`` names the source observation behind that grade, so
    the grade is attributable rather than an opinion. Both fields belong to
    motif rules only; scored models and source references carry neither.

    ``scored_endpoint`` names the :mod:`mhctools.benchmark` endpoint that a
    ``quantitative_model``'s native score answers (for example
    ``substrate_depletion`` or ``site_cleavage``), so a caller maps a scored
    site to the right measurement type from the model's own declaration
    rather than a hardcoded model name. Only quantitative models carry it.
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
    scored_endpoint: Optional[str] = None
    motif_strictness: Optional[str] = None
    strictness_basis: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "compartments", tuple(self.compartments))
        object.__setattr__(self, "references", tuple(self.references))
        if self.evidence not in _ALLOWED_SITE_STATUSES:
            raise ValueError("Unknown cleavage evidence type")
        if not self.references:
            raise ValueError("Every model must cite at least one source")
        if self.evidence == "quantitative_model":
            if not self.score_name or not self.score_units:
                raise ValueError("Quantitative models must name their native score and units")
            if not self.scored_endpoint:
                raise ValueError(
                    "Quantitative models must name the benchmark endpoint their score answers")
        elif self.score_name is not None or self.score_units is not None:
            raise ValueError("Non-quantitative evidence does not have numerical scores")
        elif self.scored_endpoint is not None:
            raise ValueError("Only quantitative models carry a benchmark endpoint for their score")
        if self.evidence == "motif_rule":
            if self.motif_strictness not in ("required", "preferred", "permissive"):
                raise ValueError(
                    "Motif rules must grade strictness as required, preferred or permissive")
            if not self.strictness_basis:
                raise ValueError("A strictness grade requires the source observation behind it")
        elif self.motif_strictness is not None or self.strictness_basis is not None:
            raise ValueError("Only motif rules carry a strictness grade")


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
        if self.status not in ("matched", "not_matched", "scored", "reported"):
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
    conditions: Tuple[Tuple[str, str], ...] = ()
    substrate_observation: Optional[str] = None

    def __post_init__(self):
        object.__setattr__(self, "sites", tuple(self.sites))
        conditions = tuple(tuple(pair) for pair in self.conditions)
        if (any(len(pair) != 2 or not all(isinstance(v, str) for v in pair)
                for pair in conditions) or len(dict(conditions)) != len(conditions)):
            raise ValueError("Conditions require distinct string key/value pairs")
        object.__setattr__(self, "conditions", conditions)
        if self.substrate_observation is not None:
            if (self.model.evidence != "substrate_reference" or self.substrate_observation not in (
                    "cleavage_reported", "no_cleavage_detected")):
                raise ValueError("Substrate observations require source-reference evidence")
            if self.substrate_observation == "no_cleavage_detected" and self.sites:
                raise ValueError("Whole-substrate non-cleavage must not create site labels")
            if self.unsupported_reason:
                raise ValueError("Unsupported inputs cannot carry substrate observations")
        if self.unsupported_reason is not None and self.sites:
            raise ValueError("Unsupported results cannot contain scored sites")
        bonds = [site.bond for site in self.sites]
        if len(set(bonds)) != len(bonds) or any(
                b >= len(self.peptide.sequence) for b in bonds):
            raise ValueError("Sites must identify distinct internal peptide bonds")
        for site in self.sites:
            # self.model is a validated CleavageModel, so this lookup cannot miss.
            if site.status not in _ALLOWED_SITE_STATUSES[self.model.evidence]:
                raise ValueError("Site values must agree with model evidence semantics")

    def to_dict(self):
        """Return JSON-compatible data including original source bond offsets."""
        result = asdict(self)
        result["sites"] = [
            dict(asdict(site), source_bond=self.peptide.source_start + site.bond)
            for site in self.sites]
        return result
