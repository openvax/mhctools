"""Per-bond peptidase evidence with explicit chemistry and native scores.

Bond ``b`` splits ``sequence[:b] | sequence[b:]``; only 1..length-1 are
peptide bonds. Motif matches are recognition rules, never probabilities.
"""

from dataclasses import asdict, dataclass, fields
import math
from typing import Optional, Tuple

from .peptide_input import PeptideInput

# The single source of truth for which CleavageModel.evidence values exist
# and which CleavageSite.status values each one legally carries. Keying the
# evidence-membership check off these same keys (below, in CleavageModel)
# means a future evidence type only needs adding here once.
_ALLOWED_SITE_STATUSES = {"quantitative_model": {"scored"}, "motif_rule": {"matched", "not_matched"},
                          "substrate_reference": {"reported"}}


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def coerce_peptide(peptide):
    """Accept a canonical sequence or exact peptide without losing fields."""
    if isinstance(peptide, str):
        peptide = PeptideInput(peptide)
    if not isinstance(peptide, PeptideInput):
        raise TypeError("Expected PeptideInput or canonical peptide string")
    if peptide.attachments:
        raise ValueError("Cleavage models do not support peptide attachments")
    return peptide


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

    ``strictness_basis`` names the source observation behind an optional
    grade, so the grade is attributable rather than an opinion. Both fields
    belong to motif rules only; scored models and source references carry
    neither.

    ``scored_endpoint`` names the :mod:`mhctools.benchmark` endpoint that a
    ``quantitative_model``'s native score answers (for example
    ``substrate_depletion`` or ``site_cleavage``), so a caller maps a scored
    site to the right measurement type from the model's own declaration
    rather than a hardcoded model name. Only quantitative models carry it.
    """

    name: str
    version: str
    enzyme: str
    uniprot: Optional[str]
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
    rule: Optional[str] = None

    def to_dict(self):
        """Keep only metadata applicable to this evidence type."""
        return {key: value for key, value in asdict(self).items() if value is not None}

    @classmethod
    def from_dict(cls, value):
        """Restore model metadata, including rules without a UniProt ID."""
        value = dict(value)
        value.setdefault("uniprot", None)
        return cls(**value)

    def __post_init__(self):
        object.__setattr__(self, "compartments", tuple(self.compartments))
        object.__setattr__(self, "references", tuple(self.references))
        if self.evidence not in _ALLOWED_SITE_STATUSES:
            raise ValueError("Unknown cleavage evidence type")
        if not self.references:
            raise ValueError("Every model must cite at least one source")
        if not all(isinstance(r, str) and r.startswith(("https://", "http://")) for r in self.references):
            raise ValueError("Every reference must be a source URL")
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
            if self.motif_strictness is not None and self.motif_strictness not in (
                    "required", "preferred", "permissive"):
                raise ValueError(
                    "Motif rules must grade strictness as required, preferred or permissive")
            if self.motif_strictness is None and self.strictness_basis is not None:
                raise ValueError("A strictness basis requires a grade")
            if self.motif_strictness is not None and not self.strictness_basis:
                raise ValueError("A strictness grade requires the source observation behind it")
            if self.rule is not None and (
                    not isinstance(self.rule, str) or not self.rule.strip()):
                raise ValueError("A motif rule description must be nonempty")
        elif self.motif_strictness is not None or self.strictness_basis is not None:
            raise ValueError("Only motif rules carry a strictness grade")
        elif self.rule is not None:
            raise ValueError("Only motif rules carry a rule description")


@dataclass(frozen=True)
class CleavageSite:
    """An assessed peptide bond, with a motif decision or a native score."""

    bond: int
    status: str
    reason: str
    score: Optional[float] = None

    @property
    def match(self):
        """Binary motif-rule indicator, distinct from a cleavage probability."""
        if self.status in ("matched", "not_matched"):
            return int(self.status == "matched")
        return None

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

    peptide: PeptideInput
    model: CleavageModel
    sites: Tuple[CleavageSite, ...] = ()
    unsupported_reason: Optional[str] = None
    conditions: Tuple[Tuple[str, str], ...] = ()
    substrate_observation: Optional[str] = None

    @property
    def pattern_mask(self):
        """Bond-aligned motif matches; ``None`` marks an unassessed bond.

        Entry ``b - 1`` describes the bond after residue ``b``. Quantitative
        and source-reference models have no motif mask and return ``None``.
        """
        if self.model.evidence != "motif_rule":
            return None
        mask = [None] * (len(self.peptide.sequence) - 1)
        for site in self.sites:
            mask[site.bond - 1] = site.match
        return tuple(mask)

    def __post_init__(self):
        if not isinstance(self.peptide, PeptideInput):
            raise TypeError("CleavageResult requires PeptideInput")
        if self.peptide.attachments and self.unsupported_reason is None:
            raise ValueError("Attached peptides require an unsupported result")
        object.__setattr__(self, "sites", tuple(self.sites))
        if self.unsupported_reason is not None and (
                not isinstance(self.unsupported_reason, str) or not self.unsupported_reason.strip()):
            raise ValueError("Unsupported results require a nonempty reason")
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
        """Return JSON-compatible evidence with compact binary motif sites."""
        result = dict(
            peptide=self.peptide.to_dict(compact=True),
            model=self.model.to_dict(), sites=[],
            unsupported_reason=self.unsupported_reason,
            conditions=self.conditions,
            substrate_observation=self.substrate_observation)
        if self.model.evidence == "motif_rule":
            for site in self.sites:
                row = dict(bond=site.bond, match=site.match,
                           source_bond=self.peptide.source_start + site.bond)
                if site.reason != self.model.rule:
                    row["reason"] = site.reason
                result["sites"].append(row)
        else:
            result["sites"] = [
                dict(asdict(site), source_bond=self.peptide.source_start + site.bond)
                for site in self.sites]
        return {key: value for key, value in result.items()
                if value is not None and value != ()}

    @classmethod
    def from_dict(cls, value):
        """Restore serialized evidence, rejecting contradictory source bonds."""
        value = dict(value)
        peptide_value = dict(value.pop("peptide"))
        if "source_id" in peptide_value:
            if "source_sequence_name" in peptide_value:
                raise ValueError("Peptide has both source_id and source_sequence_name")
            peptide_value["source_sequence_name"] = peptide_value.pop("source_id")
        extra = set(peptide_value) - {field.name for field in fields(PeptideInput)}
        if extra:
            raise ValueError("Unknown peptide fields: %s" % ", ".join(sorted(extra)))
        peptide = PeptideInput.from_dict(peptide_value)
        model = CleavageModel.from_dict(value.pop("model"))
        sites = []
        for record in value.pop("sites", ()):
            record = dict(record)
            expected = peptide.source_start + record["bond"]
            if record.pop("source_bond", expected) != expected:
                raise ValueError("Serialized source_bond contradicts peptide offset")
            if "match" in record:
                if model.evidence != "motif_rule" or "status" in record or "score" in record:
                    raise ValueError("Binary match belongs only to motif-rule sites")
                match = record.pop("match")
                if type(match) is not int or match not in (0, 1):
                    raise ValueError("Motif match must be exactly 0 or 1")
                record["status"] = "matched" if match else "not_matched"
                record.setdefault("reason", model.rule)
                if record["reason"] is None:
                    raise ValueError("Motif match needs a model rule or site reason")
            sites.append(CleavageSite(**record))
        return cls(peptide, model, tuple(sites), **value)
