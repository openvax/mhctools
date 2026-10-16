"""One enzyme entry point, with distinct source-backed prediction methods."""

from dataclasses import dataclass
import math
from typing import Optional, Tuple

from .cleavage import CleavageResult, coerce_peptide
from .enzyme_assay import assay_fraction_method, estimate_assay_parent_fraction
from .peptide_input import PeptideInput
from .peptidases import PeptidaseMotif, cleavage_models, get_cleavage_model, predict_cleavage


@dataclass(frozen=True)
class EnzymePredictions:
    """Results for one peptide and enzyme; methods keep their native meaning."""

    enzyme: str
    peptide: PeptideInput
    cleavage_results: Tuple[CleavageResult, ...]
    assay_fractions: Tuple[dict, ...] = ()
    requested_minutes: Optional[float] = None

    @property
    def by_model(self):
        """Look up an existing result by its source-specific model name."""
        return {**{result.model.name: result for result in self.cleavage_results},
                **{result["model"]: result for result in self.assay_fractions}}

    def summary_rows(self):
        """One human-readable endpoint row per method, without combining scores."""
        rows = []
        for result in self.cleavage_results:
            model = result.model
            row = {"model": model.name, "evidence_type": model.evidence,
                   "endpoint": (model.scored_endpoint if model.evidence == "quantitative_model"
                                else "recognition_pattern" if model.evidence == "motif_rule"
                                else "substrate_observation")}
            if result.unsupported_reason:
                row.update(status="unassessed", reason=result.unsupported_reason)
            elif model.evidence == "motif_rule":
                bonds = [site.bond for site in result.sites if site.match == 1]
                row.update(status="matched" if bonds else "not_matched",
                           pattern=model.rule, matched_bonds=bonds,
                           assessed_bonds=len(result.sites))
            elif model.evidence == "quantitative_model":
                row.update(status="scored" if result.sites else "unassessed",
                           score_name=model.score_name, score_units=model.score_units,
                           scores=[{"bond": site.bond, "score": site.score}
                                   for site in result.sites])
            else:
                row.update(status=result.substrate_observation or "unassessed",
                           reported_bonds=[site.bond for site in result.sites])
            rows.append(row)
        for result in self.assay_fractions:
            row = {"model": result["model"], "evidence_type": "assay_estimate",
                   "endpoint": "parent_peptide_loss_fraction", "status": result["status"],
                   "minutes": result["minutes"]}
            if result["status"] == "assay_estimate":
                row.update(parent_lost_fraction=result["parent_lost_fraction"],
                           parent_remaining_fraction=result["parent_remaining_fraction"],
                           estimated_half_life_minutes=result["estimated_half_life_minutes"])
            else:
                row["reason"] = result["reason"]
            rows.append(row)
        return tuple(rows)

    def summary(self):
        """Read the separate outputs in one compact, deliberately unranked report."""
        lines = ["%s on %s" % (self.enzyme, self.peptide.sequence)]
        for row in self.summary_rows():
            status = row["status"]
            if status == "matched":
                detail = "%s matched at bond(s) %s" % (
                    row["pattern"], ", ".join(str(bond) for bond in row["matched_bonds"]))
            elif status == "not_matched":
                detail = "%s did not match any assessed bond" % row["pattern"]
            elif status == "scored":
                detail = ", ".join("bond %d: %.4g %s" % (
                    item["bond"], item["score"], row["score_units"])
                    for item in row["scores"])
            elif status == "assay_estimate":
                detail = "%.1f%% parent lost, %.1f%% remaining after %g min in source assay" % (
                    100 * row["parent_lost_fraction"],
                    100 * row["parent_remaining_fraction"], row["minutes"])
            elif status == "unassessed":
                detail = "unassessed: " + row.get("reason", "no applicable site")
            elif status == "cleavage_reported":
                bonds = row["reported_bonds"]
                detail = ("source reported cleavage at bond(s) %s" % ", ".join(
                    str(bond) for bond in bonds) if bonds
                    else "source reported cleavage without a localized bond")
            else:
                detail = status.replace("_", " ")
            lines.append("- %s [%s]: %s" % (row["model"], row["evidence_type"], detail))
        if self.requested_minutes is not None and not self.assay_fractions:
            lines.append("- No assay-calibrated parent fraction is available for this enzyme.")
        return "\n".join(lines)

    def to_dict(self):
        """Serialize each existing result under its own evidence type."""
        return {
            "schema_version": 1,
            "enzyme": self.enzyme, "peptide": self.peptide.to_dict(compact=True),
            "cleavage_results": [result.to_dict() for result in self.cleavage_results],
            "assay_fractions": list(self.assay_fractions),
        }


class EnzymePredictor:
    """Discover and run available methods for one exact human enzyme.

    Methods are reported independently: a motif is a recognition pattern,
    a scored model keeps its native units, and a time-dependent fraction is
    available only where a specified source assay has been calibrated.
    """

    def __init__(self, enzyme, *, include_optional=False):
        if not isinstance(enzyme, str) or not enzyme.strip():
            raise ValueError("enzyme must be a nonempty name")
        models = cleavage_models(include_optional=include_optional)
        available = {model.enzyme for model in models}
        matches = [name for name in available if name.casefold() == enzyme.casefold()]
        if len(matches) != 1:
            raise ValueError("Unknown enzyme %r; choices: %s" % (
                enzyme, ", ".join(sorted(available))))
        self.enzyme = matches[0]
        self.models = tuple(model for model in models if model.enzyme == self.enzyme)

    @staticmethod
    def available(*, include_optional=False):
        """List enzymes with at least one curated cleavage method."""
        return tuple(sorted({model.enzyme for model in cleavage_models(
            include_optional=include_optional)}))

    def describe(self):
        """Show what each named method assesses, with its assay and sources."""
        methods = []
        for model in self.models:
            row = {
                "name": model.name, "evidence_type": model.evidence,
                "endpoint": (model.scored_endpoint if model.evidence == "quantitative_model"
                             else "recognition_pattern" if model.evidence == "motif_rule"
                             else "substrate_observation"),
                "description": model.rule or model.score_name or "Exact source substrate observation",
                "compartments": list(model.compartments),
                "assay": model.assay, "references": list(model.references),
                "limitations": model.limitations,
            }
            if model.evidence == "motif_rule":
                row["strictness"] = model.motif_strictness
                motif = get_cleavage_model(model.name)
                if isinstance(motif, PeptidaseMotif):
                    if motif.topology == "internal":
                        row["activity"] = "Cuts an internal peptide bond"
                    else:
                        terminus = "N" if motif.topology == "n_terminal" else "C"
                        row["activity"] = "Removes %d %s-terminal residue(s)" % (
                            motif.removed, terminus)
            if model.evidence == "quantitative_model":
                row["score_units"] = model.score_units
            methods.append(row)
        assay = assay_fraction_method(self.enzyme)
        if assay:
            methods.append(assay)
        return {"enzyme": self.enzyme, "methods": methods}

    def predict(self, peptide, *, minutes=None, enzyme_state=None):
        """Run every selected method; ``minutes`` requests assay fractions."""
        peptide = coerce_peptide(peptide)
        if minutes is not None and (
                isinstance(minutes, bool) or not isinstance(minutes, (int, float)) or
                not math.isfinite(minutes) or minutes < 0):
            raise ValueError("minutes must be finite and nonnegative")
        states = {self.enzyme: enzyme_state} if enzyme_state is not None else None
        cleavage = predict_cleavage(peptide, models=[m.name for m in self.models],
                                    enzyme_states=states)
        assays = ()
        if minutes is not None and assay_fraction_method(self.enzyme):
            assays = (estimate_assay_parent_fraction(
                peptide, enzyme=self.enzyme, minutes=minutes),)
        return EnzymePredictions(self.enzyme, peptide, cleavage, assays, minutes)
