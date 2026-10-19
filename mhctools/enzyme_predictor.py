"""One enzyme entry point, with distinct source-backed prediction methods."""

from dataclasses import dataclass
from typing import Optional, Tuple

from .cleavage import CleavageResult, coerce_peptide
from .enzyme_assay import (
    DEFAULT_EXPOSURE_MINUTES, AssayConditions, assay_fraction_method,
    estimate_assay_parent_fraction)
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

    def __repr__(self):
        """Keep interactive displays readable; details remain in the fields."""
        return ("EnzymePredictions(enzyme=%r, peptide=%r, methods=%d)" % (
            self.enzyme, self.peptide.sequence,
            len(self.cleavage_results) + len(self.assay_fractions)))

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
                   "minutes": result["minutes"],
                   "enzyme_concentration_um": result["enzyme_concentration_um"],
                   "reference_enzyme_concentration_um": result[
                       "reference_enzyme_concentration_um"]}
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
                detail = "%.1f%% parent lost, %.1f%% remaining after %g min at %g µM enzyme; estimated parent half-life %.1f min" % (
                    100 * row["parent_lost_fraction"],
                    100 * row["parent_remaining_fraction"], row["minutes"],
                    row["enzyme_concentration_um"], row["estimated_half_life_minutes"])
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
        """Return concise endpoint rows for display or routine JSON output.

        Model descriptions and source provenance are available through
        ``EnzymePredictor.describe()`` or ``to_evidence_dict()``.
        """
        rows = []
        for row in self.summary_rows():
            rows.append({key: value for key, value in row.items()
                         if key not in ("pattern", "score_name")})
        return {
            "schema_version": 2,
            "enzyme": self.enzyme, "peptide": self.peptide.to_dict(compact=True),
            "results": rows,
        }

    def to_evidence_dict(self):
        """Return the lossless, source-rich schema v1 evidence export."""
        return {
            "schema_version": 1,
            "enzyme": self.enzyme, "peptide": self.peptide.to_dict(compact=True),
            "cleavage_results": [result.to_dict() for result in self.cleavage_results],
            "assay_fractions": list(self.assay_fractions),
        }


class EnzymePredictor:
    """Discover and run methods for an enzyme class or cataloged name.

    Methods are reported independently: a motif is a recognition pattern,
    a scored model keeps its native units, and a time-dependent fraction is
    available only where a specified source assay has been calibrated.

    Use a concrete class such as ``DPP4()`` or ``FAP()`` in Python. The
    string constructor remains for existing callers; ``from_name()`` is the
    explicit registry path for configuration-driven selection.
    """

    enzyme_name = None

    def __init__(self, enzyme=None, *, include_optional=False):
        if self.enzyme_name is not None:
            if enzyme is not None:
                raise TypeError("%s already selects %s; pass no enzyme" % (
                    type(self).__name__, self.enzyme_name))
            enzyme = self.enzyme_name
        elif isinstance(enzyme, type) and issubclass(enzyme, EnzymePredictor):
            enzyme = enzyme.enzyme_name
        elif isinstance(enzyme, EnzymePredictor):
            enzyme = enzyme.enzyme
        if not isinstance(enzyme, str) or not enzyme.strip():
            raise ValueError("Use a concrete enzyme class or a nonempty cataloged name")
        models = cleavage_models(include_optional=include_optional)
        available = {model.enzyme for model in models}
        matches = [name for name in available if name.casefold() == enzyme.casefold()]
        if len(matches) != 1:
            raise ValueError("Unknown enzyme %r; choices: %s" % (
                enzyme, ", ".join(sorted(available))))
        self.enzyme = matches[0]
        self.models = tuple(model for model in models if model.enzyme == self.enzyme)

    @staticmethod
    def from_name(name, *, include_optional=False):
        """Resolve a string from configuration, using concrete classes when registered."""
        if not isinstance(name, str):
            raise TypeError("enzyme name must be a string")
        predictor_class = enzyme_predictors.get(name.casefold())
        if predictor_class is None:
            return EnzymePredictor(name, include_optional=include_optional)
        return predictor_class(include_optional=include_optional)

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

    def predict(self, peptide, *, conditions=None, minutes=None, enzyme_state=None):
        """Run each method; include a 60-minute source-assay scenario by default.

        ``conditions`` changes the assay exposure time and/or enzyme
        concentration. ``minutes`` remains a legacy exposure-time shortcut.
        Motifs and native cleavage scores do not depend on either parameter.
        """
        peptide = coerce_peptide(peptide)
        if minutes is not None and conditions is not None:
            raise ValueError("Specify either minutes or conditions, not both")
        if conditions is not None and not isinstance(conditions, AssayConditions):
            raise TypeError("conditions must be AssayConditions")
        if conditions is None:
            conditions = AssayConditions(
                exposure_minutes=DEFAULT_EXPOSURE_MINUTES if minutes is None else minutes)
        states = {self.enzyme: enzyme_state} if enzyme_state is not None else None
        cleavage = predict_cleavage(peptide, models=[m.name for m in self.models],
                                    enzyme_states=states)
        assays = ()
        if assay_fraction_method(self.enzyme):
            assays = (estimate_assay_parent_fraction(
                peptide, enzyme=self.enzyme, conditions=conditions),)
        return EnzymePredictions(self.enzyme, peptide, cleavage, assays,
                                 conditions.exposure_minutes)


class DPP4(EnzymePredictor):
    """Human DPP4: sourced motifs, qPISA score, and Keane assay parent loss."""

    enzyme_name = "DPP4"


class FAP(EnzymePredictor):
    """Human FAP: dipeptidyl and endopeptidyl motifs, plus Keane assay loss."""

    enzyme_name = "FAP"


enzyme_predictors = {"dpp4": DPP4, "fap": FAP}
"""Named enzyme classes for configuration-driven selection."""


class _EnzymeNamespace:
    """Attribute access to each default-cataloged enzyme as a predictor object."""

    @staticmethod
    def _names():
        return {"".join(word[:1].upper() + word[1:] for word in name.split()): name
                for name in EnzymePredictor.available()}

    def __getattr__(self, name):
        canonical = self._names().get(name)
        if canonical is None:
            raise AttributeError("No cataloged enzyme %r" % name)
        return EnzymePredictor.from_name(canonical)

    def __dir__(self):
        return sorted(self._names())


enzymes = _EnzymeNamespace()
"""Discoverable predictor objects, e.g. ``enzymes.ACE`` or ``enzymes.CathepsinS``."""
