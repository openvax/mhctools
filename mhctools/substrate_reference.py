"""Exact human substrate observations, explicitly separate from predictions."""

from functools import lru_cache

from ._resources import load_json_resource
from .cleavage import (
    CleavageInput, CleavageModel, CleavageResult, CleavageSite, coerce_peptide)


#: Fields every entry in a reference catalog's ``cases`` list must carry.
_REQUIRED_CASE_FIELDS = ("sequence", "n_term", "c_term", "conditions", "source",
                         "source_measurement_id", "bonds", "interpretation",
                         "substrate_observation")


class PeptidaseSubstrateReference:
    """Look up an exact chemical form; unseen sequences are never extrapolated.

    Results describe the source experiment's conditions, not an assumption
    that the same outcome occurs in the user's biological setting.
    """

    def __init__(self, metadata, cases):
        self.model = CleavageModel(**metadata)
        if self.model.evidence != "substrate_reference":
            raise ValueError("Reference catalog requires substrate_reference evidence")
        reserved = {"source", "source_measurement_id"}
        self._cases = {}
        for index, case in enumerate(cases):
            missing = [field for field in _REQUIRED_CASE_FIELDS if field not in case]
            if missing:
                raise ValueError("Reference catalog case %d for %r is missing required "
                                 "fields: %s" % (index, self.model.name, ", ".join(missing)))
            peptide = CleavageInput(case["sequence"], case["n_term"], case["c_term"])
            key = (peptide.sequence, peptide.n_term, peptide.c_term)
            if key in self._cases:
                raise ValueError("Conflicting or repeated chemical form in reference catalog")
            if reserved & set(case["conditions"]):
                raise ValueError(
                    "Case conditions must not use the reserved keys 'source' or 'source_measurement_id'")
            conditions = tuple(sorted(case["conditions"].items())) + (
                ("source", case["source"]), ("source_measurement_id", case["source_measurement_id"]))
            sites = tuple(CleavageSite(b, "reported", case["interpretation"]) for b in case["bonds"])
            result = CleavageResult(peptide, self.model, sites, conditions=conditions,
                                    substrate_observation=case["substrate_observation"])
            self._cases[key] = result

    def predict(self, peptide):
        """Return source observations for an exact match, preserving source offsets."""
        peptide = coerce_peptide(peptide)
        result = self._cases.get((peptide.sequence, peptide.n_term, peptide.c_term))
        if result is None:
            return CleavageResult(peptide, self.model, unsupported_reason=
                "No exact sequence/chemical-form observation in this source reference; no extrapolation")
        return CleavageResult(peptide, self.model, result.sites, conditions=result.conditions,
                              substrate_observation=result.substrate_observation)


@lru_cache(maxsize=None)
def substrate_references():
    """Load the small packaged factual catalog; no external weights or runtime."""
    data = load_json_resource("intracellular_substrate_evidence.json")
    models = data["models"]
    for index, item in enumerate(models):
        missing = [field for field in ("model", "cases") if field not in item]
        if missing:
            raise ValueError(
                "Reference catalog entry %d is missing required fields: %s" % (
                    index, ", ".join(missing)))
    return tuple(PeptidaseSubstrateReference(item["model"], item["cases"]) for item in models)
