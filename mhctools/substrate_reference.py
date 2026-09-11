"""Exact human substrate observations, explicitly separate from predictions."""

from functools import lru_cache
from importlib.resources import files
import json

from .cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite


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
        for case in cases:
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
        if isinstance(peptide, str):
            peptide = CleavageInput(peptide)
        if not isinstance(peptide, CleavageInput):
            raise TypeError("Expected CleavageInput or canonical peptide string")
        result = self._cases.get((peptide.sequence, peptide.n_term, peptide.c_term))
        if result is None:
            return CleavageResult(peptide, self.model, unsupported_reason=
                "No exact sequence/chemical-form observation in this source reference; no extrapolation")
        return CleavageResult(peptide, self.model, result.sites, conditions=result.conditions,
                              substrate_observation=result.substrate_observation)


@lru_cache(maxsize=None)
def substrate_references():
    """Load the small packaged factual catalog; no external weights or runtime."""
    data = json.loads(files("mhctools").joinpath("data/intracellular_substrate_evidence.json").read_text())
    return tuple(PeptidaseSubstrateReference(item["model"], item["cases"]) for item in data["models"])
