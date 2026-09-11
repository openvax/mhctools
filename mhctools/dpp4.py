"""Independent inference implementation of Gudipati et al. (2024) qPISA.

Uses the published rearranged human DPP4 coefficients in Dataset EV2;
no upstream R code or C. elegans DPF-3 coefficients are incorporated.
"""

from functools import lru_cache

from ._resources import load_json_resource
from .cleavage import (
    CleavageModel, CleavageResult, CleavageSite, coerce_peptide)


@lru_cache(maxsize=1)
def _parameters():
    return load_json_resource("dpp4_qpisa.json")["parameters"]


class DPP4qPISA:
    """Human DPP4 substrate score for the exposed N-terminal triplet.

    Higher scores indicate more substrate depletion in the source assay.
    They are neither cleavage probabilities nor half-lives. Missing
    published coefficients and unvalidated chemistry cause abstention.
    """

    model = CleavageModel(
        name="dpp4-qpisa", version="2024-EV2", enzyme="DPP4",
        uniprot="P27487", species="Homo sapiens",
        compartments=("extracellular", "plasma", "serum"),
        evidence="quantitative_model",
        references=("https://doi.org/10.1038/s44320-024-00071-4",
                    "https://www.uniprot.org/uniprotkb/P27487/entry"),
        assay="Purified human DPP4; tryptic HeLa peptides; HEPES pH 7.4, 21 C, 4 h",
        limitations=("In vitro substrate depletion model, not serum stability. "
                     "Unmodified linear peptides; structure and exposure are unmodeled. "
                     "A supported triplet is not necessarily experimentally observed. "
                     "No kinetics or successive-cleavage prediction."),
        score_name="predicted_log2_substrate_depletion",
        score_units="log2 fold change (buffer control / DPP4-treated)",
        scored_endpoint="substrate_depletion")

    def predict(self, peptide):
        """Assess only bond 2 of a :class:`CleavageInput` or canonical string."""
        peptide = coerce_peptide(peptide)
        reason = None
        if peptide.n_term != "free":
            reason = "DPP4 requires an exposed free N terminus"
        elif peptide.c_term != "free":
            reason = "Modified or unknown C-terminal chemistry is outside this model"
        elif len(peptide.sequence) < 3:
            reason = "The model requires P2, P1 and P1-prime (at least 3 residues)"
        if reason:
            return CleavageResult(peptide, self.model, unsupported_reason=reason)
        p2, p1, prime = peptide.sequence[:3]
        terms = _parameters()[p1]
        names = ("P1", "P2:" + p2, "P1':" + prime)
        missing = [name for name in names if terms[name] is None]
        if missing:
            return CleavageResult(
                peptide, self.model, unsupported_reason=(
                    "Dataset EV2 has missing coefficients for P1=%s: %s" %
                    (p1, ", ".join(missing))))
        score = sum(terms[name] for name in names)
        return CleavageResult(peptide, self.model, (CleavageSite(
            bond=2, status="scored", score=score,
            reason="Published P1 + P2:P1 + P1:P1-prime terms for " +
                   peptide.sequence[:3]),))
