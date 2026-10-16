"""Independent inference implementation of Gudipati et al. (2024) qPISA.

Uses the published rearranged human DPP4 coefficients in Dataset EV2;
no upstream R code or C. elegans DPF-3 coefficients are incorporated.
"""

from functools import lru_cache
import math

from ._resources import load_json_resource
from .cleavage import (
    CleavageModel, CleavageResult, CleavageSite, coerce_peptide)


@lru_cache(maxsize=1)
def _parameters():
    return load_json_resource("dpp4_qpisa.json")["parameters"]


def _line(xy):
    """Least-squares line for (score, log half-life) observations."""
    mean_x = sum(x for x, _ in xy) / len(xy)
    mean_y = sum(y for _, y in xy) / len(xy)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in xy) / sum(
        (x - mean_x) ** 2 for x, _ in xy)
    return mean_y - slope * mean_x, slope


@lru_cache(maxsize=1)
def _keane_calibration():
    """Fit and leave-one-triplet-out audit of the Keane 2011 DPP4 assay."""
    data = load_json_resource("dpp4_keane2011.json")
    observations = []
    for row in data["rows"]:
        result = DPP4qPISA().predict(row["triplet"])
        if not result.sites:
            raise ValueError("A calibration triplet lacks qPISA coefficients")
        observations.append((result.sites[0].score,
                             math.log(row["half_life_minutes"]), row["triplet"]))
    intercept, slope = _line([(x, y) for x, y, _ in observations])
    errors = []
    for score, observed, triplet in observations:
        held_out = [(x, y) for x, y, group in observations if group != triplet]
        held_intercept, held_slope = _line(held_out)
        errors.append(abs(held_intercept + held_slope * score - observed))
    scores = [score for score, _, _ in observations]
    return {
        "source": data["source"], "assay": data["assay"],
        "intercept": intercept, "slope": slope,
        "min_score": min(scores), "max_score": max(scores),
        "substrates": len(observations),
        "triplets": len({triplet for _, _, triplet in observations}),
        "held_out_log_mae": sum(errors) / len(errors),
    }


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
        """Assess only bond 2 of a PeptideInput or canonical string."""
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

    def estimate_keane_fraction(self, peptide, *, minutes):
        """Estimate parent loss in Keane's purified-DPP4 assay only.

        Fits log half-life to the qPISA triplet score for 11 short peptides.
        Conversion from fitted half-life to a time-specific fraction assumes
        first-order parent loss. This is not a serum or vaccine kinetic rate.

        Parameters
        ----------
        peptide : str or PeptideInput
            Canonical peptide with exposed, free termini.
        minutes : float
            Nonnegative assay exposure time.

        Returns
        -------
        dict
            Parent fractions and validation scope, or an unassessed reason.
        """
        if isinstance(minutes, bool) or not isinstance(minutes, (int, float)) or (
                not math.isfinite(minutes) or minutes < 0):
            raise ValueError("minutes must be finite and nonnegative")
        peptide = coerce_peptide(peptide)
        calibration = _keane_calibration()
        output = {
            "status": "unassessed", "bond": 2, "minutes": minutes,
            "source": calibration["source"], "assay": calibration["assay"],
        }
        if not 27 <= len(peptide.sequence) <= 42:
            output["reason"] = "Calibration includes only 27-42-residue peptides"
            return output
        result = self.predict(peptide)
        if not result.sites:
            output["reason"] = result.unsupported_reason
            return output
        score = result.sites[0].score
        if not calibration["min_score"] <= score <= calibration["max_score"]:
            output["reason"] = "qPISA score is outside the calibrated score range"
            return output
        half_life = math.exp(calibration["intercept"] + calibration["slope"] * score)
        remaining = math.exp(-math.log(2) * minutes / half_life)
        output.update(
            status="assay_estimate", qpisa_score=score,
            estimated_half_life_minutes=half_life,
            parent_remaining_fraction=remaining,
            parent_lost_fraction=-math.expm1(-math.log(2) * minutes / half_life),
            validation={
                "substrates": calibration["substrates"],
                "unique_triplets": calibration["triplets"],
                "leave_one_triplet_out_typical_factor_error": math.exp(
                    calibration["held_out_log_mae"]),
            })
        return output
