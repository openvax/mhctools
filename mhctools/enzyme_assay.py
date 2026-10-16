"""Bounded parent-loss fractions in specified purified-human-enzyme assays."""

from functools import lru_cache
import math

from ._resources import load_json_resource
from .cleavage import coerce_peptide
from .dpp4 import DPP4qPISA


def _first_order_fraction(half_life_minutes, minutes):
    exponent = -math.log(2) * minutes / half_life_minutes
    return {
        "parent_remaining_fraction": math.exp(exponent),
        "parent_lost_fraction": -math.expm1(exponent),
    }


@lru_cache(maxsize=1)
def _fap_calibration():
    """Geometric group means and substrate-held-out errors from Keane 2011."""
    data = load_json_resource("fap_keane2011.json")
    rows = data["rows"]
    groups = {aa: [row for row in rows if row["triplet"][1] == aa] for aa in "AP"}
    log_means = {
        aa: sum(math.log(row["half_life_minutes"]) for row in group) / len(group)
        for aa, group in groups.items()
    }
    errors = {}
    for aa, group in groups.items():
        errors[aa] = []
        for held_out in group:
            other = [row for row in group if row is not held_out]
            fitted = sum(math.log(row["half_life_minutes"]) for row in other) / len(other)
            errors[aa].append(abs(fitted - math.log(held_out["half_life_minutes"])))
    return data, log_means, errors


def estimate_assay_parent_fraction(peptide, *, enzyme, minutes):
    """Estimate parent loss in Keane 2011 purified-human-enzyme assays.

    DPP4 uses the independently published qPISA score and an eleven-peptide
    fit. FAP uses the geometric mean for N-terminal Ala or Pro among nine
    source substrates. Neither result estimates serum, IM-vaccine or APC
    exposure. Unsupported sequence context returns ``status='unassessed'``.

    Parameters
    ----------
    peptide : str or PeptideInput
        Canonical linear peptide, with explicit termini when modified.
    enzyme : str
        ``"DPP4"`` or ``"FAP"``; both use Keane's separate enzyme assays.
    minutes : float
        Finite, nonnegative exposure time in the source assay.

    Returns
    -------
    dict
        Source and model identity, estimated parent fractions and validation
        scope, or an unassessed status with a reason and no fraction.
    """
    if isinstance(minutes, bool) or not isinstance(minutes, (int, float)) or (
            not math.isfinite(minutes) or minutes < 0):
        raise ValueError("minutes must be finite and nonnegative")
    if enzyme == "DPP4":
        return DPP4qPISA().estimate_keane_fraction(peptide, minutes=minutes)
    if enzyme != "FAP":
        raise ValueError("Supported enzymes are DPP4 and FAP")
    peptide = coerce_peptide(peptide)
    data, log_means, errors = _fap_calibration()
    output = {
        "status": "unassessed", "enzyme": "FAP", "model": "fap-keane2011-parent-loss",
        "bond": 2, "minutes": minutes, "source": data["source"], "assay": data["assay"],
    }
    if peptide.n_term != "free" or peptide.c_term not in ("free", "amidated"):
        output["reason"] = "Requires an exposed free N terminus and free or amidated C terminus"
        return output
    if not 11 <= len(peptide.sequence) <= 42:
        output["reason"] = "Source includes only 11-42-residue peptides"
        return output
    triplet = peptide.sequence[:3]
    if triplet not in {row["triplet"] for row in data["rows"]}:
        output["reason"] = "N-terminal triplet is not represented in the FAP source cohort"
        return output
    p1 = triplet[1]
    half_life = math.exp(log_means[p1])
    group = [row for row in data["rows"] if row["triplet"][1] == p1]
    output.update(
        status="assay_estimate", n_terminal_triplet=triplet,
        estimated_half_life_minutes=half_life,
        validation={
            "substrates": len(group),
            "leave_one_substrate_out_typical_factor_error": math.exp(
                sum(errors[p1]) / len(errors[p1])),
        },
        **_first_order_fraction(half_life, minutes))
    return output
