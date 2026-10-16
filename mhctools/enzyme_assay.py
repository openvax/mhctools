"""Bounded parent-survival scenarios in specified purified-enzyme assays."""

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Optional

from ._resources import load_json_resource
from .cleavage import coerce_peptide
from .dpp4 import DPP4qPISA


_ASSAY_METHODS = {
    "DPP4": ("dpp4-qpisa-keane2011-parent-loss", "dpp4_keane2011.json"),
    "FAP": ("fap-keane2011-parent-loss", "fap_keane2011.json"),
}

DEFAULT_EXPOSURE_MINUTES = 60


@dataclass(frozen=True)
class AssayConditions:
    """Exposure scenario at the source assay's substrate loading and matrix.

    ``enzyme_concentration_um=None`` uses that enzyme's published assay
    concentration. Changing it assumes parent-loss rate scales linearly with
    active enzyme concentration; the source data do not calibrate changes in
    substrate loading, temperature, pH, or biological compartment.
    """

    exposure_minutes: float = DEFAULT_EXPOSURE_MINUTES
    enzyme_concentration_um: Optional[float] = None

    def __post_init__(self):
        if (isinstance(self.exposure_minutes, bool) or
                not isinstance(self.exposure_minutes, (int, float)) or
                not math.isfinite(self.exposure_minutes) or self.exposure_minutes < 0):
            raise ValueError("exposure_minutes must be finite and nonnegative")
        value = self.enzyme_concentration_um
        if value is not None and (isinstance(value, bool) or
                not isinstance(value, (int, float)) or
                not math.isfinite(value) or value <= 0):
            raise ValueError("enzyme_concentration_um must be finite and positive")


def assay_fraction_method(enzyme):
    """Catalog a supported Keane assay endpoint without estimating a peptide."""
    if enzyme not in _ASSAY_METHODS:
        return None
    name, resource = _ASSAY_METHODS[enzyme]
    data = load_json_resource(resource)
    return {
        "name": name, "evidence_type": "assay_estimate",
        "endpoint": "parent_peptide_loss_fraction",
        "description": "Estimated parent half-life and fraction remaining after a 60-minute reference exposure",
        "default_exposure_minutes": DEFAULT_EXPOSURE_MINUTES,
        "reference_enzyme_concentration_um": data["enzyme_concentration_um"],
        "assay": data["assay"], "references": [data["source"]],
        "limitations": data["inclusion"] + "; " + data["exclusions"],
    }


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


def _condition_estimate(result, *, enzyme, conditions):
    """Apply one shared exposure scenario to either source half-life estimate."""
    _, resource = _ASSAY_METHODS[enzyme]
    reference_concentration = load_json_resource(resource)["enzyme_concentration_um"]
    concentration = (reference_concentration if conditions.enzyme_concentration_um is None
                     else conditions.enzyme_concentration_um)
    result = dict(result)
    result.update(
        minutes=conditions.exposure_minutes,
        exposure_minutes=conditions.exposure_minutes,
        reference_enzyme_concentration_um=reference_concentration,
        enzyme_concentration_um=concentration,
        kinetic_assumption=(
            "First-order parent loss; concentration scenarios scale the loss rate "
            "linearly with active enzyme relative to the source assay."),
    )
    if result["status"] == "assay_estimate":
        reference_half_life = result["estimated_half_life_minutes"]
        half_life = reference_half_life * reference_concentration / concentration
        if not math.isfinite(half_life) or half_life <= 0:
            raise ValueError("enzyme_concentration_um yields an unrepresentable half-life")
        result.update(
            reference_half_life_minutes=reference_half_life,
            estimated_half_life_minutes=half_life,
            **_first_order_fraction(half_life, conditions.exposure_minutes),
        )
    return result


def estimate_assay_parent_fraction(peptide, *, enzyme, minutes=None, conditions=None):
    """Estimate parent survival in Keane 2011 purified-human-enzyme assays.

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
    minutes : float, optional
        Legacy exposure-time shortcut. Use ``conditions`` for new code.
    conditions : AssayConditions, optional
        Defaults to 60 minutes at the published enzyme concentration. Only
        time and active-enzyme concentration can vary; concentration scaling
        is an explicit sensitivity assumption, not a new assay calibration.

    Returns
    -------
    dict
        Source and model identity, estimated parent fractions and validation
        scope, or an unassessed status with a reason and no fraction.
    """
    if minutes is not None and conditions is not None:
        raise ValueError("Specify either minutes or conditions, not both")
    if conditions is not None and not isinstance(conditions, AssayConditions):
        raise TypeError("conditions must be AssayConditions")
    if conditions is None:
        conditions = AssayConditions(
            exposure_minutes=DEFAULT_EXPOSURE_MINUTES if minutes is None else minutes)
    if enzyme == "DPP4":
        result = DPP4qPISA().estimate_keane_fraction(
            peptide, minutes=conditions.exposure_minutes)
        return _condition_estimate(result, enzyme=enzyme, conditions=conditions)
    if enzyme != "FAP":
        raise ValueError("Supported enzymes are DPP4 and FAP")
    peptide = coerce_peptide(peptide)
    data, log_means, errors = _fap_calibration()
    output = {
        "status": "unassessed", "enzyme": "FAP", "model": "fap-keane2011-parent-loss",
        "bond": 2, "minutes": conditions.exposure_minutes,
        "source": data["source"], "assay": data["assay"],
    }
    if peptide.n_term != "free" or peptide.c_term not in ("free", "amidated"):
        output["reason"] = "Requires an exposed free N terminus and free or amidated C terminus"
        return _condition_estimate(output, enzyme=enzyme, conditions=conditions)
    if not 11 <= len(peptide.sequence) <= 42:
        output["reason"] = "Source includes only 11-42-residue peptides"
        return _condition_estimate(output, enzyme=enzyme, conditions=conditions)
    triplet = peptide.sequence[:3]
    if triplet not in {row["triplet"] for row in data["rows"]}:
        output["reason"] = "N-terminal triplet is not represented in the FAP source cohort"
        return _condition_estimate(output, enzyme=enzyme, conditions=conditions)
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
        **_first_order_fraction(half_life, conditions.exposure_minutes))
    return _condition_estimate(output, enzyme=enzyme, conditions=conditions)
