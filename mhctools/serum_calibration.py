"""Reported assay measurements and explicitly bounded reference calculations."""

import math

from ._resources import load_json_resource
from .serum_degradation import _nonnegative


def serum_calibration_evidence():
    """Return the audited perturbation inventory and reported measurements.

    This offline inventory distinguishes matched human matrix experiments
    from purified enzymes, tissue exposure, animal data and substrate discovery.
    Reported half-lives do not identify an individual enzyme's cut hazards.
    """
    return load_json_resource("serum_calibration.json")


def serum_assay_parent_reference(measurement_id, times_hours):
    """Reproduce a reported parent-loss clock for one exact assay record.

    Parameters
    ----------
    measurement_id : str
        Measurement identity from ``serum_calibration_evidence``. No sequence
        matching or transfer to another substrate, matrix or temperature occurs.
    times_hours : iterable of float
        Finite nonnegative evaluation times, in hours.

    Returns
    -------
    dict
        Parent-signal survival under the paper's first-order assumption. Source
        ranges and censored half-lives produce bounds without a midpoint.
        Replicate standard deviation is preserved separately; these bounds are
        not confidence intervals. Neither target survival nor enzyme rates are
        identified by this reference calculation.
    """
    evidence = serum_calibration_evidence()
    record = next((r for r in evidence["measurements"]
                   if r["measurement_id"] == measurement_id), None)
    if record is None:
        raise ValueError("Unknown serum assay measurement: %s" % measurement_id)
    times = tuple(_nonnegative(t, "time_hours") for t in times_hours)
    if not times:
        raise ValueError("At least one evaluation time is required")
    lower = record["half_life_lower_hours"]
    upper = record["half_life_upper_hours"]
    point = record["half_life_point_hours"]

    def survival(time, half_life):
        if time == 0 or half_life is None:
            return 1.0
        if half_life == 0:
            return 0.0
        return math.exp(-math.log(2) * time / half_life)

    return {
        "measurement": record,
        "endpoint": "reported parent signal under first-order reference assumption",
        "source": evidence["studies"][record["study_id"]]["url"],
        "curves": [{
            "time_hours": time,
            "parent_remaining": survival(time, point) if point is not None else None,
            "parent_remaining_lower": survival(time, lower),
            "parent_remaining_upper": survival(time, upper),
            "beyond_study_maximum_horizon": time > record['study_maximum_horizon_hours'],
        } for time in times],
        "enzyme_cut_rates": None,
        "target_survival": None,
        "vaccine_transfer_validated": False,
        "uncertainty": record["uncertainty"],
    }
