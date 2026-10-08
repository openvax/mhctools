"""Assay-scoped peptidase kinetics and conditional enzyme-removal effects.

Published reference substrates are not universal vaccine rate calibrators.
The plasma library supplies a boundary spectrum, not enzyme turnover rates.
"""

from dataclasses import dataclass
import json
import math
from numbers import Integral
from pathlib import Path

from .serum_degradation import (
    _nonnegative, _simulate_target_paths, degradation_curve,
)


@dataclass(frozen=True)
class EnzymeCutRate:
    """One current-fragment bond's absolute hazard, in inverse hours.

    ``enzyme`` can name an unresolved activity group. ``basis`` must explain
    whether the rate is measured, fitted or assumed; ``source`` identifies its
    evidence or supplied scenario. A native recognition score is not a rate.
    Multiple enzymes may compete for the same bond as separate channels.
    """

    bond: int
    enzyme: str
    rate_per_hour: float
    basis: str
    source: str

    def __post_init__(self):
        if isinstance(self.bond, bool) or not isinstance(self.bond, Integral) or self.bond < 1:
            raise ValueError("bond must be a positive integer")
        object.__setattr__(self, "rate_per_hour", _nonnegative(self.rate_per_hour, "cut rate"))
        for name in ("enzyme", "basis", "source"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(name + " must be nonempty")


def serum_contribution_evidence():
    """Load the source-linked factual inventory; no network or model download."""
    return json.loads((Path(__file__).parent / "data" / "serum_contributions.json").read_text())


def serum_reference_kinetics(substrate_id, *, concentration_um, serum_fraction):
    """Calculate local hazards from measured apparent Km and serum Vmax.

    Parameters
    ----------
    substrate_id : str
        Exact chemical assay substrate ID from the evidence inventory. This
        intentionally does not match arbitrary sequences by motif or family.
    concentration_um : float
        Current substrate concentration in micromolar; positive. At finite
        concentration these are LOCAL hazards, not constant decay coefficients.
    serum_fraction : float
        Serum volume / reaction volume, in (0, 1]. Linear scaling of activity
        with dilution is an explicit assumption, not a validated matrix transfer.

    Returns
    -------
    dict
        For each measured enzyme, k = 60 * serum_fraction * Vmax / (Km + S).
        pmol / microliter equals micromolar, so k has units h^-1. Shares are
        conditional on these reported channels only and are never pooled across
        different substrates. The result is not a whole-serum survival estimate.
    """
    concentration = _nonnegative(concentration_um, "concentration_um")
    fraction = _nonnegative(serum_fraction, "serum_fraction")
    if concentration == 0 or not 0 < fraction <= 1:
        raise ValueError("Positive concentration and serum_fraction in (0, 1] required")
    evidence = serum_contribution_evidence()
    substrates = evidence["kinetic_reference"]["substrates"]
    if substrate_id not in substrates:
        raise ValueError("No measured reference kinetics for this exact assay substrate")
    substrate = substrates[substrate_id]
    rows = []
    for observation in substrate["enzymes"]:
        km, vmax = observation["km_um"], observation["vmax_pmol_per_min_per_ul_serum"]
        rate = 60 * fraction * vmax / (km + concentration)
        rows.append(dict(**observation, local_rate_per_hour=rate))
    total = math.fsum(row["local_rate_per_hour"] for row in rows)
    for row in rows:
        row["share_of_reported_channels"] = row["local_rate_per_hour"] / total
    lower, upper = evidence["kinetic_reference"]["measured_concentration_range_um"]
    return dict(
        substrate_id=substrate_id, substrate=substrate,
        context=evidence["kinetic_reference"]["context"],
        source=evidence["kinetic_reference"]["source"],
        concentration_um=concentration, serum_fraction=fraction,
        extrapolated_concentration=not lower <= concentration <= upper,
        extrapolated_serum_fraction=not math.isclose(
            fraction, evidence["kinetic_reference"]["context"]["assay_serum_fraction"]),
        total_reported_local_rate_per_hour=total, enzymes=rows,
        interpretation="Local assay-reference hazards; incomplete pathways; no vaccine or systemic PK calibration")


def empirical_terminal_cut_rates(sequence, *, total_rate_per_hour, allow_transfer=False):
    """Construct an explicitly terminal-only scenario from library boundaries.

    The four immediate terminal boundary counts (46, 30, 7, 7) supply relative
    scenario weights. They are unique observed sequence contexts, not event
    frequencies. More distant boundaries may result from successive cuts and
    are NOT treated as direct jumps. Individual enzyme names remain unresolved.

    An absolute total hazard must be supplied separately, for example from ONE
    half-life estimator. This is an exploratory closed terminal-trimming model,
    not an empirical rate fit. No internal-cut hazard is modeled. Transfer from
    14-mers in pig plasma to another length requires explicit opt-in. Fewer than
    five residues are unassessed because the four terminal routes overlap.
    """
    rate = _nonnegative(total_rate_per_hour, "total_rate_per_hour")
    if not isinstance(allow_transfer, bool):
        raise ValueError("allow_transfer must be boolean")
    if len(sequence) < 5 or (len(sequence) != 14 and not allow_transfer):
        return None
    study = serum_contribution_evidence()["plasma_library"]
    counts = study["immediate_terminal_context_counts"]
    bonds = dict(n_mono=1, n_di=2, c_mono=len(sequence) - 1, c_di=len(sequence) - 2)
    total = sum(counts.values())
    return tuple(EnzymeCutRate(
        bonds[group], group, rate * count / total,
        "Assumed terminal-only hazard allocation using pig-plasma unique boundary counts; enzyme unresolved",
        study["source"]) for group, count in counts.items())


def simulate_enzyme_degradation(
        peptide, target, fragment_rates, *, scenario, estimator="enzyme-rate scenario",
        disabled_enzymes=(), n_paths=1000, horizon_hours=24.0, seed=0,
        clearance_rate_per_hour=0.0, uptake_rate_per_hour=0.0):
    """Sample successive target-bearing products with enzyme-labelled hazards.

    Parameters
    ----------
    fragment_rates : callable
        Receives each CURRENT sequence and returns EnzymeCutRate records.
        Return None when fragment kinetics are missing. An empty tuple means
        explicitly assessed zero cleavage hazard, not unavailable inference.
        Do not supply partial channels as a complete calibrated serum model.
    disabled_enzymes : iterable of str
        Remove these channels without renormalizing any other enzyme's rate.
        Biological compensation and inhibitor off-target effects are unmodeled.

    Other parameters and chemistry follow simulate_target_degradation. The
    enzyme-rate sum supplies each fragment's next-cut clock; aggregate half-life
    estimators must not be added again as an independent degradation hazard.
    Rate basis/source are retained in each sampled cut; mechanism names can be
    unresolved groups. Exact target retention does not establish presentation.
    """
    if isinstance(disabled_enzymes, str):
        raise ValueError("disabled_enzymes must be an iterable of enzyme names, not a string")
    disabled = frozenset(disabled_enzymes)
    if any(not isinstance(name, str) or not name.strip() for name in disabled):
        raise ValueError("Disabled enzyme names must be nonempty strings")

    def kinetics(sequence):
        values = fragment_rates(sequence)
        if values is None:
            return None, None
        channels = []
        for value in values:
            if not isinstance(value, EnzymeCutRate):
                raise TypeError("fragment_rates must return EnzymeCutRate records or None")
            if value.bond >= len(sequence):
                raise ValueError("Cut bond must be internal to the current fragment")
            if value.enzyme not in disabled and value.rate_per_hour > 0:
                channels.append((value.bond, value.rate_per_hour, value.enzyme,
                                 value.basis + "; source: " + value.source))
        total = math.fsum(channel[1] for channel in channels)
        if not math.isfinite(total):
            raise ValueError("Total cut rate must be finite")
        return (math.log(2) / total if total else None), tuple(channels)

    return _simulate_target_paths(
        peptide, target, kinetics, estimator=estimator, scenario=scenario,
        n_paths=n_paths, horizon_hours=horizon_hours, seed=seed,
        clearance_rate_per_hour=clearance_rate_per_hour,
        uptake_rate_per_hour=uptake_rate_per_hour)


def enzyme_removal_effects(
        peptide, target, fragment_rates, *, enzymes, times_hours, scenario,
        estimator="enzyme-rate scenario", n_paths=10000, seed=0,
        clearance_rate_per_hour=0.0, uptake_rate_per_hour=0.0):
    """Compare baseline with separate single-enzyme-removal counterfactuals.

    Returns target availability differences in percentage points and parent
    survival alongside terminal destructive-cut attribution. Differences can
    be negative when removing a protective trimming route. Sensitivities need
    not add to 100% and are NOT an enzyme's fraction of serum degradation.

    Unknown fragment kinetics give identification bounds, not a point estimate
    or a statistical confidence interval. Monte Carlo sampling uncertainty is
    separate: 100/sqrt(n_paths) percentage points is a worst-case upper bound
    on the difference's sampling SE under any coupling of two binomial arms.
    It does not collapse to zero when a finite sample observes no events. Each arm
    is simulated with a reproducible common seed. Unmodeled
    biological pathways remain outside the supplied scenario.
    """
    times = tuple(_nonnegative(t, "curve time") for t in times_hours)
    if not times:
        raise ValueError("At least one curve time is required")
    if isinstance(enzymes, str):
        raise ValueError("enzymes must be an iterable of names, not a string")
    enzymes = tuple(enzymes)
    if (len(set(enzymes)) != len(enzymes) or
            any(not isinstance(e, str) or not e.strip() for e in enzymes)):
        raise ValueError("Enzyme names must be unique and nonempty")
    kwargs = dict(estimator=estimator, scenario=scenario, n_paths=n_paths,
                  horizon_hours=max(times), seed=seed,
                  clearance_rate_per_hour=clearance_rate_per_hour,
                  uptake_rate_per_hour=uptake_rate_per_hour)
    baseline_paths = simulate_enzyme_degradation(peptide, target, fragment_rates, **kwargs)
    baseline = degradation_curve(baseline_paths, times)
    parent_unknown = sum(p.steps[0].event == "unknown" for p in baseline_paths) / n_paths
    effects = []
    for enzyme in enzymes:
        paths = simulate_enzyme_degradation(
            peptide, target, fragment_rates, disabled_enzymes=(enzyme,), **kwargs)
        removed_parent_unknown = sum(p.steps[0].event == "unknown" for p in paths) / n_paths
        for base, removed in zip(baseline, degradation_curve(paths, times)):
            difference = removed["target_in_circulation"] - base["target_in_circulation"]
            parent_difference = removed["parent_remaining"] - base["parent_remaining"]
            unknown = base["unknown"] + removed["unknown"]
            effects.append(dict(
                enzyme=enzyme, time_hours=base["time_hours"],
                target_gain_percentage_points=None if unknown else 100 * difference,
                target_gain_mc_se_upper_bound_pp=None if unknown else 100 / math.sqrt(n_paths),
                target_gain_bounds_percentage_points=(
                    100 * (difference - base["unknown"]),
                    100 * (difference + removed["unknown"])),
                parent_gain_percentage_points=(None if parent_unknown or removed_parent_unknown
                                               else 100 * parent_difference),
                parent_remaining_without_enzyme=removed["parent_remaining"],
                target_remaining_without_enzyme=removed["target_in_circulation"],
                unknown_without_enzyme=removed["unknown"]))
    causes = []
    for time in times:
        counts, first_counts = {}, {}
        for path in baseline_paths:
            first = path.steps[0]
            if first.event in ("cut", "target_destroyed") and first.event_hours <= time:
                first_counts[first.mechanism] = first_counts.get(first.mechanism, 0) + 1
            last = path.steps[-1]
            if last.event == "target_destroyed" and last.event_hours <= time:
                counts[last.mechanism] = counts.get(last.mechanism, 0) + 1
        causes.append(dict(time_hours=time,
                           first_cut_fraction_by_enzyme={e: n / n_paths for e, n in sorted(first_counts.items())},
                           destroyed_fraction_by_enzyme={e: n / n_paths for e, n in sorted(counts.items())}))
    return dict(
        estimator=estimator, scenario=scenario, n_paths=n_paths, baseline=baseline,
        enzyme_removal_effects=effects, destructive_cut_attribution=causes,
        interpretation="Conditional on supplied fragment rates; effects are nonadditive; no universal enzyme weights")
