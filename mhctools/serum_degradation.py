"""Exploratory target-preserving fragment paths, not calibrated circulation PK.

Half-life estimates set a *scenario's* exponential waiting time. Explicit
relative cut weights allocate that hazard; they are not native enzyme rates.
The target-containing product is rescored after each flank or boundary cut.
Each estimator must be simulated separately. Exact target retention does not
establish MHC presentation, cell uptake, or an immune response.
"""

from dataclasses import dataclass
import math
from numbers import Integral, Real
import random
from typing import Optional

from .peptide_input import PeptideInput, sequence_only_chemistry_error


@dataclass(frozen=True)
class DegradationTarget:
    """An exact target occurrence; coordinates are zero-based, half-open."""

    label: str
    start: int
    end: int

    def validate(self, sequence):
        if (not self.label or any(isinstance(x, bool) or not isinstance(x, Integral)
                                  for x in (self.start, self.end)) or
                not 0 <= self.start < self.end <= len(sequence)):
            raise ValueError("Target requires a label and valid [start, end) coordinates")


@dataclass(frozen=True)
class TargetFragment:
    """A target-containing fragment in original-parent coordinates."""

    sequence: str
    start: int
    end: int


@dataclass(frozen=True)
class DegradationStep:
    """One fragment residence and its terminating event.

    ``cut_bond`` is in original-parent coordinates: bond b splits parent[:b]
    from parent[b:]. A target-boundary cut retains the entire exact target.
    """

    start: int
    end: int
    entered_hours: float
    event_hours: float
    half_life_hours: Optional[float]
    event: str
    cut_bond: Optional[int] = None
    reason: str = ""


@dataclass(frozen=True)
class DegradationPath:
    """A sampled lineage of one exact target, with explicit scenario labels."""

    estimator: str
    scenario: str
    peptide_input: PeptideInput
    target: DegradationTarget
    parent_length: int
    horizon_hours: float
    steps: tuple[DegradationStep, ...]


def _input(peptide):
    item = peptide if isinstance(peptide, PeptideInput) else PeptideInput(peptide)
    error = sequence_only_chemistry_error(item)
    if error:
        raise ValueError("Fragment simulation requires canonical L-peptides with free termini: " + error)
    return item


def target_fragments(peptide, target):
    """Enumerate all contiguous products retaining the exact target occurrence.

    This permits batch half-life inference before sampling, including newly
    exposed terminal contexts. Duplicate sequences at different coordinates
    retain separate fragment records; inference can be cached by sequence.
    """
    sequence = _input(peptide).sequence
    target.validate(sequence)
    return tuple(TargetFragment(sequence[start:end], start, end)
                 for start in range(target.start + 1)
                 for end in range(target.end, len(sequence) + 1))


def _nonnegative(value, name):
    if (isinstance(value, bool) or not isinstance(value, Real) or
            not math.isfinite(value) or value < 0):
        raise ValueError(name + " must be finite and nonnegative")
    return float(value)


def simulate_target_degradation(
        peptide, target, half_lives, cut_weights, *, estimator, scenario,
        n_paths=1000, horizon_hours=24.0, seed=0,
        clearance_rate_per_hour=0.0, uptake_rate_per_hour=0.0):
    """Sample target lineages under explicitly supplied kinetic assumptions.

    Parameters
    ----------
    peptide : str or PeptideInput
        Canonical L-peptide with free termini. Chemistry is never stripped.
    target : DegradationTarget
        Exact coordinate-defined target; an internal cut destroys it.
    half_lives : mapping
        Fragment sequence to positive half-life in hours for ONE estimator.
        Missing, None, or zero estimates yield an explicit unknown outcome.
        Invalid negative/nonfinite values raise. No estimator averaging occurs.
    cut_weights : callable
        Receives the CURRENT fragment sequence and returns one finite,
        nonnegative relative weight per internal bond, in order. Total weight
        must be positive. These are scenario assumptions, not serum rates or
        enzyme score probabilities. Include a background for unmodeled cuts.
    estimator, scenario : str
        Required labels carried by every path.
    n_paths : int
        Number of independent target lineages.
    horizon_hours : float
        Finite observation horizon; surviving paths are right-censored.
    seed : int
        Local random seed; no global RNG state is changed.
    clearance_rate_per_hour, uptake_rate_per_hour : float
        Separately supplied constant competing hazards. Defaults of zero are
        a serum-digestion-only scenario, not estimates of human clearance or
        antigen-presenting-cell uptake. Uptake does not establish presentation.

    Returns
    -------
    tuple of DegradationPath
        Waiting time is exponential with total rate ln(2)/fragment_half_life
        plus the supplied clearance and uptake rates. Cleavage location is
        sampled from relative cut weights. Only the target-containing product
        is followed; products lacking this exact target need not be simulated.
        This marginal lineage cannot give joint survival of multiple targets.
    """
    peptide_input = _input(peptide)
    sequence = peptide_input.sequence
    target.validate(sequence)
    if not isinstance(estimator, str) or not estimator or not isinstance(scenario, str) or not scenario:
        raise ValueError("Estimator and scenario must be nonempty strings")
    if isinstance(n_paths, bool) or not isinstance(n_paths, Integral) or n_paths < 1:
        raise ValueError("n_paths must be a positive integer")
    horizon = _nonnegative(horizon_hours, "horizon_hours")
    clearance = _nonnegative(clearance_rate_per_hour, "clearance_rate_per_hour")
    uptake = _nonnegative(uptake_rate_per_hour, "uptake_rate_per_hour")
    if isinstance(seed, bool) or not isinstance(seed, Integral):
        raise ValueError("seed must be an integer")
    rng = random.Random(int(seed))
    cache = {}

    def kinetics(fragment):
        if fragment not in cache:
            half_life = half_lives.get(fragment)
            if half_life is not None:
                half_life = _nonnegative(half_life, "fragment half-life")
            if half_life and len(fragment) > 1:
                weights = tuple(_nonnegative(w, "cut weight") for w in cut_weights(fragment))
                if len(weights) != len(fragment) - 1:
                    raise ValueError("cut_weights must provide one weight per internal bond")
                total = math.fsum(weights)
                if not math.isfinite(total) or total <= 0:
                    raise ValueError("Positive half-life requires a positive, finite total cut weight")
                cumulative, running = [], 0.0
                for weight in weights:
                    running += weight / total
                    cumulative.append(running)
                cumulative[next(i for i in range(len(weights) - 1, -1, -1) if weights[i] > 0)] = 1.0
                cache[fragment] = (half_life, math.log(2) / half_life, tuple(cumulative))
                if not math.isfinite(cache[fragment][1]):
                    raise ValueError("Fragment half-life yields a nonfinite rate")
            else:
                cache[fragment] = (half_life, None, ())
        return cache[fragment]

    paths = []
    for _ in range(n_paths):
        start, end, time, steps = 0, len(sequence), 0.0, []
        while True:
            half_life, cleavage_rate, cumulative = kinetics(sequence[start:end])
            if cleavage_rate is None:
                steps.append(DegradationStep(start, end, time, time, half_life, "unknown",
                                             reason="Missing, zero, or unmodeled fragment kinetics; no survival inferred"))
                break
            total_rate = cleavage_rate + clearance + uptake
            if not math.isfinite(total_rate):
                raise ValueError("Total event rate must be finite")
            event_time = time + rng.expovariate(total_rate)
            if event_time > horizon:
                steps.append(DegradationStep(start, end, time, horizon, half_life, "censored"))
                break
            draw = rng.random() * total_rate
            if draw < clearance:
                steps.append(DegradationStep(start, end, time, event_time, half_life, "cleared"))
                break
            if draw < clearance + uptake:
                steps.append(DegradationStep(start, end, time, event_time, half_life, "taken_up"))
                break
            draw = rng.random()
            # Last cumulative value can be microscopically below one due to
            # rounding. It still represents the final supported internal bond.
            local_bond = next((b for b, value in enumerate(cumulative, 1) if draw < value),
                              len(cumulative))
            bond = start + local_bond
            destroyed = target.start < bond < target.end
            steps.append(DegradationStep(start, end, time, event_time, half_life,
                                         "target_destroyed" if destroyed else "cut", bond))
            if destroyed:
                break
            if bond <= target.start:
                start = bond
            else:
                end = bond
            time = event_time
        paths.append(DegradationPath(estimator, scenario, peptide_input, target, len(sequence), horizon, tuple(steps)))
    return tuple(paths)


def degradation_curve(paths, times_hours):
    """Return mutually exclusive outcome fractions at each requested time.

    ``target_in_circulation`` is intact parent + intact target-bearing fragments.
    ``unknown`` is never counted as survival or loss. Its addition to target
    retention gives an uncertainty upper bound for unassessed paths, rather
    than a statistical confidence interval. Censoring at the horizon preserves
    the last assessed state. Different estimators/scenarios must not be pooled.
    """
    paths = tuple(paths)
    if not paths:
        raise ValueError("At least one path is required")
    first = paths[0]
    def identity(path):
        return (path.estimator, path.scenario, path.peptide_input, path.target, path.parent_length, path.horizon_hours)
    if any(identity(p) != identity(first) for p in paths):
        raise ValueError("Do not pool different estimators, scenarios, or target occurrences")
    rows = []
    for value in times_hours:
        time = _nonnegative(value, "curve time")
        if time > first.horizon_hours:
            raise ValueError("Curve times cannot exceed the observation horizon")
        counts = dict(parent_remaining=0, target_in_fragment=0, target_destroyed=0,
                      cleared=0, taken_up=0, unknown=0)
        for path in paths:
            state = "parent_remaining"
            for step in path.steps:
                if step.event_hours > time or step.event == "censored":
                    break
                if step.event == "cut":
                    state = "target_in_fragment"
                else:
                    state = step.event
                    break
            counts[state] += 1
        row = {key: count / len(paths) for key, count in counts.items()}
        row.update(time_hours=time, estimator=first.estimator, scenario=first.scenario,
                   n_paths=len(paths))
        row["target_in_circulation"] = row["parent_remaining"] + row["target_in_fragment"]
        rows.append(row)
    return rows


def summarize_target_degradation(paths):
    """Summarize conditional target retention without calling it circulation PK.

    Parameters
    ----------
    paths : iterable of DegradationPath
        Lineages for one target occurrence, estimator and scenario.

    Returns
    -------
    dict
        ``retention_median_hours`` is the first time at which at least half
        the sampled targets have been split, cleared or taken up. Those causes
        remain separate in ``outcome_fractions``. This is conditional on the
        supplied kinetics and cut weights, not an experimentally calibrated
        epitope half-life. Unknown paths cause abstention. If fewer than half
        have left the tracked compartment by the horizon, report a lower
        bound rather than replacing their lifetimes with the horizon.
    """
    paths = tuple(paths)
    if not paths:
        raise ValueError("At least one path is required")
    first = paths[0]
    final = degradation_curve(paths, [first.horizon_hours])[0]
    endings = [p.steps[-1] for p in paths]
    losses = sorted(step.event_hours for step in endings
                    if step.event in ("target_destroyed", "cleared", "taken_up"))
    # Empirical first passage, not a median of truncated observation times.
    threshold = (len(paths) + 1) // 2
    median, lower_bound = None, None
    if final["unknown"]:
        status = "unassessed_paths"
    elif len(losses) >= threshold:
        status, median = "conditional_estimate", losses[threshold - 1]
    else:
        status, lower_bound = "beyond_horizon", first.horizon_hours
    return dict(
        estimator=first.estimator, scenario=first.scenario,
        target_label=first.target.label, target_start=first.target.start,
        target_end=first.target.end, n_paths=len(paths),
        horizon_hours=first.horizon_hours,
        parent_half_life_hours=first.steps[0].half_life_hours,
        retention_median_hours=median, median_status=status,
        retention_median_lower_bound_hours=lower_bound,
        outcome_fractions={key: final[key] for key in (
            "parent_remaining", "target_in_fragment", "target_destroyed",
            "cleared", "taken_up", "unknown")},
        calibrated_epitope_half_life_hours=None,
        interpretation="Conditional target retention; cut locations and event clocks are assumptions")
