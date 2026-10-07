# Following an epitope through successive cuts

An intact vaccine peptide can disappear while its target survives in a shorter
fragment. `simulate_target_degradation` follows one exact, coordinate-defined
target through successive cuts. A cut in a flank or exactly at a target boundary
retains the target-containing product. A cut strictly inside the target destroys
that exact sequence. The resulting product gets its own half-life estimate and
cut-location assessment, including its newly exposed termini.

This is an **exploratory serum-digestion scenario**, not a validated prediction
of human circulation half-life. The [half-life estimators](predictors/peptide-pk.md)
describe different assay settings. Each estimator supplies a separate time scale;
run and display them side by side. Never average their half-lives or their sampled
trajectories. Agreement is informative; disagreement is unresolved model or assay
uncertainty. PeptiVerse's sequence endpoint has only 130 training examples, and
neither estimator establishes vaccine-fragment accuracy.

Stable cleavage products are biologically possible: the
[RNase 3 peptide study](https://doi.org/10.1021/acs.jmedchem.1c00795)
measured stable byproducts after serum digestion. That experiment supports the
distinction between parent disappearance and product persistence; it does not
validate this simulation or establish antigen presentation by those products.

## What is assumed

The simulation assumes a well-mixed, dilute, first-order degradation process.
For a fragment with supplied half-life `h` hours, its next cleavage waiting time
is exponential with rate `ln(2) / h`. These kinetics are an assumption even when
`h` comes from a predictor. Relative cut weights allocate that single cleavage
hazard among the fragment's internal bonds. Enzyme scores do **not** provide
those rates: DPP4 depletion, PhageScout enrichment, MMP substrate scores and
recognition rules have different endpoints and scales.

Supply cut weights as an explicit scenario, such as uniform weights versus
recognition-weighted sensitivity scenarios. A positive background weight avoids
treating unmodeled enzymes or missing coefficients as protection. Recompute
recognition on every fragment: an N-terminal peptidase must inspect the current
N-terminus, rather than reuse its original-parent score. Weighting a bond shared
by several rules does not establish independent enzyme evidence. Do not insert
proteasome or lysosomal scores into extracellular rate allocation.

Separately measured or explicitly assumed clearance and uptake rates can compete
with degradation. A clearance event removes an intact target from circulation;
an uptake event transfers it out of this compartment. Neither establishes target
destruction or productive presentation. Defaults of zero model serum digestion
alone, omitting renal clearance, tissue distribution, uptake, binding, formulation,
injection-site release and changing enzyme activity. Do not interpret these
defaults as a patient's physiological rates.

## Run separate estimators

The sampler accepts a sequence-to-hours table so expensive inference can be
batched and audited before simulation. All possible contiguous fragments
retaining a given target are enumerated in original-parent coordinates.

```python
from mhctools import (
    CavacoHalfLife, DegradationTarget, PeptiVerse, degradation_curve,
    simulate_target_degradation, target_fragments,
)

parent = "ACDEFGHIKLMNPQRSTVWY"
target = DegradationTarget("exact target", start=5, end=14)  # [5, 14)
sequences = sorted({f.sequence for f in target_fragments(parent, target)})

def uniform_cuts(fragment):
    return [1.0] * (len(fragment) - 1)

curves = {}
for name, predictor in [
    ("PeptiVerse human serum", PeptiVerse(device="cpu")),
    ("Cavaco published baseline", CavacoHalfLife()),
]:
    results = predictor.predict(sequences, on_unsupported="record")
    # Keep predictor versions, artifact identities and measurement contexts
    # alongside this numerical table in the caller's provenance record.
    hours = {s: r.peptide_half_life.value for s, r in zip(sequences, results)}
    paths = simulate_target_degradation(
        parent, target, hours, uniform_cuts,
        estimator=name, scenario="uniform cuts; no clearance or uptake",
        n_paths=5000, horizon_hours=24, seed=42,
    )
    curves[name] = degradation_curve(paths, [0, 0.5, 1, 2, 4, 8, 24])
```

Plot `parent_remaining` and `target_in_circulation` together for each estimator.
The latter includes the original parent and every tracked fragment retaining the
entire target. Keep `target_destroyed`, `cleared`, `taken_up` and `unknown`
separate. Missing or zero half-life estimates stop in `unknown`, never survival.
There are no events after the observation horizon; the last assessed state is
right-censored. Scenario ranges are sensitivity ranges, not confidence intervals.
Sampling error decreases with more paths, but extra paths cannot resolve model
or biological uncertainty.

An exact class-I candidate and a class-II binding core are different targets.
For class II, a flank cut in a longer ligand may leave its binding core intact;
map and assess that core separately. If the mutation lies in a class-II flank,
track the core and mutant residue together: retaining only the wild-type core
does not retain the mutant target. Flanking residues can affect TCR recognition
([primary experiment](https://www.nature.com/articles/ncomms1665)), so even core
and mutation retention does not establish unchanged T-cell recognition.
Validate mutant-position inclusion before
labeling a predicted binder as mutant. An unknown target or ambiguous alignment
cannot be labeled protected. A class label on a whole vaccine construct does not
locate a binding core and must not turn all flank cuts into epitope destruction.
Individual-target paths provide marginal retention,
not the probability that **all** alternative targets are lost together.

## Validation status

Synthetic tests check original-parent coordinates, boundary cuts, sequential
trimming and re-recognition, fragment-specific half-lives, exponential survival
and competing-hazard analytical limits, missing estimates, horizon censoring,
chemistry rejection and deterministic seeds. These verify the sampler's math,
not its serum biology. Calibrating its biological inputs requires time-resolved
parent/product measurements, measured cut frequencies, assay-specific fragment
stability and applicable clearance or uptake data.
