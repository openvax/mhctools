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

## Lead with the selected epitope

Report the best qualifying mutant candidate per MHC class and vaccine construct,
with its exact interval, allele, rank and mutation mapping. Class-I and class-II
ranks have different endpoints; do not use them to invent a single winner.
Keep disclosed source epitopes as references. If the mutation or target cannot
be localized, retain that gap instead of selecting a guessed target.

The main question is how long the selected target remains intact **in any
target-bearing fragment**. A parent-peptide half-life answers when the original
sequence disappears; one harmless flank cut counts as parent loss. A separate
estimate for the released target answers how long that exact sequence would
last if already present as a free peptide. It does not account for when it is
released, protection inside its parent, or presentation after uptake. Show both
sequence estimates for each estimator without averaging them.

```python
from mhctools import summarize_target_degradation

summary = summarize_target_degradation(paths)
summary["retention_median_hours"]  # conditional first-passage time; may be None
summary["median_status"]           # conditional_estimate / beyond_horizon / unassessed_paths
summary["outcome_fractions"]       # target split, clearance and uptake remain distinct
```

Never label this simulation median as a calibrated epitope half-life. Unknown
paths cause abstention; surviving past the observation horizon gives a lower
bound rather than a lifetime capped at 24 hours. Non-cleavage removal is not
target destruction. In a digestion-only scenario, removal rates are zero and
all assessed losses are cuts inside the selected target.

## Annotate the bond and its mechanism

```python
from mhctools import DegradationTarget, annotate_target_cleavage, predict_cleavage

target = DegradationTarget("synthetic target", 2, 6)
evidence = predict_cleavage("APACDEFG", models=("dpp4-qpisa", "fap-dipeptidyl"))
annotations = annotate_target_cleavage(target, evidence)
```

Each annotation retains the enzyme, model version, native score and endpoint,
assay, limitations and compartment. A cut strictly inside the target splits it;
a cut exactly at its boundary can release it intact; a flank cut trims its
precursor. These are consequences **if the cut occurs**, not occurrence
probabilities. A native substrate score, purified-enzyme depletion estimate and
partial motif match must retain different labels. No universal score threshold
or cross-enzyme sum establishes a serum rate. For conditional fragments,
`CleavageInput.source_start` preserves original coordinates and new termini.

Use short solid red marks between residues for candidate bonds, highlight the target,
and list the responsible model and evidence type. Separate extracellular
digestion from lysosomal, proteasomal and MHC-ligand processing. Prefer the human
Pepsickle model for human input; the all-mammal model is a separate sensitivity
comparison, not a second independent vote. Display unavailable coefficients or
model input limits explicitly. A non-match is not proof of protection.

Arbitrary 3x/10x weighting belongs in supporting sensitivity detail. It gives a
matched bond weight 3 or 10 versus weight 1 at unflagged bonds; enzyme scores do
not set that multiplier and matches do not stack. An unflagged sampled cut must
read **assumed cut; no supporting enzyme prediction**, not an unnamed protease.

## Assay and chemistry matter

In a [direct matrix comparison](https://doi.org/10.1371/journal.pone.0178943),
tested peptides had different degradation profiles in mouse blood, serum and
anticoagulated plasma. Serum preparation and anticoagulants can change protease
activity. Therefore, a serum estimate is not a human circulation lifetime.

Constrained structures, terminal chemistry and attachments can change peptide
persistence. A [cyclotide study](https://doi.org/10.1021/ja405108p) measured
a human-serum half-life of 55 hours for native MCoTI-I; an engineered cyclotide
and its linearized/reduced form had very different stabilities. The
[semaglutide label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/209637s038lbl.pdf)
attributes prolonged circulation principally to albumin binding and separately
documents DPP4 resistance. These examples do not establish the lifetime of a
linear vaccine construct. Unknown formulation, attachment, binding and terminal
chemistry must remain unknown; do not silently score a modified peptide as free.

The focused local exporter is `analyses/osteosarc_vaccine_cleavage/focused_stability_report.py`.
It reads a frozen target-cascade directory and writes a new report directory:

```sh
python analyses/osteosarc_vaccine_cleavage/focused_stability_report.py \
  --input /path/to/frozen-target-cascade --output /path/to/new-focused-report
```

It reuses native inference without changing predictor outputs, keeps all target
annotations in the audit, and separates the primary mutant candidates from
disclosed source references. Install the `vaccine-report` extra for PDF rendering.

The enzyme guide distinguishes soluble blood activity, cell-surface exposure,
intracellular processing and inflammation. There is no universal DPP4/FAP/ACE
strength ranking: a [human-plasma experiment](https://doi.org/10.1152/ajpheart.2000.278.4.H1069)
found ACE dominated bradykinin degradation at low concentrations while CPN
dominated at high concentrations. [Blood-specimen experiments](https://doi.org/10.1371/journal.pone.0134427)
identified DPP4-mediated trimming of particular hormones. Neither observation
sets relative vaccine cut rates. Fragment-specific inhibitor/LC-MS time courses,
or applicable catalytic efficiencies plus active enzyme exposure, are needed
to calibrate contributions. Enzyme abundance alone does not establish activity.

Report curves as the fraction of starting copies retaining the entire displayed
epitope/span, including in shorter fragments, versus the full vaccine peptide
remaining intact. A class-II span extended to retain a mutant flank is explicitly
distinguished from its minimal binding core. Released-sequence estimates start
with that sequence alone; illustrative epitope-loss times start with the full
vaccine peptide. All frozen antigen-processing models remain visible separately,
including cathepsin assay matrices and ERAP1. Matrix time labels are source assay
conditions, not predicted fragment lifetimes.

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
