# Calibration and identifiability of vaccine core exposure

This audit addresses [#549](https://github.com/openvax/mhctools/issues/549):
can published measurements identify, for **1 µg of a specified IM vaccine
peptide window**, the fraction of its core-bearing lineage that first reaches
draining-node APCs, or the core's intrinsic blood half-life? For an mRNA
construct, can they identify antigen production, transfer, processing and
surface peptide–MHC (pMHC) kinetics? **The reviewed experiments do not identify
transferable rates for those endpoints.** The [window workflow](../vaccine-windows.md)
therefore computes them only under supplied, labeled sensitivity scenarios.

## The endpoint and its defensible bounds

For one core occurrence in one free SLP molecule, the model counts **first
arrival of a still core-bearing antigen lineage** into a draining-node APC. A
fragment shorter than the dosed peptide can qualify; carrier or label alone
cannot. The cumulative fraction is between **0 and 1** because each dosed
lineage is counted once. Those are physical bounds, not a useful calibrated
range: with no matched fragment and APC data, both very low and very high
arrival remain compatible with the reviewed evidence. If several cores occur
in a window, these are separate per-core fractions; adding them does not give
the fraction of unique dosed molecules reaching APCs.

For mRNA, the analogous output is *antigen copies first arriving per input
transcript*. Repeated translation can make it greater than 1, so it is not a
fraction of RNA dose. Neither a finite upper bound on that copy ratio nor a
source-specific value follows from the studies below without translation and
degradation measurements for the actual construct. The model can also divide
arrivals by antigen copies produced over the observation horizon when that
production is supplied.

An observed blood or pMHC decline is not an intrinsic molecular half-life
while material continues to enter the compartment. Even an ideal absolute
curve `C(t)` in a one-compartment model satisfies `dC/dt = input(t) - k C(t)`.
For any sufficiently large loss rate `k`, a nonnegative
`input(t) = dC/dt + k C(t)` reproduces the same curve. Thus the loss rate and
`ln(2)/k` cannot be uniquely recovered from the curve without measuring or
shutting off input. A falling curve could constrain `k` only under that
specific model and known absolute analyte; a fluorescence ratio across
changing APC populations cannot supply even that bound. The reported
[post-peak half-decline](../vaccine-core-exposure.md#blood-persistence-of-a-core)
is a property of the observed curve, not a blood-clearance or surface-pMHC
half-life.

## What the primary studies actually observe

| Study and intervention | Useful observation | What it does not identify |
| --- | --- | --- |
| [Hassett et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10755037/): IM reporter mRNA–LNP in cynomolgus macaques, 300 µg RNA per site | Node reporter RNA was 155 ng/g at the first 8-hour tissue sample; in a separate cohort, 9.1% of conventional DCs were EGFP-positive at 24 hours. | Reporter RNA or protein does not count active transcripts, target-bearing peptide, APC-arrival fraction or surface pMHC. The two cohorts cannot form a paired conversion rate. |
| [Kent et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11447892/): 19 people after IM bivalent mRNA-1273 booster | Linked plasma vaccine-RNA signal had a reported post-peak slope of 0.608/day and half-decline of 1.14 days. | Linked short RNA amplicons are not full-length translation-competent RNA or antigen; ongoing input prevents assigning this decline to peptide, RNA or pMHC clearance. |
| [Muro et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378541/): mice given IM OVA mRNA–LNP plus reporter mRNA–LNP, 1 µg RNA of each | Antigen-specific SIINFEKL/H-2Kᵇ staining relative to total Kᵇ on GFP-positive draining-node cDC2 peaked at day 1, remained above control through day 5 and was not detected at day 7. | This is a relative cell-surface pMHC signal, not absolute complexes or source shutoff. It establishes a presentation time course in that mouse formulation, not an intrinsic pMHC lifetime or human construct rate. |
| [Jo et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13089314/) and [Rood et al., 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC12873236/): mouse IM mRNA–LNP pathway perturbations | Preformed pMHC-I transfer can support CD8 priming, while direct APC expression can be important for MHC-II/CD4 priming. | Relative pathway interventions do not specify human route fractions or effective transfer, secretion and loading rates for a new construct. |
| [Moynihan et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6247902/): subcutaneous albumin-binding amphiphile peptide versus soluble peptide in mice | Conjugation increased node-associated signal and supported antigen presentation for at least a week; soluble-peptide presentation declined faster. | Formulation-specific signal and T-cell readouts do not quantify intact core-bearing molecules in IM lymph or APCs, nor isolate the lifetime of a single pMHC complex. |
| [Backlund et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9371699/): subcutaneous gp100–cell-penetrating-peptide conjugates in mice | With 25 nmol labeled peptide, node fluorescence was compared at 48 hours and peptide-positive APC subsets at days 2 and 6. | Fluorescence or cell positivity does not distinguish intact peptide, core-bearing fragments and released label or count a dose-normalized APC core fraction. |
| [Rabu et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6422379/): artificial SLP linkers tested in human DCs and T-cell assays | Some linker choices shifted class-I cross-presentation potency by up to about 100-fold relative to the corresponding native long peptide. | This is a within-assay concentration-response effect, not an in-vivo per-hour cleavage, uptake or transport rate. |
| [Mouse lymph DPP4 study, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4054822/) | Reported lymph activity is 4.7 ± 0.3 (abstract, n=23) or 7.2 ± 0.8 (results, n=9), versus plasma 22.3 ± 0.9 (n=8), in pmol/min/µL. | Different paper summaries are retained as reported. Assay activity on its substrate is not an inverse-hour hazard for an arbitrary SLP; mouse mesenteric free fluid is not human IM afferent lymph or APC surface exposure. |

The [mRNA measurement extract](../evidence/mrna-antigen-presentation-measurements.json)
retains the available values, units, source figures and assay context. These
observations can check an **assay-matched** model or motivate measurements. In
particular, the Muro time points can reject a scenario that predicts no mouse
surface pMHC at day 1, but cannot determine which production and clearance
rates generated the observed curve.

## Measurements that would identify the missing steps

For each proposed sequence, route and formulation, collect the starting
chemical dose and a timed, mass-balanced panel of the intact parent and
**sequenced core-bearing fragments** in injection-site fluid, afferent lymph,
node fluid, sorted APCs and blood. Measure absolute internalized copies per
APC subset and cell counts, distinguishing cell-associated surface material
from uptake. Pair those samples with allele-specific surface pMHC copies per
cell and the same core's T-cell readout; a labeled carrier or bulk antigen is
not a substitute.

For mRNA, pair active intracellular RNA with completed antigen, secretion or
release, and target-bearing fragments in producer and APC populations. Compare
direct APC expression, extracellular antigen transfer and donor-allele-matched
preformed-pMHC transfer. A pulse–chase or controlled stop of production,
uptake and transfer is needed to separate new pMHC loading from surface loss.
When testing enzymes, use the *current fragment* in the relevant lymph,
blood or APC compartment, and distinguish soluble activity from cell-surface
exposure and target-preserving from target-destroying cuts. Test uptake and
MHC loading for each relevant fragment across concentrations and competing
ligands; a single first-order rate for all fragments is a scenario assumption.

Until those matched measurements exist, the [trajectory rate contract](../vaccine-trajectory.md#kinetic-contract)
correctly leaves missing reachable rates **unassessed**. Declared assumptions
can produce transparent sensitivity scenarios, but their spread is not a
measured confidence interval or a recommended physiological dose fraction.
The model has no validated threshold converting pMHC exposure to a T-cell
response or an “effective dose.”

The remaining empirical work is tracked separately for
[APC uptake versus presentation](https://github.com/openvax/mhctools/issues/303),
[node-exposure measurements](https://github.com/openvax/mhctools/issues/306)
and [intracellular protease coverage](https://github.com/openvax/mhctools/issues/470).
