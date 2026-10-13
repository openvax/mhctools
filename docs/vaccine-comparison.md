# Interpreting SLP and mRNA exposure

The workflow answers a conditional question: **given this starting material and
these rates, how much of this target reaches APC surface MHC, where, and for how
long?** It supports free SLP, secretion-tagged mRNA-LNP and cytosolic-antigen
mRNA-LNP after local injection. See [core delivery and blood exposure](vaccine-core-exposure.md)
for first-arrival fractions, endogenous expression and preformed pMHC transfer.
It does not infer an effective immune dose from a peptide sequence.

## Three quantities to read together

| Question | Quantity | Interpretation |
|---|---|---|
| What remains available? | Target-bearing antigen copies over time | Includes shorter fragments retaining the chosen target; differs from intact-parent persistence |
| What is displayed now? | Surface pMHC copies in each APC compartment | Current complexes, after loading, trafficking and turnover; not copies per individual APC |
| How much display occurred? | Surface pMHC exposure in complex-hours over `[0,T]` | Integral of the curve; retain timing, location and horizon alongside the area |

Use draining-node presentation as one clearly named endpoint, while inspecting
local and systemic display too. Equal areas can have different peaks, timing and
cellular distributions. The current model has no T-cell or antibody response
equation, and no validated threshold defining an effective dose.

## Run and inspect an explicit scenario

Start with the [trajectory input contract](vaccine-trajectory.md#kinetic-contract).
Select the exact target occurrence and allele, the delivery route and actual
starting material. Supply measured rates only where the study and material apply;
label sensitivity assumptions as `assumed` and intentional exclusions as
`disabled`. A recognition score or binding rank is not an inverse-hour rate.

```sh
mhctools vaccine-trajectory --input scenario.json \
  --output trajectory.json --report trajectory.md
```

The report starts with the material, source normalization, context and horizon.
It then shows presentation, antigen accounting and provenance. An incomplete
scenario instead shows the missing steps and fragment assessments. Resolve those
against evidence or declare explicit scenario assumptions; do not fill unknowns
with zeros merely to obtain a curve.

The amount is material entering the model: released free peptide or carrier-bound
RNA. [Document copies, pmol or mass](vaccine-trajectory.md#starting-amount-and-readable-reports).
For RNA mass, provide the full transcript's molecular weight. For an SLP, provide
the actual peptide's molecular weight. Equal micrograms imply different molecule
counts; equal molecule counts do not imply equal antigen production.

## Compare matched endpoints

```sh
mhctools vaccine-compare --input slp-trajectory.json mrna-trajectory.json \
  --output comparison.json --report comparison.md
```

Inputs are **result JSON files**, not scenario input files. The comparison
uses results from mhctools 3.50.0 or newer; rerun older scenario inputs to obtain
integrated exposure and explicit dose accounting. The comparison
requires matching target sequence, allele, MHC class and observation horizon.
Original target coordinates can differ between constructs. Known organism
differences are rejected; unknown context and route differences remain visible.
Different MHC-II flanks or construct sequences can still affect processing and
loading, so inspect each scenario's full fragment and rate audit.

The table preserves input order and shows absolute draining-node exposure next
to exposure per input peptide/transcript. Its dose statement distinguishes equal
documented starting mass, equal starting molecule counts, unequal amounts and
one-input normalization. It does not sort by efficacy, invent a common probability
or hide an unassessed scenario behind a zero. The comparison JSON retains rate
values/sources and hashes of its input files; individual results retain all curves
and cleavage assessments.

For sensitivity analysis, change a documented rate or starting amount in a copy
of the scenario, rerun, and compare at the same horizon. Label the change in the
scenario name and rate source. A spread across assumptions is not a confidence
interval; no uncertainty distribution has been fitted.

## Reproduce the complete demonstration

```sh
python scripts/vaccine_trajectory_example.py --output /tmp/vaccine-demo --plot
mhctools vaccine-trajectory --input /tmp/vaccine-demo/slp-input.json \
  --output /tmp/slp-check.json --report /tmp/slp-check.md
mhctools vaccine-compare \
  --input /tmp/vaccine-demo/slp-trajectory.json /tmp/vaccine-demo/mrna-trajectory.json \
  --output /tmp/comparison-check.json --report /tmp/comparison-check.md
```

The script writes input files, full results, individual readable reports, a
comparison and optionally a plot (`--plot` requires the `vaccine-report` extra).
**Every rate is synthetic.** Each panel labels its peptide/transcript denominator;
the panels do not establish equal-dose efficacy. Use new output paths on each run.

## What is grounded, and what still needs calibration?

| Primary evidence | What it supports | What it does not identify |
|---|---|---|
| [Moynihan et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6247902/) | Formulation-dependent lymphatic distribution and presentation persistence in peptide-vaccine experiments | Universal free-SLP drainage or destruction rates |
| [Liang et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28958578/) | Separate LNP uptake and expression measurements in rhesus-macaque injection-site and draining-node APCs | Absolute target-specific pMHC kinetics for arbitrary human constructs |
| [Dasti et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12256567/) | An existing multiscale mRNA-vaccine model connecting modeled MHC-II exposure to cellular dynamics, with vaccine-specific calibration | Direct measurement of all latent presentation rates, or transferable SLP/secreted-mRNA defaults |

The Dasti model uses macaque APC data and human antibody data. Its calibration
is relevant prior work, but antibody agreement does not by itself validate every
intermediate rate. This implementation does not import those fitted coefficients.

To calibrate the present model, pair time-resolved measurements of actual starting
material, intact parent and target-bearing fragments, active intracellular RNA,
completed/secreted antigen, and allele-specific surface pMHC in relevant APC
locations. Retain cell counts when converting per-cell assays to compartment
totals. Establish which rates are identifiable before fitting them, and evaluate
independent constructs/donors/studies after checking source overlap. The
[serum audit](serum-calibration.md) supplies a bounded example of that evidence
discipline; it does not calibrate lymph or APC kinetics.

The machine-readable `vaccine_trajectory_evidence()` inventory preserves these
endpoint limits and calibration requirements. The
[calibration audit](research/vaccine-calibration.md) establishes the current
non-identifiability and minimal matched measurements. Future biological
calibration remains supported by
[independent benchmarks](https://github.com/openvax/mhctools/issues/291),
[uptake/presentation validation](https://github.com/openvax/mhctools/issues/303)
and [node-exposure evidence](https://github.com/openvax/mhctools/issues/306).
For primary mRNA mechanism studies and reported calibration observations, see
the [mRNA antigen-presentation bibliography](evidence/mrna-antigen-presentation.md).
