# Research and evidence

Research pages trace what was measured, in which assay and population, and
which model claims those observations can support. They do not turn a paper's
reported signal or half-life into a default vaccine rate.

Use the rest of the documentation according to the question at hand:

| Section | Use it for |
| --- | --- |
| [Getting started](../getting-started.md) | Install mhctools, run a first prediction and read its result. |
| [User guide](../predictors/index.md) | Choose a predictor or workflow, supply inputs and interpret its output. |
| [Reference](../api.md) | Look up exact APIs, result fields, units, capabilities and limits. |
| Research | Inspect primary sources, measurement context, validation and what remains unidentifiable. |

## Literature reviews already available

- [mRNA vaccine antigen presentation](../evidence/mrna-antigen-presentation.md)
  compares direct APC expression, antigen transfer and cross-dressing. Its
  [measurement extract](../evidence/mrna-antigen-presentation-measurements.json)
  retains the assay context of 23 macaque and human quantitative observations
  and a mouse pMHC time course.
- [Dose-normalized lymph-node delivery](vaccine-node-delivery.md) collects
  primary SC and IM antigen, carrier and RNA measurements with their route,
  sampled time and denominator. It explains why whole-node percent of dose
  cannot directly calibrate core-bearing arrival inside APCs.
- [Epitope position and peptide flanks](epitope-position.md) asks whether a
  core near a protein or SLP end is presented more often. It separates
  controlled terminal-placement experiments from a universal position rule.
- [Antigen processing and T-cell responses](processing-t-cell-responses.md)
  traces measured changes in processing, surface presentation and T-cell
  outcomes, which are related but distinct endpoints.
- [Human serum peptidase calibration](../serum-calibration.md) audits 12 studies
  and 27 reported hormone-stability measurements. It explains why those
  experiments do not calibrate enzyme hazards for a vaccine peptide.

The [cleavage validation audit](../cleavage/validation.md) checks model
conformance, source observations and missing independent validation. The
[protease candidate audit](../cleavage/candidates.md) compares released assets
and experimental recognition evidence. These are focused evidence and model
audits rather than standalone reviews of all antigen processing. The
[assay-aware benchmark guide](../benchmarks.md) explains how to evaluate new
measurements without conflating source reproduction with independent testing.
The [vaccine calibration audit](vaccine-calibration.md) gives physical bounds
and shows why published RNA, label and pMHC time courses do not identify
core-specific IM transport or intrinsic clearance rates.

## Reviews needed for quantitative vaccine exposure

1. **IM SLP fragment transport:** find time-resolved measurements of the
   intact parent and core-bearing fragments at the injection site, in
   afferent lymph, draining nodes and blood. Preserve route, formulation,
   species and whether a signal is peptide, label or carrier
   ([exposure issue](https://github.com/openvax/mhctools/issues/306)).
2. **APC uptake through presentation:** seek matched measurements of
   internalization, productive localization, target-bearing fragments and
   surface peptide–MHC in the same antigen and APC system. Uptake or reporter
   fluorescence alone cannot establish an epitope dose fraction
   ([endpoint issue](https://github.com/openvax/mhctools/issues/303)).
3. **Surface pMHC persistence and transfer:** separate newly loaded complexes
   from acquired preformed complexes, with allele, APC type, time course and
   loss mechanism recorded. The [mRNA review](../evidence/mrna-antigen-presentation.md)
   establishes the route but does not supply a transferable lifetime
   ([calibration issue](https://github.com/openvax/mhctools/issues/549)).
4. **APC endolysosomal cleavage:** extend the current source audit to
   quantitative, condition-matched CatS, CatL, CatB and AEP product maps that
   identify both epitope generation and destruction
   ([coverage issue](https://github.com/openvax/mhctools/issues/470)).

The [SLP/mRNA calibration audit](vaccine-calibration.md) resolves the current
identifiability question for the reviewed sources. Source-matched experiments
remain necessary before fitted human IM rates can replace explicit scenarios.
