# Interpretable SLP and mRNA exposure workflow

Implement three independently releasable PRs on the existing conditional
single-target simulator. Keep the model small, the units visible and the
evidence traceable. Calibration remains tracked in #549, with validation in
#291 and paired uptake/presentation and node-exposure evidence in #303/#306.

1. **Exposure integrals.** Integrate surface pMHC over time separately for
   local, draining-node and systemic APCs, and extracellular target availability
   separately by compartment. Integrate within the linear solver from injection
   to each observation time, without dependence on plot sampling. Keep antigen
   counts and copy-hour observables separate. Verify against analytic solutions,
   dose scaling, migration and missing-kinetics behavior. Explain that an AUC is
   an exposure descriptor whose relation to immune response is uncalibrated.
2. **Dose accounting and readable reports.** Make supplied peptide versus
   carrier-bound transcript amounts explicit, with optional documented molar or
   mass conversion requiring the actual material's molecular weight. Preserve
   normalization and administration context. Produce a compact report covering
   target, dose, observation horizon, presentation exposure, antigen fates,
   missing inputs and measured/assumed/disabled rates. Never infer RNA mass from
   the encoded peptide or turn transport/translation rates into half-lives.
3. **Grounded comparisons.** Supply an end-to-end comparison of explicit
   scenarios with compatible targets, alleles, horizons and dose bases. Show
   absolute outputs and source normalization separately, preserve assumptions,
   and avoid efficacy rankings. Document the primary evidence behind endpoints
   and the measurements still needed to identify rates. Reproduce synthetic
   examples through API and CLI; separate them visibly from study observations.

For each PR: bump version on the feature branch; lint, full tests, strict docs
and link checks; green GitHub CI; merge; deploy from clean master and verify
the published version. Use the repository release environment, preserving the
shared environment and external predictor installations. File discovered bugs
on the relevant repository. Review related open issues after each release.

Acceptance: a reader can tell what was administered, what the model follows,
what its numbers mean, which values are assumptions, why a missing-input run
has no prediction, and what evidence would support a real SLP/mRNA comparison.
This series does not claim to finish biological calibration without data.
