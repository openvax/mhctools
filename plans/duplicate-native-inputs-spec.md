# Deduplicate complete native scan inputs (#626)

## Problem

The flanked `BasePredictor.predict_proteins` path submits each source occurrence
separately. MHCflurry can assign slightly different floating-point results to
identical inputs within a batch, making otherwise identical cache records
contradict one another. The peptide-only scanner already deduplicates inputs.

## Change

- Submit each exact `(peptide, n_flank, c_flank)` tuple once per flanked scan,
  retaining first-occurrence input order. Alleles, genotype mode and model
  configuration remain fixed by the predictor instance for that invocation.
- Expand every returned prediction kind and allele back to the original
  occurrences, preserving source names, offsets, ordering and provenance.
- Keep differing, empty and boundary-shortened flank strings distinct.
- Retain backend result-count validation, applied to distinct inputs.
- Preserve the peptide-only scan path, direct prediction APIs, scoring models,
  numerical precision and downstream contradictory-cache checks.

## Verification and release

First reproduce the bug using the reported native MHCflurry scan. Add offline
regressions with an occurrence-sensitive backend, covering repeated inputs
within/across proteins, different contexts, ordering, empty scans and malformed
backend output. Exercise all MHCflurry kinds in haplotype and per-allele modes,
then verify exact duplicate agreement with the real model.

Run lint, the complete suite with `--require-all`, and documentation checks.
Bump to 4.5.1, require green CI, merge via PR, deploy from clean master and
verify the published artifacts. Use the repository-local environment.

The reported numerical discrepancy did not reproduce locally with either
MHCflurry 2.2.1/release 2.2.0 or 2.3.16/release 2.3.0. The deterministic
occurrence-sensitive regressions do fail before the scanner change; use them
to establish the duplicate-input defect independently of floating-point backend
behavior, and retain real-model coverage of duplicate agreement.
