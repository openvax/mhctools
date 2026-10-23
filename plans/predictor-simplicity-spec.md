# Predictor and test cleanup (#353, #443)

## Scope

- Remove the unused singular protein-flank helper. The live plural helper
  continues to compute separate N/C context lengths.
- Remove the MHCflurry import-error fallback from artifact discovery:
  MHCflurry is a required dependency, and broken imports should surface as
  errors rather than being reported as missing model downloads.
- The unused IEDB method lists have already been removed. Retain the public
  `BindingPrediction.elution_score` compatibility alias; lack of internal
  callers does not establish that external applications no longer use it.
  Document its supported replacement, `score`, and verify alias behavior.
- Replace signature inspection with constructor calls and observed batching,
  peptide splitting, and DataFrame/flank-validation behavior. Replace
  inheritance/attribute-presence checks with observable predictor behavior
  where coverage is needed; retain the intentional removed-hook regression.
- Release as 4.4.1 after the pandas 3 baseline ships.
- Replace the owned BigMHC `set_index(inplace=True)` call deprecated by pandas
  3.1 with assignment. Track the upstream `reset_index` warning in #615.

## Verification

- Focused offline tests must exercise the actual constructors and public
  DataFrame methods, mocking only external predictor execution/discovery.
- A regression must demonstrate that import failures propagate from artifact
  discovery while missing model files still report a missing artifact.
- Run lint, complete tests with all configured backends required, and GitHub
  CI before merge. Deploy from clean master and verify PyPI plus the tag.
