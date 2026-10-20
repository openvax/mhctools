# Optional predictor selection for embedded CLIs (#444)

## Scope

Allow a report application to reuse mhctools argument definitions when replaying
saved predictions. The standalone CLI and existing callers keep requiring an
explicit predictor. Build on the Python 3.11 baseline shipped in #605 / 4.1.5;
release the backward-compatible embedding API as 4.2.0.

## Changes

- Add keyword-only `predictor_required=True` to `add_mhc_args`.
- Optional parsing preserves all predictor, allele and runtime options without
  constructing predictors or reading allele files.
- Both predictor construction entry points reject a missing or empty selection
  with a concise ValueError. Existing allele and predictor-name validation stays
  strict when callers request fresh prediction.
- Document an embedding example and remove the stale #444 dispatch comment.
- Complete verification cleanup from the superseded #606: use stdlib tomllib,
  configure headless plotting tests (#607), and check CLI startup in the existing
  Python/pandas CI matrix. Preserve that matrix and all predictor runtime pins.

## Verification

First reproduce the missing API with behavioral regression tests. Check cached
parsing, explicit prediction, required defaults, empty/unknown selections and
missing alleles. Run focused CLI tests, ./lint.sh and ./test.sh. Merge only after
all CI passes; deploy with ./deploy.sh from clean master and verify PyPI artifacts.
