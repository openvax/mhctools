# Cleanup sweep

## Scope

- Remove the unreferenced `mhctools.common` utilities and unused route-label
  constant. Retain exported compatibility aliases and provenance metadata.
- Use one CLI table formatter for artifacts, integration status and cleavage
  catalogs; preserve column order and explanatory notes, including empty tables.
  The reproduced empty-catalog crash is tracked in #623.
- Give the immutable subcommand record its own defaults, removing the forwarding
  factory. Dispatch continues to own aliases; remove the unused alias entry point.
- Share the identical legacy-Keras subprocess environment between the runtime
  probe and both adapters, preserving caller overrides and host environment.
- Group the documentation sidebar by task, shorten duplicated homepage navigation,
  and put the historical Python migration plan alongside the other plans.
- Preserve predictor algorithms, scientific semantics, public predictor APIs,
  command names, documentation URLs and existing linked anchors.

## Verification and release

Check observed CLI output, aliases, empty tables and environment isolation.
Run `./lint.sh`, `./test.sh`, documentation links/matrix checks and a strict
MkDocs build in the repository-local environment. Refresh its stale dependencies
without changing the shared environment. Bump the patch version, require green
GitHub CI, merge via PR, then deploy from clean master and verify PyPI.
