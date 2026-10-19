# Backend download and runtime reliability

Address #568, #594, and the reproducible parts of #599 without changing predictor
scores or hiding backend failures.

- Retry PhageScout downloads a bounded number of times for transport timeouts
  and temporary HTTP failures. Restart an interrupted file from byte zero,
  preserve checksum/size validation, and publish only a complete installation.
- Reproduce the configured legacy NetMHC launch probes and multi-allele test.
  Repair a demonstrated launcher defect if it belongs here; keep unresolved
  emulation failures explicit and tracked.
- Make installed legacy backends that fail their probes fail integration tests,
  while genuinely absent optional executables may still skip. Print skip reasons
  in the default test script.

Verify retry success, exhaustion, permanent failures, and partial-transfer
cleanup with deterministic tests. Exercise real legacy predictions with Docker,
then run lint, the full test suite with no skips, and CI. Bump the patch version,
land through a PR, and deploy from clean master after CI passes.

## Findings

The #594 Python 2 error reproduces with the native bundle launcher on the
ordinary shell PATH, where Python 2 is missing. Sourcing the generated activation
file selects the Docker launcher and passes all ten legacy tests. Document this
required activation for direct commands; a container patch is not justified by
the observed error. Keep #599 open because a passing run does not resolve an
intermittent QEMU crash.
