# Python 3.11 baseline (#601)

## Scope

Raise the mhctools host minimum to Python 3.11, with no upper bound. Release
as 4.2.0 because this changes interpreter support. Keep separately provisioned
upstream runtimes at their tested versions.

## Changes

- Declare Python >=3.11 and classifiers for 3.11–3.14; target Ruff at py311.
- Remove the obsolete NumPy dependency split and use standard-library tomllib.
- Update active installation and testing documentation.
- Run installation, dependency checks, CLI startup, builds, lint and all offline
  tests on Python 3.11, 3.12, 3.13 and 3.14 in CI. Retain real-model CI jobs.

## Verification

Use fresh repository-local environments to investigate compatibility failures;
do not change the shared environment or force incompatible dependencies. File
actual dependency blockers with the responsible project and link them in the PR.
Run ./lint.sh and ./test.sh before merging green CI. Deploy from clean master
with the required-backend release gate, then verify the PyPI artifacts.
