# pandas 3 baseline (#611)

## Scope

- Require `pandas>=3,<4` in the mhctools host and test Python 3.11–3.14.
- Remove the `mixmhcpred` extra, whose pandas 2 requirement conflicts with the
  host. Document the existing isolated MixMHCpred/PRIME runtime as the install
  path; keep its backend dependencies in that environment.
- Remove obsolete host compatibility code if present. Preserve explicit dtypes
  and missing-value handling that implement result semantics under pandas 3.
- Update installation and CI documentation and release as 4.4.0.

## Verification

- Resolve a fresh repository-local install and run `pip check`, import and CLI
  smoke checks. Inspect built distribution metadata for the pandas requirement.
- Run lint, docs links/build, focused table/runtime tests, and the complete
  suite with every configured backend required.
- Require green GitHub CI for Python 3.11–3.14 and real backend integration,
  including pandas 3 in the host with pandas 2 isolated for MixMHCpred/PRIME.
- Merge, deploy from clean master using `./deploy.sh`, and verify PyPI and tag.
