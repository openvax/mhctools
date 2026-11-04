# Prediction correctness and sequence utilities release stack

## Scope and order

1. **#629 / 4.5.4:** verify native NetMHCpan 4.1 positions with a known
   protein, use the existing zero-based 4.1 parser, remove the duplicate
   misnamed parser, and cover BA/EL parsing and real protein scans.
2. **#627 / 4.6.0:** audit Biotite peptide mass conventions against primary
   sources and independent values; migrate only if average mass, terminal
   water and unsupported chemistry remain explicit. Consolidate compatible
   routine dependencies, documenting retained runtime boundaries.
3. **#634 / 4.6.1:** use Biotite FASTA iteration for MixMHCpred input/output,
   preserving full headers, ordering, alignment gaps and malformed-input
   rejection. This depends on the Biotite requirement in the second PR.

## Verification and shipping

Each PR includes its version and migration notes, targeted regressions,
`./lint.sh`, `./test.sh --require-all` with two workers, dependency checks,
and green GitHub CI. Review the complete diff before merging in order;
run `./deploy.sh` from clean master after each merge and verify the PyPI
release before merging the next PR. Use the repository-local Python 3.12
release environment; do not modify the shared environment or isolated
predictor runtimes. Retarget successive PRs after their parent merges.

Existing PRs #628 and #630 overlap other cleanup and documentation; this
stack does not absorb those unmerged changes. Related upstream issues were
triaged in topiary, mhcflurry and mhcgnomes; scientific evidence gaps and
upstream runtime warnings are not silently resolved here.

## Dependency audit adjustment

Biotite 2.0 requires Python >=3.12. Keep mhctools' Python 3.11 floor and use
`biotite>=1.6,<3`: 1.6.0 is the newest stable 3.11-compatible release, while
3.12+ resolves 2.0. Verify both mass and FASTA APIs; do not downgrade the
Python 3.12 audit environment.
