# Troubleshooting

## `mhctools ls` says `ready` but prediction fails

`ls` is an inventory: `ready` means the required path was located, not that the
tool works. Run `mhctools predictors` for a capability report. Its `LOCATED`,
`RUNNABLE` and `REPRODUCED` columns are independent, and `not checked` is never
promoted to success. See [inventory is not capability](artifacts.md#inventory-is-not-capability).

```sh
mhctools predictors netmhcpan calis --check reproduced --json
```

## `UnsupportedAllele`

The predictor does not support that allele, or the allele is the wrong class
(a DR allele passed to a class I predictor). For the NetMHC family the error
names the command that lists supported alleles, for example
`netMHCpan-4.2 -listMHC`. See [allele names](alleles.md).

## A prediction table is shorter than my input, or `ValueError` about missing pairs

`predict()` returns one result per input peptide in input order, and raises if a
backend omits a requested peptide-allele pair rather than returning a shorter
list. If you scan proteins, remember the window is narrow by default; see
[peptide lengths](predictors/index.md#peptide-lengths).

## A cleavage score is exactly 0.0

`Pepsickle` and `NetChop` score the C-terminal bond using the residues that
follow the peptide. Without `c_flanks=` there is no downstream context, so the
score is 0.0. Pass flanks or use `predict_proteins()`.

## `env: python2: No such file or directory`

NetMHC 3.4 and NetMHCcons need Python 2 and a Linux x86 executable. On macOS or
ARM Linux use the pinned container setup in [legacy NetMHC on Apple
Silicon](backends.md#legacy-netmhc-on-apple-silicon).

## NetChop fails with a Docker error

NetChop 3.1 ships 32-bit x86 Linux binaries, so on macOS and ARM Linux it runs
in a network-disabled Docker container. Docker must be running and the pinned
image must be pulled once beforehand; see [NetChop](predictors/processing.md#netchop).
Use `NetChop(execution="native")` on x86 Linux.

## Torch or TensorFlow errors from a sidecar predictor

DeepTAP, DeepImmuno, TLimmuno2, MixTCRpred, Tulip, PeptiVerse and PlifePred2
run in a separate interpreter chosen by a `*_PYTHON` variable; see
[environment variables](env-vars.md). Common causes:

- the interpreter lacks the dependency (the error names the missing module);
- `tensorflow` and `tf-keras` versions differ, which imports `tensorflow` fine
  and then raises `AttributeError` on `tensorflow.keras`; install a matching
  pair;
- the package was installed with `pip install --user`, which sidecars cannot
  see.

## Crash from duplicate OpenMP runtimes on macOS (Pepsickle)

Construct it with `Pepsickle(isolate_subprocess=True)` to run inference in a
short-lived subprocess.

## TLimmuno2 is very slow

Its percentile rank is computed against about 90,000 background peptides per
distinct allele, so a call costs about a minute per allele regardless of how
many peptides it scores. Batch peptides by allele.

## MixMHCpred refuses my output directory

`predict_detailed(..., output_dir=...)` and `predict_allele_sequences(...)`
require a path that does not exist yet, because MixMHCpred deletes and recreates
its output directory.

## Tests skip, or `--require-all` fails

Prediction tests that need an installed tool skip when it is absent. The release
gate `--require-all` turns skips into failures. See [testing](testing.md) and
[installing optional backends](backends.md).
