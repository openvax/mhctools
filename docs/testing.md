# Running all tests

`./test.sh` runs the complete suite, with a memory-aware number of pytest
workers. Some integration tests require separately installed predictors and
models. Use the release gate to ensure every collected test actually passes:

```sh
TEST_SH_MAX=2 ./test.sh --require-all -ra
```

`--require-all` fails on skips, expected failures, and unexpected passes,
including module-level skips and parallel workers. Ordinary development runs
can still use availability-based skips when optional predictors are absent.
Prediction tests do not call IEDB or other public prediction services.
Internet access is needed only when installing model assets and runtimes.

## Recorded vaccine fixtures

The [osteosarc fixtures](https://github.com/openvax/mhctools/blob/master/tests/data/osteosarc/README.md) contain source-linked
vaccine sequences and native output captures from seven real predictors.
Their parser and wrapper regression tests run offline in all four public
Python CI jobs, without downloading models or installing osteosarc:

```sh
python -m pytest tests/test_osteosarc_fixtures.py --require-all
```

The fixture README documents source provenance, experimental-label caveats,
predictor versions, and the explicit regeneration commands.

## Optional backends in tests

`./test.sh` sources `env/test-backends/activate.sh` when it exists, so the
isolated runtimes set up by `scripts/setup_test_backends.py` are found
automatically. For direct `pytest` commands, source it yourself. Point
`MHCTOOLS_TEST_ENV` at a different `activate.sh`, or set
`MHCTOOLS_TEST_ENV=/dev/null` to run without this configuration. Installing the
backends is covered in [installing optional backends](backends.md); this page
is about running the tests that exercise them.

Each backend has focused tests. They probe the runtime and skip when it is
missing, so a skip means an unprovisioned runtime, not a broken wrapper. A
`*_PYTHON` path that does not exist raises instead of skipping, because that is a
misconfiguration.

| Backend | Setup group | Command |
|---|---|---|
| [MixMHCpred](predictors/binding.md#mixmhcpred), [PRIME](predictors/immunogenicity.md#prime) | `gfeller` | `pytest tests/test_mixmhcpred_runtime_integration.py tests/test_mixmhcpred.py tests/test_prime.py --require-all` |
| [PhageScout](cleavage/phagescout.md) full DESeq2 tables | `phagescout` | `pytest tests/test_phagescout_full_profiles.py --require-all` |
| [Pepsickle](predictors/processing.md#pepsickle) gradient-boosted | `pepsickle` | `pytest tests/test_pepsickle_legacy.py tests/test_pepsickle_runtime.py --require-all` |
| [SMM](predictors/binding.md#smm-and-smm-pmbec), [SMM-PMBEC](predictors/binding.md#smm-and-smm-pmbec) | `smm` | `python -m pytest tests/test_smm.py tests/test_smm_integration.py --require-all` |
| [NetTCR](predictors/tcr.md#nettcr) (LiteRT) | `nettcr` | `python -m pytest tests/test_nettcr.py --require-all -W error` |
| [DeepImmuno](predictors/immunogenicity.md#deepimmuno), [TLimmuno2](predictors/immunogenicity.md#tlimmuno2), [NetCleave](predictors/processing.md#netcleave) | `keras` | `python -m pytest tests/test_deepimmuno.py tests/test_tlimmuno2.py tests/test_netcleave.py --require-all -W error` |
| Legacy [NetMHC](predictors/binding.md#netmhc) | `legacy` | `TEST_SH_MAX=2 ./test.sh --require-all -ra` |

Notes on what CI asserts:

- **MixMHCpred and PRIME:** the host uses pandas 3 and the isolated backend
  uses pandas below 3. Unmodified native launchers provide reference scores
  for multi-allele rows, duplicate peptides, sequence alignment, motifs/HTML,
  CLI output and nested PRIME execution. The host interpreter, pandas
  version and environment are checked before and after inference.
- **NetTCR** runs with TensorFlow absent from the host and without skips or
  warnings.
- **DeepImmuno, TLimmuno2 and NetCleave** run through the isolated interpreter
  with the host free of TensorFlow. DeepImmuno and TLimmuno2 tests probe Keras-2
  availability; NetCleave tests require their runtime whenever its model assets
  are installed.
- **SMM:** the public Python jobs replay the
  [recorded SMM outputs](https://github.com/openvax/mhctools/blob/master/tests/data/osteosarc/smm/README.md) offline, and a
  dedicated integration job installs the pinned bundle and compares real
  inference for both methods against every recorded peptide/allele pair.
- **Pepsickle:** the CI job runs the constitutive and immunoproteasome profiles
  against direct upstream inference, checks runtime provenance, and exercises
  batch save/reload. This is implementation conformance, not held-out
  biological validation; see [cleavage validation](cleavage/validation.md).
- The legacy NetMHC runtime exists only for the old NetMHC 3.4 and [NetMHCcons](predictors/binding.md#netmhccons)
  integration tests.

Maintaining the pinned SMM subset is covered in [maintaining the IEDB SMM
subset](dev/smm-subset.md).

## Documentation checks

The predictor matrix is generated and tested:

```sh
python scripts/predictor_matrix.py            # rewrite docs/predictor-matrix.md
python scripts/predictor_matrix.py --check    # exit 1 if it is stale
python -m pytest tests/test_docs_predictor_matrix.py
```

The test fails when an exported predictor class or command-line name has no row
in the matrix. To build the site locally with broken links treated as errors:

```sh
python -m pip install -e '.[docs]'
mkdocs build --strict
```

The header version is read from `mhctools/__init__.py` during each docs build.
It describes the checked-out package version and is independent of GitHub
release metadata.

## CI and release verification

CI discovers all tests on Python 3.9–3.12 and runs every test without the
`requires_external_tool` marker. New offline test files join that matrix
automatically. Mark only tests that require separately installed predictor
binaries, weights, or isolated upstream runtimes; keep parser, mocked-adapter,
and validation tests in mixed modules unmarked. Run the same selection locally:

```sh
python -m pytest tests/ -m "not requires_external_tool" --strict-markers -ra
```

CI also runs the licensed [NetMHC](predictors/binding.md#netmhc) integration
suite, and separate real-model jobs for TULIP, [CapHLA](predictors/binding.md#caphla), [MixTCRpred](predictors/tcr.md#mixtcrpred), the two
half-life predictors, the three Gfeller MHC predictors,
[DeepImmuno](predictors/immunogenicity.md#deepimmuno)/[TLimmuno2](predictors/immunogenicity.md#tlimmuno2)/[NetCleave](predictors/processing.md#netcleave), [NetTCR](predictors/tcr.md#nettcr), and local [SMM](predictors/binding.md#smm-and-smm-pmbec)/[SMM-PMBEC](predictors/binding.md#smm-and-smm-pmbec). Each focused
model job uses `--require-all`, so a missing installation cannot silently turn it green.
The complete release run requires all installed backends and zero skips:

```sh
./lint.sh
TEST_SH_MAX=2 ./test.sh --require-all -ra
# After merging, from clean master:
PYTEST_ADDOPTS=--require-all TEST_SH_MAX=2 ./deploy.sh
```

## CleaveNet source conformance

Provision with `python scripts/setup_test_backends.py cleavenet --python python3.11`,
source the generated activation file, then run
`python -m pytest tests/test_cleavenet.py --require-all`. The separate CI job runs
the real five-model ensemble; the quality matrix runs validation and mocked
output-contract regressions without TensorFlow. The recorded fixture comes
from official source inference, not an independent experimental benchmark.
