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

The [osteosarc fixtures](../tests/data/osteosarc/README.md) contain source-linked
vaccine sequences and native output captures from seven real predictors.
Their parser and wrapper regression tests run offline in all four public
Python CI jobs, without downloading models or installing osteosarc:

```sh
python -m pytest tests/test_osteosarc_fixtures.py --require-all
```

The fixture README documents source provenance, experimental-label caveats,
predictor versions, and the explicit regeneration commands.

## Optional model setup

Install mhctools in editable mode with its development dependencies. The setup
below covers CapHLA, TULIP, MixTCRpred, PeptiVerse, PlifePred2/Pfeature,
MixMHCpred 3, MixMHC2pred, PRIME, DeepImmuno, and TLimmuno2. It supplements the
other predictors listed in the [predictor reference](predictors.md);
`mhctools ls` reports their availability. This is several GB of downloads,
including ESM2's 2.6 GB weights.

Prerequisites: Python 3.11, git, Perl, and MAFFT on macOS or Linux x86-64.
The official macOS MixMHC2pred binary needs Rosetta on Apple Silicon.
Linux also needs g++ to build PRIME's supplied C++ source.
CapHLA runs in the main interpreter and requires `pip install -e '.[caphla]'`.
On Linux a CPU-only torch wheel can be installed first from
`https://download.pytorch.org/whl/cpu`.

```sh
python scripts/setup_test_backends.py half-life recognition gfeller keras smm --accept-license
```

The recognition, Gfeller, keras and SMM groups fetch separately licensed code
and weights.
`--accept-license` accepts the upstream academic/non-commercial terms: review
[MixTCRpred](https://github.com/GfellerLab/MixTCRpred),
[MixMHCpred](https://github.com/GfellerLab/MixMHCpred),
[MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred), and
[PRIME](https://github.com/GfellerLab/PRIME) before using this option.
Upstream sources and models are never bundled in mhctools distributions.

The script pins model/source revisions, installs the complete MixMHC2pred
official release (including PWM assets), and installs incompatible Python
runtimes separately under ignored `env/test-backends/`. Existing environments
are reused. The wrappers verify their pinned half-life model hashes before
inference. Re-run individual groups to repair or update the local setup.

`./test.sh` automatically sources the generated `env/test-backends/activate.sh`.
For direct pytest commands, source it yourself. An alternate setup root is
supported with `--root PATH`; point `MHCTOOLS_TEST_ENV` at its `activate.sh`.
Set `MHCTOOLS_TEST_ENV=/dev/null` to run without this configuration.

## Local SMM and SMM-PMBEC

The `smm` setup group downloads the official [IEDB MHC-I 3.1.7 bundle](https://downloads.iedb.org/tools/mhci/3.1.7/README),
verifies its pinned SHA-256 (including cache hits), and installs its Python code, allele metadata,
and model/percentile data. It does not install or execute the bundled DTU
binaries. SMM 1.0 and SMM-PMBEC 1.0 run with the current Python interpreter,
including on Apple Silicon; no extra Python dependencies are needed.
Review the archive's `LIAI_license.txt` (Non-Profit Open Software License 3.0)
before accepting the license. Upstream code and models stay outside the package.
Setup requires `curl`; cold downloads use bounded transient-error retries and
are promoted from a temporary file only after checksum verification. CI caches
the immutable archive, while model inference needs no network access.

```sh
python scripts/setup_test_backends.py smm --accept-license
source env/test-backends/activate.sh
python -m pytest tests/test_smm.py tests/test_smm_integration.py --require-all
mhctools ls smm
```

For an existing configured standalone installation, set `IEDB_MHCI_EXECUTABLE`
to an executable launcher that runs `python /absolute/path/mhc_i/src/predict_binding.py "$@"`.
Alternatively put that launcher on PATH as `iedb-mhci`, or pass it with
`SMM(program_name=...)` / `--mhc-predictor-path`.
The official CLI runs locally; errors and incomplete output fail explicitly.
`mhctools ls smm` locates the launcher, while the integration tests verify
that its models actually reproduce the recorded vaccine predictions.

The public Python jobs replay the [recorded SMM outputs](../tests/data/osteosarc/smm/README.md)
offline. A dedicated integration job installs the pinned bundle and compares
real inference for both methods against every recorded peptide/allele pair.

## DeepImmuno and TLimmuno2 (Keras 2 weights)

Both ship weights from the Keras 2 era. Modern TensorFlow reaches that API
through the `tf-keras` shim with `TF_USE_LEGACY_KERAS=1`, which the wrappers
set for their subprocess, so one runtime serves both:

```sh
python scripts/setup_test_backends.py keras --accept-license
source env/test-backends/activate.sh
python -m pytest tests/test_deepimmuno.py tests/test_tlimmuno2.py --require-all
```

The group pins `tensorflow==2.17.0` with the matching `tf-keras==2.17.0` and
sets `DEEPIMMUNO_PYTHON` and `TLIMMUNO2_PYTHON` to it. Pinning the pair
matters: a mismatched pair imports `tensorflow` successfully and then raises
`AttributeError` on `tensorflow.keras`.

Without this group both wrappers fall back to the interpreter running the
tests. The end-to-end tests probe that interpreter and skip when it cannot
load Keras 2, so an unprovisioned checkout reports a skip rather than a
failure.

## Legacy NetMHC on Apple Silicon

NetMHC 3.4 and NetMHCcons need Python 2 and Linux x86 executables. With Docker
running and an existing licensed netmhc-bundle installation:

```sh
export NETMHC_BUNDLE_HOME=/absolute/path/to/netmhc-bundle
python scripts/setup_test_backends.py legacy
TEST_SH_MAX=2 ./test.sh --require-all -ra
```

This builds a pinned Linux runtime and generates test launchers. Each execution
mounts the licensed bundle and input files read-only, disables networking, and
keeps scratch files in a disposable tmpfs. No licensed tool is copied into the
image. The runtime is only for these old integration tests. Native Linux
installations with Python 2 can continue using their existing launchers.

## CI and release verification

CI runs the public suite on Python 3.9–3.12, the licensed NetMHC integration
suite, and separate real-model jobs for TULIP, CapHLA, MixTCRpred, the two
half-life predictors, the three Gfeller MHC predictors, and local SMM/SMM-PMBEC (11 CI jobs total). Each focused model
job uses `--require-all`, so a missing installation cannot silently turn it green.
The complete release run requires all installed backends and zero skips:

```sh
./lint.sh
TEST_SH_MAX=2 ./test.sh --require-all -ra
# After merging, from clean master:
PYTEST_ADDOPTS=--require-all TEST_SH_MAX=2 ./deploy.sh
```
