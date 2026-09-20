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
The public IEDB integration tests require internet access.

## Optional model setup

Install mhctools in editable mode with its development dependencies. The setup
below covers CapHLA, TULIP, MixTCRpred, PeptiVerse, PlifePred2/Pfeature,
MixMHCpred 3, MixMHC2pred, and PRIME. It supplements the other predictors
listed in the [installation guide](../README.md); `mhctools ls` reports their
availability. This is several GB of downloads, including ESM2's 2.6 GB weights.

Prerequisites: Python 3.11, git, Perl, and MAFFT on macOS or Linux x86-64.
The official macOS MixMHC2pred binary needs Rosetta on Apple Silicon.
Linux also needs g++ to build PRIME's supplied C++ source.
CapHLA runs in the main interpreter and requires `pip install -e '.[caphla]'`.
On Linux a CPU-only torch wheel can be installed first from
`https://download.pytorch.org/whl/cpu`.

```sh
python scripts/setup_test_backends.py half-life recognition gfeller --accept-license
```

The recognition and Gfeller groups fetch separately licensed code and weights.
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
half-life predictors, and the three Gfeller MHC predictors. Each focused model
job uses `--require-all`, so a missing installation cannot silently turn it green.
The complete release run requires all installed backends and zero skips:

```sh
./lint.sh
TEST_SH_MAX=2 ./test.sh --require-all -ra
# After merging, from clean master:
PYTEST_ADDOPTS=--require-all TEST_SH_MAX=2 ./deploy.sh
```
