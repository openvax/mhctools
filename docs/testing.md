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
MixMHCpred 3, MixMHC2pred, PRIME, DeepImmuno, TLimmuno2, NetCleave, and NetTCR.
It supplements the other predictors listed in the [predictor reference](predictors.md);
`mhctools ls` reports their availability. This is several GB of downloads,
including ESM2's 2.6 GB weights.

Prerequisites: Python 3.11, git, Perl, and MAFFT on macOS or Linux x86-64.
The official macOS MixMHC2pred binary needs Rosetta on Apple Silicon.
Linux also needs g++ to build PRIME's supplied C++ source.
CapHLA runs in the main interpreter and requires `pip install -e '.[caphla]'`.
NetTCR also runs in the main interpreter and requires `pip install -e '.[nettcr]'`
for LiteRT. Neither extra requires TensorFlow; the `keras` group below keeps
TensorFlow in a separate environment for the predictors that use it.
On Linux a CPU-only torch wheel can be installed first from
`https://download.pytorch.org/whl/cpu`.

```sh
python -m pip install -e '.[dev,caphla,nettcr]'
python scripts/setup_test_backends.py half-life recognition gfeller keras nettcr smm --accept-license
```

The recognition, Gfeller, keras, NetTCR and SMM groups fetch separately licensed code
and weights, so they require `--accept-license`. Review the terms before using
that option, because they are not all the same kind of term:

- academic / non-commercial —
  [MixTCRpred](https://github.com/GfellerLab/MixTCRpred),
  [MixMHCpred](https://github.com/GfellerLab/MixMHCpred),
  [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred),
  [PRIME](https://github.com/GfellerLab/PRIME);
- [NetTCR academic software license](https://github.com/mnielLab/NetTCR-2.2/blob/7cead3fe6dcb539ff8e2d9121586dafca1e059c2/academic_software_license_agreement.pdf)
  for the `nettcr` group;
- Non-Profit Open Software License 3.0 — the IEDB MHC-I bundle used by `smm`
  (see [below](#local-smm-and-smm-pmbec));
- **no published license** —
  [TLimmuno2](https://github.com/XSLiuLab/TLimmuno2) and
  [NetCleave](https://github.com/BSC-CNS-EAPM/NetCleave), which is why the `keras`
  group is gated. Here `--accept-license` records that you have
  confirmed your own use is authorized; it cannot accept terms upstream never
  stated, and it grants no rights mhctools does not have.

DeepImmuno, also in the `keras` group, is MIT and needs no
acceptance of its own. Upstream sources and models are never bundled in
mhctools distributions.

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

The `smm` setup group installs a 2.4 MB subset of the official
[IEDB MHC-I 3.1.7 bundle](https://downloads.iedb.org/tools/mhci/3.1.7/README),
verifies its pinned SHA-256, and installs its Python code, allele metadata, and
model/percentile data. It does not install or execute the bundled DTU binaries.
SMM 1.0 and SMM-PMBEC 1.0 run with the current Python interpreter, including on
Apple Silicon; no extra Python dependencies are needed.
Review `LIAI_license.txt` (Non-Profit Open Software License 3.0) before
accepting the license. Upstream code and models stay outside the package.
Setup requires `curl`; cold downloads use bounded transient-error retries and
are promoted from a temporary file only after checksum verification. Model
inference needs no network access.

The subset is served from this repo's
[`iedb-smm-subset-3.1.7`](https://github.com/openvax/mhctools/releases/tag/iedb-smm-subset-3.1.7)
release rather than fetched from `downloads.iedb.org`. That host became
unreachable from GitHub runners on 2026-09-28 (`curl: (28) Connection timeout`,
four attempts) with the Actions cache evicted, which blocked merges on PRs that
had nothing to do with SMM. The release unpacks to 1031 MB across 38,236
members, dominated by bundled DTU executables under `method/` (netmhc-4.0 is
210 MB, netmhc-3.4 192 MB, netmhcpan-4.1 114 MB) plus 192 MB of per-method
training data under `data/`. mhctools does not run those executables from this
bundle; the netMHC family is wrapped through its own licensed distribution.
The paths SMM needs come to 9.4 MB.

The subset is a verbatim copy of `LIAI_license.txt`, `Copenhagen_license.txt`,
the upstream `README`, `src/`, `method/allele-info/`,
`method/iedbtools-utilities/` and the `smm/`, `smmpmbec/` and `consensus/`
training data. `consensus/` is required despite
the consensus method being unused: percentile ranks read
`distribution_consensus_bin.cpickle` from it. `scripts/build_iedb_smm_subset.py`
derives it from the official archive, verifying that archive's own SHA-256 and
normalizing entry order and metadata. Rebuilding with the same compression
runtime reproduces the checksum; different zlib versions may produce different
compressed bytes even when the extracted files match upstream exactly:

```sh
curl -o env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    https://downloads.iedb.org/tools/mhci/3.1.7/IEDB_MHC_I-3.1.7.tar.gz
python scripts/build_iedb_smm_subset.py \
    --archive env/test-backends/IEDB_MHC_I-3.1.7.tar.gz \
    --output dist/IEDB_MHC_I-3.1.7-smm-subset.tar.gz
```

Both paths are already git-ignored. The build exits with an error if the output
digest differs from `EXPECTED_SUBSET_SHA256`; do not publish that output. It
also fails if any allowlisted prefix contains no regular files or the archive
grows a member type it does not copy.

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

## DeepTAP and MixTCRpred (torch sidecars)

Both run out-of-process under an interpreter that defaults to the one running
mhctools, which need not have their dependencies. The end-to-end tests probe
that interpreter and skip when it cannot import them, so a skip here means the
runtime is missing rather than the wrapper being broken.

`MIXTCRPRED_PYTHON` is provisioned by the `recognition` group above. Its
sidecar additionally runs with `PYTHONNOUSERSITE=1`, so packages installed with
`pip install --user` are not visible to it — install into the interpreter
itself.

DeepTAP has no setup group; point `DEEPTAP_PYTHON` at any interpreter with
`torch` and `pytorch-lightning`, and `DEEPTAP_HOME` at a checkout
(`mhctools fetch deeptap`). The current interpreter is used when
`DEEPTAP_PYTHON` is unset, which is enough when it already has torch.

| Variable | Needs |
|---|---|
| `DEEPTAP_PYTHON` | `torch`, `pytorch_lightning`, `numpy`, `pandas` |
| `MIXTCRPRED_PYTHON` | `torch`, `torchvision`, `pytorch_lightning`, `numpy`, `pandas`, `scipy`, `sklearn` |

A `*_PYTHON` path that does not exist raises instead of skipping: that is a
misconfiguration to fix, not a backend to step over.

## NetTCR with LiteRT

NetTCR's bundled inference models run under LiteRT without TensorFlow:

```sh
python -m pip install -e '.[dev,nettcr]'
python scripts/setup_test_backends.py nettcr --accept-license
source env/test-backends/activate.sh
python -m pytest tests/test_nettcr.py --require-all -W error
```

The setup group fetches the pinned pan-model weights and sets `NETTCR_DIR`.
The existing `nettcr` extra supplies `ai-edge-litert` in the host interpreter.
The adapter prefers it over its compatibility TensorFlow Lite fallback, whose
deprecated interpreter emitted the warning in [#479](https://github.com/openvax/mhctools/issues/479).
Installing TensorFlow to
run another predictor is not needed to enable NetTCR. CI checks that TensorFlow
is absent from the host and runs the real prediction regressions without skips
or warnings.

## DeepImmuno, TLimmuno2, and NetCleave (isolated TensorFlow)

DeepImmuno and TLimmuno2 ship weights from the Keras 2 era. Modern TensorFlow
reaches that API through the `tf-keras` shim with `TF_USE_LEGACY_KERAS=1`, which the wrappers
set for their subprocess. NetCleave uses modern Keras in the same isolated
runtime, without that per-process setting:

```sh
python scripts/setup_test_backends.py keras --accept-license
source env/test-backends/activate.sh
python -m pytest tests/test_deepimmuno.py tests/test_tlimmuno2.py tests/test_netcleave.py --require-all -W error
```

The group pins `tensorflow==2.17.0` with the matching `tf-keras==2.17.0` and
sets `DEEPIMMUNO_PYTHON`, `TLIMMUNO2_PYTHON`, and `NETCLEAVE_PYTHON` to it.
It also installs NetCleave's scikit-learn, Biopython, and matplotlib dependencies. Pinning
the pair matters: a mismatched pair imports `tensorflow` successfully and then raises
`AttributeError` on `tensorflow.keras`.

An explicit NetCleave `python_executable` argument overrides `NETCLEAVE_PYTHON`;
without either it uses the current interpreter. The same precedence applies to
the other two wrappers with their respective environment variables. DeepImmuno
and TLimmuno2 tests probe Keras-2 availability; NetCleave tests require its
runtime whenever its model assets are installed. Use the setup group to make
both models and runtimes available. CI checks that the host has no TensorFlow
and executes all three predictors' regressions through the isolated interpreter.

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
half-life predictors, the three Gfeller MHC predictors,
DeepImmuno/TLimmuno2/NetCleave, NetTCR, and local SMM/SMM-PMBEC. Each focused
model job uses `--require-all`, so a missing installation cannot silently turn it green.
The complete release run requires all installed backends and zero skips:

```sh
./lint.sh
TEST_SH_MAX=2 ./test.sh --require-all -ra
# After merging, from clean master:
PYTEST_ADDOPTS=--require-all TEST_SH_MAX=2 ./deploy.sh
```
