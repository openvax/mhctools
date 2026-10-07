# Installing optional backends

Most predictors need something downloaded or installed first. There are two
routes:

- `mhctools fetch <name>` installs what mhctools can legally and safely
  fetch. This is the route for most users; see [getting models](artifacts.md).
- **The repository's `scripts/setup_test_backends.py`** provisions the heavier
  backends (isolated Python runtimes, pinned snapshots) from a source checkout.
  It is the route CI uses and what you want if you are running the full test
  suite. If you installed from PyPI, use `fetch` and the
  [environment variables](env-vars.md) instead.

Check what is installed and whether it runs with `mhctools ls` and
`mhctools predictors`. This page covers the backends that need more than a
single `fetch`.

## Source-checkout setup

Install mhctools in editable mode with its development dependencies. The setup
below covers [CapHLA](predictors/binding.md#caphla), TULIP, [MixTCRpred](predictors/tcr.md#mixtcrpred), [PeptiVerse](predictors/peptide-pk.md#peptiverse), [PlifePred2](predictors/peptide-pk.md#plifepred2)/Pfeature,
[MixMHCpred](predictors/binding.md#mixmhcpred) 3, [MixMHC2pred](predictors/binding.md#mixmhc2pred), [PRIME](predictors/immunogenicity.md#prime), [DeepImmuno](predictors/immunogenicity.md#deepimmuno), [TLimmuno2](predictors/immunogenicity.md#tlimmuno2), [NetCleave](predictors/processing.md#netcleave), and [NetTCR](predictors/tcr.md#nettcr).
It supplements the other predictors listed in the [predictor reference](predictors/index.md);
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

The recognition, Gfeller, keras, NetTCR and [SMM](predictors/binding.md#smm-and-smm-pmbec) groups fetch separately licensed code
and weights, so they require `--accept-license`. Review the terms before using
that option, because they are not all the same kind of term:

- academic / non-commercial:
  [MixTCRpred](https://github.com/GfellerLab/MixTCRpred),
  [MixMHCpred](https://github.com/GfellerLab/MixMHCpred),
  [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred),
  [PRIME](https://github.com/GfellerLab/PRIME);
- [NetTCR academic software license](https://github.com/mnielLab/NetTCR-2.2/blob/7cead3fe6dcb539ff8e2d9121586dafca1e059c2/academic_software_license_agreement.pdf)
  for the `nettcr` group;
- Non-Profit Open Software License 3.0: the IEDB MHC-I bundle used by `smm`
  (see [SMM setup](#local-smm-and-smm-pmbec));
- no published license:
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
are reused. The `half-life` group also provisions [PeptiVerse CPP](predictors/uptake.md#peptiverse-cpp)
with its required scikit-learn 1.7.2. The wrappers verify their pinned model hashes before
inference. Re-run individual groups to repair or update the local setup.

The script writes `env/test-backends/activate.sh`, which exports the variables
the wrappers read (see [environment variables](env-vars.md)); `source` it in
your shell. An alternate setup root is supported with `--root PATH`. For how the
test suite uses it, see [testing](testing.md).

## Local SMM and SMM-PMBEC

The `smm` setup group installs a 2.4 MB subset of the official
[IEDB MHC-I 3.1.7 bundle](https://downloads.iedb.org/tools/mhci/3.1.7/README),
verifies its pinned SHA-256, and installs its Python code, allele metadata, and
model/percentile data. It does not install or execute the bundled DTU binaries.
[SMM](predictors/binding.md#smm-and-smm-pmbec) 1.0 and [SMM-PMBEC](predictors/binding.md#smm-and-smm-pmbec) 1.0 run with the current Python interpreter, including on
Apple Silicon; no extra Python dependencies are needed.
Review `LIAI_license.txt` (Non-Profit Open Software License 3.0) before
accepting the license. Upstream code and models stay outside the package.
Setup requires `curl`; cold downloads use bounded transient-error retries and
are promoted from a temporary file only after checksum verification. Model
inference needs no network access.

Install it with:

```sh
python scripts/setup_test_backends.py smm --accept-license
source env/test-backends/activate.sh
mhctools ls smm
```

For an existing configured standalone installation, set `IEDB_MHCI_EXECUTABLE`
to an executable launcher that runs `python /absolute/path/mhc_i/src/predict_binding.py "$@"`.
Alternatively put that launcher on PATH as `iedb-mhci`, or pass it with
`SMM(program_name=...)` / `--mhc-predictor-path`.
The official CLI runs locally; errors and incomplete output fail explicitly.
`mhctools ls smm` locates the launcher.

## MixMHCpred and PRIME

The `gfeller` setup group installs the official pinned predictors and a
separate Python environment, then configures `MIXMHCPRED_PYTHON`. It uses the
supported API runtime selection, with no test-only launcher. The host can
keep pandas 3; the backend needs numpy, pandas below 3, scipy, logomaker and
matplotlib. See [manual setup and runtime diagnostics](predictors/binding.md#mixmhcpred).
MAFFT remains on the inherited PATH for sequence alignment. PRIME's nested
MixMHCpred call uses the same selected interpreter.

## DeepTAP and MixTCRpred (torch sidecars)

Both run out-of-process under an interpreter that defaults to the one running
mhctools, which need not have their dependencies. 
`MIXTCRPRED_PYTHON` is provisioned by the `recognition` group above. Its
sidecar additionally runs with `PYTHONNOUSERSITE=1`, so packages installed with
`pip install --user` are not visible to it, so install into the interpreter
itself.

[DeepTAP](predictors/processing.md#deeptap) has no setup group; point `DEEPTAP_PYTHON` at any interpreter with
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

[NetTCR](predictors/tcr.md#nettcr)'s bundled inference models run under LiteRT without TensorFlow:

```sh
python -m pip install -e '.[dev,nettcr]'
python scripts/setup_test_backends.py nettcr --accept-license
source env/test-backends/activate.sh
```

The setup group fetches the pinned pan-model weights and sets `NETTCR_DIR`.
The existing `nettcr` extra supplies `ai-edge-litert` in the host interpreter.
The adapter prefers it over its compatibility TensorFlow Lite fallback, whose
deprecated interpreter emitted the warning in [#479](https://github.com/openvax/mhctools/issues/479).
Installing TensorFlow to
run another predictor is not needed to enable NetTCR. 

## DeepImmuno, TLimmuno2, and NetCleave (isolated TensorFlow)

[DeepImmuno](predictors/immunogenicity.md#deepimmuno) and [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) ship weights from the Keras 2 era. Modern TensorFlow
reaches that API through the `tf-keras` shim with `TF_USE_LEGACY_KERAS=1`, which the wrappers
set for their subprocess. [NetCleave](predictors/processing.md#netcleave) uses modern Keras in the same isolated
runtime, without that per-process setting:

```sh
python scripts/setup_test_backends.py keras --accept-license
source env/test-backends/activate.sh
```

The group pins `tensorflow==2.17.0` with the matching `tf-keras==2.17.0` and
sets `DEEPIMMUNO_PYTHON`, `TLIMMUNO2_PYTHON`, and `NETCLEAVE_PYTHON` to it.
It also installs NetCleave's scikit-learn, Biopython, and matplotlib dependencies. Pinning
the pair matters: a mismatched pair imports `tensorflow` successfully and then raises
`AttributeError` on `tensorflow.keras`.

An explicit NetCleave `python_executable` argument overrides `NETCLEAVE_PYTHON`;
without either it uses the current interpreter. The same precedence applies to
the other two wrappers with their respective environment variables. 

## Pepsickle gradient-boosted digestion models

The upstream artifact was trained with scikit-learn 0.23.2 and cannot be loaded
by current scikit-learn. Provision its pinned Python 3.8.20/Linux x86-64 runtime
with Docker (emulation is used on Apple Silicon):

```sh
python scripts/setup_test_backends.py pepsickle
source env/test-backends/activate.sh
```

The generated launcher binds to the built image ID, disables networking, and
does not mount host files. CI runs both constitutive and immunoproteasome profiles against direct upstream
inference and checks runtime provenance. Host neural [Pepsickle](predictors/processing.md#pepsickle) dependencies are
unchanged. This is implementation conformance, not held-out biological
validation; see [cleavage validation](cleavage/validation.md).

## Legacy NetMHC on Apple Silicon

[NetMHC](predictors/binding.md#netmhc) 3.4 and [NetMHCcons](predictors/binding.md#netmhccons) need Python 2 and Linux x86 executables. With Docker
running and an existing licensed netmhc-bundle installation:

```sh
export NETMHC_BUNDLE_HOME=/absolute/path/to/netmhc-bundle
python scripts/setup_test_backends.py legacy
```

This builds a pinned Linux runtime and generates test launchers. Each execution
mounts the licensed bundle and input files read-only, disables networking, and
keeps scratch files in a disposable tmpfs. No licensed tool is copied into the
image. Native Linux installations with Python 2 can continue using their existing
launchers.
