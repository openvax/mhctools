# Environment variables

Most wrappers find their tool from an explicit constructor argument first, then
an environment variable, then conventional locations such as `~/<Tool>`, then a
`mhctools fetch` snapshot. The exact order varies a little by wrapper, and
`mhctools ls` shows which one was found. "Interpreter" variables name a separate Python that has the model's
dependencies, so heavy stacks (torch, TensorFlow) stay out of the environment
mhctools runs in.

## Locating tools and models

| Variable | Used by | Meaning |
|---|---|---|
| `MHCTOOLS_DATA_DIR` | `fetch`, `ls`, every fetched snapshot | Root for mhctools-managed snapshots. Default is the platform user-data directory. |
| `NETMHC_BUNDLE_HOME` | NetChop | A bundle directory that contains your licensed DTU tools. |
| `NETCHOP_HOME` | NetChop | The directory containing `bin/netChop`. |
| `IEDB_MHCI_EXECUTABLE` | SMM, SMM-PMBEC | Launcher that runs IEDB's `predict_binding.py`. |
| `MIXMHCPRED_PATH` | MixMHCpred | Your MixMHCpred release. |
| `PRIME_EXECUTABLE` | PRIME | The PRIME executable. |
| `BIGMHC_DIR` | BigMHC | A BigMHC clone. |
| `CAPHLA_HOME` | CapHLA | A CapHLA snapshot. |
| `DEEPIMMUNO_HOME` | DeepImmuno | A DeepImmuno checkout. |
| `DEEPTAP_HOME` | DeepTAP | A DeepTAP checkout. |
| `ERAMER_HOME`, `ERAMER_PWM` | ERAMER, `eramer-step` | An ERAMER checkout, or the path to its `PWM.xlsx`. |
| `MIXTCRPRED_HOME` | MixTCRpred | A MixTCRpred checkout. |
| `NETCLEAVE_DIR` | NetCleave | A NetCleave clone. |
| `NETTCR_DIR` | NetTCR | A NetTCR checkout. |
| `PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME` | PeptiVerse | The pinned PeptiVerse snapshot and local ESM2 weights. |
| `PLIFEPRED2_HOME`, `PFEATURE_HOME` | PlifePred2 | The `plifepred2` package and the pinned Pfeature checkout. |
| `TLIMMUNO2_HOME` | TLimmuno2 | A TLimmuno2 clone. |
| `TULIP_HOME` | Tulip | A TULIP-TCR checkout. |

## Choosing an interpreter

| Variable | Used by | Needs |
|---|---|---|
| `DEEPIMMUNO_PYTHON` | DeepImmuno | TensorFlow with Keras 2, or newer TensorFlow plus `tf-keras` |
| `TLIMMUNO2_PYTHON` | TLimmuno2 | Same as DeepImmuno |
| `DEEPTAP_PYTHON` | DeepTAP | `torch`, `pytorch_lightning`, `numpy`, `pandas` |
| `MIXTCRPRED_PYTHON` | MixTCRpred | `torch`, `torchvision`, `pytorch_lightning`, `numpy`, `pandas`, `scipy`, `sklearn` |
| `TULIP_PYTHON` | Tulip | An isolated Python 3.11 with `torch` and `transformers==4.32.1` |
| `PEPTIVERSE_PYTHON` | PeptiVerse | torch, `transformers==4.46.0`, xgboost, lightning |
| `PLIFEPRED2_PYTHON` | PlifePred2 | The `plifepred2` runtime |
| `NETCLEAVE_PYTHON` | NetCleave | The NetCleave runtime |
| `PEPSICKLE_PYTHON` | Pepsickle | Interpreter for subprocess-isolated inference |
| `PEPSICKLE_GB_PYTHON` | Pepsickle gradient-boosted models | A Python 3.8 runtime with scikit-learn 0.23.2; see [cleavage models](cleavage/models.md#proteasome-models) |

A `*_PYTHON` path that does not exist raises an error rather than falling back.
Sidecar interpreters run with `PYTHONNOUSERSITE=1`, so packages installed with
`pip install --user` are invisible to them; install into the interpreter itself.

`TF_USE_LEGACY_KERAS=1` is set by the DeepImmuno and TLimmuno2 wrappers for
their subprocess; you do not need to set it yourself.

## Testing

| Variable | Meaning |
|---|---|
| `MHCTOOLS_TEST_ENV` | Path to an `activate.sh` for `./test.sh`; `/dev/null` disables it. See [testing](testing.md). |
| `TEST_SH_MAX`, `TEST_SH_MIN`, `PER_WORKER_GB` | Worker limits for `./test.sh`. |
