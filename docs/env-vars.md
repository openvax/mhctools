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
| `NETMHC_BUNDLE_HOME` | [NetChop](predictors/processing.md#netchop) | A bundle directory that contains your licensed DTU tools. |
| `NETCHOP_HOME` | NetChop | The directory containing `bin/netChop`. |
| `IEDB_MHCI_EXECUTABLE` | [SMM](predictors/binding.md#smm-and-smm-pmbec), [SMM-PMBEC](predictors/binding.md#smm-and-smm-pmbec) | Launcher that runs IEDB's `predict_binding.py`. |
| `MIXMHCPRED_PATH` | [MixMHCpred](predictors/binding.md#mixmhcpred) | Your MixMHCpred release. |
| `PRIME_EXECUTABLE` | [PRIME](predictors/immunogenicity.md#prime) | The PRIME executable. |
| `BIGMHC_DIR` | [BigMHC](predictors/binding.md#bigmhc) | A BigMHC clone. |
| `CAPHLA_HOME` | [CapHLA](predictors/binding.md#caphla) | A CapHLA snapshot. |
| `DEEPIMMUNO_HOME` | [DeepImmuno](predictors/immunogenicity.md#deepimmuno) | A DeepImmuno checkout. |
| `DEEPTAP_HOME` | [DeepTAP](predictors/processing.md#deeptap) | A DeepTAP checkout. |
| `ERAMER_HOME`, `ERAMER_PWM` | [ERAMER](predictors/processing.md#eramer), `eramer-step` | An ERAMER checkout, or the path to its `PWM.xlsx`. |
| `MIXTCRPRED_HOME` | [MixTCRpred](predictors/tcr.md#mixtcrpred) | A MixTCRpred checkout. |
| `NETCLEAVE_DIR` | [NetCleave](predictors/processing.md#netcleave) | A NetCleave clone. |
| `NETTCR_DIR` | [NetTCR](predictors/tcr.md#nettcr) | A NetTCR checkout. |
| `PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME` | [PeptiVerse](predictors/peptide-pk.md#peptiverse) | The pinned PeptiVerse snapshot and local ESM2 weights. |
| `PLIFEPRED2_HOME`, `PFEATURE_HOME` | [PlifePred2](predictors/peptide-pk.md#plifepred2) | The `plifepred2` package and the pinned Pfeature checkout. |
| `TLIMMUNO2_HOME` | [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) | A TLimmuno2 clone. |
| `TULIP_HOME` | [Tulip](predictors/tcr.md#tulip) | A TULIP-TCR checkout. |

## Choosing an interpreter

| Variable | Used by | Needs |
|---|---|---|
| `DEEPIMMUNO_PYTHON` | [DeepImmuno](predictors/immunogenicity.md#deepimmuno) | TensorFlow with Keras 2, or newer TensorFlow plus `tf-keras` |
| `TLIMMUNO2_PYTHON` | [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) | Same as DeepImmuno |
| `DEEPTAP_PYTHON` | [DeepTAP](predictors/processing.md#deeptap) | `torch`, `pytorch_lightning`, `numpy`, `pandas` |
| `MIXTCRPRED_PYTHON` | [MixTCRpred](predictors/tcr.md#mixtcrpred) | `torch`, `torchvision`, `pytorch_lightning`, `numpy`, `pandas`, `scipy`, `sklearn` |
| `TULIP_PYTHON` | [Tulip](predictors/tcr.md#tulip) | An isolated Python 3.11 with `torch` and `transformers==4.32.1` |
| `PEPTIVERSE_PYTHON` | [PeptiVerse](predictors/peptide-pk.md#peptiverse) | torch, `transformers==4.46.0`, xgboost, lightning |
| `PLIFEPRED2_PYTHON` | [PlifePred2](predictors/peptide-pk.md#plifepred2) | The `plifepred2` runtime |
| `NETCLEAVE_PYTHON` | [NetCleave](predictors/processing.md#netcleave) | The NetCleave runtime |
| `PEPSICKLE_PYTHON` | [Pepsickle](predictors/processing.md#pepsickle) | Interpreter for subprocess-isolated inference |
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
