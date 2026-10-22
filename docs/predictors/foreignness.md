# Foreignness

`WEIRDO` exposes an allele-independent **category score**, using the verified model
bundled with WEIRDO 2.2.0. Install it explicitly:

```bash
pip install 'mhctools[weirdo]'
mhctools predict-table --input peptides.csv --out scored.csv \
  --predictor weirdo:weirdo_foreignness:foreignness
```

`peptides.csv` must have a `peptide` column. The table command preserves input
rows and adds the requested score column; use the Python prediction frame below
when retaining full model provenance for serialization or cache replay.

```python
from mhctools import WEIRDO

predictor = WEIRDO()
predictions = predictor.predict(["ACDEFGHI", "MTMDKSELV", "ACDEFGHI"])
frame = predictor.predict_dataframe(["ACDEFGHI"])
print(predictor.predictor_version)
```

The released `human-self-complement-v1` policy returns one minus the model's
human-category output. Higher scores mean lower predicted human-category membership
under the training distribution. This is distinct from MHC presentation and immune
recognition. No independent immune-response ranking benefit has been established;
use in ranking is explicitly experimental. See [WEIRDO's model and evaluation
contract](https://github.com/pirl-unc/weirdo/blob/main/docs/category-predictor.md).

The model accepts canonical ASCII 8–11-mers, normalizes lowercase, rejects whitespace
and ambiguous residues, preserves duplicate order and returns an empty result for an
empty batch. Longer supported peptides average the category outputs of their 8-mer
windows. Other lengths and reference hosts are unsupported. In particular, the model's
rodent category is not a mouse-specific self reference.

Every result has `kind="foreignness"`, an empty allele, policy semantics in its
measurement context, and a `predictor_version` containing the adapter/package versions,
immutable artifact, score policy and reference host. Keep the entire version in caches.

Prediction never downloads weights. `artifact_status("weirdo")` reports locally
available files; `integration_status("weirdo", check="reproduced")` verifies the actual
artifact and reference output. `models_path` accepts an identical explicit local copy.
