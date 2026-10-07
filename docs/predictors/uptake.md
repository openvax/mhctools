# Cell penetration classification

These models classify cell-penetrating peptides (CPPs). They do not estimate
the fraction entering cells, uptake by antigen-presenting cells, endosomal
escape, intracellular localization or productive antigen presentation.
See [exposure results](../exposure-results.md) for the distinct measurement kinds.

## PeptiVerse CPP

`PeptiVerseCPP` wraps the released canonical-sequence SVM from
[PeptiVerse](https://doi.org/10.1038/s41467-026-74167-w), using the same pinned
source and ESM2 embedding snapshot as the [half-life adapter](peptide-pk.md#peptiverse).
Despite its upstream directory name `svm_gpu_wt`, this checkpoint is a CPU
`sklearn.svm.SVC`. No GPU or retraining is needed. Prediction runs offline in
`PEPTIVERSE_PYTHON`; the host environment does not need its ML dependencies.

Provision the isolated runtime and public assets:

```sh
python3.11 -m venv /models/peptiverse-env
/models/peptiverse-env/bin/python -m pip install \
  'torch>=2.1' 'transformers==4.46.0' 'lightning==2.5.5' \
  xgboost 'scikit-learn==1.7.2' joblib mapie pandas SmilesPE rdkit
/models/peptiverse-env/bin/huggingface-cli download ChatterjeeLab/PeptiVerse \
  --revision 8cf0b21dae356278ae96b414a088e4360357d16c \
  --include inference.py best_models.txt 'tokenizer/*.py' \
    training_classifiers/permeability_penetrance/svm_gpu_wt/best_model.joblib \
  --local-dir /models/PeptiVerse
/models/peptiverse-env/bin/huggingface-cli download facebook/esm2_t33_650M_UR50D \
  --revision 08e4846e537177426273712802403f7ba8261b6c \
  --include config.json tokenizer_config.json special_tokens_map.json vocab.txt model.safetensors \
  --local-dir /models/esm2_t33_650M_UR50D
export PEPTIVERSE_HOME=/models/PeptiVerse
export PEPTIVERSE_ESM_HOME=/models/esm2_t33_650M_UR50D
export PEPTIVERSE_PYTHON=/models/peptiverse-env/bin/python
```

For development, `python scripts/setup_test_backends.py half-life` provisions
both PeptiVerse endpoints and PlifePred2. CPP-only installations need no
half-life weights. The serialized SVC requires **scikit-learn 1.7.2**; another
version fails before source import or checkpoint loading. The adapter verifies
the exact source, model, original threshold manifest and ESM2 resources. Only
the sequence CPP head is loaded; unused SMILES models are disabled.

```python
from mhctools import PeptideInput, PeptiVerseCPP

predictor = PeptiVerseCPP(device="cpu", uncertainty=True)
results = predictor.predict([
    PeptideInput("GRKKRRQRRRPQ", occurrence_id="sample-1"),
    "SIINFEKL",
])
prediction = results[0].preds[0]
prediction.kind                              # cpp_classification
prediction.score                             # native positive-class P(CPP), 0-1
prediction.value                             # None: no physical uptake quantity
prediction.measurement_context.class_label   # CPP or non-CPP
predictor.last_qc                             # threshold, entropy and applicability
predictor.artifact_inventory.to_dict()        # exact hashes and runtime contract
```

On the command line:

```sh
mhctools predict-table --input peptides.csv --out cpp.csv \
  --predictor peptiverse-cpp:cpp_probability:score
```

Scores at or above the native threshold **0.5493** are classified as CPP.
The score is the SVC's positive-class output, not a calibrated percentage of
peptide delivered. No cell type, matrix, delivery route or exposure time is
inferred. With `uncertainty=True`, `last_qc.predictive_entropy_nats` contains
the native single-model binary predictive entropy. It is not ensemble
variation, an accuracy estimate or a physical error bar; it stays outside
the prediction's score/value fields.

The pinned source metadata contains 1,859 training and 465 validation rows;
training sequences span **3-61 residues**. Those limits describe the observed
training domain, not absolute biological cutoffs. Each result flags lengths
outside it. The tokenizer's independent computational capacity is 1,020
residues; longer inputs fail rather than being truncated. Accuracy on long
vaccine peptides remains unestablished. The five fixed native controls verify
source conformance, not external biological accuracy or freedom from training
overlap. Broader uptake and benchmark work remains in
[#302](https://github.com/openvax/mhctools/issues/302),
[#303](https://github.com/openvax/mhctools/issues/303) and
[#291](https://github.com/openvax/mhctools/issues/291).

Only canonical L-peptides with free termini are accepted. Exact chemical form,
context and occurrences are preserved; unsupported modifications are rejected.
Use `on_unsupported="record"` to keep unavailable entries in a mixed batch.
Upstream declares Apache-2.0 on the model card and MIT in the README; the
separate UI Space is not used. The SVC uses joblib pickle serialization, so
load only source/model snapshots you trust. A matching hash establishes
identity, not serialization safety.
