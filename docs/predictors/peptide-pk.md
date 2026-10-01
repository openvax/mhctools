# Peptide half-life predictors

Predictors of how long a free peptide survives. They emit `peptide_half_life`,
read with `result.peptide_half_life`; the sample matrix (serum, plasma, whole
blood) lives in the measurement context. For the wider set of pharmacokinetic,
uptake and tissue-exposure kinds, see [peptide PK, uptake, and tissue
exposure](../exposure-results.md).

| Predictor | What it estimates | Needs |
|---|---|---|
| [PeptiVerse](#peptiverse) | Parent-peptide half-life in human serum | pinned PeptiVerse + ESM2 snapshots (`PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME`) + a torch/transformers Python |
| [PlifePred2](#plifepred2) | Parent-peptide half-life, **matrix unknown** | `plifepred2==1.0` (`PLIFEPRED2_HOME`) + pinned Pfeature (`PFEATURE_HOME`) |

`peptide_half_life` records how long the parent peptide persists, in hours. Its
context distinguishes a defined solution, serum/plasma/whole blood, cellular
compartments, and systemic in-vivo PK without multiplying kind strings.

It is deliberately a separate kind from `pMHC_stability`, which is the
dissociation half-life of an assembled peptide-MHC complex — a different
molecule in a different assay — and from the cleavage kinds, which are
site-resolved and intracellular. `MeasurementContext` preserves the matrix,
compartment, analyte, and systemic scope when they are known.

## PeptiVerse

`PeptiVerse` wraps one endpoint of the upstream multi-property platform. Its
dependencies (torch, `transformers==4.46.0`, xgboost, lightning, and ESM2) stay
out of the mhctools environment: inference runs offline in a subprocess under
`PEPTIVERSE_PYTHON`. Provision the exact snapshots before prediction:

```bash
git clone https://huggingface.co/ChatterjeeLab/PeptiVerse
git -C PeptiVerse checkout 8cf0b21dae356278ae96b414a088e4360357d16c
huggingface-cli download facebook/esm2_t33_650M_UR50D \
  --revision 08e4846e537177426273712802403f7ba8261b6c \
  --include config.json tokenizer_config.json special_tokens_map.json vocab.txt model.safetensors \
  --local-dir /models/esm2_t33_650M_UR50D
export PEPTIVERSE_HOME="$PWD/PeptiVerse"
export PEPTIVERSE_ESM_HOME=/models/esm2_t33_650M_UR50D
```

```python
from mhctools import PeptideContext, PeptideInput, PeptiVerse

predictor = PeptiVerse(device="cpu")       # resolves PEPTIVERSE_HOME / ~/PeptiVerse
exact_input = PeptideInput(
    "SIINFEKL",
    occurrence_id="sample-1:occurrence-2",
    context=PeptideContext(matrix="serum", assay_species="Homo sapiens"),
)
results = predictor.predict([exact_input, "KLGGALQAK"])
results[0].peptide_half_life.value         # hours, higher = longer-lived
results[0].serum_half_life.peptide_input   # exact chemistry + context
results[0].serum_half_life.cache_key       # input + assets + settings
predictor.artifact_inventory.to_dict()     # exact files, hashes, capability
```

Sequence input only. Upstream's SMILES models return a number that is *not* on
the hours scale — the `expm1` inverse transform is applied only to the sequence
model. Although mhctools records exact chemical form, this adapter does not
consume it, so terminal modifications, attachments, and non-standard residues
are rejected rather than scored as their unmodified sequence. Pass
`on_unsupported="record"` to retain unsupported entries in a mixed batch.

> ⚠️ The sequence half-life model was fit on **130 examples** and evaluated by
> cross-validation only, from a preprint, with no external test set and no
> evaluation on long vaccine peptides. Upstream declares Apache-2.0 on its model
> card and MIT in its README. mhctools verifies the exact inference source,
> model, calibration, ESM2 weights, configuration, and tokenizer files before
> launch. The PeptiVerse checkpoint and calibration still use unsafe pickle
> serialization; matching a checksum establishes identity, not safety. Use
> only snapshots you trust.

## PlifePred2

> ⚠️ **This endpoint's semantics are not established.** PlifePred2 ships no
> publication, no training data and no target definition, so its units,
> transform, species and assay matrix are all inferred from the artifacts. By
> default the wrapper reports only the model's native output and claims no
> duration at all.

```python
from mhctools import PlifePred2

predictor = PlifePred2()                       # PLIFEPRED2_HOME + PFEATURE_HOME
results = predictor.predict(["SIINFEKLGGALQAKKY"])
results[0].peptide_half_life.score             # native output, higher = longer-lived
results[0].peptide_half_life.value             # None by default
predictor.artifact_inventory.to_dict()         # exact files, hashes, capability
predictor.last_qc["log10_seconds"]             # the same value, named

# Opt in to a duration, accepting the inference below:
opted_in = PlifePred2(assume_log10_seconds=True)
opted_in.predict(["SIINFEKLGGALQAKKY"])[0].peptide_half_life.value # hours
```

**What is known.** Both shipped models are `RandomForestRegressor`, verified by
loading them. So the output is not a class probability — upstream's docs
("Halflife … Predicted probability") and its CLI's `predict_proba` branch are
both wrong, and the branch is dead code. Being monotone in half-life, the
output ranks correctly whatever the transform turns out to be.

**What is inferred.** `log10(half-life in seconds)` is the strongest reading:
inverting the forests' extreme leaf values under it gives round durations —
exactly 7.000 days for the natural model, 95.0 days for the modified one — to
about seven significant figures, where log2 and ln both invert the whole
training range to a few seconds up to a couple of minutes. The minimum also
lands on 20.2 s, matching the 20-second floor in the lineage paper. That last
point is corroboration rather than proof: the same forests hold targets past
that paper's 24-hour ceiling, so PlifePred2 was trained on a different dataset
and the old filter cannot establish the new target. Note also that the lineage
paper states log2, not log10.

**What is not established.** The species and assay matrix. The result therefore
uses generic `peptide_half_life` with `matrix=None`; do not report it as a
measured whole-blood property or treat it as interchangeable with PeptiVerse's
human-serum endpoint.

Natural peptides only, 12–100 residues. Upstream's CLI silently drops
out-of-range and modified sequences into an `eliminated_sequences.csv` and
returns a shorter result set; mhctools rejects them instead so a caller never
gets a quietly truncated answer.

The Linux-only `pfeature_comp` binary that `plifepred2` bundles is **not** used
— it is a PyInstaller freeze of Pfeature's `pfeature_comp.py`, and that plain
Python source computes the same descriptor on any platform. Both upstreams are
GPLv3, so neither is vendored and neither is imported into the mhctools
interpreter.

> ⚠️ PlifePred2 cites no publication of its own, so its training set is
> unverified beyond what the artifacts reveal. In the lineage paper the
> composition-based natural model was the weaker of the pair (r = 0.643 against
> 0.743), and sequences up to 90% similar were deliberately kept in the data, so
> reported accuracy is optimistic for novel peptides.
