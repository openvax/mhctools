# Peptide half-life predictors

Predictors of how long a free peptide survives. They emit `peptide_half_life`,
read with `result.peptide_half_life`; the sample matrix (serum, plasma, whole
blood) lives in the measurement context. For the wider set of pharmacokinetic,
uptake and tissue-exposure kinds, see [peptide PK, uptake, and tissue
exposure](../exposure-results.md).

| Predictor | What it estimates | Needs |
|---|---|---|
| [Cavaco](#cavaco) | Published-equation baseline; mixed/unknown assay matrix | built in, no extra runtime |
| [PeptiVerse](#peptiverse) | Parent-peptide half-life in human serum | pinned PeptiVerse + ESM2 snapshots (`PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME`) + a torch/transformers Python |
| [PlifePred2](#plifepred2) | Parent-peptide half-life, matrix unknown | `plifepred2==1.1.0` (`PLIFEPRED2_HOME`; bundles Pfeature under `comp/`) |

`peptide_half_life` records how long the parent peptide persists, in hours. Its
context distinguishes a defined solution, serum/plasma/whole blood, cellular
compartments, and systemic in-vivo PK without multiplying kind strings.

It is deliberately a separate kind from `pMHC_stability`, which is the
dissociation half-life of an assembled peptide-MHC complex (a different
molecule in a different assay), and from the cleavage kinds, which are
site-resolved and intracellular. `MeasurementContext` preserves the matrix,
compartment, analyte, and systemic scope when they are known.

## Cavaco

`CavacoHalfLife` independently implements [Cavaco et al. (2021)](https://doi.org/10.1111/cts.12985)
Equation 1 / Table S20 with explicit descriptor provenance. No download,
GUI or additional model runtime is needed. This is a baseline, with
unestablished accuracy on long vaccine peptides.

```python
from mhctools import CavacoHalfLife

predictor = CavacoHalfLife()
prediction = predictor.predict(["SIINFEKL"])[0].peptide_half_life
prediction.score   # ln(half-life in minutes), higher = longer-lived
prediction.value   # exp(score) / 60, in hours
predictor.last_qc  # per-input pI, nonpolar percentage, W/Y counts and status

# Use a separately established pI definition when appropriate; None falls back
# to the author-app descriptor. The choice and value are recorded per result.
predictor.predict(["SIINFEKL"], isoelectric_points=[6.2])
```

On the command line:

```sh
mhctools predict-table --input peptides.csv --out half-life.csv \
  --predictor cavaco:half_life_hours:peptide_half_life
```

The formula is `ln(minutes) = 2.226 + 0.053 * nonpolar_percent - 1.515 *
I(W >= 1) + 1.290 * I(Y >= 2) - 1.052 * I(pI >= 10)`.
Nonpolar residues are `ACFILMPVWY`, following the deposited app's descriptor
data. Default pI reproduces its free-terminal pKa algorithm, including
first-occurrence summation order and JavaScript rounding. Each prediction
records the bundled data hash, descriptor choice and exact peptide identity;
provided pI values also enter the cache fingerprint.

The training set contains 129 peptides from heterogeneous stability and
proteolysis assays, so `measurement_context.matrix` remains `None`. The
sixteen related 15-residue validation peptides were C-terminal carboxamides
tested in 50% human serum at 37 C. The sequence-only API accepts canonical
L-peptides with free termini and explicitly rejects amidation, capping and
attachments; `on_unsupported="record"` preserves unavailable inputs in place.
The equation has no length cutoff, which does not establish long-peptide
accuracy or chemical-form interchangeability. This is neither a per-enzyme
cleavage model nor a whole-serum survival percentage.

### What the adapter reproduces

Two source inconsistencies are handled explicitly: the prose after Equation 1
describes an increased half-life for high pI despite its negative coefficient;
the equation's description says `>10`, whereas methods and app use `>=10`.
The adapter uses the negative coefficient from Equation 1/Table S20 and the
`>=10` boundary. The app's W/Y counting is broken (uppercase residue keys are
tested as lowercase), so its displayed half-life is not the reference
endpoint. The published equation and app descriptor algorithm are distinct
choices; this adapter does not claim native app half-life conformance.

The source oracle covers 31 sequences, using the unmodified app for pI and
an independent equation calculation. On the sixteen source-panel means,
the published equation with these descriptors gives log-scale r² = 0.7602,
log-minute RMSE = 0.9206 and median absolute fold error = 1.808; this is
source-panel numerical reproduction, not an external long-peptide benchmark.
For example, panel peptide 4 estimates 568 minutes against an observed mean
of 89 minutes. App pI also differs from some tabulated source descriptors
(ABD_CL1: 9.6193 versus 9.83); it is not a reconstruction of all training
preprocessing. Broader benchmarking remains tracked in
[#284](https://github.com/openvax/mhctools/issues/284) and
[#291](https://github.com/openvax/mhctools/issues/291).

### Bundled data and source oracle

The app's root manifest declares CC0-1.0. Only descriptor data, factual
coefficients, exact source hashes and attribution are bundled; no article
figures or Electron dependencies are redistributed. The independent adapter
is Apache-2.0. See `mhctools/data/CAVACO_LICENSE.txt` and the fixed
`tests/data/cavaco_source_reference.json` source records (verified against
the author app's `calc.js`, `amino.json`, `codes.json` and `package.json`).
Node is never needed for prediction or tests.

See [target-preserving degradation scenarios](../target-degradation.md) to follow
an exact epitope through successive fragment cuts, with separate estimators and
explicit cut-location assumptions.

## PeptiVerse

[Cell penetration classification](uptake.md#peptiverse-cpp) uses a separate
`PeptiVerseCPP` adapter and returns a class score with no physical uptake value.

[PeptiVerse](#peptiverse) wraps one endpoint of the upstream multi-property platform. Its
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
the hours scale, because the `expm1` inverse transform is applied only to the
sequence model. Although mhctools records exact chemical form, this adapter does not
consume it, so terminal modifications, attachments, and non-standard residues
are rejected rather than scored as their unmodified sequence. Pass
`on_unsupported="record"` to retain unsupported entries in a mixed batch.

PeptiVerse was [published in Nature Communications on July 16, 2026](https://doi.org/10.1038/s41467-026-74167-w);
this adapter retains its explicitly pinned source/model snapshot. The sequence
half-life model was fit on 130 examples and evaluated by cross-validation
only, with no external test set and no
evaluation on long vaccine peptides. Upstream declares Apache-2.0 on its model
card and MIT in its README. mhctools verifies the exact inference source, model,
calibration, ESM2 weights, configuration and tokenizer files before launch. The
PeptiVerse checkpoint and calibration still use unsafe pickle serialization, and
a matching checksum shows identity, not safety, so use only snapshots you trust.

## PlifePred2

This endpoint's physical unit and assay matrix are not established.
[PlifePred2](#plifepred2) 1.1.0 documents ``y = log10(t)`` in its inference
script, but ships no publication, no training data, and no confirmed target
unit or assay-matrix definition. By default the wrapper reports only the
model's native ``log10`` output and claims no duration at all.

```python
from mhctools import PlifePred2

predictor = PlifePred2()                       # PLIFEPRED2_HOME (bundles comp/ Pfeature)
results = predictor.predict(["SIINFEKLGGALQAKKY"])
results[0].peptide_half_life.score             # native output, higher = longer-lived
results[0].peptide_half_life.value             # None by default
predictor.artifact_inventory.to_dict()         # exact files, hashes, capability
predictor.last_qc["log10_seconds"]             # the same value, named

# Opt in to a duration, accepting the unit inference below:
opted_in = PlifePred2(assume_log10_seconds=True)
opted_in.predict(["SIINFEKLGGALQAKKY"])[0].peptide_half_life.value # hours
```

### What is known

The shipped sequence model (`models/Sequence_model.sav`) is a scikit-learn
`ExtraTreesRegressor` over 100 selected Pfeature sequence descriptors plus
five binary modification flags (`d`, `ct`, `nt`, `cyc`, `ptm` — all zero for
natural peptides). Upstream's `plifepred2.py` applies `10**y` to the raw
estimator output. Being monotone in half-life, the native output ranks
peptides identically regardless of the underlying time unit.

### What is inferred

`log10(half-life in seconds)` is the strongest reading. Inverting the forest's
extreme leaf values under it gives `22.4 s` at the floor (consistent with the
~20-second floor in the lineage paper) and `24.85 days` at the ceiling. Under
`log2` or `ln` the same range inverts instead to between a few seconds and
~9 minutes.

Because the forest holds targets past the 2018 lineage paper's 24-hour
ceiling, [PlifePred2](#plifepred2) was trained on a different dataset, and the
old filter cannot independently establish the target unit.

### What is not established

The species and assay matrix are not established. The result therefore
uses generic `peptide_half_life` with `matrix=None`; do not report it as a
measured whole-blood property or treat it as interchangeable with [PeptiVerse](#peptiverse)'s
human-serum endpoint.

Natural peptides only, 12–100 residues; non-standard residues and modified
`PeptideInput` values are rejected explicitly rather than scored as unmodified
sequences.

PlifePred2 1.1.0 bundles plain-Python Pfeature (`comp/pfeature_comp.py` and
`comp/Data/`) inside the wheel. Because PlifePred2 and Pfeature are GPLv3,
nothing is vendored and neither is imported into the mhctools interpreter;
mhctools runs them inside a separate isolated interpreter with verified
checksums and an isolated scratch workspace.

PlifePred2 cites no publication of its own, so its training set is unverified
beyond what the artifacts reveal.
