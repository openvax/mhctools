# Prediction kinds, units, and MHC context

Every `Prediction` says what it measures (`kind`), how confident or favourable
it is (`score`), and — when the measurement has a physical unit — how much
(`value`). This page is the reference for all three.

- [The kinds](#the-kinds)
- [score, value, and percentile_rank](#score-value-and-percentile_rank)
- [Units are the wrapper's job](#units-are-the-wrappers-job)
- [Measurement context](#measurement-context)
- [MHC dependence and class](#mhc-dependence-and-class)
- [What each predictor emits](#what-each-predictor-emits)

## The kinds

The canonical strings live in `mhctools.pred.Kind`.

| Kind | Meaning | `value` unit |
|---|---|---|
| `pMHC_affinity` | Peptide-MHC binding affinity | `nM` (IC50) |
| `pMHC_presentation` | Likelihood of surface presentation (EL/processing) | — |
| `pMHC_stability` | Peptide-MHC complex stability | `hours` (Thalf) |
| `pMHC_TCR_binding` | TCR recognition of a peptide-MHC (pMHC:TCR binding) | — |
| `immunogenicity` | T-cell immunogenicity | — |
| `antigen_processing` | Combined processing score | — |
| `proteasome_cleavage` | Proteasomal (MHC-I, cytosolic) C-terminal cleavage score | — |
| `endolysosomal_cleavage` | Endolysosomal (MHC-II, cathepsin) C-terminal cleavage score | — |
| `tap_transport` | TAP transport / binding score | `nM` |
| `erap_trimming` | ERAP1 N-terminal trimming score | — |
| `peptide_half_life` | Parent-peptide half-life; matrix and systemic scope live in context | `hours` |
| `systemic_clearance` | Systemic or apparent clearance | context-defined |
| `distribution_volume` | Systemic or apparent distribution volume | context-defined |
| `systemic_exposure` | Systemic exposure, such as AUC | context-defined |
| `cpp_classification` | CPP class label and confidence | — |
| `cellular_uptake` | Quantitative uptake in a named cellular context | context-defined |
| `tissue_concentration` | Concentration in a named tissue/compartment and timepoint | context-defined |

Two kinds that share a unit can still be different measurements.
`pMHC_stability` is the lifetime of a peptide-MHC complex; `peptide_half_life`
is the lifetime of the parent peptide. Serum, plasma, whole blood, cellular,
and systemic settings all share the latter kind and stay distinct through
[`MeasurementContext`](#measurement-context).

## score, value, and percentile_rank

Kind and unit are independent. Every prediction has a `kind`, because every
prediction measures *something*; only some kinds have a unit. A model that
emits a bare 0–1 confidence is still a prediction of a kind — it fills `score`
and leaves `value` empty. Wrappers fill both wherever the predictor supports it.

**`score`** is always present and always orders higher-is-better. Its scale is
predictor-specific: it may be a probability, an uncalibrated model output, a
transformed estimate, or a copy of `value`. Check the predictor's own notes
before comparing or thresholding it. PeptiVerse, for instance, repeats its
predicted hours in both `score` and `value`, while PlifePred2 keeps its
unresolved native output in `score` and leaves `value` empty unless its
inferred conversion is explicitly enabled.

Higher-is-better is a numerical selection convention within one documented
endpoint. A larger score is not universally better for a vaccine, and scores
from different predictors or endpoints are not interchangeable.

**`value`** appears only for the kinds marked with a unit above, and carries a
physical quantity on a **linear** scale in that unit — never a log, never a
rescaling, never whatever the upstream tool happened to print. A kind having a
unit does not oblige every predictor to fill it: a wrapper whose transform to
that unit is unresolved leaves `value` empty rather than guessing (see
[`PlifePred2`](predictors.md#plifepred2)).

**`percentile_rank`** appears when the predictor scores against a background
distribution, and is always lower-is-better.

For affinity predictions, `score` is commonly the monotone `1-log50k`
rescaling and `value` is the estimated IC50 in nM. That rescaling is useful for
ordering predictions. It is not a calibrated probability or confidence, and it
is not inherently bounded to 0–1.

Ask a kind for its unit directly:

```python
from mhctools import Kind
from mhctools.pred import value_unit

value_unit(Kind.pMHC_affinity)     # 'nM'
value_unit(Kind.peptide_half_life) # 'hours'
value_unit(Kind.immunogenicity)    # None
```

## Units are the wrapper's job

Converting to the canonical unit happens in the wrapper, and it long predates
this registry. Affinity predictors commonly work in `1-log50k` space
internally, and every affinity wrapper here inverts it to nM — so a NetMHCpan
IC50 and an MHCflurry IC50 are directly comparable. PeptiVerse's upstream
sequence model applies its `log1p(hours)` inverse and the wrapper reports
hours. PlifePred2's target transform and assay provenance remain unresolved, so
that wrapper reports only the native score by default;
`assume_log10_seconds=True` opts into the inferred conversion to hours.

A predictor's native output is kept, just not in a units-bearing field.
Wrappers park it on `last_qc`:

```python
predictor.predict(["SIINFEKLGGALQAKKY"])
predictor.last_qc["log10_seconds"]     # PlifePred2's raw model output
```

## Measurement context

Every prediction carries a small immutable `MeasurementContext`. Ordinary
predictors get a shared default automatically; assay-specific wrappers fill
only the fields they know. Equal contexts are interned and reused. PK, uptake,
tissue, and peptide-half-life results add explicit units, matrices,
compartments, scope, or time identity as needed.

Full details: [Peptide PK, uptake, and tissue-exposure results](exposure-results.md).

## MHC dependence and class

Predictors expose `kind_support()` so downstream code can tell what MHC context
is meaningful for each kind they emit:

```python
support = predictor.kind_support()
support["pMHC_affinity"]
# {"mhc_dependence": "single_allele", "mhc_class": "I"}
```

`mhc_dependence` is one of:

| Value | Meaning |
|---|---|
| `none` | The prediction is MHC-independent; `Prediction.allele` is empty. |
| `single_allele` | The prediction is for one peptide/MHC allele pair; `Prediction.allele` is part of the key. |
| `haplotype` | The prediction uses the requested MHC repertoire jointly; `Prediction.allele` may carry best-allele attribution but is not the prediction key. |

`mhc_class` is one of `none`, `I`, `II`, or `both`.

The allowed metadata values are defined in `mhctools.pred` as
`MHC_DEPENDENCE_VALUES` and `MHC_CLASS_VALUES`.

## What each predictor emits

| Predictor | Kind | `mhc_dependence` | `mhc_class` |
|---|---|---|---|
| `NetMHCpan41` | `pMHC_affinity` | `single_allele` | `I` |
| `NetMHCpan41` | `pMHC_presentation` | `single_allele` | `I` |
| `NetMHCIIpan4_EL` | `pMHC_presentation` | `single_allele` | `II` |
| `CapHLA` | `pMHC_affinity` | `single_allele` | `both` |
| `CapHLA` | `pMHC_presentation` | `single_allele` | `both` |
| `MixMHC2pred` | `pMHC_presentation` | `single_allele` | `II` |
| `NetMHCstabpan` | `pMHC_stability` | `single_allele` | `I` |
| `MHCflurry` | `pMHC_affinity` | `single_allele` | `I` |
| `MHCflurry` haplotype mode | `pMHC_presentation` | `haplotype` | `I` |
| `MHCflurry` per-allele panel mode | `pMHC_presentation` | `single_allele` | `I` |
| `MHCflurry` | `antigen_processing` | `none` | `none` |
| `Pepsickle` | `proteasome_cleavage` | `none` | `none` |
| `NetCleave_I` | `proteasome_cleavage` | `none` | `I` |
| `NetCleave_II` | `endolysosomal_cleavage` | `none` | `II` |
| `DeepTAP` | `tap_transport` | `none` | `none` |
| `ERAMER` | `erap_trimming` | `none` | `I` |
| `NetTCR` | `pMHC_TCR_binding` | `none` | `I` |
| `Tulip` | `pMHC_TCR_binding` | `single_allele` | `I` |
| `MixTCRpred` | `pMHC_TCR_binding` | `single_allele` | model-specific |
| `BigMHC_IM` | `immunogenicity` | `single_allele` | `I` |
| `PRIME` | `immunogenicity` | `single_allele` | `I` |
| `DeepImmuno` | `immunogenicity` | `single_allele` | `I` |
| `TLimmuno2` | `immunogenicity` | `single_allele` | `II` |
| `Calis` | `immunogenicity` | `none` | `I` |
| `PeptiVerse` | `peptide_half_life` | `none` | `none` |
| `PlifePred2` | `peptide_half_life` | `none` | `none` |
