# Results and DataFrames

What `predict()` returns and how to read it.

**`predict()` returns a list of `PeptideResult`.** Each one carries the peptide
string and gives you accessors for each kind of prediction. An accessor returns
`None` when the predictor doesn't produce that kind — so `r.stability` is
`None` from an affinity-only model, rather than an error or a zero.

**Results preserve input order and repeated peptides.** `results[i]` corresponds
to `peptides[i]`, so you can use `zip(peptides, results)`. A repeated peptide
gets its own `PeptideResult` for each occurrence. The NetMHC family and the
legacy binding-prediction fallback score each distinct peptide once and expand
its predictions back to every occurrence. If those backends omit a requested
peptide/allele pair, `predict()` raises `ValueError` instead of returning a
shorter, misaligned list.

```python
from mhctools import MHCflurry

results = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"]).predict(["SIINFEKL", "GILGFVFTL"])
r = results[0]

r.peptide                      # "SIINFEKL"
r.offset                       # position in source protein (if scanned)
r.kinds                        # {"pMHC_affinity", "pMHC_presentation", "antigen_processing"}
r.alleles                      # {"HLA-A*02:01", "HLA-B*07:02"}

# best prediction of each kind, or None when the kind is absent
r.affinity
r.presentation
r.stability

if r.affinity:
    r.affinity.value            # IC50 in nM
    r.affinity.percentile_rank  # 0-100, lower = better
    r.affinity.score            # predictor-specific scale, higher = better
    r.affinity.allele           # best allele for this kind

r.best_affinity_by_rank        # by lowest percentile rank instead of score

r.preds                        # tuple of every underlying Prediction
r.filter(kind="pMHC_affinity")
r.filter(allele="HLA-A*02:01")
```

**Underneath, each `PeptideResult` wraps a tuple of `Prediction` objects** —
frozen dataclasses, one per allele-kind combination, each self-contained:

```python
from mhctools import Prediction

pred = Prediction(
    kind="pMHC_affinity",
    score=0.85,           # predictor-specific scale, higher = better
    peptide="SIINFEKL",
    allele="HLA-A*02:01",
    value=120.5,          # IC50 in nM
    percentile_rank=0.8,
    source_sequence_name="TP53",
    offset=42,
    predictor_name="netMHCpan",
    predictor_version="4.1",
)
```

Three fields do most of the work, and they mean the same thing everywhere:

- **`score`** — always present, always higher-is-better, but on a
  predictor-specific scale.
- **`value`** — a physical quantity on a linear scale, in that kind's canonical
  unit (nM for affinity, hours for stability). Empty when the kind has no unit,
  or when the wrapper cannot honestly convert to it.
- **`percentile_rank`** — 0-100, lower is stronger, present when the predictor
  scores against a background distribution.

A predictor can emit more than one kind. NetMHCpan 4.1, for example, produces
both `pMHC_affinity` and `pMHC_presentation` for every peptide-allele pair.

Full reference: **[prediction kinds, units, and MHC context](kinds.md)**
— the kinds, what fills `value`, and how `MeasurementContext` works. Which
predictor emits which kind is in the [predictor matrix](predictor-matrix.md).

## Which accessor for which kind

Each kind has an accessor on `PeptideResult` that returns the best prediction of
that kind, or `None`:

| Kind | Accessor |
|---|---|
| `pMHC_affinity` | `result.affinity` |
| `pMHC_presentation` | `result.presentation` |
| `pMHC_stability` | `result.stability` |
| `antigen_processing` | `result.processing` |
| `proteasome_cleavage` | `result.cleavage` |
| `endolysosomal_cleavage` | `result.endolysosomal_cleavage` |
| `tap_transport` | `result.tap_transport` |
| `erap_trimming` | `result.erap_trimming` |
| `immunogenicity` | `result.immunogenicity` |
| `pMHC_TCR_binding` | `result.tcr_binding` |
| `peptide_half_life` | `result.peptide_half_life` (also `serum_half_life`, `plasma_half_life`, `blood_half_life` by matrix) |

The other kinds (`systemic_clearance`, `cellular_uptake`, ...) are reached with
`result.filter(kind=...)`; they have no best-of ordering, see [peptide PK,
uptake, and tissue exposure](exposure-results.md#ordering-and-downstream-projection).
`result.to_dict()` and `result.to_dataframe()` serialise every underlying
prediction.

## DataFrames

Every level has a `_dataframe` variant that flattens to a pandas DataFrame:

```python
df = predictor.predict_dataframe(["SIINFEKL"], sample_name="pat001")
df = predictor.predict_proteins_dataframe({"TP53": "MEEPQ..."}, sample_name="pat001")
```

The columns are the same for every predictor, and are defined once as
`mhctools.pred.COLUMNS`:

| Column | |
|---|---|
| `sample_name` | who this was predicted for |
| `peptide`, `n_flank`, `c_flank` | the peptide and its flanking residues |
| `source_sequence_name`, `offset` | where it came from, if scanned |
| `predictor_name`, `predictor_version` | what produced it |
| `allele`, `tcr` | the MHC allele and/or TCR, empty when not applicable |
| `kind`, `score`, `value`, `percentile_rank` | [the prediction itself](kinds.md) |
| `measurement_context`, `peptide_input`, `cache_key` | assay context and exact-input identity |
