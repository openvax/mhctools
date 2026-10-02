# Results and DataFrames

The `predict()` method returns a list of `PeptideResult` objects, one per input
peptide. Each result contains predictions across the requested alleles and
prediction kinds.

## Read a peptide result

An accessor such as `result.affinity` selects the best prediction of that kind.
It returns `None` when the predictor does not produce the kind. For example,
an affinity-only model has no stability result.

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])
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

## Read a prediction

Each peptide result contains a tuple of `Prediction` objects, one per allele
and kind. These immutable objects also carry the predictor and source-sequence
information:

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

The main numerical fields are:

- `score`: always present and higher-is-better, on a predictor-specific scale.
- `value`: a physical quantity on a linear scale, in that kind's canonical
  unit (nM for affinity, hours for stability). Empty when the kind has no unit,
  or when the wrapper cannot convert to it.
- `percentile_rank`: 0-100, lower is stronger, present when the predictor
  scores against a background distribution.

A predictor can emit more than one kind. [NetMHCpan](predictors/binding.md#netmhcpan) 4.1, for example, produces
both `pMHC_affinity` and `pMHC_presentation` for every peptide-allele pair.

See [prediction kinds](kinds.md) for units and measurement context, and the
[predictor matrix](predictor-matrix.md) for each model's output kinds.

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

## Ordering and repeated peptides

Results preserve input order and repeated peptides. `results[i]` corresponds
to `peptides[i]`, so you can use `zip(peptides, results)`. A repeated peptide
gets its own `PeptideResult` for each occurrence. The [NetMHC](predictors/binding.md#netmhc) family and the
legacy binding-prediction fallback score each distinct peptide once and expand
its predictions back to every occurrence. If those backends omit a requested
peptide/allele pair, `predict()` raises `ValueError` instead of returning a
shorter, misaligned list.
