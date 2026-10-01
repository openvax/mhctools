# Migration guide

Old spellings that still work, and what replaces them.

## Prediction types

| Old | New |
|---|---|
| `BindingPrediction`, `BindingPredictionCollection` | `Prediction` and `PeptideResult` (see [results](results.md)) |
| `Pred`, `PeptidePreds` | Aliases of `Prediction` and `PeptideResult`; use the new names |
| `predict_peptides(peptides)` | `predict(peptides)` |
| `predict_subsequences(sequences, peptide_lengths)` | `predict_proteins(sequences, peptide_lengths=...)` |
| `predict_peptides_dataframe()` | `predict_dataframe()`; now returns the canonical `mhctools.pred.COLUMNS` schema rather than the legacy columns (`prediction_method_name`, `length`, ...) |

Convert legacy results with `collection.to_preds()` (a list of `Prediction`) or
`collection.to_peptide_preds()` (a list of `PeptideResult`). The standalone
`to_peptide_preds()` conversion groups by `(peptide, offset,
source_sequence_name)` because it has no input list, so unlike `predict()` it
does not preserve repeated input peptides.

```python
collection = predictor.predict_subsequences({"1L2Y": "NLYIQWLKDGGPSSGRPPPS"},
                                            peptide_lengths=[9])
df = collection.to_dataframe()
preds = collection.to_preds()
```

## Predictor names

| Old | New |
|---|---|
| `IedbNetMHCpan`, `IedbNetMHCcons`, `IedbNetMHCIIpan`, `IedbSMM`, `IedbSMM_PMBEC` and their `*-iedb` CLI names | Still work as local wrappers. Prefer the explicit local names (`NetMHCpan41_BA`, `NetMHCcons`, `NetMHCIIpan43_BA`, `SMM`, `SMMPMBEC`) when recording provenance. See [compatibility names](predictors/binding.md#compatibility-names-for-the-old-iedb-predictors). |
| IEDB HTTP constructor arguments `url`, `request_timeout`, `raise_on_error` and the CLI option `--do-not-raise-on-error` | Removed. There is no HTTP fallback; missing installations and unsupported inputs raise errors. |

## Commands and options

| Old | New |
|---|---|
| `mhctools integrations` | `mhctools predictors` (the old spelling is an alias) |
| `--mhc-epitope-lengths` | `--mhc-peptide-lengths` (deprecated name still accepted) |

## Kinds

`serum_half_life`, `plasma_half_life`, `blood_half_life` and
`systemic_elimination_half_life` are accepted as input and canonicalized to
`peptide_half_life`, with the matrix kept in the [measurement
context](kinds.md#measurement-context). `canonical_kind()` performs the
mapping.
