# Predictors

mhctools wraps about 30 predictors behind one `predict()` call. This section
has one page per family; the [predictor matrix](../predictor-matrix.md) lists
every predictor, class, command-line name, input and license on a single page.

| Family | What it predicts | Page |
|---|---|---|
| MHC binding and presentation | Affinity, eluted-ligand presentation, complex stability (class I and II) | [binding](binding.md) |
| Antigen processing | Proteasomal and endolysosomal cleavage, TAP transport, ERAP1 trimming | [processing](processing.md) |
| Immunogenicity | Whether a peptide elicits a T-cell response | [immunogenicity](immunogenicity.md) |
| TCR specificity | Whether a given TCR recognises a peptide-MHC | [TCR](tcr.md) |
| Peptide half-life | How long a free peptide survives | [peptide half-life](peptide-pk.md) |

Not sure which to use? Start from [choosing a predictor](../choosing.md). Before
you rely on a score, read [known limits](../limitations.md).

## Input shapes

"Every predictor answers `predict()`" is true, but what you pass to it differs:

| Shape | You pass | Predictors |
|---|---|---|
| Peptides + alleles | `predict(peptides)` on a predictor built with `alleles=` | NetMHCpan, NetMHC, NetMHCcons, NetMHCIIpan, NetMHCstabpan, MHCflurry, BigMHC, CapHLA, MixMHCpred, MixMHC2pred, SMM, PRIME, DeepImmuno, TLimmuno2 |
| Peptides + flanks | `predict(peptides, n_flanks=..., c_flanks=...)` | Pepsickle, NetChop, NetCleave (class II needs a C-terminal flank of at least 3 residues) |
| Peptides only | `predict(peptides)` | Calis, DeepTAP, ERAMER, PeptiVerse, PlifePred2 |
| Peptides + TCR | `predict_pairs([(peptide, tcr)])` or `predict(peptides, tcrs)` | NetTCR, Tulip (also `mhc=`) |
| TCRs against a fixed target | `predict_tcrs(tcrs)` on a model chosen for one pMHC | MixTCRpred |
| Exact chemical form | `predict([PeptideInput(...)])` | PeptiVerse (strings still work) |

Every predictor also has `predict_proteins()` to scan protein sequences and
`predict_dataframe()` / `predict_proteins_dataframe()`. See [results and
DataFrames](../results.md).

## Peptide lengths

Scanning a protein needs window lengths, and a predictor's **default window is
much narrower than what it supports**:

| Predictors | Default lengths |
|---|---|
| NetMHCpan, NetMHC, NetMHCcons, MHCflurry, MixMHCpred, SMM, PRIME, Pepsickle, NetChop | 9 |
| NetMHCIIpan | 15-20 |
| MixMHC2pred | 15 |
| NetCleave | 9 (class I), 15 (class II) |
| Historical `Iedb*` names | 8-11 (class I), 15-20 (class II) |
| NetMHCstabpan | none; pass `peptide_lengths=` |

Override per call or per predictor:

```python
predictor.predict_proteins(proteins, peptide_lengths=[8, 9, 10, 11])
NetMHCpan42(alleles=["HLA-A*02:01"], default_peptide_lengths=[8, 9, 10, 11])
```

On the command line use `--mhc-peptide-lengths 8-11`. Passing explicit peptides
to `predict()` is not limited by these defaults, but each model still has its
own supported range (for example `MixMHCpred` 8-14, `CapHLA` 7-25, `DeepImmuno`
9-10 only, `PlifePred2` 12-100, `ERAMER` 9-16).
