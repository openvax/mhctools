# Predictors

mhctools provides a common Python interface to MHC binding, presentation,
antigen-processing, immunogenicity, TCR recognition, and peptide half-life
predictors. Each family guide explains the inputs, installation, and output
for its models.

Start with [choosing a predictor](../choosing.md) if you need help selecting a
model. For exact Python classes, command-line names, and installation routes,
see the [predictor matrix](../predictor-matrix.md).

## Available predictors

| Predict | Predictors |
|---|---|
| Binding affinity | [NetMHCpan](binding.md#netmhcpan), [NetMHC](binding.md#netmhc), [NetMHCIIpan](binding.md#netmhciipan), [NetMHCcons](binding.md#netmhccons), [MHCflurry](binding.md#mhcflurry), [CapHLA](binding.md#caphla), [SMM](binding.md#smm-and-smm-pmbec), [SMM-PMBEC](binding.md#smm-and-smm-pmbec) |
| Presentation | [NetMHCpan 4.1/4.2](binding.md#netmhcpan), [NetMHCIIpan](binding.md#netmhciipan), [MHCflurry](binding.md#mhcflurry), [CapHLA](binding.md#caphla), [MixMHCpred](binding.md#mixmhcpred) (I), [MixMHC2pred](binding.md#mixmhc2pred) (II), [BigMHC](binding.md#bigmhc) |
| Binding stability | [NetMHCstabpan](binding.md#netmhcstabpan) |
| Antigen processing | [MHCflurry](binding.md#mhcflurry) |
| Proteasomal cleavage | [Pepsickle](processing.md#pepsickle), [NetChop](processing.md#netchop), [NetCleave (class I)](processing.md#netcleave) |
| Endolysosomal cleavage | [NetCleave (class II)](processing.md#netcleave) |
| TAP transport | [DeepTAP](processing.md#deeptap) |
| ERAP1 trimming | [ERAMER](processing.md#eramer) |
| Immunogenicity | [Calis](immunogenicity.md#calis), [PRIME](immunogenicity.md#prime), [BigMHC (IM)](binding.md#bigmhc), [DeepImmuno](immunogenicity.md#deepimmuno), [TLimmuno2](immunogenicity.md#tlimmuno2) (II) |
| TCR recognition | [NetTCR](tcr.md#nettcr), [Tulip](tcr.md#tulip), [MixTCRpred](tcr.md#mixtcrpred) |
| Peptide half-life | [PeptiVerse](peptide-pk.md#peptiverse), [PlifePred2](peptide-pk.md#plifepred2) |
| Per-bond cleavage | [cleavage API](../cleavage/index.md) |

Read the [known limits](../limitations.md) before interpreting a score.

`RandomBindingPredictor` is built in and produces random affinities, which is
occasionally useful as a null baseline.

## Input shapes

The prediction method and its inputs depend on the family:

| Shape | You pass | Predictors |
|---|---|---|
| Peptides + alleles | `predict(peptides)` on a predictor built with `alleles=` | [NetMHCpan](binding.md#netmhcpan), [NetMHC](binding.md#netmhc), [NetMHCcons](binding.md#netmhccons), [NetMHCIIpan](binding.md#netmhciipan), [NetMHCstabpan](binding.md#netmhcstabpan), [MHCflurry](binding.md#mhcflurry), [BigMHC](binding.md#bigmhc), [CapHLA](binding.md#caphla), [MixMHCpred](binding.md#mixmhcpred), [MixMHC2pred](binding.md#mixmhc2pred), [SMM](binding.md#smm-and-smm-pmbec), [PRIME](immunogenicity.md#prime), [DeepImmuno](immunogenicity.md#deepimmuno), [TLimmuno2](immunogenicity.md#tlimmuno2) |
| Peptides + flanks | `predict(peptides, n_flanks=..., c_flanks=...)` | [Pepsickle](processing.md#pepsickle), [NetChop](processing.md#netchop), [NetCleave](processing.md#netcleave) (class II needs a C-terminal flank of at least 3 residues) |
| Peptides only | `predict(peptides)` | [Calis](immunogenicity.md#calis), [DeepTAP](processing.md#deeptap), [ERAMER](processing.md#eramer), [PeptiVerse](peptide-pk.md#peptiverse), [PlifePred2](peptide-pk.md#plifepred2) |
| Peptides + TCR | `predict_pairs([(peptide, tcr)])` or `predict(peptides, tcrs)` | [NetTCR](tcr.md#nettcr), [Tulip](tcr.md#tulip) (also `mhc=`) |
| TCRs against a fixed target | `predict_tcrs(tcrs)` on a model chosen for one pMHC | [MixTCRpred](tcr.md#mixtcrpred) |
| Exact chemical form | `predict([PeptideInput(...)])` | PeptiVerse (strings still work) |

Every predictor also has `predict_proteins()` to scan protein sequences and
`predict_dataframe()` / `predict_proteins_dataframe()`. See [results and
DataFrames](../results.md).

## Peptide lengths

Protein scans use default window lengths unless you specify them. These
defaults are usually narrower than the model's supported range:

| Predictors | Default lengths |
|---|---|
| [NetMHCpan](binding.md#netmhcpan), [NetMHC](binding.md#netmhc), [NetMHCcons](binding.md#netmhccons), [MHCflurry](binding.md#mhcflurry), [MixMHCpred](binding.md#mixmhcpred), [SMM](binding.md#smm-and-smm-pmbec), [PRIME](immunogenicity.md#prime), [Pepsickle](processing.md#pepsickle), [NetChop](processing.md#netchop) | 9 |
| [NetMHCIIpan](binding.md#netmhciipan) | 15-20 |
| [MixMHC2pred](binding.md#mixmhc2pred) | 15 |
| [NetCleave](processing.md#netcleave) | 9 (class I), 15 (class II) |
| Historical `Iedb*` names | 8-11 (class I), 15-20 (class II) |
| [NetMHCstabpan](binding.md#netmhcstabpan) | none; pass `peptide_lengths=` |

Override per call or per predictor:

```python
predictor.predict_proteins(proteins, peptide_lengths=[8, 9, 10, 11])
NetMHCpan42(alleles=["HLA-A*02:01"], default_peptide_lengths=[8, 9, 10, 11])
```

On the command line use `--mhc-peptide-lengths 8-11`. Passing explicit peptides
to `predict()` is not limited by these defaults, but each model still has its
own supported range (for example MixMHCpred 8-14, [CapHLA](binding.md#caphla) 7-25, [DeepImmuno](immunogenicity.md#deepimmuno)
9-10 only, [PlifePred2](peptide-pk.md#plifepred2) 12-100, [ERAMER](processing.md#eramer) 9-16).
