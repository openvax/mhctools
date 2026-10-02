# Choosing a predictor

Choose the biological question first, then check the model's inputs,
installation requirements, and license. The [predictor matrix](predictor-matrix.md)
lists every supported class and command-line name.

## By question

| Question | Predictors |
|---|---|
| Class I binding affinity | [NetMHCpan](predictors/binding.md#netmhcpan) or [MHCflurry](predictors/binding.md#mhcflurry) |
| Class I presentation | [NetMHCpan 4.1/4.2](predictors/binding.md#netmhcpan), [MHCflurry](predictors/binding.md#mhcflurry), [MixMHCpred](predictors/binding.md#mixmhcpred), [BigMHC (EL)](predictors/binding.md#bigmhc), [CapHLA](predictors/binding.md#caphla) |
| Class II binding or presentation | [NetMHCIIpan](predictors/binding.md#netmhciipan) for affinity or presentation; [MixMHC2pred](predictors/binding.md#mixmhc2pred) for presentation |
| Peptide-MHC complex stability | [NetMHCstabpan](predictors/binding.md#netmhcstabpan) |
| Proteasomal cleavage | [Pepsickle](predictors/processing.md#pepsickle) or [NetChop](predictors/processing.md#netchop) |
| Class II cleavage | [NetCleave (class II)](predictors/processing.md#netcleave) |
| TAP transport | [DeepTAP](predictors/processing.md#deeptap) |
| ERAP1 trimming | [ERAMER](predictors/processing.md#eramer) |
| T-cell immunogenicity | [Calis](predictors/immunogenicity.md#calis), [PRIME](predictors/immunogenicity.md#prime), [BigMHC (IM)](predictors/binding.md#bigmhc), [DeepImmuno](predictors/immunogenicity.md#deepimmuno); [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) for class II |
| Recognition by a specific TCR | [NetTCR](predictors/tcr.md#nettcr), [Tulip](predictors/tcr.md#tulip), [MixTCRpred](predictors/tcr.md#mixtcrpred) |
| Free-peptide half-life | [PeptiVerse](predictors/peptide-pk.md#peptiverse), [PlifePred2](predictors/peptide-pk.md#plifepred2) |
| Per-bond peptidase evidence | [Peptidase activity](cleavage/index.md) |

The family guides explain each model's output and validation limits.
[Prediction kinds](kinds.md) defines the corresponding result fields and units.

## By constraint

### License

[MHCflurry](predictors/binding.md#mhcflurry), [CapHLA](predictors/binding.md#caphla), [SMM](predictors/binding.md#smm-and-smm-pmbec)/[SMM-PMBEC](predictors/binding.md#smm-and-smm-pmbec), [Pepsickle](predictors/processing.md#pepsickle), [DeepTAP](predictors/processing.md#deeptap), [DeepImmuno](predictors/immunogenicity.md#deepimmuno), and [Calis](predictors/immunogenicity.md#calis)
are open source or built in. The DTU tools require a license from DTU. The
Gfeller lab tools, [BigMHC](predictors/binding.md#bigmhc), and [NetTCR](predictors/tcr.md#nettcr) have academic, non-commercial terms.
See [licensing](licensing.md) before installing a model.

### Downloads

[Calis](predictors/immunogenicity.md#calis) and `RandomBindingPredictor` need no download. Other models need weights,
reference data, or an external tool; see [getting models](artifacts.md).

### Runtime

[MHCflurry](predictors/binding.md#mhcflurry), [CapHLA](predictors/binding.md#caphla), [SMM](predictors/binding.md#smm-and-smm-pmbec), [Calis](predictors/immunogenicity.md#calis), [Pepsickle](predictors/processing.md#pepsickle)'s neural models, and [NetTCR](predictors/tcr.md#nettcr) run in
the current Python environment. The DTU and Gfeller tools use external
executables. [DeepTAP](predictors/processing.md#deeptap), [DeepImmuno](predictors/immunogenicity.md#deepimmuno), [TLimmuno2](predictors/immunogenicity.md#tlimmuno2), [MixTCRpred](predictors/tcr.md#mixtcrpred), [Tulip](predictors/tcr.md#tulip), [PeptiVerse](predictors/peptide-pk.md#peptiverse),
and [PlifePred2](predictors/peptide-pk.md#plifepred2) use a separate interpreter; see [optional backends](backends.md)
and [environment variables](env-vars.md).

### Inputs

[Calis](predictors/immunogenicity.md#calis), [DeepTAP](predictors/processing.md#deeptap), [ERAMER](predictors/processing.md#eramer), [PeptiVerse](predictors/peptide-pk.md#peptiverse), and [PlifePred2](predictors/peptide-pk.md#plifepred2) accept peptides without
alleles. Cleavage predictors also use flanking residues. See
[input shapes](predictors/index.md#input-shapes) for the other families.

## Good habits

- Compare predictors of the same endpoint. The [recipes](recipes.md) show how
  to combine results in a table.
- Compare physical values only when the units and measurement context agree.
  Model scores use predictor-specific scales; see
  [score, value, and percentile rank](kinds.md#score-value-and-percentile_rank).
- Read the [immunogenicity caveats](predictors/immunogenicity.md#read-this-before-trusting-a-score)
  before ranking neoepitopes. Performance falls toward chance on unseen tumor neoepitopes.
- Check the [known limits](limitations.md) of each model you select.
