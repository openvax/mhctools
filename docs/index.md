# mhctools

mhctools runs MHC binding, presentation, immunogenicity and antigen-processing
predictors through a single `predict()` call and returns the same result
objects whichever one you use. Swapping NetMHCpan for MHCflurry is a one-line
change, and comparing them gives you a DataFrame.

## Available predictors

| I want to predict… | Predictors |
|---|---|
| Binding affinity to an allele | [`NetMHCpan`](predictors/binding.md#netmhcpan), [`NetMHC`](predictors/binding.md#netmhc), [`NetMHCIIpan`](predictors/binding.md#netmhciipan), [`NetMHCcons`](predictors/binding.md#netmhccons), [`MHCflurry`](predictors/binding.md#mhcflurry), [`CapHLA`](predictors/binding.md#caphla), [`SMM`](predictors/binding.md#smm-and-smm-pmbec), [`SMMPMBEC`](predictors/binding.md#smm-and-smm-pmbec) |
| Surface presentation | [`NetMHCpan41`/`42`](predictors/binding.md#netmhcpan), [`NetMHCIIpan`](predictors/binding.md#netmhciipan), [`MHCflurry`](predictors/binding.md#mhcflurry), [`CapHLA`](predictors/binding.md#caphla), [`MixMHCpred`](predictors/binding.md#mixmhcpred) (I), [`MixMHC2pred`](predictors/binding.md#mixmhc2pred) (II), [`BigMHC`](predictors/binding.md#bigmhc) |
| How long the pMHC complex lasts | [`NetMHCstabpan`](predictors/binding.md#netmhcstabpan) |
| Combined antigen processing | [`MHCflurry`](predictors/binding.md#mhcflurry) |
| Proteasomal cleavage | [`Pepsickle`](predictors/processing.md#pepsickle), [`NetChop`](predictors/processing.md#netchop), [`NetCleave_I`](predictors/processing.md#netcleave) |
| Endolysosomal cleavage (class II) | [`NetCleave_II`](predictors/processing.md#netcleave) |
| TAP transport into the ER | [`DeepTAP`](predictors/processing.md#deeptap) |
| ERAP1 N-terminal trimming | [`ERAMER`](predictors/processing.md#eramer) |
| Whether a T cell responds | [`Calis`](predictors/immunogenicity.md#calis), [`PRIME`](predictors/immunogenicity.md#prime), [`BigMHC_IM`](predictors/binding.md#bigmhc), [`DeepImmuno`](predictors/immunogenicity.md#deepimmuno), [`TLimmuno2`](predictors/immunogenicity.md#tlimmuno2) (II) |
| Whether a specific TCR recognises it | [`NetTCR`](predictors/tcr.md#nettcr), [`Tulip`](predictors/tcr.md#tulip), [`MixTCRpred`](predictors/tcr.md#mixtcrpred) |
| How long the free peptide survives | [`PeptiVerse`](predictors/peptide-pk.md#peptiverse), [`PlifePred2`](predictors/peptide-pk.md#plifepred2) |
| Which peptidase cuts which bond | [cleavage API](cleavage/index.md) |

- [Predictor matrix](predictor-matrix.md): every predictor, class, command-line name, input, install route and license on one page.
- [Choosing a predictor](choosing.md) and [known limits](limitations.md). Several of these models are weaker than their own papers suggest; read the limits before you trust a score.

`RandomBindingPredictor` is built in and produces random affinities, which is
occasionally useful as a null baseline.

## Quickstart

```sh
pip install mhctools
mhctools fetch mhcflurry     # MHCflurry ships as a dependency; this downloads its weights
```

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

for r in results:
    if r.affinity:
        print(f"{r.peptide} -> {r.affinity.allele} IC50={r.affinity.value:.1f}nM")
```

`predict()` returns one `PeptideResult` per peptide, in input order, with an
accessor for each [kind of prediction](kinds.md) (`r.affinity`, `r.presentation`,
`r.immunogenicity`, ...) that is `None` when the predictor does not produce it.

Every predictor answers `predict()`, though what you
pass differs by family (alleles, flanks, TCRs); see [input
shapes](predictors/index.md#input-shapes). `Calis` needs no download at all.

## Find what you need

**I am choosing a predictor**

- [Choosing a predictor](choosing.md): by question and by constraint
- [Predictor matrix](predictor-matrix.md): every predictor, class, CLI name, input, install route and license
- [Known limits](limitations.md): read before trusting a score

**I am using one**

- [Predictors](predictors/index.md): one page per family, with an example for each
- [Results and DataFrames](results.md) and [recipes](recipes.md)
- [Allele names](alleles.md) and [peptide lengths](predictors/index.md#peptide-lengths)
- [Command line](cli.md)

**I am installing something**

- [Getting models](artifacts.md): `mhctools fetch`, `ls` and `predictors`
- [Installing optional backends](backends.md), [environment variables](env-vars.md) and [licensing](licensing.md)
- [Troubleshooting](troubleshooting.md)

**I need to understand the output**

- [Prediction kinds, units and MHC context](kinds.md)
- [Peptide PK, uptake and tissue exposure](exposure-results.md)
- [Optional backend conformance](optional-backends.md)

**Cleavage, vaccines and benchmarks**

- [Peptidase cleavage evidence](cleavage/index.md): [models](cleavage/models.md), [batch assessments](cleavage/batch.md), [validation](cleavage/validation.md)
- [Route-aware vaccine reports](vaccine-reports.md)
- [Assay-aware benchmarks](benchmarks.md)

**Reference and maintenance**

- [API reference](api.md)
- [Migration guide](migration.md) and [known gaps](known-gaps.md)
- [Testing](testing.md); releases are described in
  [RELEASING.md](https://github.com/openvax/mhctools/blob/master/RELEASING.md)
