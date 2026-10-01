# mhctools

One Python interface to about 30 MHC binding, presentation, immunogenicity and
antigen-processing predictors.

Each predictor has its own input format, output format, allele spelling and
installation ritual. mhctools gives them all the same `predict()` call and the
same result objects, so swapping NetMHCpan for MHCflurry is a one-line change
and comparing them is a DataFrame. Everything runs locally.

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

Every predictor is built the same way and answers `predict()`, though what you
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
