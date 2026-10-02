# Getting started

This guide shows how to predict peptide-MHC binding with [MHCflurry](predictors/binding.md#mhcflurry) and read
the results. It assumes you have Python 3.9 or later and a list of peptide
sequences and MHC alleles.

## Install

```sh
pip install mhctools
mhctools fetch mhcflurry
```

[MHCflurry](predictors/binding.md#mhcflurry) is installed with mhctools; the second command downloads its model
weights. Other predictors may need a separate executable or Python environment.
See [getting models](artifacts.md) for download commands and
[optional backends](backends.md) for installation details.

## Predict for peptides

Create a predictor with the alleles you want to evaluate, then pass the peptide
sequences to its `predict()` method:

```python
from mhctools import MHCflurry

peptides = ["SIINFEKL", "GILGFVFTL"]
predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(peptides)
```

The method returns one `PeptideResult` per input peptide, in the same order.
Each result contains the predictions for that peptide across the requested
alleles.

## Read the results

Use `result.affinity` to select the strongest affinity prediction across the
alleles. Its `value` is the predicted IC50 in nM:

```python
for result in results:
    affinity = result.affinity
    if affinity is not None:
        print(result.peptide, affinity.allele, affinity.value)
```

Other accessors include `result.presentation`, `result.processing`, and
`result.immunogenicity`. An accessor returns `None` when the predictor does
not produce that kind of prediction. To examine every allele-specific
prediction, use `result.preds` or `result.filter(allele="HLA-A*02:01")`.

For a table, use the corresponding DataFrame method:

```python
df = predictor.predict_dataframe(peptides)
print(df[["peptide", "allele", "kind", "value", "percentile_rank"]])
```

A table has one row per prediction, so a peptide may appear in several rows
for different alleles and prediction kinds. See
[results and DataFrames](results.md) for the complete output format.

## Use another predictor

Most MHC predictors use the same interface. For example, with a licensed
[NetMHCpan](predictors/binding.md#netmhcpan) 4.2 installation available on your path:

```python
from mhctools import NetMHCpan42

predictor = NetMHCpan42(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(peptides)
```

The [binding and presentation guide](predictors/binding.md) explains the
available versions and modes. Predictors in other families may need flanking
residues or TCR sequences; see [input shapes](predictors/index.md#input-shapes).

## Next steps

- [Choose a predictor](choosing.md) for your question and installation constraints.
- [Scan proteins or annotate a table](recipes.md) with an existing predictor.
- [Use the command line](cli.md) to run predictions without writing Python.
- [Understand prediction kinds](kinds.md) and [model limits](limitations.md) before comparing scores.
