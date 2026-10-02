# mhctools

mhctools is a Python library for running MHC binding, presentation,
immunogenicity, and antigen-processing predictors. It provides a common
interface to tools such as [NetMHCpan](predictors/binding.md#netmhcpan) and [MHCflurry](predictors/binding.md#mhcflurry), with results you can inspect
in Python or export as a pandas DataFrame.

## Quickstart

Install mhctools and download the [MHCflurry](predictors/binding.md#mhcflurry) model weights:

```sh
pip install mhctools
mhctools fetch mhcflurry
```

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])
df = predictor.predict_dataframe(["SIINFEKL", "GILGFVFTL"])
```

The [getting started guide](getting-started.md) explains the results and shows
how to use another predictor.

## Documentation

- [Getting started](getting-started.md): install the library and make your first prediction.
- [User guide](predictors/index.md): supported predictors, inputs, and examples.
- [Recipes](recipes.md): scan proteins, run multiple samples, and annotate tables.
- [API reference](api.md): Python classes and functions.

For model selection, read [choosing a predictor](choosing.md) and the
[known limits](limitations.md). For installation help, see
[getting models](artifacts.md) and [troubleshooting](troubleshooting.md).

## Processing and vaccine analysis

The [antigen-processing guide](predictors/processing.md) covers proteasomes,
peptidases, transport, and trimming, with [model recommendations](cleavage/choosing.md)
and [batch assessments](cleavage/batch.md). Related workflows include
[vaccine reports](vaccine-reports.md) and [assay-aware benchmarks](benchmarks.md).
