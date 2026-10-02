[![Tests](https://github.com/openvax/mhctools/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/mhctools/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/mhctools.svg?maxAge=1000)](https://pypi.org/project/mhctools/)
[![Python versions](https://img.shields.io/pypi/pyversions/mhctools.svg)](https://pypi.org/project/mhctools/)
[![License](https://img.shields.io/pypi/l/mhctools.svg)](https://github.com/openvax/mhctools/blob/master/LICENSE)
[![Docs](https://img.shields.io/badge/docs-openvax.github.io%2Fmhctools-blue)](https://openvax.github.io/mhctools/)

# mhctools

mhctools is a Python library for running MHC binding, presentation,
immunogenicity, and antigen-processing predictors. It provides a common
interface to tools such as [NetMHCpan](https://openvax.github.io/mhctools/predictors/binding/#netmhcpan) and [MHCflurry](https://openvax.github.io/mhctools/predictors/binding/#mhcflurry), with results you can inspect
in Python or export as a pandas DataFrame.

Read the [getting started guide](https://openvax.github.io/mhctools/getting-started/)
or browse the [documentation](https://openvax.github.io/mhctools/).

## Install and predict

```sh
pip install mhctools
mhctools fetch mhcflurry
```

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

for result in results:
    affinity = result.affinity
    if affinity is not None:
        print(result.peptide, affinity.allele, affinity.value)
```

Each result corresponds to one input peptide. The affinity accessor selects
the strongest prediction across alleles, with IC50 in nM. Use
`predict_dataframe()` for a pandas table or `predict_proteins()` to scan protein
sequences. See [results and DataFrames](https://openvax.github.io/mhctools/results/).

Most predictors need model weights or an external installation. Use
`mhctools ls` to locate models and `mhctools predictors` to check which can run.
The [installation guide](https://openvax.github.io/mhctools/artifacts/) explains
downloads, optional backends, and licensing.

## Guides

- [Choosing a predictor](https://openvax.github.io/mhctools/choosing/) and
  [known limits](https://openvax.github.io/mhctools/limitations/)
- [Predictor matrix](https://openvax.github.io/mhctools/predictor-matrix/): Python classes, CLI names, inputs, and installation routes
- [Recipes](https://openvax.github.io/mhctools/recipes/): protein scans, multiple samples, and table annotation
- [Command line](https://openvax.github.io/mhctools/cli/)
- [Peptidase activity](https://openvax.github.io/mhctools/cleavage/),
  [vaccine reports](https://openvax.github.io/mhctools/vaccine-reports/), and
  [benchmarks](https://openvax.github.io/mhctools/benchmarks/)

## Development

```sh
./develop.sh
./lint.sh
./test.sh
```

See the [testing guide](https://openvax.github.io/mhctools/testing/) and
[release instructions](RELEASING.md).
