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

In the openvax stack, mhctools is the predictor layer: it parses allele names
with [mhcgnomes](https://github.com/pirl-unc/mhcgnomes), runs
[mhcflurry](https://github.com/openvax/mhcflurry) and external tools, and feeds
epitope pipelines such as [topiary](https://github.com/openvax/topiary) and
[vaxrank](https://github.com/openvax/vaxrank).

## Install

Requires Python 3.11 or later and pandas 3.x.

```sh
pip install mhctools
mhctools fetch mhcflurry
```

[MHCflurry](https://openvax.github.io/mhctools/predictors/binding/#mhcflurry) is
installed with mhctools; the second command downloads its model weights. Most
other predictors need model weights or an external installation. Use
`mhctools ls` to locate models and `mhctools predictors` to check which can run.
The [installation guide](https://openvax.github.io/mhctools/artifacts/) explains
downloads, optional backends, and licensing.

## Quick start

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

for result in results:
    affinity = result.affinity
    if affinity is not None:
        print(result.peptide, affinity.allele, affinity.value)
```

```
SIINFEKL HLA-A*02:01 6249.062082358259
GILGFVFTL HLA-A*02:01 25.84136114436916
```

Each result corresponds to one input peptide. The affinity accessor selects
the strongest prediction across alleles, with IC50 in nM. Use
`predict_dataframe()` for a pandas table or `predict_proteins()` to scan protein
sequences. See [results and DataFrames](https://openvax.github.io/mhctools/results/).

## Command line

`mhctools` is one executable. The bare command predicts for peptides or FASTA
sequences and writes a tab-separated table to stdout:

```sh
mhctools --sequence SIINFEKL GILGFVFTL --mhc-predictor mhcflurry --mhc-alleles 'HLA-A*02:01'
```

```
source_sequence_name	offset	peptide	allele	score	affinity	percentile_rank	prediction_method_name	length
	0	SIINFEKL	HLA-A*02:01	0.192203	6249.06	7.96416	mhcflurry	8
	0	GILGFVFTL	HLA-A*02:01	0.699441	25.8414	0.130863	mhcflurry	9
```

Subcommands cover artifacts (`ls`, `fetch`, `predictors`), table annotation
(`predict-table`), peptidase activity (`cleavage`, `cleavenet`), TCR scoring
(`mixtcrpred`), vaccine modeling (`vaccine-report`, `vaccine-trajectory`,
`vaccine-windows`, `vaccine-compare`) and evaluation (`benchmark`). See the
[command-line guide](https://openvax.github.io/mhctools/cli/).

## Documentation

Full documentation: [openvax.github.io/mhctools](https://openvax.github.io/mhctools/).

- [Getting started](https://openvax.github.io/mhctools/getting-started/): first prediction, start here
- [Choosing a predictor](https://openvax.github.io/mhctools/choosing/) and [known limits](https://openvax.github.io/mhctools/limitations/)
- [Predictor matrix](https://openvax.github.io/mhctools/predictor-matrix/): Python classes, CLI names, inputs, and installation routes
- [Recipes](https://openvax.github.io/mhctools/recipes/): protein scans, multiple samples, and table annotation
- [Command line](https://openvax.github.io/mhctools/cli/) and [API reference](https://openvax.github.io/mhctools/api/)
- [Peptidase activity](https://openvax.github.io/mhctools/cleavage/), [vaccine reports](https://openvax.github.io/mhctools/vaccine-reports/), and [benchmarks](https://openvax.github.io/mhctools/benchmarks/)

## Development

```sh
./develop.sh
./lint.sh
./test.sh
```

See the [testing guide](https://openvax.github.io/mhctools/testing/) and
[release instructions](https://github.com/openvax/mhctools/blob/master/RELEASING.md).

## License

[Apache License 2.0](https://github.com/openvax/mhctools/blob/master/LICENSE).
