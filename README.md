[![Tests](https://github.com/openvax/mhctools/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/mhctools/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/mhctools.svg?maxAge=1000)](https://pypi.org/project/mhctools/)
[![Python versions](https://img.shields.io/pypi/pyversions/mhctools.svg)](https://pypi.org/project/mhctools/)
[![License](https://img.shields.io/pypi/l/mhctools.svg)](https://github.com/openvax/mhctools/blob/master/LICENSE)
[![Docs](https://img.shields.io/badge/docs-openvax.github.io%2Fmhctools-blue)](https://openvax.github.io/mhctools/)

# mhctools

mhctools runs MHC binding, presentation, immunogenicity and antigen-processing
predictors through a single `predict()` call and returns the same result
objects whichever one you use. Swapping NetMHCpan for MHCflurry is a one-line
change, and comparing them gives you a DataFrame.

**Documentation: <https://openvax.github.io/mhctools/>**

## Available predictors

| I want to predict… | Predictors |
|---|---|
| Binding affinity to an allele | [`NetMHCpan`](https://openvax.github.io/mhctools/predictors/binding/#netmhcpan), [`NetMHC`](https://openvax.github.io/mhctools/predictors/binding/#netmhc), [`NetMHCIIpan`](https://openvax.github.io/mhctools/predictors/binding/#netmhciipan), [`NetMHCcons`](https://openvax.github.io/mhctools/predictors/binding/#netmhccons), [`MHCflurry`](https://openvax.github.io/mhctools/predictors/binding/#mhcflurry), [`CapHLA`](https://openvax.github.io/mhctools/predictors/binding/#caphla), [`SMM`](https://openvax.github.io/mhctools/predictors/binding/#smm-and-smm-pmbec), [`SMMPMBEC`](https://openvax.github.io/mhctools/predictors/binding/#smm-and-smm-pmbec) |
| Surface presentation | [`NetMHCpan41`/`42`](https://openvax.github.io/mhctools/predictors/binding/#netmhcpan), [`NetMHCIIpan`](https://openvax.github.io/mhctools/predictors/binding/#netmhciipan), [`MHCflurry`](https://openvax.github.io/mhctools/predictors/binding/#mhcflurry), [`CapHLA`](https://openvax.github.io/mhctools/predictors/binding/#caphla), [`MixMHCpred`](https://openvax.github.io/mhctools/predictors/binding/#mixmhcpred) (I), [`MixMHC2pred`](https://openvax.github.io/mhctools/predictors/binding/#mixmhc2pred) (II), [`BigMHC`](https://openvax.github.io/mhctools/predictors/binding/#bigmhc) |
| How long the pMHC complex lasts | [`NetMHCstabpan`](https://openvax.github.io/mhctools/predictors/binding/#netmhcstabpan) |
| Combined antigen processing | [`MHCflurry`](https://openvax.github.io/mhctools/predictors/binding/#mhcflurry) |
| Proteasomal cleavage | [`Pepsickle`](https://openvax.github.io/mhctools/predictors/processing/#pepsickle), [`NetChop`](https://openvax.github.io/mhctools/predictors/processing/#netchop), [`NetCleave_I`](https://openvax.github.io/mhctools/predictors/processing/#netcleave) |
| Endolysosomal cleavage (class II) | [`NetCleave_II`](https://openvax.github.io/mhctools/predictors/processing/#netcleave) |
| TAP transport into the ER | [`DeepTAP`](https://openvax.github.io/mhctools/predictors/processing/#deeptap) |
| ERAP1 N-terminal trimming | [`ERAMER`](https://openvax.github.io/mhctools/predictors/processing/#eramer) |
| Whether a T cell responds | [`Calis`](https://openvax.github.io/mhctools/predictors/immunogenicity/#calis), [`PRIME`](https://openvax.github.io/mhctools/predictors/immunogenicity/#prime), [`BigMHC_IM`](https://openvax.github.io/mhctools/predictors/binding/#bigmhc), [`DeepImmuno`](https://openvax.github.io/mhctools/predictors/immunogenicity/#deepimmuno), [`TLimmuno2`](https://openvax.github.io/mhctools/predictors/immunogenicity/#tlimmuno2) (II) |
| Whether a specific TCR recognises it | [`NetTCR`](https://openvax.github.io/mhctools/predictors/tcr/#nettcr), [`Tulip`](https://openvax.github.io/mhctools/predictors/tcr/#tulip), [`MixTCRpred`](https://openvax.github.io/mhctools/predictors/tcr/#mixtcrpred) |
| How long the free peptide survives | [`PeptiVerse`](https://openvax.github.io/mhctools/predictors/peptide-pk/#peptiverse), [`PlifePred2`](https://openvax.github.io/mhctools/predictors/peptide-pk/#plifepred2) |
| Which peptidase cuts which bond | [cleavage API](https://openvax.github.io/mhctools/cleavage/) |

- [Predictor matrix](https://openvax.github.io/mhctools/predictor-matrix/): every predictor, class, command-line name, input, install route and license on one page.
- [Choosing a predictor](https://openvax.github.io/mhctools/choosing/) and [known limits](https://openvax.github.io/mhctools/limitations/). Several of these models are weaker than their own papers suggest; read the limits before you trust a score.

`RandomBindingPredictor` is built in and produces random affinities, which is
occasionally useful as a null baseline.

## Install and predict

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

`predict()` returns one `PeptideResult` per input
peptide, in input order. Each exposes an accessor per [kind of prediction](https://openvax.github.io/mhctools/kinds/)
(`r.affinity`, `r.presentation`, `r.immunogenicity`, ...) that is `None` when
the predictor does not produce that kind. Scan proteins with
`predict_proteins()`, and get a pandas DataFrame from any `*_dataframe()` method.
See [results and DataFrames](https://openvax.github.io/mhctools/results/).

`Calis` needs no download. Most other predictors need model weights or an
external tool first (see below). The command line does the same job:

```sh
mhctools --sequence SIINFEKL SIINFEKLQ --mhc-predictor mhcflurry --mhc-alleles A0201
```

## Getting models

Most predictors need something downloaded first, with one command for all of it:

```sh
mhctools ls                       # what exists, where it lives, who manages it
mhctools fetch mhcflurry          # get it
mhctools predictors               # can it actually run?
```

`fetch` is idempotent. Academic-licensed tools need an explicit
`--accept-license`, and the DTU NetMHC family needs a license you request from
DTU directly. See [getting models](https://openvax.github.io/mhctools/artifacts/)
and [licensing](https://openvax.github.io/mhctools/licensing/).

## Beyond peptide-MHC

- **[Per-bond peptidase evidence](https://openvax.github.io/mhctools/cleavage/)**, including
  [contextual batches](https://openvax.github.io/mhctools/cleavage/batch/) for epitopes with
  flanks and complete vaccine constructs under tumor, APC and extracellular scenarios.
- **[Route-aware vaccine reports](https://openvax.github.io/mhctools/vaccine-reports/)**:
  `mhctools vaccine-report` writes a sequence-centered PDF with route policies and checksums.
- **[Assay-aware benchmarks](https://openvax.github.io/mhctools/benchmarks/)**:
  `mhctools benchmark` evaluates source-linked observations with training-overlap reporting.
- **[Peptide PK, uptake and exposure](https://openvax.github.io/mhctools/exposure-results/)** result kinds.

## Documentation

| Start here | |
|---|---|
| [Command line](https://openvax.github.io/mhctools/cli/) | Every `mhctools` subcommand |
| [Results and DataFrames](https://openvax.github.io/mhctools/results/) | `PeptideResult`, `Prediction`, columns |
| [Recipes](https://openvax.github.io/mhctools/recipes/) | Scan proteins, many genotypes, annotate a table |
| [Allele names](https://openvax.github.io/mhctools/alleles/) | Accepted spellings and errors |
| [Troubleshooting](https://openvax.github.io/mhctools/troubleshooting/) | Common failures |
| [Migration guide](https://openvax.github.io/mhctools/migration/) | Old names and what replaced them |

## Development

```sh
./develop.sh    # editable install
./lint.sh       # ruff
./test.sh       # pytest
```

See the [testing guide](https://openvax.github.io/mhctools/testing/) for a complete run with no skipped tests.
Releases are described in [RELEASING.md](RELEASING.md).
