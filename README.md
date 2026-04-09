[![Tests](https://github.com/openvax/mhctools/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/mhctools/actions/workflows/tests.yml)
<a href="https://pypi.python.org/pypi/mhctools/">
<img src="https://img.shields.io/pypi/v/mhctools.svg?maxAge=1000" alt="PyPI" />
</a>

# mhctools

Python interface to MHC binding, presentation, immunogenicity, and antigen processing predictors.

## Installation

```sh
pip install mhctools
```

For MHCflurry support, also run:

```sh
mhcflurry-downloads fetch
```

## Quick start

```python
from mhctools import NetMHCpan41

predictor = NetMHCpan41(alleles=["HLA-A*02:01", "HLA-B*07:02"])

# predict for specific peptides
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

# results is a list of PeptidePreds — one per peptide
for pp in results:
    best = pp.best_affinity
    if best:
        print(f"{best.peptide} -> {best.allele} IC50={best.value:.1f}nM")
```

## Python API

### Predicting peptides

`predict()` takes a list of peptide sequences and returns a `list[PeptidePreds]`.
Each `PeptidePreds` contains `Pred` objects for every allele and measurement
kind the predictor supports.

```python
from mhctools import NetMHCpan41

predictor = NetMHCpan41(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

pp = results[0]
pp.best_affinity                # Pred with highest affinity score
pp.best_affinity.allele         # "HLA-A*02:01"
pp.best_affinity.value          # IC50 in nM
pp.best_affinity.score          # higher = better (~0-1)
pp.best_affinity.percentile_rank  # lower = better (0-100)

pp.best_affinity_by_rank        # Pred with lowest percentile rank
pp.best_presentation            # best EL/presentation score
pp.best_presentation_by_rank    # best EL percentile rank
pp.best_stability               # best pMHC stability (if available)
pp.best_stability_by_rank

# filter by kind or allele
pp.filter(kind=Kind.pMHC_affinity)
pp.filter(allele="HLA-A*02:01")
```

NetMHCpan 4.1 automatically emits both `pMHC_affinity` and `pMHC_presentation`
predictions per peptide-allele pair.

### Scanning proteins

`predict_proteins()` takes a dictionary of protein sequences and returns
`{sequence_name: list[PeptidePreds]}`:

```python
proteins = predictor.predict_proteins(
    {"TP53": "MEEPQSDPSVEPPLSQETFS...", "KRAS": "MTEYKLVVVGAGGVGKS..."},
    peptide_lengths=[9, 10],
)

for pp in proteins["TP53"]:
    best = pp.best_affinity
    if best and best.value < 500:
        print(f"  offset={best.offset} {best.peptide} IC50={best.value:.0f}")
```

### DataFrames

Every level has a `_dataframe` variant that flattens to a pandas DataFrame
with consistent columns:

```python
df = predictor.predict_dataframe(["SIINFEKL"], sample_name="pat001")
df = predictor.predict_proteins_dataframe({"TP53": "MEEPQ..."}, sample_name="pat001")
```

Columns: `sample_name`, `peptide`, `n_flank`, `c_flank`,
`source_sequence_name`, `offset`, `predictor_name`, `predictor_version`,
`allele`, `kind`, `score`, `value`, `percentile_rank`.

### Multi-sample predictions

`MultiSample` runs a predictor across multiple samples, each with its own
HLA genotype:

```python
from mhctools import MultiSample, NetMHCpan41

ms = MultiSample(
    samples={
        "pat001": ["HLA-A*02:01", "HLA-B*07:02"],
        "pat002": ["HLA-A*01:01", "HLA-B*08:01"],
    },
    predictor_class=NetMHCpan41,
)

# {sample_name: list[PeptidePreds]}
results = ms.predict(["SIINFEKL", "GILGFVFTL"])

# {sample_name: {seq_name: list[PeptidePreds]}}
protein_results = ms.predict_proteins({"TP53": "MEEPQ..."})

# flat DataFrames with sample_name column
df = ms.predict_dataframe(["SIINFEKL"])
df = ms.predict_proteins_dataframe({"TP53": "MEEPQ..."})
```

### Measurement kinds

The `Kind` enum describes what biological quantity a `Pred` measures:

| Kind | Meaning |
|---|---|
| `pMHC_affinity` | Peptide-MHC binding affinity |
| `pMHC_presentation` | Likelihood of surface presentation (EL) |
| `pMHC_stability` | Peptide-MHC complex stability |
| `immunogenicity` | T-cell immunogenicity |
| `cellular_presentation` | Cross-allele presentation (e.g. MHCflurry) |
| `antigen_processing` | Combined processing score |
| `proteasome_cleavage` | Proteasomal cleavage score |
| `tap_transport` | TAP transport score |
| `erap_trimming` | ERAP trimming score |

### The Pred object

Every prediction is a frozen, self-contained `Pred` dataclass:

```python
from mhctools import Pred, Kind

pred = Pred(
    kind=Kind.pMHC_affinity,
    score=0.85,           # ~0-1, higher = better
    peptide="SIINFEKL",
    allele="HLA-A*02:01",
    value=120.5,          # IC50 in nM
    percentile_rank=0.8,
    source_sequence_name="TP53",
    offset=42,
    predictor_name="netMHCpan",
    predictor_version="4.1",
)
```

`score` is always higher-is-better. `value` is in native units (nM for
affinity, hours for stability). `percentile_rank` is always optional,
0-100, lower = stronger.

## Supported predictors

### MHC binding & presentation

| Predictor | Kinds produced | Requires |
|---|---|---|
| `NetMHCpan` / `NetMHCpan41` / `NetMHCpan42` | affinity + presentation | [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) |
| `NetMHCpan4` | affinity or presentation | NetMHCpan 4.0 |
| `NetMHCpan3` / `NetMHCpan28` | affinity | older NetMHCpan |
| `NetMHC` / `NetMHC3` / `NetMHC4` | affinity | [NetMHC](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) |
| `NetMHCIIpan` / `NetMHCIIpan43` | affinity or presentation | [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) |
| `NetMHCcons` | affinity | [NetMHCcons](https://services.healthtech.dtu.dk/services/NetMHCcons-1.1/) |
| `NetMHCstabpan` | stability | [NetMHCstabpan](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/) |
| `MHCflurry` | affinity + presentation | `pip install mhcflurry` + `mhcflurry-downloads fetch` |
| `BigMHC` | presentation or immunogenicity | [BigMHC](https://github.com/KarchinLab/bigmhc) clone (set `BIGMHC_DIR`) |
| `MixMHCpred` | presentation | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) |
| `IedbNetMHCpan` / `IedbSMM` / `IedbNetMHCIIpan` | affinity | IEDB web API |
| `RandomBindingPredictor` | affinity | (built-in) |

### Antigen processing

| Predictor | Kinds produced | Requires |
|---|---|---|
| `Pepsickle` | proteasome cleavage | `pip install mhctools[pepsickle]` |
| `NetChop` | proteasome cleavage | [NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/) |

Processing predictors use configurable scoring to aggregate per-position
cleavage probabilities into peptide-level scores. See `ProcessingPredictor`
and `ProteasomePredictor` for details.

## Commandline examples

### Prediction for user-supplied peptide sequences

```sh
mhctools --sequence SIINFEKL SIINFEKLQ --mhc-predictor netmhc --mhc-alleles A0201
```

### Automatically extract peptides as subsequences of specified length

```sh
mhctools --sequence AAAQQQSIINFEKL --extract-subsequences --mhc-peptide-lengths 8-10 --mhc-predictor mhcflurry --mhc-alleles A0201
```

## Legacy API

The old `predict_peptides()` and `predict_subsequences()` methods still work
and return `BindingPredictionCollection` objects:

```python
predictor = NetMHCpan(alleles=["A*02:01"])
collection = predictor.predict_subsequences(
    {"1L2Y": "NLYIQWLKDGGPSSGRPPPS"},
    peptide_lengths=[9],
)
df = collection.to_dataframe()

for bp in collection:
    if bp.affinity < 100:
        print("Strong binder: %s" % bp)
```

To convert legacy results to the new types:

```python
preds = collection.to_preds()           # list of Pred
pp_list = collection.to_peptide_preds() # list of PeptidePreds
```
