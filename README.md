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

# predict() returns a list of PeptideResult — one per peptide
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

for r in results:
    if r.affinity:
        print(f"{r.peptide} -> {r.affinity.allele} IC50={r.affinity.value:.1f}nM")
```

## Data model

`predict()` returns a list of `PeptideResult` — one per peptide. Each
result carries the peptide string and provides accessors for each
prediction kind (affinity, presentation, stability, etc.). Accessors
return `None` when a predictor doesn't produce that kind.

```python
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])
r = results[0]

r.peptide                    # "SIINFEKL"
r.affinity.value             # IC50 in nM
r.affinity.percentile_rank   # 0-100, lower = better
r.affinity.allele            # best allele for this kind
r.presentation               # None if predictor doesn't produce it
```

Under the hood, each `PeptideResult` wraps a tuple of `Prediction` objects —
frozen dataclasses, one per allele-kind combination. Everything converts
to DataFrames with consistent column names.

## Python API

### Predicting peptides

```python
from mhctools import NetMHCpan41

predictor = NetMHCpan41(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

r = results[0]
r.peptide                      # "SIINFEKL"
r.offset                       # position in source protein (if scanned)
r.kinds                        # {"pMHC_affinity", "pMHC_presentation"}
r.alleles                      # {"HLA-A*02:01", "HLA-B*07:02"}

# best prediction by kind — None when the kind is absent
r.affinity                     # Prediction or None
r.presentation                 # Prediction or None
r.stability                    # None (predictor doesn't produce it)

if r.affinity:
    r.affinity.value            # IC50 in nM
    r.affinity.percentile_rank  # 0-100, lower = better
    r.affinity.score            # ~0-1, higher = better
    r.affinity.allele           # best allele for this kind

# by rank instead of score
r.best_affinity_by_rank        # Prediction with lowest percentile rank, or None

# all predictions
r.preds                        # tuple of all Prediction objects
r.filter(kind="pMHC_affinity")
r.filter(allele="HLA-A*02:01")
```

NetMHCpan 4.1 automatically emits both `pMHC_affinity` and `pMHC_presentation`
predictions per peptide-allele pair.

### Scanning proteins

`predict_proteins()` takes a dictionary of protein sequences and returns
`{sequence_name: list[PeptideResult]}`:

```python
proteins = predictor.predict_proteins(
    {"TP53": "MEEPQSDPSVEPPLSQETFS...", "KRAS": "MTEYKLVVVGAGGVGKS..."},
    peptide_lengths=[9, 10],
)

for r in proteins["TP53"]:
    if r.affinity and r.affinity.value < 500:
        print(f"  offset={r.offset} {r.peptide} IC50={r.affinity.value:.0f}")
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

# {sample_name: list[PeptideResult]}
results = ms.predict(["SIINFEKL", "GILGFVFTL"])

# {sample_name: {seq_name: list[PeptideResult]}}
protein_results = ms.predict_proteins({"TP53": "MEEPQ..."})

# flat DataFrames with sample_name column
df = ms.predict_dataframe(["SIINFEKL"])
df = ms.predict_proteins_dataframe({"TP53": "MEEPQ..."})
```

### Measurement kinds and MHC context

Each `Prediction` has a `kind` string describing what it measures:

The canonical prediction kind strings are defined in `mhctools.pred.Kind`.

| Kind | Meaning |
|---|---|
| `pMHC_affinity` | Peptide-MHC binding affinity |
| `pMHC_presentation` | Likelihood of surface presentation (EL/processing) |
| `pMHC_stability` | Peptide-MHC complex stability |
| `pMHC_TCR_binding` | TCR recognition of a peptide-MHC (pMHC:TCR binding) |
| `immunogenicity` | T-cell immunogenicity |
| `antigen_processing` | Combined processing score |
| `proteasome_cleavage` | Proteasomal (MHC-I, cytosolic) C-terminal cleavage score |
| `endolysosomal_cleavage` | Endolysosomal (MHC-II, cathepsin) C-terminal cleavage score |
| `tap_transport` | TAP transport score (reserved, not yet used) |
| `erap_trimming` | ERAP trimming score (reserved, not yet used) |

Predictors also expose `kind_support()` so downstream code can tell what MHC
context is meaningful for each emitted kind:

```python
support = predictor.kind_support()
support["pMHC_affinity"]
# {"mhc_dependence": "single_allele", "mhc_class": "I"}
```

`mhc_dependence` is one of:

| Value | Meaning |
|---|---|
| `none` | The prediction is MHC-independent; `Prediction.allele` is empty. |
| `single_allele` | The prediction is for one peptide/MHC allele pair; `Prediction.allele` is part of the key. |
| `haplotype` | The prediction uses the requested MHC repertoire jointly; `Prediction.allele` may carry best-allele attribution but is not the prediction key. |

`mhc_class` is one of `none`, `I`, `II`, or `both`.

The allowed metadata values are defined in `mhctools.pred` as
`MHC_DEPENDENCE_VALUES` and `MHC_CLASS_VALUES`.

Examples:

| Predictor | Kind | `mhc_dependence` | `mhc_class` |
|---|---|---|---|
| `NetMHCpan41` | `pMHC_affinity` | `single_allele` | `I` |
| `NetMHCpan41` | `pMHC_presentation` | `single_allele` | `I` |
| `NetMHCIIpan4_EL` | `pMHC_presentation` | `single_allele` | `II` |
| `NetMHCstabpan` | `pMHC_stability` | `single_allele` | `I` |
| `MHCflurry` | `pMHC_affinity` | `single_allele` | `I` |
| `MHCflurry` haplotype mode | `pMHC_presentation` | `haplotype` | `I` |
| `MHCflurry` per-allele panel mode | `pMHC_presentation` | `single_allele` | `I` |
| `MHCflurry` | `antigen_processing` | `none` | `none` |
| `Pepsickle` | `proteasome_cleavage` | `none` | `none` |
| `NetCleave_I` | `proteasome_cleavage` | `none` | `I` |
| `NetCleave_II` | `endolysosomal_cleavage` | `none` | `II` |
| `NetTCR` | `pMHC_TCR_binding` | `none` | `I` |
| `Tulip` | `pMHC_TCR_binding` | `single_allele` | `I` |
| `BigMHC_IM` | `immunogenicity` | `single_allele` | `I` |
| `PRIME` | `immunogenicity` | `single_allele` | `I` |

### TCR predictors (`NetTCR`, `Tulip`)

`NetTCR` and `Tulip` predict pMHC:TCR binding — whether a paired αβ T-cell
receptor (an `mhctools.TCR`, described by its CDR loops) recognises a peptide.
Both take `(peptide, TCR)` inputs; `Tulip` additionally takes the presenting
MHC allele.

```python
from mhctools import Tulip, TCR

tcr = TCR(cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF", name="clone1")
predictor = Tulip()                       # needs TULIP_HOME + TULIP_PYTHON
results = predictor.predict(["GILGFVFTL"], [tcr], mhc="HLA-A*02:01")
results[0].preds[0].score                 # higher = more likely binding
```

[TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR) is **GPLv3** and pinned to
`transformers==4.32.1`; mhctools is Apache-2.0 and depends on neither torch nor
transformers. The `Tulip` wrapper therefore vendors none of TULIP — it runs a
user-provided checkout out-of-process, in an isolated interpreter, via TULIP's
own `predict.py`. Set two things up first (see `scripts/setup_tulip_env.sh`,
which does both):

- `TULIP_HOME` — a clone of TULIP-TCR (provides `predict.py`, `src/`, tokenizers,
  and the released `model_weights/`);
- `TULIP_PYTHON` — an isolated **Python 3.11** interpreter with `torch` and
  `transformers==4.32.1` (3.11 so `tokenizers` installs from a prebuilt wheel and
  needs no Rust toolchain).

For MHCflurry presentation, `presentation_allele_mode="haplotype"` treats the
requested alleles as one sample genotype and emits one `pMHC_presentation`
record per peptide. The `allele` field carries MHCflurry's `best_allele`
attribution when available. `presentation_allele_mode="per_allele"` treats each
allele as a separate one-allele synthetic sample and emits one presentation
record per peptide/allele pair. The default `"auto"` mode uses haplotype mode
for up to six alleles and per-allele mode for larger allele panels.

### The Prediction object

Every prediction is a frozen, self-contained `Prediction` dataclass:

```python
from mhctools import Prediction

pred = Prediction(
    kind="pMHC_affinity",
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
| `MHCflurry` | affinity + presentation + processing | `pip install mhcflurry` + `mhcflurry-downloads fetch` |
| `MHCflurry_Affinity` | affinity | `pip install mhcflurry` + `mhcflurry-downloads fetch` |
| `BigMHC` | presentation or immunogenicity | [BigMHC](https://github.com/KarchinLab/bigmhc) clone (set `BIGMHC_DIR`) |
| `MixMHCpred` | presentation | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) |
| `IedbNetMHCpan` / `IedbSMM` / `IedbNetMHCIIpan` | affinity | IEDB web API |
| `RandomBindingPredictor` | affinity | (built-in) |

### Antigen processing

| Predictor | Kinds produced | Requires |
|---|---|---|
| `Pepsickle` | proteasome cleavage | `pip install pepsickle` ([paper](https://doi.org/10.1093/bioinformatics/btab628)) |
| `NetChop` | proteasome cleavage | [NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/) |
| `NetCleave_I` / `NetCleave_II` | proteasomal (I) / endolysosomal (II) C-terminal cleavage | [NetCleave](https://github.com/BSC-CNS-EAPM/NetCleave) clone (set `NETCLEAVE_DIR`) |

`Pepsickle` and `NetChop` use configurable scoring to aggregate per-position
cleavage probabilities into peptide-level scores (see `ProcessingPredictor`
and `ProteasomePredictor`).

`NetCleave` is different: it emits a **single C-terminal cleavage score per
peptide** and covers **both** the MHC-I proteasomal (`NetCleave_I` →
`proteasome_cleavage`) and MHC-II endolysosomal (`NetCleave_II` →
`endolysosomal_cleavage`) pathways — MHC-II processing is otherwise a gap in
the predictor set. It needs the residues downstream of the peptide to build
the cleavage site, so pass `c_flanks` (or scan proteins). Its weights ship in
the git repo; the R dependency in NetCleave's README is only for its training
pipeline, not prediction.

```python
from mhctools import NetCleave_II

predictor = NetCleave_II()                 # resolves NETCLEAVE_DIR / ~/NetCleave
# score peptides with their C-terminal flanking residues (>= 3)
results = predictor.predict(["SIINFEKL"], c_flanks=["DGH"])
results[0].endolysosomal_cleavage.score

# or scan a protein so each peptide is scored in real context
by_protein = predictor.predict_proteins({"TP53": "MEEPQ..."}, peptide_lengths=[15])
```

> ⚠️ NetCleave's own paper reports class-II C-terminal cleavage is a much
> weaker signal than class I (AUC ~0.66 vs ~0.91). Treat
> `endolysosomal_cleavage` scores accordingly.

### Immunogenicity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `BigMHC_IM` | immunogenicity | [BigMHC](https://github.com/KarchinLab/bigmhc) clone (set `BIGMHC_DIR`) |
| `PRIME` | immunogenicity | [PRIME](https://github.com/GfellerLab/PRIME) clone + MixMHCpred |

`PRIME` predicts CD8+ T-cell immunogenicity of class-I peptides by combining
MHC-I binding (via MixMHCpred, which it calls internally) with a TCR-recognition
propensity model. It emits one `immunogenicity` prediction per (peptide, allele):
`score` is the PRIME score (higher = more immunogenic) and `percentile_rank` is
the PRIME %Rank (lower = better). PRIME is academic / non-commercial licensed, so
mhctools shells out to a user-provided install rather than vendoring it.

```python
from mhctools import PRIME

predictor = PRIME(
    alleles=["HLA-A*02:01", "HLA-B*07:02"],
    program_name="PRIME",                    # or an absolute path
    mixmhcpred_path="/path/to/MixMHCpred")    # optional if MixMHCpred is on PATH
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score
```

> ⚠️ Every current CD8 immunogenicity predictor — `PRIME` and `BigMHC_IM`
> included — ranks well in the characterized regime but generalizes poorly to
> truly novel neoepitopes; independent benchmarks put the field near AUC
> 0.5–0.65, with most of the usable signal coming from expression / agretopicity
> features rather than the peptide alone. Use these scores to prioritize, not as
> ground truth.

### TCR specificity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `NetTCR` | pMHC:TCR binding | [NetTCR-2.2](https://github.com/mnielLab/NetTCR-2.2) clone (set `NETTCR_DIR`) + a TFLite runtime (`pip install mhctools[nettcr]`) |

`NetTCR` predicts whether a paired αβ T-cell receptor recognises a
(class-I) peptide. Unlike the MHC-ligand predictors, its input is a peptide
plus a `TCR` (the six CDR loops), not an allele, and it emits the
`pMHC_TCR_binding` kind. NetTCR ships its pretrained weights in its git
repository as small TFLite models; this wrapper runs the pan cross-validation
ensemble in-process and does not need NetTCR's conda environment.

```python
from mhctools import NetTCR, TCR

predictor = NetTCR()   # resolves NETTCR_DIR / ~/NetTCR-2.2
tcr = TCR(
    cdr1a="NSASQS", cdr2a="VYSSG", cdr3a="VVEGDKVI",
    cdr1b="MGHRA", cdr2b="YSYEKL", cdr3b="ASSHSGYEQF", name="clone1")

# Score explicit (peptide, TCR) pairs...
results = predictor.predict_pairs([("LLWNGPMAV", tcr)])
results[0].tcr_binding.score        # ensemble-mean recognition probability

# ...or every peptide x TCR combination.
results = predictor.predict(["LLWNGPMAV", "GILGFVFTL"], [tcr])
```

## Commandline examples

### Prediction for user-supplied peptide sequences

```sh
mhctools --sequence SIINFEKL SIINFEKLQ --mhc-predictor netmhc --mhc-alleles A0201
```

### Automatically extract peptides as subsequences of specified length

```sh
mhctools --sequence AAAQQQSIINFEKL --extract-subsequences --mhc-peptide-lengths 8-10 --mhc-predictor mhcflurry --mhc-alleles A0201
```

### Annotate an existing table with predictor scores (`predict-table`)

Downstream evaluation workflows often start from an annotated benchmark table
(with columns like `sample_id`, `hit`, `peptide`, and per-row genotype/allele
info) and just need external predictor scores appended. `mhctools
predict-table` reads a CSV, runs each requested predictor once, and appends one
score column per predictor — choosing the best allele per row — while
preserving every input column:

```sh
mhctools predict-table \
    --input benchmark.csv.bz2 \
    --peptide-column peptide \
    --alleles-column hla \
    --predictor netmhcpan42-ba:netmhcpan4.2.ba:affinity \
    --predictor netmhcpan42-el:netmhcpan4.2.el:score \
    --out benchmark.with_scores.csv.bz2
```

Each `--predictor` spec is `NAME[:OUTPUT_COLUMN[:FIELD]]`, where `FIELD` is
`affinity`, `score`, or `percentile_rank` (lower is better for `affinity` and
`percentile_rank`; higher for `score`). Rows may hold several alleles per cell
(whitespace-, comma-, or semicolon-separated); the best one per peptide is
chosen and recorded in a `<OUTPUT_COLUMN>_best_allele` provenance column.
Pass `--predictor-info info.csv` to also write a sidecar describing each
column's `score_field` and `higher_is_better`.

The same thing from Python (I/O-free, works on any `DataFrame`):

```python
from mhctools import annotate_table, AnnotationSpec, NetMHCpan42_BA

annotated = annotate_table(
    df,
    [AnnotationSpec(
        predictor=lambda alleles: NetMHCpan42_BA(alleles=alleles),
        output_column="netmhcpan4.2.ba",
        field="affinity")],
    peptide_column="peptide",
    allele_column="hla")
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
preds = collection.to_preds()           # list of Prediction
pp_list = collection.to_peptide_preds() # list of PeptideResult
```
