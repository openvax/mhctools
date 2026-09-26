[![Tests](https://github.com/openvax/mhctools/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/mhctools/actions/workflows/tests.yml)
<a href="https://pypi.python.org/pypi/mhctools/">
<img src="https://img.shields.io/pypi/v/mhctools.svg?maxAge=1000" alt="PyPI" />
</a>

# mhctools

One Python interface to ~30 MHC binding, presentation, immunogenicity, and
antigen-processing predictors.

Each predictor has its own input format, output format, allele spelling, and
installation ritual. mhctools gives them all the same `predict()` call and the
same result objects, so swapping NetMHCpan for MHCflurry is a one-line change
and comparing them is a DataFrame. Everything runs locally.

- [Install and predict](#install-and-predict)
- [Core concepts](#core-concepts)
- [Common tasks](#common-tasks)
- [Predictor catalog](#predictor-catalog)
- [Getting models](#getting-models)
- [Command line](#command-line)
- [Beyond peptide-MHC](#beyond-peptide-mhc)
- [Legacy API](#legacy-api)
- [Further reading](#further-reading)

## Install and predict

```sh
pip install mhctools
```

`Calis` works immediately with no download. Most predictors need model weights
or an external tool first — see [getting models](#getting-models). MHCflurry
ships as a dependency, so it only needs its weights:

```sh
mhctools fetch mhcflurry
```

Then predict:

```python
from mhctools import MHCflurry

predictor = MHCflurry(alleles=["HLA-A*02:01", "HLA-B*07:02"])
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

for r in results:
    if r.affinity:
        print(f"{r.peptide} -> {r.affinity.allele} IC50={r.affinity.value:.1f}nM")
```

That is the whole pattern. Every predictor in the catalog below is constructed
the same way and answers the same `predict()` call.

## Core concepts

**`predict()` returns a list of `PeptideResult`.** Each one carries the peptide
string and gives you accessors for each kind of prediction. An accessor returns
`None` when the predictor doesn't produce that kind — so `r.stability` is
`None` from an affinity-only model, rather than an error or a zero.

**Results preserve input order and repeated peptides.** `results[i]` corresponds
to `peptides[i]`, so you can use `zip(peptides, results)`. The NetMHC family and
the legacy binding-prediction fallback score each distinct peptide once and
expand its predictions into a separate `PeptideResult` for every occurrence.
If those backends omit a requested peptide/allele pair, `predict()` raises
`ValueError` instead of returning a shorter, misaligned list. The standalone
`BindingPredictionCollection.to_peptide_preds()` conversion still groups by
`(peptide, offset, source_sequence_name)` because it has no input list.

```python
r = results[0]

r.peptide                      # "SIINFEKL"
r.offset                       # position in source protein (if scanned)
r.kinds                        # {"pMHC_affinity", "pMHC_presentation", "antigen_processing"}
r.alleles                      # {"HLA-A*02:01", "HLA-B*07:02"}

# best prediction of each kind, or None when the kind is absent
r.affinity
r.presentation
r.stability

if r.affinity:
    r.affinity.value            # IC50 in nM
    r.affinity.percentile_rank  # 0-100, lower = better
    r.affinity.score            # predictor-specific scale, higher = better
    r.affinity.allele           # best allele for this kind

r.best_affinity_by_rank        # by lowest percentile rank instead of score

r.preds                        # tuple of every underlying Prediction
r.filter(kind="pMHC_affinity")
r.filter(allele="HLA-A*02:01")
```

**Underneath, each `PeptideResult` wraps a tuple of `Prediction` objects** —
frozen dataclasses, one per allele-kind combination, each self-contained:

```python
from mhctools import Prediction

pred = Prediction(
    kind="pMHC_affinity",
    score=0.85,           # predictor-specific scale, higher = better
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

Three fields do most of the work, and they mean the same thing everywhere:

- **`score`** — always present, always higher-is-better, but on a
  predictor-specific scale.
- **`value`** — a physical quantity on a linear scale, in that kind's canonical
  unit (nM for affinity, hours for stability). Empty when the kind has no unit,
  or when the wrapper cannot honestly convert to it.
- **`percentile_rank`** — 0-100, lower is stronger, present when the predictor
  scores against a background distribution.

A predictor can emit more than one kind. NetMHCpan 4.1, for example, produces
both `pMHC_affinity` and `pMHC_presentation` for every peptide-allele pair.

📖 Full reference: **[prediction kinds, units, and MHC context](docs/kinds.md)**
— the 17 kinds, what fills `value`, how `MeasurementContext` works, and which
kind each predictor emits.

## Common tasks

### Scan proteins instead of peptides

`predict_proteins()` takes a dictionary of sequences and returns
`{sequence_name: list[PeptideResult]}`, with each result's `offset` set:

```python
proteins = predictor.predict_proteins(
    {"TP53": "MEEPQSDPSVEPPLSQETFS...", "KRAS": "MTEYKLVVVGAGGVGKS..."},
    peptide_lengths=[9, 10],
)

for r in proteins["TP53"]:
    if r.affinity and r.affinity.value < 500:
        print(f"  offset={r.offset} {r.peptide} IC50={r.affinity.value:.0f}")
```

### Get a DataFrame

Every level has a `_dataframe` variant that flattens to a pandas DataFrame:

```python
df = predictor.predict_dataframe(["SIINFEKL"], sample_name="pat001")
df = predictor.predict_proteins_dataframe({"TP53": "MEEPQ..."}, sample_name="pat001")
```

The columns are the same for every predictor, and are defined once as
`mhctools.pred.COLUMNS`:

| Column | |
|---|---|
| `sample_name` | who this was predicted for |
| `peptide`, `n_flank`, `c_flank` | the peptide and its flanking residues |
| `source_sequence_name`, `offset` | where it came from, if scanned |
| `predictor_name`, `predictor_version` | what produced it |
| `allele`, `tcr` | the MHC allele and/or TCR, empty when not applicable |
| `kind`, `score`, `value`, `percentile_rank` | [the prediction itself](docs/kinds.md) |
| `measurement_context`, `peptide_input`, `cache_key` | assay context and exact-input identity |

### Run many samples with different genotypes

```python
from mhctools import MultiSample, MHCflurry

ms = MultiSample(
    samples={
        "pat001": ["HLA-A*02:01", "HLA-B*07:02"],
        "pat002": ["HLA-A*01:01", "HLA-B*08:01"],
    },
    predictor_class=MHCflurry,
)

results = ms.predict(["SIINFEKL", "GILGFVFTL"])       # {sample: [PeptideResult]}
protein_results = ms.predict_proteins({"TP53": "MEEPQ..."})  # {sample: {seq: [...]}}

df = ms.predict_dataframe(["SIINFEKL"])               # flat, with sample_name
df = ms.predict_proteins_dataframe({"TP53": "MEEPQ..."})
```

### Add predictor scores to an existing table

Evaluation workflows often start from an annotated benchmark table — columns
like `sample_id`, `hit`, `peptide`, and a per-row genotype — and just need
scores appended. `annotate_table` is I/O-free and works on any `DataFrame`:

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

There's a [CLI equivalent](#annotate-a-table-predict-table) that reads and
writes CSV.

## Predictor catalog

Start from what you want to know:

| I want to predict… | Kind | Predictors |
|---|---|---|
| Binding affinity to an allele | `pMHC_affinity` | `NetMHCpan`, `NetMHC`, `NetMHCIIpan`, `NetMHCcons`, `MHCflurry`, `CapHLA`, `SMM`, `SMMPMBEC` |
| Surface presentation | `pMHC_presentation` | `NetMHCpan41`/`42`, `NetMHCIIpan`, `MHCflurry`, `CapHLA`, `MixMHCpred` (I), `MixMHC2pred` (II), `BigMHC` |
| How long the pMHC complex lasts | `pMHC_stability` | `NetMHCstabpan` |
| Combined antigen processing | `antigen_processing` | `MHCflurry` |
| Proteasomal cleavage | `proteasome_cleavage` | `Pepsickle`, `NetChop`, `NetCleave_I` |
| Endolysosomal cleavage (class II) | `endolysosomal_cleavage` | `NetCleave_II` |
| TAP transport into the ER | `tap_transport` | `DeepTAP` |
| ERAP1 N-terminal trimming | `erap_trimming` | `ERAMER` |
| Whether a T cell responds | `immunogenicity` | `Calis`, `PRIME`, `BigMHC_IM`, `DeepImmuno`, `TLimmuno2` (II) |
| Whether a specific TCR recognises it | `pMHC_TCR_binding` | `NetTCR`, `Tulip`, `MixTCRpred` |
| How long the free peptide survives | `peptide_half_life` | `PeptiVerse`, `PlifePred2` |
| Which peptidase cuts which bond | — | [cleavage API](docs/cleavage.md) |

📖 **[Predictor reference](docs/predictors.md)** — what each model does, what it
needs installed, and a worked example for every one.

⚠️ **[Known limits](docs/limitations.md)** — several of these models are weaker
than their own papers suggest. Worth reading before you trust a score.

`RandomBindingPredictor` is built in and produces random affinities, which is
occasionally useful as a null baseline.

## Getting models

Most predictors need something downloaded first. There is one command for all
of it:

```sh
mhctools ls                       # what exists, where it lives, who manages it
mhctools fetch mhcflurry          # get it
mhctools predictors               # can it actually run?
```

`fetch` is idempotent, so a provisioning script can call it over a whole list
without special-casing tools it cannot install. Where a predictor has its own
download manager, mhctools reports that manager's path instead of making a
second copy. Academic-licensed tools need an explicit `--accept-license`, and
the DTU NetMHC family needs a license you request from DTU directly.

```python
from mhctools import MHCflurry, fetch, list_artifacts

MHCflurry.fetch()
for artifact in list_artifacts():
    print(artifact.name, artifact.manager, artifact.version, artifact.path)
```

📖 **[Getting models](docs/artifacts.md)** — the full `fetch`/`ls`/`predictors`
reference, cache locations, the four ownership tiers, and licensing.

## Command line

| Command | Does |
|---|---|
| `mhctools` | Predict for peptides or FASTA sequences |
| `mhctools ls` | List model artifacts and where they live |
| `mhctools fetch` | Download model weights and tool snapshots |
| `mhctools predictors` | Report which predictors can actually run |
| `mhctools predict-table` | Append predictor scores to a CSV |
| `mhctools cleavage` | Per-bond peptidase evidence |
| `mhctools vaccine-report` | Route-aware vaccine construct report |
| `mhctools benchmark` | Assay-aware model evaluation |
| `mhctools mixtcrpred` | Score paired TCRs against a fixed target |

### Predict for peptides you supply

```sh
mhctools --sequence SIINFEKL SIINFEKLQ --mhc-predictor netmhc --mhc-alleles A0201
```

`--sequence` may be repeated and all occurrences accumulate. Or use
`--input-peptides-file` for one peptide per line (blank lines ignored), or
`--input-fasta-file` for protein sequences. Pick exactly one of the three.

### Extract subsequences automatically

```sh
mhctools --sequence AAAQQQSIINFEKL --extract-subsequences \
    --mhc-peptide-lengths 8-10 --mhc-predictor mhcflurry --mhc-alleles A0201
```

### Annotate a table (`predict-table`)

Reads a CSV, runs each requested predictor once, and appends one score column
per predictor — choosing the best allele per row — while preserving every input
column:

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
`affinity`, `score`, or `percentile_rank`. Lower is better for `affinity` and
`percentile_rank`, higher for `score`.

A row may hold several alleles per cell (whitespace-, comma-, or
semicolon-separated); the best one per peptide is chosen and recorded in a
`<OUTPUT_COLUMN>_best_allele` provenance column. Missing or blank
peptide/allele cells stay unscored — they are never coerced into a literal
sequence or allele string and sent to a predictor.

Pass `--predictor-info info.csv` to also write a sidecar describing each
column's `score_field`, `units`, and `higher_is_better`. Empty `units` means
the field is dimensionless or predictor-specific.

### Output conventions

CLI prediction tables follow one convention across every predictor. Plain
peptide inputs get an empty `source_sequence_name` and offset `0`, while FASTA
and subsequence inputs keep their source and zero-based offset and are ordered
by those coordinates. `prediction_method_name` is the exact CLI predictor name
you selected, including version and mode. `affinity` is IC50 in nM,
`percentile_rank` is a 0–100 percentile, and `score` stays
predictor-specific. CSV floats are serialized with six significant digits.

Default stdout is streamed as tab-separated values, so an empty source name
survives as an empty field and large tables don't need a second formatted copy
in memory. A downstream closed pipe (`| head`, say) exits cleanly.

## Beyond peptide-MHC

Three features that don't fit the peptide-in, score-out shape:

**[Per-bond peptidase evidence](docs/cleavage.md).** `DPP4qPISA` evaluates the
published human DPP4 N-terminal triplet model locally. `CleavageInput` and
`CleavageResult` preserve terminal chemistry, native scores, assay provenance,
and parent-sequence bond coordinates. The panel also provides motif rules for
CPN, aminopeptidase P, FAP, aminopeptidases A/N, DPP8/9, TPP2,
puromycin-sensitive aminopeptidase, PREP, and ERAP2, each stating how strict it
is (`required`, `preferred`, or `permissive`) together with the source
observation behind that grade — so a non-match can be read for what it is
worth. Three enzymes whose published specificity does not generalize (THOP1,
neurolysin, endosomal IRAP) ship as exact-sequence source references instead:
they report what an experiment found for that precise chemical form and abstain
on anything else.

An optional `eramer-step` model exposes ERAP1's existing length-specific PWM
score at the initial trimming bond. Installed Pepsickle epitope models are
available through the same per-bond contract, with exact weight,
inference-code, and feature-code hashes; their forced endpoint zero is excluded
rather than mislabeled as an internal bond.

```sh
mhctools cleavage --list-models
mhctools cleavage --sequence RPPGFSPFR --model app2-xp --model cpn-basic
```

The [cleavage guide](docs/cleavage.md) also carries a wider candidate inventory
and prioritized follow-up issues. qPISA scores are substrate-depletion
estimates, not serum half-lives or probabilities.

**[Route-aware vaccine reports](docs/vaccine-reports.md).** `mhctools
vaccine-report` reads a structured construct manifest and writes a timestamped
directory containing a sequence-centered PDF, declared processing-route policy,
selected MHC windows, placement assessments, and checksums:

```sh
mhctools vaccine-report \
    --input docs/vaccine-report-example.json \
    --output-dir vaccine-reports
```

Injected SLP and RNA-encoded constructs use different route policies. Cytosolic
RNA emphasizes proteasome/ER/class-I processing and omits free serum-peptide
tracks; an injected SLP begins extracellularly/endolysosomally while retaining
proteasome/TAP cross-presentation as a conditional route.

**[Assay-aware benchmarks](docs/benchmarks.md).** `mhctools benchmark`
evaluates source-linked observations in separate assay, endpoint, and
native-unit strata, reporting training overlap, repeated measurements,
unsupported inputs, and missing target-domain evidence. See also `mhctools
benchmark --lineage-inventory`.

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

The historical `Iedb*` predictor names also still work, as local wrappers — see
[compatibility names](docs/predictors.md#compatibility-names-for-the-old-iedb-predictors).

## Further reading

| Guide | Covers |
|---|---|
| [Prediction kinds](docs/kinds.md) | The 17 kinds, units, `MeasurementContext`, MHC dependence |
| [Predictor reference](docs/predictors.md) | Per-predictor setup and examples |
| [Known limits](docs/limitations.md) | Every caveat, collected and indexed |
| [Getting models](docs/artifacts.md) | `fetch`, `ls`, `predictors`, caches, licensing |
| [Cleavage evidence](docs/cleavage.md) | Per-bond peptidase API and model panel |
| [Vaccine reports](docs/vaccine-reports.md) | Input schema, route policies, Vaxrank integration |
| [Benchmarks](docs/benchmarks.md) | Assay-aware evaluation and training-overlap reporting |
| [Peptide PK and uptake](docs/exposure-results.md) | Half-life, clearance, uptake, tissue exposure |
| [Optional backend conformance](docs/optional-backends.md) | Artifact verification and capability levels |
| [Testing](docs/testing.md) | Running the full suite, optional model setup, CI |

## Development

```sh
./develop.sh    # editable install
./lint.sh       # ruff
./test.sh       # pytest
```

For a complete run with no skipped tests, and for setting up the optional
predictor backends, see the [testing guide](docs/testing.md). Releases are
described in [RELEASING.md](RELEASING.md).
