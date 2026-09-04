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
mhctools fetch mhcflurry
```

## Predictor artifacts

mhctools exposes one acquisition command for model weights, reference files,
and external tool snapshots required by its wrappers. Predictors with their own
download manager keep using it; mhctools reports the manager and resolved path
instead of copying the files into a second cache.

```sh
# Show packaged and optional artifacts, where they live, and who manages them.
mhctools ls

# Fetch the upstream package's default compatible release.
mhctools fetch mhcflurry

# Reproducibility runs may request an explicit artifact release.
mhctools fetch mhcflurry --version 2.2.0

# Fetch a pinned open-source snapshot (code plus the wrapper's model files).
mhctools fetch eramer

# Academic licenses must be reviewed and accepted explicitly.
mhctools fetch nettcr --accept-license

# Machine-readable inventory, optionally rooted somewhere else.
mhctools ls --json
mhctools ls --data-dir /shared/models
```

The same operations are available in Python:

```python
from mhctools import ERAMER, MHCflurry, fetch, list_artifacts

MHCflurry.fetch()
fetch("mhcflurry-affinity")
ERAMER.fetch()
for artifact in list_artifacts():
    print(artifact.name, artifact.manager, artifact.version, artifact.path)
```

`fetch()` obtains every safely and legally downloadable artifact needed by the
named wrapper. It does not install Python packages, execute upstream setup
scripts, or duplicate caches owned by another package.

mhctools-managed snapshots default to the platform's user data directory
(`~/Library/Application Support/mhctools` on macOS,
`~/.local/share/mhctools` on Linux). Set `MHCTOOLS_DATA_DIR`, pass
`--data-dir`, or use the Python `data_dir=` argument to put them on shared or
scratch storage. Every snapshot lives under
`artifacts/<tool>/<git-commit>/` and includes `.mhctools-artifact.json` with
its source repository, exact commit, sparse paths, and license provenance.

The `MANAGER` column distinguishes four ownership models:

- `mhctools package` / `<package> package`: weights shipped in an installed
  Python package;
- `mhcflurry`: MHCflurry's own native download cache;
- `mhctools`: a pinned snapshot fetched into the directory above;
- `user` / `manual`: an existing checkout or licensed executable owned by the
  user. Manual artifacts are listed but `fetch` will not redistribute them.

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
| `tap_transport` | TAP transport / binding score |
| `erap_trimming` | ERAP1 N-terminal trimming score |

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
| `CapHLA` | `pMHC_affinity` | `single_allele` | `both` |
| `CapHLA` | `pMHC_presentation` | `single_allele` | `both` |
| `MixMHC2pred` | `pMHC_presentation` | `single_allele` | `II` |
| `NetMHCstabpan` | `pMHC_stability` | `single_allele` | `I` |
| `MHCflurry` | `pMHC_affinity` | `single_allele` | `I` |
| `MHCflurry` haplotype mode | `pMHC_presentation` | `haplotype` | `I` |
| `MHCflurry` per-allele panel mode | `pMHC_presentation` | `single_allele` | `I` |
| `MHCflurry` | `antigen_processing` | `none` | `none` |
| `Pepsickle` | `proteasome_cleavage` | `none` | `none` |
| `NetCleave_I` | `proteasome_cleavage` | `none` | `I` |
| `NetCleave_II` | `endolysosomal_cleavage` | `none` | `II` |
| `DeepTAP` | `tap_transport` | `none` | `none` |
| `ERAMER` | `erap_trimming` | `none` | `I` |
| `NetTCR` | `pMHC_TCR_binding` | `none` | `I` |
| `Tulip` | `pMHC_TCR_binding` | `single_allele` | `I` |
| `BigMHC_IM` | `immunogenicity` | `single_allele` | `I` |
| `PRIME` | `immunogenicity` | `single_allele` | `I` |
| `DeepImmuno` | `immunogenicity` | `single_allele` | `I` |
| `TLimmuno2` | `immunogenicity` | `single_allele` | `II` |
| `Calis` | `immunogenicity` | `none` | `I` |

### TCR predictors (`NetTCR`, `Tulip`)

`NetTCR` and `Tulip` predict pMHC:TCR binding — whether a paired αβ T-cell
receptor (an `mhctools.TCR`, described by its CDR loops) recognises a peptide.
Both take `(peptide, TCR)` inputs; `Tulip` additionally takes the presenting
MHC allele.

```python
from mhctools import Tulip, TCR

tcr = TCR(cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF", name="clone1")
Tulip.fetch()                             # pinned code, tokenizers, and weights
predictor = Tulip()                       # also needs a TULIP-capable Python
results = predictor.predict(["GILGFVFTL"], [tcr], mhc="HLA-A*02:01")
results[0].preds[0].score                 # higher = more likely binding
```

[TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR) is **GPLv3** and pinned to
`transformers==4.32.1`; mhctools is Apache-2.0 and depends on neither torch nor
transformers. The `Tulip` wrapper therefore vendors none of TULIP — it runs an
upstream checkout out-of-process, in an isolated interpreter, via TULIP's own
`predict.py`. `mhctools fetch tulip` obtains the tested code, tokenizers, and
weights; `scripts/setup_tulip_env.sh` can build the separate runtime. You may
instead provide your own checkout and interpreter:

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
| `MHCflurry` | affinity + presentation + processing | `pip install mhcflurry` + `mhctools fetch mhcflurry` |
| `MHCflurry_Affinity` | affinity | `pip install mhcflurry` + `mhctools fetch mhcflurry-affinity` |
| `BigMHC` | presentation or immunogenicity | `mhctools fetch bigmhc --accept-license` + PyTorch, or set `BIGMHC_DIR` |
| `CapHLA` / `CapHLA_EL` / `CapHLA_BA` | presentation + affinity (class I and II) | `pip install "mhctools[caphla]"` + `mhctools fetch caphla` |
| `MixMHCpred` | presentation (class I) | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) |
| `MixMHC2pred` | presentation (class II) | [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred) release (has `PWMdef/`) |
| `IedbNetMHCpan` / `IedbSMM` / `IedbNetMHCIIpan` | affinity | IEDB web API |
| `RandomBindingPredictor` | affinity | (built-in) |

`CapHLA` is a 2025 MIT-licensed PyTorch model family ([Chang & Wu,
*Briefings in Bioinformatics*](https://doi.org/10.1093/bib/bbae595)) covering
human and mouse MHC class I and II, with peptides from 7–25 residues. The default wrapper emits
both outputs for every peptide/allele pair: EL `presentation_score` as
`pMHC_presentation`, and the BA normalized score as `pMHC_affinity`. For BA,
mhctools also inverts CapHLA's training transform to provide predicted IC50 nM
in `value`. Upstream provides neither percentile ranks nor binder thresholds,
so the wrapper does not invent them. `CapHLA_EL` and `CapHLA_BA` load only the
five-fold ensemble they need.

```sh
pip install "mhctools[caphla]"
mhctools fetch caphla
```

```python
from mhctools import CapHLA

predictor = CapHLA(alleles=[
    "HLA-A*02:01",
    "HLA-DPA1*01:03-DPB1*04:01",
])
results = predictor.predict(["GILGFVFTL", "GELIGTLNAAKVPAD"])
results[0].presentation.score
results[0].affinity.score
results[0].affinity.value       # predicted IC50, nM

# Explicit pairs preserve order and duplicates without a cross product.
paired = predictor.predict_pairs([
    ("GILGFVFTL", "HLA-A*02:01"),
    ("GELIGTLNAAKVPAD", "HLA-DPA1*01:03-DPB1*04:01"),
])
```

The wrapper loads the pinned upstream model definitions and weights unchanged,
batches inference deterministically in-process, and preserves canonical
mhcgnomes allele identity in its outputs. CapHLA performance numbers are
author-reported; treat it as a complementary research predictor rather than a
default or independent validation.

`MixMHCpred` 3.0 predicts **class-I presentation** for peptides of length
8-14. Version 3.0 adds pan-allele inference, MHC-I sequence alignment and
sequence-driven prediction, and optional binding-motif/peptide-length plots.
mhctools exposes all per-allele scores and percentile ranks through the
canonical prediction API. `predict_detailed` also retains MixMHCpred's raw
`Score_bestAllele`, `BestAllele`, and `%Rank_bestAllele` columns plus each
allele's closest training allele, sequence distance, and pan-allele status.

MixMHCpred 3.0 is licensed for academic, non-commercial research and prohibits
redistribution without written permission, so its approximately 200 MB of
code, models, and reference data are not included in mhctools. Review the
[upstream license and installation guide](https://github.com/GfellerLab/MixMHCpred)
before downloading the official tagged release:

```sh
git clone --branch v3.0 --depth 1 \
    https://github.com/GfellerLab/MixMHCpred.git
chmod +x MixMHCpred/MixMHCpred
export MIXMHCPRED_PATH="$PWD/MixMHCpred"
pip install "mhctools[mixmhcpred]"
```

The `mixmhcpred` extra installs the upstream Python dependencies. Sequence
alignment additionally needs the `mafft` executable. The upstream
`install_packages` script is another way to install both sets of dependencies.

```python
from mhctools import MixMHCpred

predictor = MixMHCpred(
    alleles=["HLA-A*02:01", "HLA-A*01:02"],  # A*01:02 uses v3 pan inference
)

# Canonical mhctools output: one pMHC_presentation Prediction per allele.
results = predictor.predict(["SIINFEKL"])
results[0].presentation.score

# Complete native output and v3 quality/provenance metadata.
detailed = predictor.predict_detailed(["SIINFEKL"])
detailed.table[["Score_bestAllele", "BestAllele", "%Rank_bestAllele"]]
detailed.allele_info[1].closest_training_allele
detailed.allele_info[1].distance
detailed.allele_info[1].pan_allele

# Retain Binding_predictions.txt, PWM/PLD files and images, and the HTML view.
motifs = predictor.predict_detailed(
    ["SIINFEKL"], output_dir="mixmhcpred-output", output_motifs=True)
motifs.artifacts.files

# Align novel MHC-I sequences, then optionally predict and render their motifs.
sequence_result = predictor.predict_allele_sequences(
    "unaligned-mhc-i.fasta",
    peptides=["SIINFEKL"],
    output_dir="mixmhcpred-sequence-output",
    output_motifs=True,
)
sequence_result.aligned_sequences
sequence_result.table
sequence_result.allele_info[0].closest_database_allele
sequence_result.artifacts.files
```

Both artifact APIs require a new output path: the wrapper refuses an existing
path because MixMHCpred itself deletes and recreates its output directory.
`exclude_peptides_with_cysteine=True` is implemented by mhctools before the
external call, including under v3.0 where the legacy `-c` option was removed.

`MixMHC2pred` is a pan-allele **class-II** presentation predictor and a strong
complement to `NetMHCIIpan` (independently co-best in the Frontiers in
Immunology 2024 class-II benchmark). It emits one `pMHC_presentation`
prediction per (peptide, allele): `score` is the raw MixMHC2pred score (higher
= better), `percentile_rank` is its %Rank (lower = better). It's academic /
non-commercial licensed, so mhctools shells out to a user-provided install
(download a **release**, not a bare clone — the release ships the `PWMdef/`
allele definitions). Alleles may be given in the usual spellings
(`HLA-DRB1*15:01`) or MixMHC2pred's own (`DRB1_15_01`,
`DQA1_01_02__DQB1_06_02`).

```python
from mhctools import MixMHC2pred

predictor = MixMHC2pred(
    alleles=["HLA-DRB1*15:01", "HLA-DQA1*01:02-DQB1*06:02"],
    program_name="/path/to/MixMHC2pred_unix")   # MixMHC2pred on macOS
results = predictor.predict(["GELIGTLNAAKVPAD"])   # class-II length peptides
results[0].presentation.score
```

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

### TAP transport

| Predictor | Kinds produced | Requires |
|---|---|---|
| `DeepTAP` | TAP transport (`tap_transport`) | `mhctools fetch deeptap` + a DeepTAP-capable Python |

TAP (transporter associated with antigen processing) is the step that shuttles
cytosolic peptides into the ER for MHC-I loading — a distinct part of the
processing pathway from proteasomal cleavage, and otherwise a gap in the
predictor set. `DeepTAP` is a BiGRU that scores each peptide once
(**allele-independent**, like the cleavage predictors), emitting one
`tap_transport` prediction per peptide with an empty `allele`. `score` is in
0-1 (higher = stronger TAP binding); in `task_type="reg"` mode the predicted
affinity in nM is also surfaced as `value` (lower = stronger).

DeepTAP ships its weights in-repo and is Apache-2.0, but pins an old
`pytorch-lightning`, so mhctools shells out to DeepTAP's own CLI in a
separate interpreter (the checkpoints load fine under modern Lightning too).
Run `mhctools fetch deeptap`; if the current interpreter lacks torch, set
`DEEPTAP_PYTHON` to one that has it. `DEEPTAP_HOME` can still select a manual
checkout.

```python
from mhctools import DeepTAP

DeepTAP.fetch()
predictor = DeepTAP(task_type="cla")       # resolves DEEPTAP_HOME / ~/DeepTAP
results = predictor.predict(["SIINFEKL", "AEASAAAAY"])
results[1].tap_transport.score             # 0-1, higher = stronger TAP binding
```

> ⚠️ DeepTAP's evaluation is self-reported, and no independent TAP benchmark
> exists for any tool (true of the whole TAP field). Treat the score as a useful
> pathway signal for prioritization, not a validated oracle.

### ERAP1 trimming

| Predictor | Kinds produced | Requires |
|---|---|---|
| `ERAMER` | ERAP1 trimming (`erap_trimming`) | `mhctools fetch eramer` + `openpyxl` |

ERAP1 trims the N-termini of 9–16mer precursor peptides in the ER down to the
8–10mers MHC-I presents — the step between TAP transport and MHC loading, and
otherwise the last empty stage in the pathway. `ERAMER` scores a precursor by
averaging a per-length position-weight-matrix specificity over each residue
trimmed off as it is cut toward a target epitope length (allele-independent, one
`erap_trimming` prediction per peptide; `score` roughly −1…1, higher = more
likely trimmed).

ERAMER is **GPLv3** and its PWM ships in a GPL-licensed `PWM.xlsx`, so mhctools
vendors neither: this is a clean-room Python-3 reimplementation of the
(Python-2.7) tool's trimming-cascade average that loads the PWM from a
upstream ERAMER checkout at runtime. Run `mhctools fetch eramer`, or point at a
manual clone with `ERAMER_HOME`.

```python
from mhctools import ERAMER

ERAMER.fetch()
predictor = ERAMER(epitope_length=8)       # resolves ERAMER_HOME / ~/ERAMER
results = predictor.predict(["GGGGGVVVVVVAAAEE"])   # a 9-16mer precursor
results[0].erap_trimming.score
```

> ⚠️ ERAMER's evaluation is self-reported and ERAP1 trimming is an intrinsically
> noisy signal; treat the score as a pathway prior, not a validated oracle.

### Immunogenicity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `Calis` | immunogenicity | nothing — self-contained |
| `BigMHC_IM` | immunogenicity | `mhctools fetch bigmhc --accept-license` + PyTorch, or set `BIGMHC_DIR` |
| `PRIME` | immunogenicity | [PRIME](https://github.com/GfellerLab/PRIME) clone + MixMHCpred |
| `DeepImmuno` | immunogenicity | `mhctools fetch deepimmuno` + a TensorFlow/Keras-2-capable Python |
| `TLimmuno2` | immunogenicity (class II) | [TLimmuno2](https://github.com/XSLiuLab/TLimmuno2) clone (set `TLIMMUNO2_HOME`) |

`Calis` is the classic sequence-only IEDB class-I immunogenicity model (Calis et
al. 2013): a fixed per-amino-acid log-enrichment scale weighted by per-position
importance, with the anchor positions (P1/P2/C-terminus) masked out. It needs
**no external install and no downloaded weights** — the ~30 published parameters
(from the open-access CC-BY paper) are built in — so it is a fast,
dependency-free, allele-independent baseline. It emits one `immunogenicity`
prediction per peptide (empty `allele`); `score > 0` leans immunogenic.

```python
from mhctools import Calis

predictor = Calis()
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score            # 0.30484 (higher = more immunogenic)
```

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

`DeepImmuno` predicts class-I CD8+ immunogenicity from the peptide and its
HLA-A/B/C allele with a small CNN (Li et al. 2021). It scores **9- and 10-mers
only** and supports a fixed set of ~62 alleles, snapping anything else to the
nearest it knows. It emits one `immunogenicity` prediction per (peptide,
allele); `score` is in 0–1 (higher = more immunogenic). DeepImmuno ships its
weights in-repo and is MIT-licensed, but its script loads them with an old
Keras 2 / TensorFlow stack, so mhctools shells out to DeepImmuno's own CLI in a
separate checkout. Run `mhctools fetch deepimmuno`, or point at a manual clone
with `DEEPIMMUNO_HOME`, and set
`DEEPIMMUNO_PYTHON` to an interpreter that has TensorFlow (with Keras 2, or
newer TensorFlow plus the `tf-keras` shim — the wrapper sets
`TF_USE_LEGACY_KERAS=1` for the subprocess).

```python
from mhctools import DeepImmuno

DeepImmuno.fetch()
predictor = DeepImmuno(alleles=["HLA-A*02:01"])   # resolves DEEPIMMUNO_HOME / ~/DeepImmuno
results = predictor.predict(["NLVPMVATV", "GILGFVFTL"])
results[0].immunogenicity.score                   # 0.9568 (higher = more immunogenic)
```

> ⚠️ Every current CD8 immunogenicity predictor — `PRIME`, `BigMHC_IM`, and
> `DeepImmuno` included — ranks well in the characterized regime but generalizes poorly to
> truly novel neoepitopes; independent benchmarks put the field near AUC
> 0.5–0.65 on unseen tumor neoepitopes (ITSNdb ~0.52–0.60, ICERFIRE ~0.56,
> IMPROVE ~0.60). In the one neutral head-to-head that scored both (NeoaPred,
> *Bioinformatics* 2024), **`BigMHC_IM` edged `PRIME` on cancer neoepitopes**,
> while PRIME tends to do better on viral / infectious-disease epitopes — its
> training positives are mostly viral and cancer-testis antigens, with only
> ~129 (v1) / ~596 (v2) true immunogenic neoepitopes. PRIME's higher
> self-reported numbers are partly attributable to documented train/test
> overlap (IMPROVE flagged ~70% overlap with its evaluation set). Use these
> scores to prioritize, not as ground truth.

`TLimmuno2` is the odd one out: it predicts **class-II (CD4+)** immunogenicity —
the only class-II immunogenicity model here, filling a gap the class-I models
(`Calis`, `PRIME`, `BigMHC_IM`, `DeepImmuno`) leave. It scores a peptide against a class-II
allele (transfer-learned from class-II binding) and emits one `immunogenicity`
prediction per (peptide, allele): `score` in 0–1 (higher = more immunogenic)
and `percentile_rank` from its %Rank against a background set, rescaled to
0–100 (lower = more immunogenic). Native NetMHCIIpan-style keys (`DRB1_0803`,
`HLA-DPA10103-DPB10101`) pass through; common DR forms (`HLA-DRB1*08:03`) are
converted; anything TLimmuno2 does not know raises. Its upstream license is
ambiguous (an Apache-2.0 README badge, no LICENSE file), so mhctools does not
vendor it — it shells out to a user-provided checkout (`TLIMMUNO2_HOME`), with
`TLIMMUNO2_PYTHON` naming an interpreter that has TensorFlow (Keras 2, or newer
TensorFlow plus `tf-keras`).

```python
from mhctools import TLimmuno2

predictor = TLimmuno2(alleles=["DRB1_0803"])       # resolves TLIMMUNO2_HOME / ~/TLimmuno2
results = predictor.predict(["FHTMWHVTRGAVLMY"])
results[0].immunogenicity.score                    # 0.9874 (higher = more immunogenic)
```

> ⚠️ TLimmuno2's %Rank is computed against ~90,000 background peptides **per
> distinct allele**, so a call costs about a minute per allele regardless of how
> many peptides you pass — batch peptides by allele. Class-II immunogenicity is
> noisier than class-I; a prioritization aid, not ground truth.

### TCR specificity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `NetTCR` | pMHC:TCR binding | `mhctools fetch nettcr --accept-license` + a TFLite runtime (`pip install mhctools[nettcr]`) |

`NetTCR` predicts whether a paired αβ T-cell receptor recognises a
(class-I) peptide. Unlike the MHC-ligand predictors, its input is a peptide
plus a `TCR` (the six CDR loops), not an allele, and it emits the
`pMHC_TCR_binding` kind. NetTCR ships its pretrained weights in its git
repository as small TFLite models; this wrapper runs the pan cross-validation
ensemble in-process and does not need NetTCR's conda environment.

```python
from mhctools import NetTCR, TCR

NetTCR.fetch(accept_license=True)  # downloads only the ~8 MB pan ensemble
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
