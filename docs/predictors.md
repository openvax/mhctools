# Predictor reference

Per-predictor notes: what each model does, what it needs installed, and how to
call it. For the master list of every predictor in one table, see the
[README catalog](../README.md#predictor-catalog). For what the output fields
mean, see [prediction kinds](kinds.md). For the known limits of these models
collected in one place, see [limitations](limitations.md).

- [MHC binding and presentation](#mhc-binding-and-presentation)
  - [NetMHC family](#netmhc-family)
  - [Compatibility names for the old IEDB predictors](#compatibility-names-for-the-old-iedb-predictors)
  - [MHCflurry](#mhcflurry)
  - [CapHLA](#caphla)
  - [MixMHCpred](#mixmhcpred)
  - [MixMHC2pred](#mixmhc2pred)
- [Antigen processing](#antigen-processing)
  - [Pepsickle and NetChop](#pepsickle-and-netchop)
  - [NetCleave](#netcleave)
- [TAP transport](#tap-transport)
  - [DeepTAP](#deeptap)
- [ERAP1 trimming](#erap1-trimming)
  - [ERAMER](#eramer)
- [Peptide half-life](#peptide-half-life)
  - [PeptiVerse](#peptiverse)
  - [PlifePred2](#plifepred2)
- [Immunogenicity](#immunogenicity)
  - [Calis](#calis)
  - [PRIME](#prime)
  - [DeepImmuno](#deepimmuno)
  - [TLimmuno2](#tlimmuno2)
- [TCR specificity](#tcr-specificity)
  - [NetTCR](#nettcr)
  - [Tulip](#tulip)
  - [MixTCRpred](#mixtcrpred)

## MHC binding and presentation

| Predictor | Kinds produced | Requires |
|---|---|---|
| `NetMHCpan` / `NetMHCpan41` / `NetMHCpan42` | affinity + presentation | [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) |
| `NetMHCpan4` | affinity or presentation | NetMHCpan 4.0 |
| `NetMHCpan3` / `NetMHCpan28` | affinity | older NetMHCpan |
| `NetMHC` / `NetMHC3` / `NetMHC4` | affinity | [NetMHC](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) |
| `NetMHCIIpan` / `NetMHCIIpan43` | affinity or presentation | [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) |
| `NetMHCcons` | affinity | [NetMHCcons](https://services.healthtech.dtu.dk/services/NetMHCcons-1.1/) |
| `NetMHCstabpan` | stability | [NetMHCstabpan](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/) |
| `MHCflurry` | affinity + presentation + processing | `mhctools fetch mhcflurry` (the package is a dependency) |
| `MHCflurry_Affinity` | affinity | `mhctools fetch mhcflurry-affinity` |
| `BigMHC` | presentation or immunogenicity | `mhctools fetch bigmhc --accept-license` + PyTorch, or set `BIGMHC_DIR` |
| `CapHLA` / `CapHLA_EL` / `CapHLA_BA` | presentation + affinity (class I and II) | `pip install "mhctools[caphla]"` + `mhctools fetch caphla` |
| `MixMHCpred` | presentation (class I) | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) |
| `MixMHC2pred` | presentation (class II) | [MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred) release (has `PWMdef/`) |
| `SMM` / `SMMPMBEC` | affinity | [Local IEDB standalone tools](testing.md#local-smm-and-smm-pmbec) |
| `RandomBindingPredictor` | affinity | (built-in) |

### NetMHC family

NetMHCpan 4.1 emits both `pMHC_affinity` and `pMHC_presentation` for every
peptide-allele pair. The other family members are listed in the table above.

These are DTU tools under identity-bound academic licenses, so mhctools calls a
licensed installation you provide rather than fetching one. See
[licensing](artifacts.md#licensing).

### Compatibility names for the old IEDB predictors

Since 3.44.55, every predictor runs locally. The historical Python names
`IedbNetMHCpan`, `IedbNetMHCcons`, `IedbNetMHCIIpan`, `IedbSMM`, and
`IedbSMM_PMBEC` (and their `*-iedb` CLI names) remain as local compatibility
wrappers. They require installed NetMHCpan **4.1 BA**, NetMHCcons, NetMHCIIpan
**4.3 BA**, SMM, or SMM-PMBEC respectively.

A few differences are worth knowing before you rely on them:

- The class-II alias now uses 4.3 rather than the hosted service's 4.1 default.
  Local versions and percentile calibration can produce different values.
- The class-I compatibility names keep their original default windows of 8–11
  residues; the canonical local classes keep their own defaults.
- `IedbNetMHCpan` emits affinity only, while `NetMHCpan41_BA.predict()` can emit
  both affinity and presentation.
- Affinity is still IC50 in nM.

Prefer the explicit local predictor names when you are recording model
provenance.

There is no HTTP fallback. The HTTP-only `url`, `request_timeout`, and
`raise_on_error` constructor arguments and the CLI `--do-not-raise-on-error`
option have been removed. Missing installations and unsupported inputs raise
errors, so predictions are never silently dropped by an IEDB error policy. For
the standalone matrix methods, use the CLI names `smm` and `smm-pmbec`.

### MHCflurry

`presentation_allele_mode` controls how the requested alleles are interpreted:

- `"haplotype"` treats them as one sample genotype and emits one
  `pMHC_presentation` record per peptide. The `allele` field carries
  MHCflurry's `best_allele` attribution when available.
- `"per_allele"` treats each allele as a separate one-allele synthetic sample
  and emits one presentation record per peptide/allele pair.
- `"auto"` (the default) uses haplotype mode for up to six alleles and
  per-allele mode for larger panels.

MHCflurry predictions carry both the Python package and the official
model-release identity in `predictor_version` — for example
`2.2.1+release-2.2.0`. The version is captured when weights are loaded and
retained with the cached model object. `mhcflurry_composite_version()` exposes
the same rule publicly; it checks the selected directory, including environment
overrides, against the official bundle path. This is release provenance, not a
checksum of the weights.

For custom paths or injected predictors, pass `predictor_version="my-model-id"`
to `MHCflurry` or `MHCflurry_Affinity` if the predictions need a cacheable
identity. Otherwise they stay unversioned rather than being mislabeled as the
active default release. The modern prediction and DataFrame APIs retain the
version; legacy `BindingPrediction` objects keep their original unversioned
schema.

### CapHLA

`CapHLA` is a 2025 MIT-licensed PyTorch model family ([Chang & Wu, *Briefings
in Bioinformatics*](https://doi.org/10.1093/bib/bbae595)) covering human and
mouse MHC class I and II, with peptides from 7–25 residues.

The default wrapper emits both outputs for every peptide/allele pair: the EL
`presentation_score` as `pMHC_presentation`, and the BA normalized score as
`pMHC_affinity`. For BA, mhctools also inverts CapHLA's training transform to
provide predicted IC50 nM in `value`. Upstream provides neither percentile
ranks nor binder thresholds, so the wrapper does not invent them. `CapHLA_EL`
and `CapHLA_BA` load only the five-fold ensemble they need.

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
mhcgnomes allele identity in its outputs.

> ⚠️ CapHLA performance numbers are author-reported. Treat it as a
> complementary research predictor rather than a default or an independent
> validation.

### MixMHCpred

`MixMHCpred` 3.0 predicts **class-I presentation** for peptides of length 8-14.
Version 3.0 adds pan-allele inference, MHC-I sequence alignment and
sequence-driven prediction, and optional binding-motif/peptide-length plots.

mhctools exposes all per-allele scores and percentile ranks through the
canonical prediction API. `predict_detailed` additionally retains MixMHCpred's
raw `Score_bestAllele`, `BestAllele`, and `%Rank_bestAllele` columns plus each
allele's closest training allele, sequence distance, and pan-allele status.

MixMHCpred 3.0 is licensed for academic, non-commercial research and prohibits
redistribution without written permission, so its roughly 200 MB of code,
models, and reference data are not included in mhctools. Review the
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

### MixMHC2pred

`MixMHC2pred` is a pan-allele **class-II** presentation predictor and a strong
complement to `NetMHCIIpan` — the two were independently co-best in the
*Frontiers in Immunology* 2024 class-II benchmark.

It emits one `pMHC_presentation` prediction per (peptide, allele): `score` is
the raw MixMHC2pred score (higher = better) and `percentile_rank` is its %Rank
(lower = better).

It is academic / non-commercial licensed, so mhctools shells out to an install
you provide. Download a **release**, not a bare clone — the release ships the
`PWMdef/` allele definitions. Alleles may be given in the usual spellings
(`HLA-DRB1*15:01`) or in MixMHC2pred's own (`DRB1_15_01`,
`DQA1_01_02__DQB1_06_02`).

```python
from mhctools import MixMHC2pred

predictor = MixMHC2pred(
    alleles=["HLA-DRB1*15:01", "HLA-DQA1*01:02-DQB1*06:02"],
    program_name="/path/to/MixMHC2pred_unix")   # MixMHC2pred on macOS
results = predictor.predict(["GELIGTLNAAKVPAD"])   # class-II length peptides
results[0].presentation.score
```

## Antigen processing

| Predictor | Kinds produced | Requires |
|---|---|---|
| `Pepsickle` | proteasome cleavage | `pip install pepsickle` ([paper](https://doi.org/10.1093/bioinformatics/btab628)) |
| `NetChop` | proteasome cleavage | [NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/) |
| `NetCleave_I` / `NetCleave_II` | proteasomal (I) / endolysosomal (II) C-terminal cleavage | `mhctools fetch netcleave --accept-license`, or your own clone via `NETCLEAVE_DIR` |

### Pepsickle and NetChop

Both use configurable scoring to aggregate per-position cleavage probabilities
into peptide-level scores (see `ProcessingPredictor` and `ProteasomePredictor`).

NetChop 3.1 ships as 32-bit x86 Linux binaries. On macOS and ARM Linux,
`NetChop` automatically runs a user-supplied licensed installation in a
digest-pinned compatibility container. Set `NETCHOP_HOME` to the directory
containing `bin/netChop` (or set `NETMHC_BUNDLE_HOME` to its parent bundle) and
preload the runtime image once:

```sh
docker pull --platform linux/386 \
  i386/debian@sha256:75efd55b326373cf69989912388c0d50c5390638af7378d2fedc3aeb9d100e46
```

Inference runs with Docker network access disabled and image pulling forbidden.
The licensed NetChop files and input directory are mounted read-only. Use
`NetChop(execution="native")` or
`NetChop(execution="container", netchop_dir="/path/to/netchop-3.1")` to select a
backend explicitly.

### NetCleave

`NetCleave` works differently from the other two: it emits a **single
C-terminal cleavage score per peptide**, and it covers **both** the MHC-I
proteasomal (`NetCleave_I` → `proteasome_cleavage`) and MHC-II endolysosomal
(`NetCleave_II` → `endolysosomal_cleavage`) pathways. MHC-II processing is
otherwise a gap in the predictor set.

It needs the residues downstream of the peptide to build the cleavage site, so
pass `c_flanks` or scan proteins. Its weights ship in the git repo; the R
dependency mentioned in NetCleave's README is only for its training pipeline,
not for prediction.

```python
from mhctools import NetCleave_II

predictor = NetCleave_II()   # NETCLEAVE_DIR, ~/NetCleave, ~/code/NetCleave, then snapshot
# score peptides with their C-terminal flanking residues (>= 3)
results = predictor.predict(["SIINFEKL"], c_flanks=["DGH"])
results[0].endolysosomal_cleavage.score

# or scan a protein so each peptide is scored in real context
by_protein = predictor.predict_proteins({"TP53": "MEEPQ..."}, peptide_lengths=[15])
```

`mhctools fetch netcleave --accept-license` installs a pinned snapshot of about
10 MB: the entry script, `predictor/`, and `data/models/`. It skips the ~118 MB
of IEDB and UniParc databases that only upstream's `--generate`/`--train` paths
use. Upstream publishes no license file, so here the flag acknowledges that you
have confirmed your own use is authorized rather than accepting stated terms;
the recorded manifest says `"license": "none published"`. A checkout you manage
yourself still takes precedence.

> ⚠️ NetCleave's own paper reports that class-II C-terminal cleavage is a much
> weaker signal than class I (AUC ~0.66 vs ~0.91). Treat
> `endolysosomal_cleavage` scores accordingly.

## TAP transport

| Predictor | Kinds produced | Requires |
|---|---|---|
| `DeepTAP` | TAP transport (`tap_transport`) | `mhctools fetch deeptap` + a DeepTAP-capable Python |

### DeepTAP

TAP (transporter associated with antigen processing) shuttles cytosolic
peptides into the ER for MHC-I loading. It is a distinct step from proteasomal
cleavage, and otherwise a gap in the predictor set.

`DeepTAP` is a BiGRU that scores each peptide once — **allele-independent**,
like the cleavage predictors — emitting one `tap_transport` prediction per
peptide with an empty `allele`. `score` is in 0-1 (higher = stronger TAP
binding); in `task_type="reg"` mode the predicted affinity in nM is also
surfaced as `value` (lower = stronger).

DeepTAP ships its weights in-repo and is Apache-2.0, but pins an old
`pytorch-lightning`, so mhctools shells out to DeepTAP's own CLI in a separate
interpreter. (The checkpoints load fine under modern Lightning too.) Run
`mhctools fetch deeptap`; if the current interpreter lacks torch, set
`DEEPTAP_PYTHON` to one that has it. `DEEPTAP_HOME` selects a manual checkout.

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

## ERAP1 trimming

| Predictor | Kinds produced | Requires |
|---|---|---|
| `ERAMER` | ERAP1 trimming (`erap_trimming`) | `mhctools fetch eramer` + `openpyxl` |

### ERAMER

ERAP1 trims the N-termini of 9–16mer precursor peptides in the ER down to the
8–10mers MHC-I presents — the step between TAP transport and MHC loading, and
otherwise the last empty stage in the pathway.

`ERAMER` scores a precursor by averaging a per-length position-weight-matrix
specificity over each residue trimmed off as it is cut toward a target epitope
length. It is allele-independent, emitting one `erap_trimming` prediction per
peptide, with `score` roughly −1…1 (higher = more likely trimmed).

ERAMER is **GPLv3** and its PWM ships in a GPL-licensed `PWM.xlsx`, so mhctools
vendors neither. This is a clean-room Python-3 reimplementation of the
(Python-2.7) tool's trimming-cascade average, loading the PWM from an upstream
ERAMER checkout at runtime. Run `mhctools fetch eramer`, or point at a manual
clone with `ERAMER_HOME`.

```python
from mhctools import ERAMER

ERAMER.fetch()
predictor = ERAMER(epitope_length=8)       # resolves ERAMER_HOME / ~/ERAMER
results = predictor.predict(["GGGGGVVVVVVAAAEE"])   # a 9-16mer precursor
results[0].erap_trimming.score
```

> ⚠️ ERAMER's evaluation is self-reported and ERAP1 trimming is an intrinsically
> noisy signal; treat the score as a pathway prior, not a validated oracle.

## Peptide half-life

| Predictor | Kinds produced | Requires |
|---|---|---|
| `PeptiVerse` | Parent-peptide half-life in human serum | pinned PeptiVerse + ESM2 snapshots (`PEPTIVERSE_HOME`, `PEPTIVERSE_ESM_HOME`) + a torch/transformers Python |
| `PlifePred2` | Parent-peptide half-life, matrix unknown ⚠️ | `plifepred2==1.0` (`PLIFEPRED2_HOME`) + pinned Pfeature (`PFEATURE_HOME`) |

`peptide_half_life` records how long the parent peptide persists, in hours. Its
context distinguishes a defined solution, serum/plasma/whole blood, cellular
compartments, and systemic in-vivo PK without multiplying kind strings.

It is deliberately a separate kind from `pMHC_stability`, which is the
dissociation half-life of an assembled peptide-MHC complex — a different
molecule in a different assay — and from the cleavage kinds, which are
site-resolved and intracellular. `MeasurementContext` preserves the matrix,
compartment, analyte, and systemic scope when they are known.

### PeptiVerse

`PeptiVerse` wraps one endpoint of the upstream multi-property platform. Its
dependencies (torch, `transformers==4.46.0`, xgboost, lightning, and ESM2) stay
out of the mhctools environment: inference runs offline in a subprocess under
`PEPTIVERSE_PYTHON`. Provision the exact snapshots before prediction:

```bash
git clone https://huggingface.co/ChatterjeeLab/PeptiVerse
git -C PeptiVerse checkout 8cf0b21dae356278ae96b414a088e4360357d16c
huggingface-cli download facebook/esm2_t33_650M_UR50D \
  --revision 08e4846e537177426273712802403f7ba8261b6c \
  --include config.json tokenizer_config.json special_tokens_map.json vocab.txt model.safetensors \
  --local-dir /models/esm2_t33_650M_UR50D
export PEPTIVERSE_HOME="$PWD/PeptiVerse"
export PEPTIVERSE_ESM_HOME=/models/esm2_t33_650M_UR50D
```

```python
from mhctools import PeptideContext, PeptideInput, PeptiVerse

predictor = PeptiVerse(device="cpu")       # resolves PEPTIVERSE_HOME / ~/PeptiVerse
exact_input = PeptideInput(
    "SIINFEKL",
    occurrence_id="sample-1:occurrence-2",
    context=PeptideContext(matrix="serum", assay_species="Homo sapiens"),
)
results = predictor.predict([exact_input, "KLGGALQAK"])
results[0].peptide_half_life.value         # hours, higher = longer-lived
results[0].serum_half_life.peptide_input   # exact chemistry + context
results[0].serum_half_life.cache_key       # input + assets + settings
predictor.artifact_inventory.to_dict()     # exact files, hashes, capability
```

Sequence input only. Upstream's SMILES models return a number that is *not* on
the hours scale — the `expm1` inverse transform is applied only to the sequence
model. Although mhctools records exact chemical form, this adapter does not
consume it, so terminal modifications, attachments, and non-standard residues
are rejected rather than scored as their unmodified sequence. Pass
`on_unsupported="record"` to retain unsupported entries in a mixed batch.

> ⚠️ The sequence half-life model was fit on **130 examples** and evaluated by
> cross-validation only, from a preprint, with no external test set and no
> evaluation on long vaccine peptides. Upstream declares Apache-2.0 on its model
> card and MIT in its README. mhctools verifies the exact inference source,
> model, calibration, ESM2 weights, configuration, and tokenizer files before
> launch. The PeptiVerse checkpoint and calibration still use unsafe pickle
> serialization; matching a checksum establishes identity, not safety. Use
> only snapshots you trust.

### PlifePred2

> ⚠️ **This endpoint's semantics are not established.** PlifePred2 ships no
> publication, no training data and no target definition, so its units,
> transform, species and assay matrix are all inferred from the artifacts. By
> default the wrapper reports only the model's native output and claims no
> duration at all.

```python
from mhctools import PlifePred2

predictor = PlifePred2()                       # PLIFEPRED2_HOME + PFEATURE_HOME
results = predictor.predict(["SIINFEKLGGALQAKKY"])
results[0].peptide_half_life.score             # native output, higher = longer-lived
results[0].peptide_half_life.value             # None by default
predictor.artifact_inventory.to_dict()         # exact files, hashes, capability
predictor.last_qc["log10_seconds"]             # the same value, named

# Opt in to a duration, accepting the inference below:
opted_in = PlifePred2(assume_log10_seconds=True)
opted_in.predict(["SIINFEKLGGALQAKKY"])[0].peptide_half_life.value # hours
```

**What is known.** Both shipped models are `RandomForestRegressor`, verified by
loading them. So the output is not a class probability — upstream's docs
("Halflife … Predicted probability") and its CLI's `predict_proba` branch are
both wrong, and the branch is dead code. Being monotone in half-life, the
output ranks correctly whatever the transform turns out to be.

**What is inferred.** `log10(half-life in seconds)` is the strongest reading:
inverting the forests' extreme leaf values under it gives round durations —
exactly 7.000 days for the natural model, 95.0 days for the modified one — to
about seven significant figures, where log2 and ln both invert the whole
training range to a few seconds up to a couple of minutes. The minimum also
lands on 20.2 s, matching the 20-second floor in the lineage paper. That last
point is corroboration rather than proof: the same forests hold targets past
that paper's 24-hour ceiling, so PlifePred2 was trained on a different dataset
and the old filter cannot establish the new target. Note also that the lineage
paper states log2, not log10.

**What is not established.** The species and assay matrix. The result therefore
uses generic `peptide_half_life` with `matrix=None`; do not report it as a
measured whole-blood property or treat it as interchangeable with PeptiVerse's
human-serum endpoint.

Natural peptides only, 12–100 residues. Upstream's CLI silently drops
out-of-range and modified sequences into an `eliminated_sequences.csv` and
returns a shorter result set; mhctools rejects them instead so a caller never
gets a quietly truncated answer.

The Linux-only `pfeature_comp` binary that `plifepred2` bundles is **not** used
— it is a PyInstaller freeze of Pfeature's `pfeature_comp.py`, and that plain
Python source computes the same descriptor on any platform. Both upstreams are
GPLv3, so neither is vendored and neither is imported into the mhctools
interpreter.

> ⚠️ PlifePred2 cites no publication of its own, so its training set is
> unverified beyond what the artifacts reveal. In the lineage paper the
> composition-based natural model was the weaker of the pair (r = 0.643 against
> 0.743), and sequences up to 90% similar were deliberately kept in the data, so
> reported accuracy is optimistic for novel peptides.

## Immunogenicity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `Calis` | immunogenicity | nothing — self-contained |
| `BigMHC_IM` | immunogenicity | `mhctools fetch bigmhc --accept-license` + PyTorch, or set `BIGMHC_DIR` |
| `PRIME` | immunogenicity | [PRIME](https://github.com/GfellerLab/PRIME) clone + MixMHCpred |
| `DeepImmuno` | immunogenicity | `mhctools fetch deepimmuno` + a TensorFlow/Keras-2-capable Python |
| `TLimmuno2` | immunogenicity (class II) | `mhctools fetch tlimmuno2 --accept-license`, or your own clone via `TLIMMUNO2_HOME` |

### Calis

`Calis` is the classic sequence-only IEDB class-I immunogenicity model (Calis
et al. 2013): a fixed per-amino-acid log-enrichment scale weighted by
per-position importance, with the anchor positions (P1/P2/C-terminus) masked
out.

It needs **no external install and no downloaded weights** — the ~30 published
parameters (from the open-access CC-BY paper) are built in — so it is a fast,
dependency-free, allele-independent baseline. It emits one `immunogenicity`
prediction per peptide (empty `allele`); `score > 0` leans immunogenic.

```python
from mhctools import Calis

predictor = Calis()
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score            # 0.30484 (higher = more immunogenic)
```

### PRIME

`PRIME` predicts CD8+ T-cell immunogenicity of class-I peptides by combining
MHC-I binding (via MixMHCpred, which it calls internally) with a
TCR-recognition propensity model. It emits one `immunogenicity` prediction per
(peptide, allele): `score` is the PRIME score (higher = more immunogenic) and
`percentile_rank` is the PRIME %Rank (lower = better).

PRIME is academic / non-commercial licensed, so mhctools shells out to an
install you provide rather than vendoring it.

```python
from mhctools import PRIME

predictor = PRIME(
    alleles=["HLA-A*02:01", "HLA-B*07:02"],
    program_name="PRIME",                    # or an absolute path
    mixmhcpred_path="/path/to/MixMHCpred",   # v3.0+, optional if on PATH
    timeout=300)
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score
```

mhctools verifies MixMHCpred's reported version before PRIME inference and
rejects versions older than 3.0 or an unparseable version. The timeout covers
the PRIME process tree, including its nested MixMHCpred call.

### DeepImmuno

`DeepImmuno` predicts class-I CD8+ immunogenicity from the peptide and its
HLA-A/B/C allele with a small CNN (Li et al. 2021). It scores **9- and 10-mers
only** and supports a fixed set of ~62 alleles, snapping anything else to the
nearest it knows. It emits one `immunogenicity` prediction per (peptide,
allele); `score` is in 0–1 (higher = more immunogenic).

DeepImmuno ships its weights in-repo and is MIT-licensed, but its script loads
them with an old Keras 2 / TensorFlow stack, so mhctools shells out to
DeepImmuno's own CLI in a separate checkout. Run `mhctools fetch deepimmuno`,
or point at a manual clone with `DEEPIMMUNO_HOME`, and set `DEEPIMMUNO_PYTHON`
to an interpreter that has TensorFlow — with Keras 2, or newer TensorFlow plus
the `tf-keras` shim, since the wrapper sets `TF_USE_LEGACY_KERAS=1` for the
subprocess.

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

### TLimmuno2

`TLimmuno2` is the odd one out: it predicts **class-II (CD4+)** immunogenicity
— the only class-II immunogenicity model here, filling a gap the class-I models
(`Calis`, `PRIME`, `BigMHC_IM`, `DeepImmuno`) leave.

It scores a peptide against a class-II allele (transfer-learned from class-II
binding) and emits one `immunogenicity` prediction per (peptide, allele):
`score` in 0–1 (higher = more immunogenic) and `percentile_rank` from its %Rank
against a background set, rescaled to 0–100 (lower = more immunogenic).

Native NetMHCIIpan-style keys (`DRB1_0803`, `HLA-DPA10103-DPB10101`) pass
through; common DR forms (`HLA-DRB1*08:03`) are converted; anything TLimmuno2
does not know raises. Its upstream license is ambiguous (an Apache-2.0 README
badge, no LICENSE file), which mhctools treats the same way as NetCleave: it
can fetch a pinned snapshot, but only when you confirm your own use is
authorized, so the first fetch requires `--accept-license`. `TLIMMUNO2_PYTHON`
names an interpreter that has TensorFlow (Keras 2, or newer TensorFlow plus
`tf-keras`).

```python
from mhctools import TLimmuno2

predictor = TLimmuno2(alleles=["DRB1_0803"])  # TLIMMUNO2_HOME, ~/TLimmuno2, then snapshot
results = predictor.predict(["FHTMWHVTRGAVLMY"])
results[0].immunogenicity.score                    # 0.9874 (higher = more immunogenic)
```

> ⚠️ TLimmuno2's %Rank is computed against ~90,000 background peptides **per
> distinct allele**, so a call costs about a minute per allele regardless of how
> many peptides you pass — batch peptides by allele. Class-II immunogenicity is
> noisier than class-I; a prioritization aid, not ground truth.

## TCR specificity

| Predictor | Kinds produced | Requires |
|---|---|---|
| `NetTCR` | pMHC:TCR binding | `mhctools fetch nettcr --accept-license` + a TFLite runtime (`pip install mhctools[nettcr]`) |
| `Tulip` | pMHC:TCR binding | `mhctools fetch tulip` + a TULIP-capable Python (`TULIP_HOME`, `TULIP_PYTHON`) |
| `MixTCRpred` | fixed-pMHC:TCR binding | `pip install "mhctools[mixtcrpred]"` + `mhctools fetch mixtcrpred --accept-license` |

`NetTCR` and `Tulip` predict pMHC:TCR binding — whether a paired αβ T-cell
receptor (an `mhctools.TCR`, described by its CDR loops) recognises a peptide.
Both take `(peptide, TCR)` inputs; `Tulip` additionally takes the presenting
MHC allele.

### NetTCR

`NetTCR` predicts whether a paired αβ T-cell receptor recognises a (class-I)
peptide. Unlike the MHC-ligand predictors, its input is a peptide plus a `TCR`
(the six CDR loops), not an allele, and it emits the `pMHC_TCR_binding` kind.

NetTCR ships its pretrained weights in its git repository as small TFLite
models. This wrapper runs the pan cross-validation ensemble in-process and does
not need NetTCR's conda environment.

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

### Tulip

```python
from mhctools import Tulip, TCR

tcr = TCR(cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF", name="clone1")
Tulip.fetch()                             # pinned code, tokenizers, and weights
predictor = Tulip()                       # also needs a TULIP-capable Python
results = predictor.predict(["GILGFVFTL"], [tcr], mhc="HLA-A*02:01")
results[0].preds[0].score                 # higher = more likely binding
```

[TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR) is **GPLv3** and pinned
to `transformers==4.32.1`; mhctools is Apache-2.0 and depends on neither torch
nor transformers. The `Tulip` wrapper therefore vendors none of TULIP — it runs
an upstream checkout out-of-process, in an isolated interpreter, via TULIP's own
`predict.py`. `mhctools fetch tulip` obtains the tested code, tokenizers, and
weights; `scripts/setup_tulip_env.sh` can build the separate runtime.

You may instead provide your own checkout and interpreter:

- `TULIP_HOME` — a clone of TULIP-TCR (provides `predict.py`, `src/`,
  tokenizers, and the released `model_weights/`);
- `TULIP_PYTHON` — an isolated **Python 3.11** interpreter with `torch` and
  `transformers==4.32.1`. Python 3.11 specifically, so `tokenizers` installs
  from a prebuilt wheel and needs no Rust toolchain.

### MixTCRpred

MixTCRpred has a different, deliberately explicit shape: each checkpoint is
trained for one fixed peptide/MHC target. Its catalog currently contains 146
models (43 marked high-confidence by upstream), spanning human/mouse class I
and II targets.

The input `TCR` stores paired CDR3s and optional `trav`, `traj`, `trbv`, and
`trbj` assignments. Released models use CDR3alpha/beta plus CDR1/2 derived from
the V genes; J assignments are accepted and QC-reported but are not network
inputs.

Each checkpoint is about 31.8 MB. The immutable Zenodo record contains 146
checkpoints totaling 4.64 GB (4.33 GiB); the 43 high-confidence checkpoints
total 1.37 GB (1.27 GiB). To avoid turning every mhctools installation into a
multi-gigabyte download, the pinned upstream artifact includes its two bundled
reference checkpoints, and additional models are fetched individually with
bounded retries, atomic installation, and checksum verification.

```sh
pip install "mhctools[mixtcrpred]"
mhctools fetch mixtcrpred --accept-license
mhctools ls mixtcrpred --models --high-confidence
mhctools fetch mixtcrpred --model A0201_GILGFVFTL
mhctools mixtcrpred --model A0201_GILGFVFTL \
  --input paired-tcrs.csv --out scored-tcrs.csv
```

```python
from mhctools import MixTCRpred, TCR

models = MixTCRpred.catalog()
target = MixTCRpred.resolve_model("GILGFVFTL", "HLA-A*02:01")
predictor = MixTCRpred(target.name)
tcr = TCR(
    cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF",
    trav="TRAV27", traj="TRAJ42", trbv="TRBV19", trbj="TRBJ2-6",
)
prediction = predictor.predict_tcrs([tcr])[0]
prediction.score
prediction.percentile_rank
```

`annotate_dataframe()` and the CSV command retain the original table and add
the raw score, percentile rank, fixed target metadata, corrected V/J names,
V-derived CDR1/2 loops, and upstream-equivalent QC warning.

MixTCRpred code is academic/non-commercial and fetched directly from its
original repository only after explicit acceptance. Optional checkpoints come
from the authors' immutable CC-BY-4.0 Zenodo record and are checksum-verified
before use.

> ⚠️ As with any PyTorch checkpoint, explicit overrides and user-managed model
> files must come from a trusted source, because loading can execute serialized
> code.
