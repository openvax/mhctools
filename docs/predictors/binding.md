# MHC binding and presentation predictors

Class I and class II predictors of peptide-MHC binding affinity, eluted-ligand
presentation, and complex stability. Every row of the
[predictor matrix](../predictor-matrix.md#mhc-binding-and-presentation) links
to the section here that explains it.

All of these take **peptides plus alleles** and return one prediction per
peptide-allele pair (or one per peptide for MHCflurry's haplotype mode). Spell
alleles any way [mhcgnomes](../alleles.md) understands, for example
`HLA-A*02:01`. The default scanning window is narrow (9 residues for class I,
15-20 for class II); see [peptide lengths](index.md#peptide-lengths).

The DTU tools (`NetMHCpan`, `NetMHC`, `NetMHCcons`, `NetMHCIIpan`,
`NetMHCstabpan`) are identity-bound academic licenses, so mhctools calls an
installation you provide rather than fetching one: put the executable on
`PATH` or pass `program_name="/path/to/netMHCpan"`. See
[licensing](../licensing.md).

## NetMHCpan

`NetMHCpan(...)` runs `netMHCpan --version` and returns the class for whatever
is installed. Use a version-specific class when you want the version fixed in
your code rather than discovered: `NetMHCpan42`, `NetMHCpan41`, `NetMHCpan4`,
`NetMHCpan3`, `NetMHCpan28`. The `_BA` and `_EL` variants restrict a class to
binding affinity or eluted-ligand presentation.

NetMHCpan 4.1 and 4.2 emit both `pMHC_affinity` and `pMHC_presentation` for
every peptide-allele pair from a single run.

```python
from mhctools import NetMHCpan42

predictor = NetMHCpan42(alleles=["HLA-A*02:01"])   # program_name="netMHCpan-4.2"
results = predictor.predict(["SIINFEKL", "GILGFVFTL"])

results[0].affinity.value                 # predicted IC50, nM
results[0].presentation.percentile_rank   # %Rank of the eluted-ligand score
```

The default window is 9 residues. Pass `default_peptide_lengths=[8, 9, 10, 11]`
to the constructor or `peptide_lengths=` to `predict_proteins()` for more.
`process_limit`, `max_peptides_per_file` and `max_alleles_per_command` control
how large batches are split across `netMHCpan` processes.

## NetMHC

`NetMHC(...)` inspects the installed tool and returns `NetMHC3` (3.4) or
`NetMHC4` (4.0). Both emit class I `pMHC_affinity`.

```python
from mhctools import NetMHC

results = NetMHC(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
results[0].affinity.value                 # IC50, nM
```

NetMHC 3.4 needs Python 2 and a Linux x86 executable. On Apple Silicon, see
[legacy NetMHC on Apple Silicon](../backends.md#legacy-netmhc-on-apple-silicon).

## NetMHCcons

`NetMHCcons` 1.1 is the consensus of several NetMHC-family methods and emits
class I `pMHC_affinity`.

```python
from mhctools import NetMHCcons

results = NetMHCcons(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
```

Like NetMHC 3.4 it needs Python 2 and a Linux x86 executable; see
[legacy NetMHC on Apple Silicon](../backends.md#legacy-netmhc-on-apple-silicon).

## NetMHCIIpan

`NetMHCIIpan(...)` returns `NetMHCIIpan43` for a 4.3 install and the 4.x or 3.x
class otherwise. The `_BA` and `_EL` variants select affinity or eluted-ligand
presentation: `NetMHCIIpan43` and `NetMHCIIpan4` default to presentation, and
`NetMHCIIpan43_BA` / `NetMHCIIpan4_BA` / `NetMHCIIpan3` emit affinity.

```python
from mhctools import NetMHCIIpan

predictor = NetMHCIIpan(alleles=["HLA-DRB1*15:01", "HLA-DPA1*01:03-DPB1*04:01"])
results = predictor.predict(["GELIGTLNAAKVPAD"])
results[0].presentation.score
```

The default window is 15-20 residues. Write DP and DQ alleles as alpha-beta
pairs (`HLA-DPA1*01:03-DPB1*04:01`); see [allele names](../alleles.md).

## NetMHCstabpan

`NetMHCstabpan` predicts the half-life of the assembled peptide-MHC complex,
emitting class I `pMHC_stability` with `value` in hours. That is a different
quantity from `peptide_half_life`; see [prediction kinds](../kinds.md#the-kinds).

```python
from mhctools import NetMHCstabpan

results = NetMHCstabpan(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
results[0].stability.value                # complex half-life, hours
```

It has no default scanning window, so pass `peptide_lengths=` to
`predict_proteins()`.

## MHCflurry

MHCflurry ships as a dependency of mhctools, so only its model weights need
downloading:

```sh
mhctools fetch mhcflurry
```

`MHCflurry` uses the modern presentation API and emits three kinds: per-allele
`pMHC_affinity`, `pMHC_presentation`, and the allele-independent
`antigen_processing` score (read it with `result.processing`).
`MHCflurry_Affinity` uses the older affinity-only API and emits
`pMHC_affinity` alone; fetch its weights with `mhctools fetch mhcflurry-affinity`.

```python
from mhctools import MHCflurry

results = MHCflurry(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
results[0].affinity.value        # IC50, nM
results[0].processing.score      # antigen-processing score, not allele-specific
```

`presentation_allele_mode` controls how the requested alleles are interpreted:

- `"haplotype"` treats them as one sample genotype and emits one
  `pMHC_presentation` record per peptide. The `allele` field carries
  MHCflurry's `best_allele` attribution when available.
- `"per_allele"` treats each allele as a separate one-allele synthetic sample
  and emits one presentation record per peptide/allele pair.
- `"auto"` (the default) uses haplotype mode for up to six alleles and
  per-allele mode for larger panels.

MHCflurry predictions carry both the Python package and the official
model-release identity in `predictor_version`, for example
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

## BigMHC

`BigMHC` wraps two class I models behind one constructor: `BigMHC_EL`
(eluted-ligand `pMHC_presentation`) and `BigMHC_IM` (`immunogenicity`). The
generic `BigMHC(alleles, mode="el" | "im")` selects between them. Models load
on the first `predict()` call and stay in memory.

```python
from mhctools import BigMHC_EL, BigMHC_IM

presentation = BigMHC_EL(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
presentation[0].presentation.score

immunogenicity = BigMHC_IM(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
immunogenicity[0].immunogenicity.score
```

```sh
mhctools fetch bigmhc --accept-license        # academic license; review it first
```

Set `BIGMHC_DIR` (or pass `bigmhc_path=`) to use your own clone. PyTorch is
required; `device="cpu"` is the default. Read the
[immunogenicity caveats](immunogenicity.md#read-this-before-trusting-a-score)
before using `BigMHC_IM` to rank neoepitopes.

## CapHLA

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

CapHLA's performance numbers are author-reported. Treat it as a complementary
research predictor, not a default and not an independent validation.

## MixMHCpred

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

## MixMHC2pred

`MixMHC2pred` is a pan-allele **class-II** presentation predictor and a strong
complement to `NetMHCIIpan`. The two were independently co-best in the
*Frontiers in Immunology* 2024 class-II benchmark.

It emits one `pMHC_presentation` prediction per (peptide, allele): `score` is
the raw MixMHC2pred score (higher = better) and `percentile_rank` is its %Rank
(lower = better).

It is academic / non-commercial licensed, so mhctools shells out to an install
you provide. Download a release, not a bare clone: the release ships the
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

## SMM and SMM-PMBEC

`SMM` and `SMMPMBEC` run IEDB's official standalone matrix methods locally and
emit class I `pMHC_affinity` as IC50 in nM. Unsupported allele and length
pairs fail explicitly.

```python
from mhctools import SMM

results = SMM(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
results[0].affinity.value
```

The launcher is found through `program_name=`, then `IEDB_MHCI_EXECUTABLE`,
then `iedb-mhci` on `PATH`. To install the pinned IEDB subset, see
[SMM and SMM-PMBEC setup](../backends.md#local-smm-and-smm-pmbec).

## RandomBindingPredictor

A built-in predictor that returns random class I affinities. It needs no
install and exists as a null baseline for evaluations and for testing code
that consumes predictions. Never use its output as a prediction.

```python
from mhctools import RandomBindingPredictor

results = RandomBindingPredictor(alleles=["HLA-A*02:01"]).predict(["SIINFEKL"])
```

## Compatibility names for the old IEDB predictors

Every predictor runs locally. The historical Python names
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
provenance. The full list of renamed and removed spellings is in the
[migration guide](../migration.md).

There is no HTTP fallback. The HTTP-only `url`, `request_timeout`, and
`raise_on_error` constructor arguments and the CLI `--do-not-raise-on-error`
option have been removed. Missing installations and unsupported inputs raise
errors, so predictions are never silently dropped by an IEDB error policy. For
the standalone matrix methods, use the CLI names `smm` and `smm-pmbec`.
