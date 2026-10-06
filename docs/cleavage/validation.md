# Cleavage coverage and validation status

This page accompanies the [batch API](batch.md). Software conformance,
source-observation reproduction and independent prediction validation are
different claims. No new held-out biological performance number is claimed.

## Proteasome variants

The [Pepsickle author implementation](https://github.com/pdxgx/pepsickle)
provides epitope, gradient-boosted digestion and neural digestion families.
The [paper](https://doi.org/10.1093/bioinformatics/btab628) describes its
independent digestion evaluation. The batch adapter's live checks compare
human/all-mammal neural C/I profiles with the upstream functions on the same
sequence and assert that C/I can change the result. These are conformance
checks, not independent accuracy measurements.

The [paper repository's MASTER.sh](https://github.com/pdxgx/pepsickle-paper/blob/master/MASTER.sh)
expects held-out inputs under `data/validation_data/digestion_data/raw/` and
processed, overlap-filtered validation FASTAs. Those paths are absent from
the inspected repository tree at `c448c4db81925afad78477e74a7d25e0209d3bce`. Its available `data/raw/digestion_map_files`
instead feed the training-data extraction path. Reusing those maps and
calling them held-out would be incorrect. The Wada case study below recovers
source products from one of the paper's held-out studies and audits the raw
training maps. Reconstructing the full processed validation partition and
establishing family-level independence are open work ([known gaps](../known-gaps.md#benchmarks-and-validation));
the case study does not claim to reproduce the paper's overall metrics.

The gradient-boosted artifact records scikit-learn **0.23.2**. It fails under
current scikit-learn (`sklearn.ensemble._gb_losses` is missing). The isolated
Python 3.8.20/0.23.2 runtime for the gradient-boosted models
executes both C/I routes and matches direct upstream inference with networking
disabled. It records the actual subprocess's package, code and weight identity.
The MAGE-A3 sequence in this conformance test occurs in upstream training data;
it is explicitly **not** held-out validation. Neural and gradient-boosted
predictions are not silently substituted.

## APC endolysosomal enzymes

The built-in panel now includes [ITCell](itcell.md) sequence-specific human
cathepsin B/S internal profiles and H initial N-terminal trimming, with native
scores verified against the released author script. These are pH-6.5 assay
specificity profiles, not independently validated cellular processing or
serum-stability models. Cathepsin L and AEP predictors remain unavailable.
IRAP is an exact-substrate source catalog, and [NetCleave](../predictors/processing.md#netcleave)-II is a class-II
C-terminal processing proxy, not an enzyme-specific cathepsin predictor.
The remaining gaps are explicit in batch coverage reports and tracked in
[known gaps](../known-gaps.md#cleavage-validation-and-coverage).

Primary sources for the next validation block:

| Source | Useful evidence | Applicability constraints |
|---|---|---|
| [CatS/AEP processing of MBP](https://pubmed.ncbi.nlm.nih.gov/11745393/) | Product mapping, epitope generation/destruction and protection by HLA-DR | Purified/experimental processing conditions; not a general peptide survival model |
| [DIPPS specificity profiling](https://pubmed.ncbi.nlm.nih.gov/28733325/) | Cathepsin profiles and pH-dependent legumain specificity | Denatured in-gel substrates, condition-specific preferences; cannot turn sequence logos into calibrated probabilities |
| [Legumain activation](https://pubmed.ncbi.nlm.nih.gov/22232165/) | Activation and substrate-dependent pH response | Both activation state and P1 residue affect the result; unconditional cut-after-Asn rules omit this context |
| [Cathepsin MSP-MS comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC10399199/) | Human B/L/S and other cathepsins on 228 14-mers at pH 4.6/7.2, 15/60 minutes, 25 C | Distinct CatB dipeptidyl-carboxypeptidase and endopeptidase activity; source-specific product detection thresholds |

The MSP-MS publication deposits mass-spectrometry data as **MSV000090043 /
PXD035641**. Its workflow uses quadruplicates and significance/fold-change
criteria for detected products. Missing products are not automatically
non-cleaved bonds. The available source supports assay-aware curation; it
does not supply an already validated predictor of long-vaccine processing.
The [CatL example](https://github.com/openvax/mhctools/blob/master/tests/data/tusar2023/README.md) now includes seven
experimentally identified products from two protected synthetic peptides in
[Tusar et al. 2023](https://doi.org/10.1038/s42003-023-04772-8), Supplementary
Data 4 page 1. The five observed internal substrate/bond pairs agree with
Supplementary Table 17. Each original product, terminal modification, source
measurement ID and assay condition survives batch save/reload. The assay is
purified CatL at pH 5.5, 37 C for 2 hours; changing sequence or chemical form
causes abstention. This is a selected source panel, not a transferable model
or validation of whole-cell processing.

The batch `reference_panels` input can preserve curated experiments from
these studies now, with each condition in its own named panel. Actual
novel-sequence inference needs a separately verified adapter/data model.
The model audit below records why novel-sequence coverage remains open.

### Open-model availability audit (updated 2026-10-06)

Only openly licensed models are eligible for this work. Publicly downloadable
files without a project license do not satisfy that constraint.

| Candidate | Verified artifact availability | Remaining requirement |
|---|---|---|
| [Tusar et al. cathepsin S/L/B SVMs](https://doi.org/10.1038/s42003-023-04772-8) | CC BY 4.0 Supplementary Data 3 supplies six SVM-light files; B/L/S each declare 192 input features | Reproduce the original sequence/structure features and independent reference scores before adapting them to vaccine inputs |
| [PCSS backend](https://github.com/salilab/pcss/tree/ea4c3ec81ef30b7a30f3c03508ee2bf1dd78ce34) | LGPL-2.1 code at `ea4c3ec81ef30b7a30f3c03508ee2bf1dd78ce34` | Author README explicitly reports that hard-coded databases/programs prevent operation outside the Sali lab; weights alone do not supply this pipeline |
| [ProsperousPlus](https://github.com/lifuyi774/ProsperousPlus/tree/66a9d08cd5a44febf64950caf9684c81aa0e8807) | Directories C01.060 (CatB), C01.032 (CatL), C01.034 (CatS), C13.004 (animal legumain) are present | GitHub license metadata is null and the root has no project license; excluded under the open-model requirement ([known gaps](../known-gaps.md#cleavage-validation-and-coverage)) |
| [panCleave](https://gitlab.com/machine-biology-group-public/pancleave) | Author describes a pooled, protease-agnostic random forest | Its output cannot supply enzyme-specific CatS/L/B/AEP coverage |
| [DIPPS legumain study](https://doi.org/10.15252/embj.201796750) | Experimental pH-dependent specificity evidence | It is not a published fitted AEP predictor; an unconditional cut-after-Asn rule would discard the reported context |
| [PhageScout](https://doi.org/10.3390/ijms27177593) | Ten bundled [sequence profiles](phagescout.md) reproduce [Zenodo 21387981](https://zenodo.org/records/21387981) native values and missing scores; all four phage-only boosters load in R XGBoost 3.2.1.1 | Classifiers lack saved training-imputation medians and require complete-mature-protein normalization; structural inputs and independent free-peptide accuracy remain unverified |
| [CatS/L/B octamer neural ensembles](https://doi.org/10.3390/ijms20194843) | Published sequence-property/JMP methodology and [CC BY 4.0 observed-product datasets](https://doi.org/10.6084/m9.figshare.9777725) | No trained ensemble parameters or licensed portable inference runtime verified; experimental CSV files are not model weights |

The CatB/L/S files were extracted from the Europe PMC open-access supplement
archive for PMC10124925 and inspected. Their SHA-256 digests are:

- `SuppData3_CatB_SVMmodel.txt`: `3acfa1cf903657694a6f5d68ba5151e6ea1ca191b2d7c212d0723d9fa1b5caac`
- `SuppData3_CatL_SVMmodel.txt`: `dc013a06685c727384ec64cf9e0dc6166cf514fb5e849f66e0a07872cfba5951`
- `SuppData3_CatS_SVMmodel.txt`: `6c2c857c37f56d916c8706f43d008cefe6b47ee2feea43a35d66225fca6df60b`

The article describes protein secondary-structure and solvent-exposure inputs.
No unverified feature values, replacement model, or independent accuracy claim
are supplied here. Completing this open-model runtime and obtaining an openly
licensed, verified AEP predictor remain concrete blockers ([known gaps](../known-gaps.md#cleavage-validation-and-coverage)). Requests
for novel-sequence CatL/AEP predictions must continue to report unsupported
coverage. ITCell B/S/H profiles and experimental source imports are available now.

The [additional candidate audit](candidates.md) distinguishes published
recognition evidence, released trained artifacts and missing inference inputs.
The PhageScout Zenodo assets correct an earlier GitHub-only availability audit.
The native sequence adapter reproduces 3,822 PWM and 36 peptide-profile values
on two source mature proteins, including their complete missing-score masks.
This establishes implementation conformance, not independent accuracy or
calibrated long-peptide degradation. Learned-classifier inference still needs
the source preprocessing state.

## Source-backed long-peptide case study

[Wada et al. 2018](https://doi.org/10.1371/journal.pone.0199249) is explicitly
assigned to validation in [Pepsickle](../predictors/processing.md#pepsickle)'s Table 1. The fixture curates all 47
detected products from Figure 2A and 2C: two 31-residue vaccine constructs
containing the same epitopes in different orders, joined by RR linkers.
The dataset retains first-detection times, figure row IDs, parent endpoints
and the purified murine-i20S assay conditions.

The real gradient-boosted immunoproteasome model scores the complete constructs.
Its native scores are paired with 24 observed internal construct/bond pairs;
parent sequence ends never become cleavage labels. The pinned raw training-map
audit covers 79 files and 58 distinct source sequences. It finds no exact
construct or observed seven-residue cleavage-context overlap. Wada's DOI is
absent from the reported training DOI fields; unresolved source identifiers,
homology-family overlap and the original fitted partition remain explicit
limitations. This is a small study-held-out case study, not a new calibrated
performance estimate or a reconstruction of all 225 author validation windows.

The [fixture README](https://github.com/openvax/mhctools/blob/master/tests/data/wada2018/README.md) gives the primary source,
license, curation scope and reproducible command. The report preserves observed
products through batch save/reload, with endpoint/boundary/internal overlays.
No missing product is relabeled as a negative, and no product detection time is
converted into a predicted cleavage rate or presentation outcome.

Additional extracellular coverage is tracked in [known gaps](../known-gaps.md#cleavage-validation-and-coverage):
CleaveNet's MMP substrate scores need their own native endpoint and assay
validation rather than conversion to per-bond probabilities.

## Required validation before biological ranking

Use exact source substrates/products, enzymes, species, chemical forms,
activation, pH, temperature, time and detection limits. Separate positive
product observations from verified non-cleavage measurements. Audit training
overlap by study and sequence/family before claiming held-out performance.
Preserve unsupported inputs and failures in the denominator.

Evaluate purified-enzyme cleavage, whole-matrix degradation and antigen
presentation as distinct endpoints. In particular, [human-DC long-peptide
cross-presentation](https://pubmed.ncbi.nlm.nih.gov/24587108/) supports a
possible proteasome/TAP route; it does not calibrate every construct or imply
that cleavage scores alone predict presentation. [IFN-gamma tumor-organoid
experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC10655636/) likewise show
why a single sequence-only profile cannot encode the tumor's cellular state.
