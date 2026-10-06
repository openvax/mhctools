# Additional protease candidates

Primary-source and author-artifact audit, checked 2026-10-06. Bundled
[PhageScout sequence profiles](phagescout.md) now reproduce native scores;
the other routes below remain integration candidates.
[Known gaps](../known-gaps.md#cleavage-validation-and-coverage) tracks their
implementation and shared validation dependencies.

The remaining gaps have different causes: some predictors publish weights but
require unreproduced features, some describe neural models without supplying
verified inference parameters, and many experiments supply useful recognition
preferences without a transferable quantitative predictor. Recognition flags
can be added under their experimental scope without calling them calibrated
probabilities or serum-survival estimates.

## Released predictor assets

**PhageScout sequence profiles are available now.** Its
[2026 paper](https://doi.org/10.3390/ijms27177593) describes neutrophil elastase
and cathepsin G specificity scoring. [Zenodo 21387981](https://zenodo.org/records/21387981)
declares CC BY 4.0 and supplies PWMs, peptide profiles, reference scores and
trained R/XGBoost objects, including phage-only models. The earlier audit
missed this separate release. The
[GitHub revision](https://github.com/YuEnoch/PhageScout/tree/663e5189d31b6ce74493f26bec27a4ec2b7d1963)
has no project license; its code and the licensed deposited assets require
separate treatment.

The independently implemented sequence-only route reproduces deposited
native scores. Five-mer substrate enrichment, inferred aligned-nine-mer
anchors and known MEROPS cuts have different evidentiary status. The author
notebook contains raw, normalized and structural features; exact anchor,
terminal handling, normalization population and trained-model feature order
remain requirements for new-sequence classifier inference. Structural classifiers
require actual applicable features. Neither a normalized landscape score nor
a classifier trained on balanced sites is a serum-loss percentage.

Downloaded files matched the author's MD5 checksums. SHA-256 identities:

| Asset | SHA-256 |
|---|---|
| `elastase_relaxed_aligned_pwm.txt` | `5063e9e26da2c0bf21e3c4b2b668838e7c0bc3856e6a6ab6d05f753b453d2993` |
| `cathepsin G_relaxed_aligned_pwm.txt` | `16302980d91388fe826ca11eb66d9043f09976e7123db788eba7c81a2efd7481` |
| `elastase_phage_balanced_XGB.rds` | `713d4b42cfb7443061b93dba51bb905b7afefac78a910fd69256ffdd91df2e28` |
| `cathepsin_G_phage_balanced_XGB.rds` | `f98247ce7b8669f90e4ac7c168091e9b18e7ad90e009c2d5daf700701a3d9a28` |

Ten native sequence profiles reproduce deposited values and missing scores.
All four phage-only boosters load in the official R XGBoost 3.2.1.1 runtime;
their 18 feature names are retained, but no training-imputation medians are
saved in their attributes. Automatic new-sequence classifier inference and
independent accuracy remain unverified.

**Other trained-model routes remain distinct.**
[Tusar et al. CatS/L/B SVMs](https://doi.org/10.1038/s42003-023-04772-8)
supply model files but require a sequence/structure feature pipeline.
[ProsperousPlus](https://doi.org/10.1093/bib/bbad372) supplies pretrained
protease directories but still has no project license at the audited
[revision](https://github.com/lifuyi774/ProsperousPlus/tree/66a9d08cd5a44febf64950caf9684c81aa0e8807).
[CatS/L/B octamer neural ensembles](https://doi.org/10.3390/ijms20194843)
describe sequence-property/JMP models; no released trained parameters or
portable inference runtime were verified. Their
[CC BY 4.0 Figshare deposit](https://doi.org/10.6084/m9.figshare.9777725)
contains three observed-product CSVs, not the neural weights.

[UniZyme](https://github.com/Ao-LiChen/UniZyme/tree/8f932ccdd82bd9551c9046095856dc20f445d88e)
has MIT code and a [CC BY 4.0 data/weight archive](https://zenodo.org/records/20673344).
Its structure/energy input pipeline and application to free peptide
conformations still need an independent applicability audit. Availability of
that archive alone does not establish runnable predictions for every enzyme.

## Experimental recognition evidence

| Enzyme | Primary evidence worth curating | Scope to preserve |
|---|---|---|
| Human CatL | [PICS, 845 mapped sites](https://doi.org/10.1021/pr200621z), [DIPPS](https://doi.org/10.15252/embj.201796750), [cathepsin MSP-MS comparisons](https://pmc.ncbi.nlm.nih.gov/articles/PMC10399199/) | Aromatic P2 is a preference, not a universal requirement. Retain library, pH, time, reducing conditions and source erratum. |
| Human AEP/legumain | [pH-specific DIPPS datasets](https://doi.org/10.15252/embj.201796750), [counter-selection libraries](https://pmc.ncbi.nlm.nih.gov/articles/PMC5939596/), [peptidase/ligase experiments](https://doi.org/10.1021/acscatal.1c02057) | Asn/Asp acceptance changes with pH and reaction context. Keep hydrolysis, ligation, activation and human/plant enzymes distinct. |
| Human thrombin | [extended phage/linker specificity](https://doi.org/10.1371/journal.pone.0031756), [motif/loop validation](https://pubmed.ncbi.nlm.nih.gov/22050556/), [high-throughput profiling](https://pmc.ncbi.nlm.nih.gov/articles/PMC5809430/) | Include cooperative alternatives: P2 Pro is preferred but not universally required. Actual thrombin activity, access and exosite/cofactor context matter. |
| Human plasmin | [P1-Lys-fixed positional library](https://doi.org/10.1038/nbt0200_187), [subsite cooperativity](https://pubmed.ncbi.nlm.nih.gov/21877690/), [P1-Arg-fixed plasma-protease profiling](https://www.seas.upenn.edu/~diamond/Pubs/2006_Gosalia_Biotech_Bioeng.pdf) | Keep fixed-P1 assumptions, reporter chemistry, prime-side gaps and noncanonical residues explicit. Inhibitor affinities are not substrate cleavage rates. |
| Human plasma kallikrein KLKB1 | [fluorogenic microarrays](https://doi.org/10.1074/mcp.M500004-MCP200), [human-plasma activation comparisons](https://www.seas.upenn.edu/~diamond/Pubs/2006_Gosalia_Biotech_Bioeng.pdf), [separately identified human MSP-MS comparator](https://doi.org/10.1371/journal.pntd.0006446) | Distinguish KLKB1 from [tissue KLK1/KLK6](https://pmc.ncbi.nlm.nih.gov/articles/PMC2271166/) and the parasite enzyme in the comparator study. Preserve contact activation, inhibitors and sample preparation. |
| Human ELANE/PRTN3/CTSG | [ELANE/PRTN3 extended preferences](https://doi.org/10.3389/fimmu.2018.02387), [CTSG dual specificity](https://pmc.ncbi.nlm.nih.gov/articles/PMC5898719/), [protein-substrate tests](https://pmc.ncbi.nlm.nih.gov/articles/PMC7014372/) | Keep each enzyme distinct, including CTSG secondary Lys specificity, tentative phage anchors, detergent effects and protein accessibility. |

The practical first additions are conditional recognition flags and exact
observed-bond imports, followed by published profile scores where the native
calculation can be reproduced. This ordering is an implementation assessment,
not an accuracy ranking. Model and source-sequence overlap must be audited
before calling any comparison independent validation. Enzyme activation and
exposure are separate from sequence specificity; purified-enzyme preferences
cannot determine whole-serum survival or productive antigen presentation.
