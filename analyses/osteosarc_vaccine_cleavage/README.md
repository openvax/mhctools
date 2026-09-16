# Osteosarc vaccine-sequence cleavage analysis

This directory contains a reproducible, local-only analysis of the vaccine
sequences disclosed on [osteosarc.com](https://osteosarc.com/vaccines/).
It answers a narrow question: which bonds in the disclosed synthetic long
peptides (SLPs) are flagged by several proteasome/endolysosomal predictors or
by curated human-peptidase recognition rules?

The generated [report](results/REPORT.md), [tables](results/tables), and
[figures](results/figures) are checked in. The complete model/code/weight
inventory and checksums are in [provenance.json](results/provenance.json) and
[SHA256SUMS.json](results/SHA256SUMS.json).
All four figures are also collected in
[`all-figures.pdf`](results/all-figures.pdf), one full page per figure.

For a compact answer to “which SLP is flagged by which model,” start with
[`slp_predictor_matrix.csv`](results/tables/slp_predictor_matrix.csv). Its
column suffixes distinguish within-model fractions above 0.5, native scores,
and motif-match counts; those unlike quantities must not be combined or
ranked as though they shared a scale. Exact bond-level outputs remain in the
long-form tables.

## Scope

- Inventory every disclosed vaccine sequence, separating mRNA encoded
  contexts, displayed mRNA minimal epitopes, and JLF/CeGaT SLPs.
- Keep vaccine targets without a disclosed sequence in a separate table.
- Keep identical sequences assigned to different variants visible as a
  provenance warning; never silently choose one label.
- Score SLP internal bonds with Pepsickle (all-mammal and human-only in-vivo
  models), NetChop 3.1 (Cterm and 20S), and NetCleave (class I and II).
- Evaluate the intact, free-terminal SLP against mhctools' quantitative and
  motif-based human peptidase panel, including the optional ERAMER ERAP1 model.
- Preserve native scores and model-specific applicability. No cross-model
aggregate or biological stability score is calculated.

For NetCleave-I, each candidate bond is represented by the eight residues
ending at that bond plus the three downstream residues required by the model.
For NetCleave-II the corresponding ending peptide is 13 residues. Bonds
without the required upstream or downstream context are not assessed.

Proteasome scores are relevant only after cytosolic access. NetCleave-II is a
C-terminal MHC-II-processing model, not a named-cathepsin assay. Peptidase
motifs are partial recognition rules, not turnover probabilities. None of the
models includes formulation, uptake, abundance, enzyme activation, ordered
digestion, peptide structure, kinetics, or MHC protection.

## Reproduce

The recorded run used the exact commits and file hashes in `provenance.json`.
With those repositories checked out locally and Docker running:

```sh
python analyses/osteosarc_vaccine_cleavage/analyze.py \
  --osteosarc-repo /path/to/osteosarc.com \
  --netchop-dir /path/to/netchop-3.1 \
  --netcleave-dir /path/to/NetCleave \
  --eramer-dir /path/to/ERAMER
```

NetChop is licensed software and is not redistributed here. The script runs a
user-supplied installation inside a pinned 32-bit Debian container because the
available executable is a 32-bit Linux binary. No peptide sequence is uploaded.
NetCleave and ERAMER are run from local checkouts. Pepsickle and the curated
motif panel run in the current Python environment.

## Primary model references

- Pepsickle: Weeder et al., 2021, <https://doi.org/10.1093/bioinformatics/btab628>
- NetChop 3.1: Nielsen et al., 2005, <https://doi.org/10.1007/s00251-005-0781-7>
- NetCleave: Amengual-Rigo and Guallar, 2021,
  <https://doi.org/10.1038/s41598-021-92632-y>
- ERAMER: Al-okaily et al., 2024,
  <https://doi.org/10.1016/j.jim.2024.113713>
- Individual peptidase references and applicability limits are retained in
  `results/tables/model_catalog.csv` and documented in
  [the cleavage guide](../../docs/cleavage.md).
