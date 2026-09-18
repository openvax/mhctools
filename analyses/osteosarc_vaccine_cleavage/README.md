# Osteosarc vaccine-sequence cleavage analysis

This directory contains a reproducible, local-only analysis of the vaccine
sequences disclosed on [osteosarc.com](https://osteosarc.com/vaccines/).
It answers a narrow question: which bonds in the disclosed synthetic long
peptides (SLPs) are flagged by several proteasome/endolysosomal predictors or
by curated human-peptidase recognition rules?

Every invocation writes a new date/time-stamped directory under `results/` so
an earlier analysis is never overwritten. The [generated-run index](results/README.md)
links the checked-in report, tables, and
[figure atlas](results/2026-09-17T012843-301909-0400/mhctools-all-figures.pdf).
The complete model/code/weight inventory and checksums travel with that run.
The PDF begins with a clustered predictor/SLP agreement overview and then
uses each disclosed SLP sequence as the central visual axis. Large residue
letters carry bond-aligned cleavage profiles and motif flags, disclosed minimal
epitope windows, and candidate MHC-I/MHC-II ligand spans. Cleavage tracks are
continuous piecewise-linear profiles through the exact per-bond scores: there
is no smoothing, averaging, or interpolation across unassessed gaps. Candidate
ligand windows hug the sequence from above (MHC-I) and below (MHC-II), with the
corresponding processing tracks farther outward. A pale blue or purple tint
directly behind the residue letters marks coverage by at least one displayed
window; it is binary coverage, not prediction count or probability. Red ticks inside a
ligand span show relevant pre-binding internal cleavage evidence; they do not
claim post-binding cleavage or protection. The 80-aa outlier is split across three
continuation pages rather than compressed. Use
[`atlas_sequence_order.csv`](results/2026-09-17T012843-301909-0400/tables/atlas_sequence_order.csv) to jump
from a sequence record to its PDF page. Every map page is also exported as a
vector PDF and 300 dpi PNG; `slp_map_exports.csv` indexes those files.

The atlas selects up to ten windows per MHC class and sequence segment, using
five dedicated non-overlapping lanes. It first
retains the strongest eligible window overlapping a disclosed intended epitope,
then tries to represent distinct alleles, and finally fills free, non-overlapping
lane capacity by native percentile rank. This is a display-selection rule, not
an ensemble score. Each map header states “shown / eligible” instead of an
ambiguous “+N more” footnote. `slp_mhc_display_selection.csv` records every
displayed row and its reason; `slp_mhc_ligand_predictions.csv` retains every raw prediction.
The run also includes a four-page, full-size manuscript subset selected by four
declared criteria rather than visual preference. Its exact choices and metrics
are in `manuscript_figure_selection.csv`.

For a compact answer to “which SLP is flagged by which model,” start with
[`slp_predictor_matrix.csv`](results/2026-09-17T012843-301909-0400/tables/slp_predictor_matrix.csv). Its
column suffixes distinguish within-model fractions above 0.5, native scores,
and motif-match counts; those unlike quantities must not be combined or
ranked as though they shared a scale. Exact bond-level outputs remain in the
long-form tables.

The sequence pages treat the injected SLP route explicitly. Four complementary
cytosolic cross-presentation tracks sit above the sequence: human-only
Pepsickle, NetChop Cterm, NetChop 20S, and NetCleave-I. This is a conditional
route after dendritic-cell uptake/export, not direct exposure of an injected
peptide to the cytosol. A
red bond mark is drawn only when all four assess the bond and at least three
reach the common 0.5 display threshold. Four small beads on the mark encode the
support count (filled beads are hits and an open bead is a miss); this is not a
probability. A faint guide connects that exact bond to the four intracellular
tracks. NetChop Cterm is retained because these maps
emphasize candidate MHC-I ligand boundaries, but it is ligand-trained and not a
pure proteasome assay. NetChop 20S is shown separately as an in-vitro
proteasome view. The human-only Pepsickle model is species-matched but
experimental and trained on less data than the all-mammal model; the
near-redundant all-mammal output remains in the tables. The NetCleave-II track
has its own primary endolysosomal/class-II section below the sequence, while
DPP4 and matched MME, FAP, ANPEP, and ENPEP rules each receive a readable
enzyme-specific extracellular track. FAP is labeled tumor-stroma conditional.
Plasma-oriented ACE, CPB2, and CPN, XPNPEP2, intact-SLP cytosolic aminopeptidase
motifs, and ERAP1-on-the-intact-SLP are omitted from the map because injection
does not establish their exposure or substrate state; their raw assessments
remain in the tables. The context-separated
`slp_vulnerable_bonds.csv` retains its separate, conservative three-of-five
rule and never combines biological contexts.

For a reusable input contract rather than the osteosarc-specific adapter, use
`mhctools vaccine-report`; see [`docs/vaccine-reports.md`](../../docs/vaccine-reports.md).
Its manifest distinguishes `synthetic_long_peptide` from `rna_encoded` and
requires an explicit RNA routing policy. Cytosolic RNA emphasizes endogenous
proteasome/ER/class-I processing and omits free serum-peptide tracks unless
secretion or extracellular exposure is declared.

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
- Scan 8–11-mers with local MHCflurry class-I presentation models across the
  five disclosed classical class-I alleles.
- Scan 13–21-mers with local NetMHCIIpan 4.3 EL models. Only class-II
  combinations already named in the source candidate fields are used; no
  alpha/beta phase is guessed from the unphased HLA table.
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
  --eramer-dir /path/to/ERAMER \
  --netmhciipan-path /path/to/netMHCIIpan-4.3
```

NetChop is licensed software and is not redistributed here. The script runs a
user-supplied installation inside a pinned 32-bit Debian container because the
available executable is a 32-bit Linux binary. No peptide sequence is uploaded.
NetCleave and ERAMER are run from local checkouts. Pepsickle and the curated
motif panel run in the current Python environment. MHCflurry and NetMHCIIpan
also run locally; the complete file-level model inventories are recorded in
the timestamped run's `tables/mhc_model_file_inventory.csv`.

For a visualization-only revision, reuse the exact checksummed prediction
tables instead of silently rerunning models:

```sh
python analyses/osteosarc_vaccine_cleavage/rerender.py \
  --source-run analyses/osteosarc_vaccine_cleavage/results/<timestamp>
```

The command verifies every source-run checksum, creates another timestamped
directory, and records that inference was not repeated in `provenance.json`.

## Primary model references

- Pepsickle: Weeder et al., 2021, <https://doi.org/10.1093/bioinformatics/btab628>
- NetChop 3.1: Nielsen et al., 2005, <https://doi.org/10.1007/s00251-005-0781-7>
- NetCleave: Amengual-Rigo and Guallar, 2021,
  <https://doi.org/10.1038/s41598-021-92632-y>
- ERAMER: Al-okaily et al., 2024,
  <https://doi.org/10.1016/j.jim.2024.113713>
- MHCflurry 2.0: O'Donnell et al., 2020,
  <https://doi.org/10.1016/j.cels.2020.06.010>
- NetMHCIIpan 4.3: Nilsson et al., 2023,
  <https://doi.org/10.1126/sciadv.adj6367>
- Individual peptidase references and applicability limits are retained in
  the timestamped run's `tables/model_catalog.csv` and documented in
  [the cleavage guide](../../docs/cleavage.md).
