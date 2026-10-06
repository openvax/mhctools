# Cleavage site and mechanism scorecards

The [PDF scorecards](results/2026-10-06T142004-779643-0400-mechanisms/cleavage-mechanisms.pdf)
show every internal bond of 40 disclosed Sid vaccine SLP records (39 unique
sequences). Every processing row names its training endpoint and every named
peptidase retains its internal or terminal topology. The PDF has bookmarks
and a [page index](results/2026-10-06T142004-779643-0400-mechanisms/page_index.csv).

The [site table](results/2026-10-06T142004-779643-0400-mechanisms/cleavage_sites.csv)
retains all 7,248 frozen assessments, including native score strings,
categorical motif states, abstention reasons, vaccine/source links, original
context annotations, corrected mechanism descriptions, and core coordinates.
The [mechanism catalog](results/2026-10-06T142004-779643-0400-mechanisms/mechanism_catalog.csv)
names each model's enzyme attribution, topology, assay/training endpoint,
native units, location and exposure requirements, limitations, and primary
references. The [provenance receipt](results/2026-10-06T142004-779643-0400-mechanisms/provenance.json)
retains the source run's complete model/runtime provenance and records that
inference was not repeated. All 107 source artifacts were checksum-verified.
The [output checksums](results/2026-10-06T142004-779643-0400-mechanisms/SHA256SUMS.json)
cover the generated PDF and tables.

## Interpretation

- **Pepsickle:** human-only is preferred for a declared human context;
  all-mammal is preferred for nonhuman, mixed, or unknown context. Both recorded
  scores remain separate. This display policy does not assert that human-only
  is more accurate: upstream calls it experimental, with a smaller training
  set. The epitope family is proteasome-type agnostic, not constitutive-only.
  All-mammal training does not validate arbitrary non-mammalian organisms.
- **NetChop Cterm versus 20S:** Cterm is trained on MHC-I ligand C-termini and
  is a ligand-boundary processing proxy. 20S is trained on in-vitro proteasome
  digests. DTU reports Cterm performs best for CTL epitope boundaries; that is
  not a claim about extracellular vaccine degradation. Neither assigns a
  particular proteasome catalytic subunit.
- **DPP4 qPISA:** only the exposed N-terminal dipeptide bond (SLP bond 2),
  displayed as whole-percent predicted loss in the source four-hour assay;
  native log2 substrate depletion remains unchanged in CSV. It is not a 0-1 processing score,
  serum half-life, or internal-proline cleavage rule. Repeated trimming and
  hypothetical newly exposed fragments are not simulated.
- **ERAMER:** ERAP1's initial N-terminal trimming step (bond 1), using native
  PWM specificity and a 9-16-residue input domain; ER access is conditional.
- **PREP, MME, FAP and other motifs:** categorical recognition evidence with
  its strictness grade. A match is not a measured cut, numerical probability,
  or rate. Terminal rules only assess the currently exposed end. CPB2's
  frozen producer explicitly assumed activated enzyme.

Gold headers mark bonds strictly inside the source-listed minimal epitope
only when its recorded offset matches the SLP. The source annotation is an
mRNA minimal/candidate epitope; it does not establish the experimentally used
minimal target for every peptide vaccine. Missing or nonmatching offsets are
left unlocated. In particular, TECPR1's CeGaT peptide contains a different
candidate core; its internal bonds are shown, without borrowing the mRNA
core annotation.

Numerical grid cells are rounded to three decimals; CSV retains native
precision. DPP4 uses the same simple percent-loss display as the atlas,
converted as `100 * (1 - 2**(-score))`; negative estimates say `No predicted loss`
and small positive estimates say `<1% predicted loss`. The display describes
predicted relative peptide-signal loss in the source assay, not observed
vaccine-peptide loss or in-vivo degradation.
Bold cells meet the raw score's within-model 0.5 display threshold.
No threshold is applied to DPP4 or ERAMER. `NA` is unassessed/unsupported,
never zero. Large triangles mark motif matches and hollow circles mark
assessed non-matches; motifs remain qualitative evidence.
The models and their biological contexts are not combined into a consensus,
whole-peptide survival probability, uptake estimate, or clinical conclusion.
The original immutable predictions and their annotations remain available.
There are no CleaveNet predictions in this view.

## Reproduce offline

Install the rendering extra, then choose a new output directory:

```sh
python -m pip install -e '.[vaccine-report]'
python analyses/osteosarc_vaccine_cleavage/mechanism_scorecards.py \
  --source-run analyses/osteosarc_vaccine_cleavage/results/2026-09-18T175125-855754-0400 \
  --output-dir /tmp/sid-cleavage-mechanisms-new \
  --organism human
```

`--organism` accepts `human`, `nonhuman`, `mixed`, or `unknown`. It changes
Pepsickle display preference and provenance, not prediction values. The
renderer requires a complete valid source manifest, rejects duplicate or
inconsistent assessments, and refuses an existing output directory. It does
not require predictor binaries, model downloads, Docker, or network access.

Validation:

```sh
python -m pytest tests/test_cleavage_mechanism_scorecards.py \
  tests/test_osteosarc_visualization.py
```

Primary sources:
[NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/),
[Pepsickle](https://github.com/pdxgx/pepsickle),
[NetCleave](https://doi.org/10.1038/s41598-021-92632-y),
[DPP4 qPISA](https://doi.org/10.1038/s44320-024-00071-4),
[PREP/POP profiling](https://pubmed.ncbi.nlm.nih.gov/22750443/).
Every model's full reference list is in the mechanism catalog.
