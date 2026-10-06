# PhageScout sequence profiles

PhageScout supplies separate native recognition scores for **human neutrophil
elastase (ELANE)** and **human cathepsin G (CTSG)**. These are extracellular
and inflammatory enzyme candidates under active-enzyme exposure. A sequence
score does not establish enzyme presence, activation, concentration,
accessibility, inhibition, or degradation in serum.

The adapter independently implements the scoring described in the
[2026 paper](https://doi.org/10.3390/ijms27177593), using ten unmodified
matrices and significant-peptide profiles from
[Zenodo 21387981](https://zenodo.org/records/21387981). The deposited data are
CC BY 4.0. Attribution, the license, exact source-file SHA-256 identities and
a checksum for the bundled JSON container accompany the package. No author
notebook code is included; there is no external inference runtime to install.

## Select a native model

```python
from mhctools import PhageScout

result = PhageScout(enzyme="ELANE", profile="pwm-relaxed-aligned").predict("MAACVYTLA")
for site in result.sites:
    print(site.bond, site.score, site.reason)
```

```sh
mhctools cleavage --sequence MAACVYTLA --model phagescout-elane-pwm-relaxed-aligned
mhctools cleavage --sequence AAFLDAAFF --model phagescout-ctsg-peptide-relaxed-unaligned
```

| ELANE model | CTSG model | Native calculation |
|---|---|---|
| `phagescout-elane-pwm-deseq2` | `phagescout-ctsg-pwm-deseq2` | Mean of complete five-mer sums of DESeq2-derived positional weights |
| `phagescout-elane-pwm-relaxed-unaligned` | `phagescout-ctsg-pwm-relaxed-unaligned` | Mean of complete five-mer sums of significant-peptide log2 enrichment weights |
| `phagescout-elane-pwm-relaxed-aligned` | `phagescout-ctsg-pwm-relaxed-aligned` | Sum of P5–P4′ positional log2 enrichment weights |
| `phagescout-elane-peptide-relaxed-unaligned` | `phagescout-ctsg-peptide-relaxed-unaligned` | Mean log2 fold change of matching released five-mers spanning a bond |
| `phagescout-elane-peptide-relaxed-aligned` | `phagescout-ctsg-peptide-relaxed-aligned` | Log2 fold change of the first matching released aligned wildcard profile |
| `phagescout-elane-peptide-deseq2` | `phagescout-ctsg-peptide-deseq2` | Mean DESeq2 log2 fold change of matching released five-mers; optional full tables |

### Optional full DESeq2 peptide profiles

```sh
mhctools fetch phagescout
mhctools ls phagescout --json
mhctools cleavage --sequence AAACD --model phagescout-elane-peptide-deseq2
```

Each enzyme's full table contains 652,509 distinct five-mers, including
negative DESeq2 log2 fold changes. These are peptide lookups, separate from
the DESeq2 PWM summaries. The two unmodified files total about 32.6 MB and
are excluded from the package. Fetching streams the fixed Zenodo record
21387981, checks file sizes and SHA-256, and atomically installs both tables
with source provenance, attribution and the CC BY 4.0 license. A failed
download leaves no partial installation. An invalid existing destination
must be moved aside before fetching again.

The [managed data directory](../artifacts.md#where-snapshots-live) can be
selected with `MHCTOOLS_DATA_DIR` or `mhctools fetch phagescout --data-dir PATH`.
For existing copies of the exact deposited files, set `PHAGESCOUT_HOME` to
their directory, or pass `PhageScout(enzyme="CTSG", profile="peptide-deseq2",
profile_dir=PATH)`. Explicit `--data-dir` overrides `PHAGESCOUT_HOME` for
fetch/inventory; inference uses `PHAGESCOUT_HOME` or `MHCTOOLS_DATA_DIR`.
The selected file is verified before inference. Invalid explicit paths raise
an actionable error. `mhctools cleavage --list-models` discovers both models
without reading their tables; only explicit selection loads the lookup.

Bond `b` splits `sequence[:b] | sequence[b:]`. The aligned ninth-position
window puts residue `b` at column 5, the inferred P1 anchor; the whole window
is P5–P4′. DECIPHER alignment of enriched phage five-mers infers this anchor.
It does not experimentally locate a scissile bond. An unaligned five-mer can
contribute to four candidate bonds; the adapter retains that ambiguity.

Five-mer modes average up to four complete windows spanning each bond.
Aligned PWM scoring omits missing terminal flanks, following the native
calculation, and explains that omission in the site's reason. Aligned
peptide-profile matching requires all nine positions. Dashes in the released
profiles are wildcards, including internal gaps; when multiple rows match,
the first **source-ordered** row supplies the score. Row order is preserved.

Only canonical uppercase L-peptide sequences with assumed free termini are
supported. Modified or unknown terminal chemistry abstains. The computational
context requirements are not absolute enzyme length cutoffs. Enzyme exposure,
subsequent fragments and repeated cleavage are not simulated.

## Missing scores and interpretation

Every available native score is retained, including negative and small PWM
scores. No threshold, min-max transformation, cross-model average or peptide
loss percentage is added. Significant-peptide profiles are sparse: no match
means **unassessed**, not zero or resistance. The result's
`conditions["unassessed_bonds"]` JSON object gives a reason for each omitted
bond; an entirely unassessed input also has `unsupported_reason`. API, CLI
and batch epitope overlays use the same canonical bond coordinates.

The assay screened immobilized randomized five-mer phages against active
elastase (5 nM, 1.5 hours) or cathepsin G (10 nM, 6 hours), with FLAG-eluted
controls. These times describe the source assay. They do not predict a new
peptide's half-life. Transfer from phage display to free vaccine peptides,
cellular processing and whole-serum survival remains unvalidated.

## Source reproduction and remaining models

Offline fixtures retain the authors' deposited scores on mature alpha-1
antitrypsin (P01009, 393 bonds) and 14-3-3 theta (P27348, 244 bonds), with
matching UniProt sequence identities. All ten models reproduce their source
values and missing scores: 3,822 numeric PWM scores and 36 numeric
significant-peptide-profile scores, plus every unassessed source position.
The optional full tables additionally reproduce both enzymes' deposited
`pep_deseq2` values and missing masks on these proteins. CI fetches the real
files and runs that comparison without skips. These are
implementation-conformance examples, not independent accuracy validation or
experimentally measured non-cleavages.

The optional distance-weighted PWM variant is excluded: the source notebook
multiplies the accumulated sum at successive positions, rather than weighting
each residue independently. Substituting an ordinary distance-weighted sum
would change its released scores.

All four deposited phage-only XGBoost boosters load in the official R
XGBoost 3.2.1.1 runtime. They require 18 features: raw, within-mature-protein
min-max, and sample-standard-deviation Z scores for six sequence features.
Their saved attributes do not contain the training median vectors used to
impute sparse profile matches. New-sequence classifier inference needs that
preprocessing state and the exact mature-protein normalization context.
The adapter exposes the reproduced native sequence scores while that work
remains tracked in [known gaps](../known-gaps.md#cleavage-validation-and-coverage).
Structural classifiers additionally require actual applicable structure
features. Balanced-classifier output and normalized landscape scores are
not calibrated serum-loss probabilities.
