# Peptidase models

The model panel behind `predict_cleavage()` and `mhctools cleavage`. Which
enzyme, where it acts, what pattern it assesses and how strongly a match should
be read is below; what the result states mean is in [reading the
evidence](evidence.md).

For model selection by biological endpoint and compartment, see
[choosing processing and peptidase models](choosing.md).

## Running the panel

```sh
mhctools cleavage --list-models
mhctools cleavage --list-models --json
mhctools cleavage --sequence RPPGFSPFR --model app2-xp --model cpn-basic
mhctools cleavage --sequence VPYGSFKHV --compartment cytosol --out cleavage.json
mhctools cleavage --sequence HAEGTFTSD --model dpp4-qpisa --n-term acetylated
mhctools cleavage --sequence SIINFEKL --model pepsickle-in-vivo-human-only
```

```python
from mhctools import predict_cleavage
results = predict_cleavage("TSGPNQ", models=["fap-endo-gp", "prep-pro"])
```

### Define a sourced motif rule

Use `motif_rule()` to make a rule outside the built-in catalog. It returns a
predictor, so no registry edit is needed:

```python
from mhctools import PeptideInput, motif_rule

nep = motif_rule(
    name="study-nep", enzyme="MME",
    references="https://pubmed.ncbi.nlm.nih.gov/2417254/",
    topology="internal", right="[FILY]",
    strictness="preferred",
    basis="Reported substance P cuts precede Phe or Leu",
    assay="Purified human kidney NEP acting on substance P",
    limitations="Other bonds and substrates can be cleaved",
    c_chemistry=("free", "amidated"),
)
peptide = PeptideInput("RPKPQQFFGLM", c_term="amidated")
result = nep.predict(peptide)
mask = result.pattern_mask  # (0, 0, 0, 0, 0, 1, 1, 0, 1, 0)
matches = [bond for bond, match in enumerate(mask, 1) if match]  # [6, 7, 9]
```

`left` and `right` are Python regexes for the residues before and after each
candidate bond. Omit a flank or pass `None` when it is unconstrained. In the
example, `[FILY]` means one of those four residues immediately after the bond;
`FILY` without brackets would require all four in order. For a terminal
enzyme, use `topology="n_terminal"` or `"c_terminal"` and state how many
residues it removes with `removed=`.
`motif_rule()` also accepts optional UniProt accession, species, compartment,
length, terminal chemistry and activation constraints. A custom rule can omit
`strictness` and `basis`; if it assigns a grade, it must give the source
observation behind that grade. Grades do not change pattern matches. A custom
rule reports binary motif matches, never numerical cleavage probabilities.
The built-in `mme-hydrophobic` rule already covers the example above; custom rules should
retain the assay and limits supported by their own references.

The default panel evaluates all 39 built-in models and returns separate
results. `--model` and `--sequence` can be repeated. `--list-models` prints a
compact discovery table; add `--json` for its full machine-readable catalog.
Nine bundled [ITCell cathepsin profiles](itcell.md) add human B/S internal
specificity and H initial N-terminal trimming. Select their 15/60/240-minute
source profiles explicitly; these times are assay-training scope, not a new
peptide's predicted degradation time.

| ITCell models | Mechanism | Native score |
| --- | --- | --- |
| `itcell-catb-15`, `itcell-catb-60`, `itcell-catb-240` | Human cathepsin B internal specificity | Sum of log2 profile/background ratios |
| `itcell-cats-15`, `itcell-cats-60`, `itcell-cats-240` | Human cathepsin S internal specificity | Sum of log2 profile/background ratios |
| `itcell-cath-15`, `itcell-cath-60`, `itcell-cath-240` | Human cathepsin H initial N-terminal trimming | Sum of log2 profile/background ratios |

Ten bundled [PhageScout sequence profiles](phagescout.md) add separate human
ELANE and CTSG native scores. The exact names are
`phagescout-elane-pwm-deseq2`, `phagescout-elane-pwm-relaxed-unaligned`,
`phagescout-elane-pwm-relaxed-aligned`, `phagescout-elane-peptide-relaxed-unaligned`,
`phagescout-elane-peptide-relaxed-aligned`, `phagescout-ctsg-pwm-deseq2`,
`phagescout-ctsg-pwm-relaxed-unaligned`, `phagescout-ctsg-pwm-relaxed-aligned`,
`phagescout-ctsg-peptide-relaxed-unaligned` and `phagescout-ctsg-peptide-relaxed-aligned`.
These are phage-derived recognition features with inferred aligned P1 anchors.
Sparse profile non-matches are unassessed; neither their native scores nor
negative PWM values are serum-loss percentages.

Two optional full DESeq2 peptide lookups, `phagescout-elane-peptide-deseq2`
and `phagescout-ctsg-peptide-deseq2`, retain the mean native log2 fold change
over matching complete five-mers. Install their verified data separately
with `mhctools fetch phagescout`; catalog discovery does not load the tables.

The eight optional [Pepsickle](../predictors/processing.md#pepsickle) models cover epitope and C/I digestion families;
see [proteasome models](#proteasome-models).
Missing assets and uninspected external runtimes are listed as unresolved.
Predictions carry SHA-256 identities for the actual runtime's weights,
inference code, feature code and dependency metadata.
Prediction JSON uses schema version 3: its top-level `models` object stores each
full provenance record once, keyed by model name, and each item in `results`
references that name in its `model` field. Motif sites carry `match: 0|1`,
with the rule documented once in model metadata; this is a binary pattern
indicator, not a cleavage probability. The output retains non-matches and
unsupported results, source coordinates and chemistry. Numerical scores are
written with six significant digits. `--source-start` is a zero-based offset
shared by the supplied inputs; use separate calls when fragments have different
offsets.
The `--list-models --json` catalog uses schema version 2 and omits metadata
fields that do not apply to a model. Schema-v2 prediction JSON remains readable
by the vaccine-report command.

Compartment filtering uses exact, conservative enzyme-location annotations.
`serum`, `plasma`, `extracellular`, `cytosol`, `endosome` and `er` are
distinct. The `extracellular` filter is useful for broader candidate
screening; `serum` is not an exhaustive inventory of everything potentially
present in a serum sample. Presence, concentration, activation, inhibitors and
exposure are not inferred. No combination is converted into overall stability.

## The models

Strictness grades are explained in [reading the evidence](evidence.md#how-strict-is-each-motif).

| Model | Assessed recognition pattern | Strictness | Main scope and primary evidence |
| --- | --- | --- | --- |
| `dpp4-qpisa` | N-terminal P2-P1\|P1′ score | scored | Human DPP4; quantitative model described above |
| `ace-dipeptidyl` | C-terminal \|non-Pro–non-Asp/Glu | required | Ordinary human ACE dipeptide activity; [angiotensin assays](https://doi.org/10.1042/BJ20040634) |
| `mme-hydrophobic` | Selected P1′ residues Phe/Ile/Leu/Tyr | preferred | Human neprilysin; [kidney peptide assays](https://pubmed.ncbi.nlm.nih.gov/6349683/); incomplete whole-sequence specificity |
| `cpb2-basic` | C-terminal \|Lys/Arg, explicitly active enzyme | required | Human TAFIa; [chemerin cleavage](https://pmc.ncbi.nlm.nih.gov/articles/PMC2613638/); unknown/zymogen/inactive states abstain |
| `cpn-basic` | C-terminal \|Lys/Arg | required | Human plasma CPN; [Oshima et al. 1975](https://doi.org/10.1016/0003-9861(75)90104-6) |
| `app1-xp` | N-terminal X\|Pro | required | Cytosolic human XPNPEP1; [Cottrell et al. 2000](https://pubmed.ncbi.nlm.nih.gov/11106490/); manganese dependent, distinct gene product from XPNPEP2 |
| `app2-xp` | N-terminal X\|Pro | required | Human XPNPEP2; [Molinaro et al.](https://pubmed.ncbi.nlm.nih.gov/15361070/) |
| `fap-dipeptidyl` | N-terminal X-Ala/Pro\|non-Pro | preferred | Human FAP; [substrate profiling](https://pubmed.ncbi.nlm.nih.gov/16410248/) and [direct Ala/Pro cleavage](https://doi.org/10.1111/j.1742-4658.2011.08051.x). Ala was much slower in the tested cohort. |
| `fap-endo-gp` | Gly-Pro\|non-Pro | required | FAP endopeptidase; [substrate profiling](https://pubmed.ncbi.nlm.nih.gov/16480718/), [prime-side constraint](https://pubmed.ncbi.nlm.nih.gov/22750443/) |
| `enpep-acidic` | N-terminal Asp/Glu\|X | preferred | Aminopeptidase A; [human specificity study](https://pubmed.ncbi.nlm.nih.gov/23888046/); calcium and sequence affect activity |
| `anpep-ala` | N-terminal Ala\|X preference | permissive | Aminopeptidase N; [human structure/biochemistry](https://pubmed.ncbi.nlm.nih.gov/22932899/); many other substrates omitted |
| `dpp8-xp-xa` | N-terminal X-(Pro/Ala)\|non-Pro | required | Cytosolic DPP8; [characterization](https://pubmed.ncbi.nlm.nih.gov/11012666/), [degradomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC3656252/) |
| `dpp9-xp-xa` | N-terminal X-(Pro/Ala)\|non-Pro | required | Cytosolic DPP9; [antigen processing](https://pubmed.ncbi.nlm.nih.gov/19667070/), same degradomics study |
| `tpp2-tripeptidyl` | N-terminal tripeptide, non-Pro P1 and P1′ | permissive | Cytosolic TPP2; [RU1 precursor processing](https://doi.org/10.4049/jimmunol.169.8.4161); topology, not selectivity |
| `npepps-n-terminal` | First bond lacking Gly/Pro context | permissive | Puromycin-sensitive aminopeptidase; same RU1 study; broad enzyme, weak flag |
| `prep-pro` | Internal X-Pro\|X | required | PREP/POP, conservative 4–30-residue domain; [human profiling](https://pubmed.ncbi.nlm.nih.gov/22750443/); flanking preferences omitted |
| `erap2-basic` | N-terminal Arg/Lys\|X preference | preferred | ERAP2 in ER; [biochemistry](https://pubmed.ncbi.nlm.nih.gov/12799365/), [peptide structures](https://pubmed.ncbi.nlm.nih.gov/26381406/); not a full-context predictor |
| `thop1-observed` | Exact-sequence source lookup | source observations | Cytosolic THOP1; [Knight et al. 1995](https://pubmed.ncbi.nlm.nih.gov/7755557/) |
| `nln-observed` | Exact-sequence source lookup | source observations | Cytosolic NLN; [human neurolysin structures and LC-MS](https://pubmed.ncbi.nlm.nih.gov/39117724/) |
| `lnpep-observed` | Exact-sequence source lookup | source observations | Endosomal LNPEP/IRAP; [Georgiadou et al. 2010](https://pubmed.ncbi.nlm.nih.gov/20592285/) |
| `eramer-step` (optional) | Initial N-terminal bond; length-specific PWM score | scored | ERAP1 in ER; [ERAMER](https://pubmed.ncbi.nlm.nih.gov/38925438/), 9–16 residues |

The motif models return decisions without numerical scores. For example,
APP removes the first residue of `RPPGFSPFR` at `R|PPGFSPFR`; DPP-like
activity would remove two residues and is a separate assessment. FAP's
endopeptidase rule can also assess an N-acetylated input, supported by its
blocked-substrate activity. ACE also accepts N-acetylation, and MME accepts
C-amidation. Other rules conservatively accept free termini only;
an unsupported modified input does not establish that
its bonds resist enzymatic cleavage.

## Proteasome models

Eight optional [Pepsickle](../predictors/processing.md#pepsickle) models cover the epitope and digestion families, with
explicit constitutive (C) or immunoproteasome (I) selection where the family has
it. They are excluded from the default panel and selected by exact name:

| Model | Family | Context | Proteasome |
|---|---|---|---|
| `pepsickle-in-vivo-human-only` | epitope-trained neural ensemble | 8 residues before and after P1 | agnostic |
| `pepsickle-in-vivo-all-mammal` | epitope-trained neural ensemble | 8 residues before and after P1 | agnostic |
| `pepsickle-in-vitro-2-human-only-constitutive` | digestion-trained neural ensemble | 3 residues before and after P1 | C |
| `pepsickle-in-vitro-2-human-only-immunoproteasome` | digestion-trained neural ensemble | 3 residues before and after P1 | I |
| `pepsickle-in-vitro-2-all-mammal-constitutive` | digestion-trained neural ensemble | 3 residues before and after P1 | C |
| `pepsickle-in-vitro-2-all-mammal-immunoproteasome` | digestion-trained neural ensemble | 3 residues before and after P1 | I |
| `pepsickle-in-vitro-all-mammal-constitutive` | gradient-boosted digestion model | needs the isolated scikit-learn 0.23.2 runtime (see [below](#the-legacy-gradient-boosted-runtime)) | C |
| `pepsickle-in-vitro-all-mammal-immunoproteasome` | gradient-boosted digestion model | needs the isolated scikit-learn 0.23.2 runtime | I |

The Python Pepsickle wrapper takes `model_type` and `proteasome_type` (`C`
or `I`) for digestion models. It rejects a proteasome type for the epitope
model and `human_only=True` for gradient boosting, which upstream ignores.
Model metadata identifies actual weights, code, population and model family.
The upstream endpoint sentinel is excluded from canonical bond results.

### The legacy gradient-boosted runtime

Provision it with Docker:

```sh
python scripts/setup_test_backends.py pepsickle
source env/test-backends/activate.sh
```

This pins Python 3.8.20, scikit-learn 0.23.2 and its companion packages in a
separate container. Its generated `PEPSICKLE_GB_PYTHON` launcher selects that
runtime only for gradient boosting and runs inference with networking
disabled. Neural models keep their existing runtime. A separately managed
interpreter can be selected with `Pepsickle(python_executable=...)` or
`PEPSICKLE_PYTHON` for every model family.

Prediction provenance comes from the selected interpreter, including package
versions and actual code/weight hashes. A changed identity between inspection
and inference causes failure. Catalog listing does not start external
runtimes; their identities remain unresolved until the predictor is selected.
No model is silently substituted when a runtime is missing or incompatible.

Adapter agreement with upstream inference establishes implementation
conformance, not accuracy for tumor/APC processing or vaccine-peptide survival.

## Human DPP4 qPISA

The [Gudipati et al. 2024 paper](https://doi.org/10.1038/s44320-024-00071-4)
reports a model fitted to substrate depletion by purified human DPP4. The
source assay used tryptic HeLa peptides in HEPES pH 7.4, at 21 C for 4 hours.
The implementation independently evaluates the three rearranged terms in
Dataset EV2: `P1 + P2:P1 + P1:P1-prime`, where the first three peptide
residues are P2, P1 and P1-prime. Only bond 2 is assessed.

Higher scores indicate greater predicted log2 depletion relative to buffer
control in that assay. Negative values are retained. Scores are **not
cleavage probabilities, serum half-lives or stability ranks across enzymes**.
Compartment metadata records where the enzyme can act, not where the model
has been calibrated. Physiological exposure, structure, concentration and
competing enzymes are not modeled.

All triplets with complete coefficients can be evaluated, including those
without Pro or Ala at P1. Of 8,000 canonical triplets, 6,420 have complete
coefficients; the other 1,580 return an explicit missing-coefficient reason.
Evaluability does not establish that an individual triplet occurred in the
training data. The related C. elegans DPF-3 model is not used for human DPP8/9.

### Parameter provenance

`mhctools/data/dpp4_qpisa.json` contains the numeric cells from
`44320_2024_71_MOESM3_ESM.xlsx`, sheet `dpp4_modelParams`. Published `NA`
cells become JSON `null`; no numerical imputation or refitting is performed.
The source workbook SHA-256 is
`ee449da13b5ec66fd6fb08c16203da8c44f2e9c10572a355abdcb985694c6ddc`.
It was retrieved via the [Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11612144/supplementaryFiles).
The article assigns associated data to [CC0](https://creativecommons.org/publicdomain/zero/1.0/),
unless otherwise credited; Dataset EV2 has no separate credit restriction.
Attribution: Rajani Kanth Gudipati and colleagues, 2024, DOI above. No figures
or upstream R source code are redistributed.

### A bounded purified-DPP4 fraction estimate

`DPP4qPISA().estimate_keane_fraction(peptide, minutes=...)` maps the qPISA
triplet score to a **parent-peptide loss fraction in one purified-enzyme assay**.
It fits log half-life to the eleven 27–42-residue peptide measurements in
[Keane et al. 2011, Tables 1–2](https://doi.org/10.1111/j.1742-4658.2011.08051.x):
0.1 µM recombinant human DPP4 at 37 °C, with MALDI-TOF parent-peptide time
courses. The fit gives each substrate one observation; the published mean
half-lives, SDs, replicate counts, N-terminal
triplets and lengths are transcribed in `mhctools/data/dpp4_keane2011.json`.
BNP and substance P lack a qPISA coefficient; censored non-cuts are never
assigned a zero rate. The fit treats the reported amide-terminated source
peptides as part of the training set, but estimates only free-terminus inputs
accepted by `DPP4qPISA.predict()`.

```python
from mhctools import DPP4qPISA

# Synthetic 30-residue illustration, not a measured vaccine peptide.
estimate = DPP4qPISA().estimate_keane_fraction("HAE" + "G" * 27, minutes=10)
print(estimate["parent_lost_fraction"])  # about 0.33 in this assay model
print(estimate["parent_remaining_fraction"])  # about 0.67
```

The method evaluates only the first, exposed DPP4 bond (bond 2), lengths
27–42, supported unmodified qPISA inputs, and scores within the training
range. Otherwise it returns `status="unassessed"` with a reason. It retains
the native qPISA score separately. A first-order loss assumption converts
the fitted half-life into a fraction at the requested time:
`parent_remaining_fraction = 2 ** (-minutes / estimated_half_life_minutes)`.
This is a modeling assumption; Keane's half-lives are experimental parent-loss
summaries, not proof that every new peptide follows one exponential.

There are only eight unique triplets among the eleven training substrates.
Leaving out each complete triplet group in turn gives a mean absolute error
in log half-life corresponding to a **2.8-fold typical error**. That is a
validation error summary, not a
confidence interval or a guarantee for new substrates. qPISA was trained in
a separate, median-normalized 21 °C assay; its score cannot itself be read
as an absolute fraction. Other sequence positions, substrate structure,
enzyme amount and matrix all limit transfer. In particular, this estimate
is **not** a serum, lymph, IM-vaccine or APC-arrival fraction and is not
automatically passed to `assess_vaccine_windows()`.

### The corresponding purified-FAP estimate

The same [Keane et al. 2011 Table 2](https://doi.org/10.1111/j.1742-4658.2011.08051.x)
reports nine numerical parent half-lives for 11–42-residue hormones and
neuropeptides incubated with **0.2 µM purified human FAP at 37 °C**. Five
substrates with Ala at the second residue had measured half-lives of about
15.5–39.1 hours; four with Pro ranged from 5.78 to 60.2 minutes. These are
different enzyme concentrations from the DPP4 assay and must not be compared
as intrinsic in-vivo clearance rates.

`estimate_assay_parent_fraction()` gives a common entry point for the two
bounded estimates. For FAP it uses the geometric mean measured half-life of
the Ala or Pro group, then assumes first-order parent loss. It accepts only
11–42-residue peptides with one of the **eight measured N-terminal triplets**,
an exposed free N terminus, and a free or amidated C terminus. The rest of
the sequence is unmodeled; this is a rough cohort prior, not a fit to full
peptide sequence or a claim that an arbitrary motif match has this rate.

```python
from mhctools import estimate_assay_parent_fraction

peptide = "HAE" + "G" * 27  # synthetic illustration, not a measured substrate
for enzyme in ("DPP4", "FAP"):
    result = estimate_assay_parent_fraction(peptide, enzyme=enzyme, minutes=10)
    print(enzyme, result["status"], result.get("parent_lost_fraction"))
```

The FAP Ala group has 5 source substrates and 1.4-fold typical
leave-one-substrate-out error; the Pro group has 4 and **3.0-fold** typical
error, with PYY a conspicuous outlier. Those checks reuse N-terminal groups
already seen during fitting. They do **not** validate transfer to a novel
sequence with the same triplet, and the source includes censored non-cuts
which are excluded rather than assigned zero rates. The output names the
source assay, model, estimated half-life, lost and remaining **parent**
fractions, and group validation. Unassessed inputs contain a reason and no
fabricated fraction. Neither estimate is core-epitope survival, serum
half-life, IM exposure or APC arrival; broader calibration remains [#587](https://github.com/openvax/mhctools/issues/587).

For names, keep three distinct things visible: `fap-dipeptidyl` is a
**binary recognition pattern**, `dpp4-qpisa` is a **native numerical
substrate-depletion score**, and `fap-keane2011-parent-loss` or
`dpp4-qpisa-keane2011-parent-loss` is an **assay-specific fraction at a
specified time**. A bare name like `FAP_HalfLife` would hide the assay and
suggest a transferable biological half-life. The stable Python call is
`estimate_assay_parent_fraction(peptide, enzyme=..., minutes=...)`; the
source and endpoint remain explicit in its returned `model` and `assay`.

Other enzymes need different evidence before they get a time-specific
parent-loss estimate for an arbitrary canonical peptide. A [human CatS
substitution series](https://pmc.ncbi.nlm.nih.gov/articles/PMC2662293/)
reports digestion and kinetic constants for internally quenched peptides
with dye, D-amino-acid and amide modifications around one scaffold; treating
its fitted constants as rates for unmodified vaccine peptides would silently
change the substrate chemistry. A [human neprilysin kinetic
panel](https://pmc.ncbi.nlm.nih.gov/articles/PMC4121237/) reports constants
for a small set of named substrates, but no validated sequence-to-rate
mapping for novel windows. Their motif or native scores can still be used for
the endpoints each measures; neither is an absolute per-minute fraction.

## Serum and extracellular candidates

```sh
mhctools cleavage --sequence DRVYIHPFHL --model ace-dipeptidyl --model mme-hydrophobic
mhctools cleavage --sequence YFPGQFAFSK --model cpb2-basic --enzyme-state CPB2=active
mhctools benchmark --reference-cleavage serum --out serum-reference.json
```

In Python, supply `enzyme_states={"CPB2": "active"}` to `predict_cleavage`,
or `get_cleavage_model("cpb2-basic", enzyme_state="active")`. The result
records that assumption in `conditions`. CPB2 requires proteolytic activation;
its presence as a zymogen does not establish activity. `unknown`, `zymogen`
and `inactive` return no assessed sites. [Activation experiments](https://doi.org/10.1074/jbc.274.49.35046)
also show that the surrounding coagulation environment matters.

ACE's ordinary rule assesses only the bond before the last two residues.
It recognizes angiotensin I cleavage at bond 8 and the [N-acetyl-SDKP](https://doi.org/10.1038/srep13742)
bond 2. It abstains on amidated substance P even though human ACE is known
to cleave that peptide at bonds 8 and 9 through exceptional processing.
MME flags selected hydrophobic residues after a bond, including substance P's
reported bonds 6, 7 and 9; it accepts that peptide's C-terminal amide.
These [human enzyme observations](https://pubmed.ncbi.nlm.nih.gov/2417254/)
do not establish a transferable model of all substrates. The MME rule's
2–30-residue domain is a conservative implementation scope; longer substrates
can exist. Use the `extracellular` filter to include MME: location annotations
do not assert that purified-kidney specificity was calibrated in serum.

The packaged serum-candidate reference contains 17 source observations,
including two explicitly reported non-cleavages. Fifteen are assessable by
the corresponding rules; the two exceptional ACE/substance-P observations
abstain. Cross-enzyme combinations are unassessed. The source IDs, chemical
forms and known or unspecified experimental conditions are retained.
These small reproduction controls do not establish specificity, serum
half-life performance, or validity on long peptides. The report explicitly
shows missing serum-half-life and long-peptide evidence.

## Cytosolic and endosomal candidates

```sh
mhctools cleavage --sequence VPYGSFKHV --compartment cytosol
mhctools cleavage --sequence KSLYNTVATL --compartment endosome
mhctools cleavage --sequence RPPGFSPFR --model thop1-observed --model app1-xp
mhctools benchmark --reference-cleavage intracellular --out intracellular-reference.json
```

Antigen-processing peptidases are the reason this panel exists, but three of
them publish their specificity as whole-substrate outcomes rather than as a
transferable pattern. Inventing a motif from those papers would misrepresent
them, so `thop1-observed`, `nln-observed` and `lnpep-observed` are **source
references** instead: an exact sequence and chemical form returns what the
experiment reported, and anything else abstains. There is no nearest-neighbour
matching and no extrapolation. Results carry a `substrate_observation` of
`cleavage_reported` or `no_cleavage_detected` alongside any bonds the source's
product identities actually pin down. Where a source saw degradation but no
intermediate, cleavage is recorded with no bond rather than guessing one.

A reported non-cleavage never becomes a per-bond label. `no_cleavage_detected`
means that assay saw no loss of that peptide under its own conditions and
detection limits, which is not the same as a bond that cannot be cleaved.

The remaining cytosolic enzymes are ordinary motif rules, and their grades
matter. `app1-xp` is `required`: aminopeptidase P is defined by hydrolysing
the X-Pro bond, so a non-match is real evidence. `tpp2-tripeptidyl` and
`npepps-n-terminal` are `permissive`. Removing three residues describes TPP2's
topology, not which peptides it turns over, and puromycin-sensitive
aminopeptidase is broad enough that its Gly/Pro exclusion is only a hint.
Do not read a TPP2 match as a prediction that the peptide is consumed.

THOP1 and NLN are closely related and are deliberately kept apart. They cleave
neurotensin at different bonds and their specificities can be swapped by
mutating two active-site residues, so neither model's observations transfer to
the other. Only human-enzyme observations are curated: the widely cited
bradykinin, enkephalin and neurotensin results for neurolysin come from rat or
species-unspecified preparations and are excluded rather than relabelled human.

`lnpep-observed` covers IRAP as an endosomal cross-presentation candidate, not
a second ER enzyme. Its source digested peptides at pH 8.0 with purified
enzyme, so the records describe that experiment, not an acidified endosome.

The packaged intracellular reference contains 42 source observations across
four models, including six reported non-cleavages. Forty-one reproduce; the
amidated substance P record abstains because C-amidation is outside the
aminopeptidase P rule's documented input domain. One resistant IRAP precursor
from Georgiadou 2010 is excluded entirely: the paper prints `DIRSSVQNKL` in
its results and Table I but `DIRSSQVNKL` in the Figure 3F caption, and an
exact-sequence catalog cannot silently pick one. Both strings are retained in
the dataset notice, and the discrepancy is tracked in
[known gaps](../known-gaps.md#cleavage-validation-and-coverage).

Thimet oligopeptidase contributes no non-cleavage record at all. Every
resistant peptide in its source is a hydroxyproline analogue or carries an
N-terminal pyroglutamate, and both are outside the canonical-peptide input
type. Its absence from the negatives is a curation limit, not a finding.

The report shows no evidence for antigen presentation in primary dendritic
cells and none for cleavage measured in cytosol rather than purified enzyme.
Those gaps are the point: nothing here calibrates how long a peptide survives
in the cytosol of an antigen-presenting cell, where the proteasome, competing
aminopeptidases and TAP transport all act at once.

## ERAP1 and existing processing models

```sh
mhctools fetch eramer
# Install openpyxl in the environment if it is not already available.
mhctools cleavage --sequence LAAAFGAAA --model eramer-step
```

Alternatively, construct `ERAMERCleavage(pwm_path="/path/to/PWM.xlsx")` or set
`ERAMER_HOME`. This optional model is listed without loading assets and is
excluded from the default panel. Explicit selection fails clearly if the
external asset or its runtime is missing.

`eramer-step` computes the existing [ERAMER](../predictors/processing.md#eramer) intermediate PWM specificity for
one exposed 9–16-residue precursor, assigning it to bond 1. It does not report
the average of later trimming intermediates. Its version contains the SHA-256
of the actual workbook snapshot used for inference. The GPL-licensed workbook
is loaded at runtime and is not included in the mhctools distribution.

The existing ERAMER cascade API, [NetChop](../predictors/processing.md#netchop), [Pepsickle](../predictors/processing.md#pepsickle) and other proteasome
predictors remain available through their existing interfaces. Pepsickle also
has a canonical `PepsickleCleavage.predict()` facade and the two CLI model names
shown above. Array index `i` maps to internal bond `i + 1`; the upstream final
zero is an endpoint sentinel and is never emitted as a bond. Canonical results
default to subprocess isolation because the upstream package deserializes a
pickle, while still requiring users to trust the installed model artifact.
Processing-model output scales must not be mixed with qPISA scores or motif
decisions.
