# Antigen processing

Antigen processing includes the generation, trimming, destruction, and
transport of peptides before MHC loading. Peptidases are the enzymes that
hydrolyze peptide bonds within these pathways. The relevant compartment depends
on the pathway:

| Location | Role | Models and evidence |
|---|---|---|
| Cytosol | Proteasomal peptide generation and further trimming or degradation | [Pepsickle](#pepsickle), [NetChop](#netchop), [cytosolic peptidases](../cleavage/choosing.md#cytosolic-trimming-and-degradation) |
| ER | N-terminal trimming of class I precursors after TAP transport | [ERAMER](#eramer), [ERAP2 motif evidence](../cleavage/choosing.md#er-trimming) |
| Endosomes and lysosomes | Class II antigen processing and some class I cross-presentation routes | [NetCleave](#netcleave), [IRAP source observations](../cleavage/choosing.md#endosomal-cross-presentation), [coverage limits](../cleavage/validation.md#apc-endolysosomal-enzymes) |
| Cell surface and extracellular fluids | Peptide turnover before or outside cellular uptake | [Extracellular peptidase models](../cleavage/choosing.md#extracellular-peptide-degradation) |

Experimental work establishes roles for
[proteasomes in class I peptide generation](https://pubmed.ncbi.nlm.nih.gov/8087844/),
[ERAP1 in ER trimming](https://pubmed.ncbi.nlm.nih.gov/12436110/), and
[cathepsin S in class II processing](https://pubmed.ncbi.nlm.nih.gov/9616206/).
[IRAP-mediated cross-presentation](https://pubmed.ncbi.nlm.nih.gov/19498108/)
also illustrates why an endosomal location does not imply class II alone.
Extracellular turnover can affect which peptides reach a cell; its prediction
does not establish uptake or presentation.

## Choose the output you need

Use the peptide-level predictors below for a processing score, TAP transport,
or precursor trimming. Use the [peptidase activity API](../cleavage/index.md)
when you need individual bonds, enzyme identity, and experimental or motif
evidence. Both interfaces include [Pepsickle](#pepsickle) and
[ERAMER](#eramer); the distinction is the result format and endpoint.

The [model-selection guide](../cleavage/choosing.md) recommends models for
specific biological questions. For peptide stability and delivery endpoints,
see [peptide half-life](peptide-pk.md) and [uptake and tissue exposure](../exposure-results.md).

## Peptide-level results

The following predictors are allele-independent: `Prediction.allele` is empty.
Each emits one kind, available through a dedicated result accessor:

| Kind | Accessor | Predictors |
|---|---|---|
| `proteasome_cleavage` | `result.cleavage` | [Pepsickle](#pepsickle), [NetChop](#netchop), [NetCleave](#netcleave) (class I) |
| `endolysosomal_cleavage` | `result.endolysosomal_cleavage` | [NetCleave](#netcleave) (class II) |
| `tap_transport` | `result.tap_transport` | [DeepTAP](#deeptap) |
| `erap_trimming` | `result.erap_trimming` | [ERAMER](#eramer) |

The `antigen_processing` kind is specifically [MHCflurry](binding.md#mhcflurry)'s combined processing
score (`result.processing`); it is not the name of a compartment or a complete
simulation of antigen processing. See [MHCflurry](binding.md#mhcflurry).

Proteasome predictors summarize cleavage evidence into a peptide-level score.
C-terminal scores depend on the residues after the peptide. Pass `c_flanks=`
and `n_flanks=` where supported, or scan a protein with `predict_proteins()`.

## Pepsickle

[Pepsickle](#pepsickle) is the in-vivo epitope proteasome model of [Weeder et al.
(2021)](https://doi.org/10.1093/bioinformatics/btab628). Install the package
with `pip install pepsickle`; there is nothing to fetch.

```python
from mhctools import Pepsickle

predictor = Pepsickle()
results = predictor.predict(["SIINFEKL"], c_flanks=["GGG"])
results[0].cleavage.score        # C-terminal cleavage score

by_protein = predictor.predict_proteins({"TP53": "MEEPQSDPSVEPPLSQETFS"},
                                        peptide_lengths=[9])
```

Without `c_flanks` the C-terminal position has no downstream context and
scores 0.0, so a flank-free call is not a prediction that the peptide is not
cleaved. For a declared human organism context, prefer
`Pepsickle(human_only=True)`; use `human_only=False` for nonhuman, mixed, or
uncertain context. This is a population-selection policy, not an accuracy
claim: [upstream](https://github.com/pdxgx/pepsickle) labels human-only
experimental and notes its smaller training set. All-mammal does not establish
validity for arbitrary non-mammalian organisms. Use
`isolate_subprocess=True` to run inference in a subprocess (this avoids macOS
OpenMP crashes). The epitope models are proteasome-type agnostic. The separate
digestion-trained families support explicit constitutive/immunoproteasome
selection. These families are exposed per bond through the
[cleavage API](../cleavage/models.md#proteasome-models).

## NetChop

[NetChop](#netchop) 3.1 ships as 32-bit x86 Linux binaries. On macOS and ARM Linux,
NetChop automatically runs a user-supplied licensed installation in a
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

```python
from mhctools import NetChop

predictor = NetChop()                       # NETCHOP_HOME / NETMHC_BUNDLE_HOME / PATH
results = predictor.predict(["SIINFEKL"], c_flanks=["GGG"])
results[0].cleavage.score
```

Choose the NetChop model for the endpoint being interpreted:

| Model | Training endpoint | Interpretation |
|---|---|---|
| `NetChop(model_variant=0)` / Cterm 3.0 | MHC-I ligand C-termini | Ligand-boundary processing proxy; not an isolated protease assay |
| `NetChop(model_variant=1)` / 20S 3.0 | In-vitro proteasome degradation | More direct proteasome-digestion model; no individual catalytic-subunit assignment |

[DTU's model documentation](https://services.healthtech.dtu.dk/services/NetChop-3.1/)
reports that Cterm performs best for CTL epitope boundaries. That performance
statement does not establish better prediction of extracellular SLP turnover.
The two scores can disagree because their training endpoints differ. Neither
predicts DPP4 activity: DPP4's [qPISA model](../cleavage/models.md#human-dpp4-qpisa)
assesses the exposed N-terminal dipeptide bond and reports log2 substrate
depletion, rather than a NetChop score. Keep model, native units, bond topology,
and exposure requirements together in site-level reports.

## NetCleave

[NetCleave](#netcleave) works differently from [Pepsickle](#pepsickle) and [NetChop](#netchop): it emits a **single
C-terminal cleavage score per peptide**, and it covers **both** the MHC-I
proteasomal (`NetCleave_I` → `proteasome_cleavage`) and MHC-II endolysosomal
(`NetCleave_II` → `endolysosomal_cleavage`) pathways.

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

[NetCleave](#netcleave)'s own paper reports that class-II C-terminal cleavage is a much weaker
signal than class I (AUC ~0.66 vs ~0.91), so weigh `endolysosomal_cleavage`
scores accordingly.

## DeepTAP

TAP (transporter associated with antigen processing) shuttles cytosolic
peptides into the ER for MHC-I loading. It is a distinct step from proteasomal
cleavage.

[DeepTAP](#deeptap) is a BiGRU that scores each peptide once, independent of allele,
like the cleavage predictors. It emits one `tap_transport` prediction per
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

DeepTAP's evaluation is self-reported, and no independent TAP benchmark exists
for any tool. Treat the score as a pathway signal for prioritizing, not a
validated one.

## ERAMER

ERAP1 trims the N-termini of 9–16mer precursor peptides in the ER down to the
8–10mers MHC-I presents, the step between TAP transport and MHC loading.

[ERAMER](#eramer) scores a precursor by averaging a per-length position-weight-matrix
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

ERAMER's evaluation is self-reported and ERAP1 trimming is inherently noisy.
Treat the score as a pathway prior, not a validated one.
