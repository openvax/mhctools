# Antigen-processing predictors

Predictors for the steps between a protein and an MHC molecule: proteasomal and
endolysosomal cleavage, TAP transport, and ERAP1 trimming. All of them are
**allele-independent**, so `Prediction.allele` is empty, and each emits a
single kind that you read through a dedicated accessor:

| Kind | Accessor | Predictors |
|---|---|---|
| `proteasome_cleavage` | `result.cleavage` | [Pepsickle](#pepsickle), [NetChop](#netchop), [NetCleave](#netcleave) (class I) |
| `endolysosomal_cleavage` | `result.endolysosomal_cleavage` | [NetCleave](#netcleave) (class II) |
| `tap_transport` | `result.tap_transport` | [DeepTAP](#deeptap) |
| `erap_trimming` | `result.erap_trimming` | [ERAMER](#eramer) |

MHCflurry additionally emits an `antigen_processing` score (`result.processing`);
see [MHCflurry](binding.md#mhcflurry). For per-bond, per-enzyme evidence
rather than a per-peptide score, use the [cleavage API](../cleavage/index.md).

Proteasome predictors score each cleavage position and aggregate them into one
peptide-level number (see `ProcessingPredictor` and `ProteasomePredictor`). A
C-terminal cleavage score depends on the residues that follow the peptide, so
**pass flanks** (`c_flanks=`, and `n_flanks=` where supported) or scan a
protein with `predict_proteins()`.

## Pepsickle

`Pepsickle` is the in-vivo epitope proteasome model of [Weeder et al.
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
cleaved. Use `human_only=True` for the human-trained model, and
`isolate_subprocess=True` to run inference in a subprocess (this avoids macOS
OpenMP crashes). The same models, with explicit constitutive/immunoproteasome
selection, are exposed per bond through the
[cleavage API](../cleavage/models.md#proteasome-models).

## NetChop

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

```python
from mhctools import NetChop

predictor = NetChop()                       # NETCHOP_HOME / NETMHC_BUNDLE_HOME / PATH
results = predictor.predict(["SIINFEKL"], c_flanks=["GGG"])
results[0].cleavage.score
```

## NetCleave

`NetCleave` works differently from `Pepsickle` and `NetChop`: it emits a **single
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

NetCleave's own paper reports that class-II C-terminal cleavage is a much weaker
signal than class I (AUC ~0.66 vs ~0.91), so weigh `endolysosomal_cleavage`
scores accordingly.

## DeepTAP

TAP (transporter associated with antigen processing) shuttles cytosolic
peptides into the ER for MHC-I loading. It is a distinct step from proteasomal
cleavage.

`DeepTAP` is a BiGRU that scores each peptide once, independent of allele,
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

ERAMER's evaluation is self-reported and ERAP1 trimming is inherently noisy.
Treat the score as a pathway prior, not a validated one.
