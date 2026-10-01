# TCR specificity predictors

Predictors of whether a T-cell receptor recognises a peptide-MHC complex. Their
input is a **peptide plus a `TCR`** (the CDR loops, optionally V/J genes), not
an allele, and they emit `pMHC_TCR_binding`, read with `result.tcr_binding`.

| Predictor | Input | Needs |
|---|---|---|
| [NetTCR](#nettcr) | `(peptide, TCR)` pairs | `mhctools fetch nettcr --accept-license` + a TFLite runtime (`pip install mhctools[nettcr]`) |
| [Tulip](#tulip) | `(peptide, TCR)` pairs plus the presenting allele | `mhctools fetch tulip` + a TULIP-capable Python (`TULIP_HOME`, `TULIP_PYTHON`) |
| [MixTCRpred](#mixtcrpred) | TCRs, scored against one fixed pMHC target per model | `pip install "mhctools[mixtcrpred]"` + `mhctools fetch mixtcrpred --accept-license` |

`NetTCR` and `Tulip` take explicit pairs through `predict_pairs()` or every
peptide-by-TCR combination through `predict(peptides, tcrs)`. They have no
`--mhc-predictor` name on the command line, because the input does not fit that
interface; `MixTCRpred` has its own `mhctools mixtcrpred` subcommand.

## NetTCR

`NetTCR` predicts whether a paired αβ T-cell receptor recognises a (class-I)
peptide. Unlike the MHC-ligand predictors, its input is a peptide plus a `TCR`
(the six CDR loops), not an allele, and it emits the `pMHC_TCR_binding` kind.

NetTCR ships its pretrained weights in its git repository as small TFLite
models. This wrapper runs the pan cross-validation ensemble in-process and does
not need NetTCR's conda environment.

```python
from mhctools import NetTCR, TCR

NetTCR.fetch(accept_license=True)  # downloads only the ~8 MB pan ensemble
predictor = NetTCR()   # resolves NETTCR_DIR / ~/NetTCR-2.2
tcr = TCR(
    cdr1a="NSASQS", cdr2a="VYSSG", cdr3a="VVEGDKVI",
    cdr1b="MGHRA", cdr2b="YSYEKL", cdr3b="ASSHSGYEQF", name="clone1")

# Score explicit (peptide, TCR) pairs...
results = predictor.predict_pairs([("LLWNGPMAV", tcr)])
results[0].tcr_binding.score        # ensemble-mean recognition probability

# ...or every peptide x TCR combination.
results = predictor.predict(["LLWNGPMAV", "GILGFVFTL"], [tcr])
```

## Tulip

```python
from mhctools import Tulip, TCR

tcr = TCR(cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF", name="clone1")
Tulip.fetch()                             # pinned code, tokenizers, and weights
predictor = Tulip()                       # also needs a TULIP-capable Python
results = predictor.predict(["GILGFVFTL"], [tcr], mhc="HLA-A*02:01")
results[0].preds[0].score                 # higher = more likely binding
```

[TULIP-TCR](https://github.com/barthelemymp/TULIP-TCR) is **GPLv3** and pinned
to `transformers==4.32.1`; mhctools is Apache-2.0 and depends on neither torch
nor transformers. The `Tulip` wrapper therefore vendors none of TULIP. It runs
an upstream checkout out-of-process, in an isolated interpreter, via TULIP's own
`predict.py`. `mhctools fetch tulip` obtains the tested code, tokenizers, and
weights; `scripts/setup_tulip_env.sh` can build the separate runtime.

You may instead provide your own checkout and interpreter:

- `TULIP_HOME`: a clone of TULIP-TCR (provides `predict.py`, `src/`,
  tokenizers, and the released `model_weights/`);
- `TULIP_PYTHON`: an isolated Python 3.11 interpreter with `torch` and
  `transformers==4.32.1`. Python 3.11 specifically, so `tokenizers` installs
  from a prebuilt wheel and needs no Rust toolchain.

## MixTCRpred

MixTCRpred has a different, deliberately explicit shape: each checkpoint is
trained for one fixed peptide/MHC target. Its catalog currently contains 146
models (43 marked high-confidence by upstream), spanning human/mouse class I
and II targets.

The input `TCR` stores paired CDR3s and optional `trav`, `traj`, `trbv`, and
`trbj` assignments. Released models use CDR3alpha/beta plus CDR1/2 derived from
the V genes; J assignments are accepted and QC-reported but are not network
inputs.

Each checkpoint is about 31.8 MB. The immutable Zenodo record contains 146
checkpoints totaling 4.64 GB (4.33 GiB); the 43 high-confidence checkpoints
total 1.37 GB (1.27 GiB). To avoid turning every mhctools installation into a
multi-gigabyte download, the pinned upstream artifact includes its two bundled
reference checkpoints, and additional models are fetched individually with
bounded retries, atomic installation, and checksum verification.

```sh
pip install "mhctools[mixtcrpred]"
mhctools fetch mixtcrpred --accept-license
mhctools ls mixtcrpred --models --high-confidence
mhctools fetch mixtcrpred --model A0201_GILGFVFTL
mhctools mixtcrpred --model A0201_GILGFVFTL \
  --input paired-tcrs.csv --out scored-tcrs.csv
```

```python
from mhctools import MixTCRpred, TCR

models = MixTCRpred.catalog()
target = MixTCRpred.resolve_model("GILGFVFTL", "HLA-A*02:01")
predictor = MixTCRpred(target.name)
tcr = TCR(
    cdr3a="CAGASGNTGKLIF", cdr3b="CASSIRASYEQYF",
    trav="TRAV27", traj="TRAJ42", trbv="TRBV19", trbj="TRBJ2-6",
)
prediction = predictor.predict_tcrs([tcr])[0]
prediction.score
prediction.percentile_rank
```

`annotate_dataframe()` and the CSV command retain the original table and add
the raw score, percentile rank, fixed target metadata, corrected V/J names,
V-derived CDR1/2 loops, and upstream-equivalent QC warning.

MixTCRpred code is academic/non-commercial and fetched directly from its
original repository only after explicit acceptance. Optional checkpoints come
from the authors' immutable CC-BY-4.0 Zenodo record and are checksum-verified
before use.

As with any PyTorch checkpoint, explicit overrides and user-managed model files
must come from a trusted source, because loading can execute serialized code.
