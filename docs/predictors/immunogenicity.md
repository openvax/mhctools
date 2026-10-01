# Immunogenicity predictors

Predictors of whether a peptide elicits a T-cell response. They emit
`immunogenicity`, read with `result.immunogenicity`. `Calis` needs only the
peptide; the others also need an allele.

| Predictor | Class | Needs | Notes |
|---|---|---|---|
| [Calis](#calis) | I | peptides only | Built in; the baseline |
| [PRIME](#prime) | I | peptides + alleles | Calls MixMHCpred |
| [DeepImmuno](#deepimmuno) | I | peptides + alleles | 9- and 10-mers only |
| [TLimmuno2](#tlimmuno2) | II | peptides + class II alleles | Slow percentile rank |
| [BigMHC_IM](binding.md#bigmhc) | I | peptides + alleles | Described under BigMHC |

## Read this before trusting a score

> ⚠️ Every current CD8 immunogenicity predictor — `PRIME`, `BigMHC_IM`, and
> `DeepImmuno` included — ranks well in the characterized regime but generalizes poorly to
> truly novel neoepitopes; independent benchmarks put the field near AUC
> 0.5–0.65 on unseen tumor neoepitopes (ITSNdb ~0.52–0.60, ICERFIRE ~0.56,
> IMPROVE ~0.60). In the one neutral head-to-head that scored both (NeoaPred,
> *Bioinformatics* 2024), **`BigMHC_IM` edged `PRIME` on cancer neoepitopes**,
> while PRIME tends to do better on viral / infectious-disease epitopes — its
> training positives are mostly viral and cancer-testis antigens, with only
> ~129 (v1) / ~596 (v2) true immunogenic neoepitopes. PRIME's higher
> self-reported numbers are partly attributable to documented train/test
> overlap (IMPROVE flagged ~70% overlap with its evaluation set). Use these
> scores to prioritize, not as ground truth.


## Calis

`Calis` is the classic sequence-only IEDB class-I immunogenicity model (Calis
et al. 2013): a fixed per-amino-acid log-enrichment scale weighted by
per-position importance, with the anchor positions (P1/P2/C-terminus) masked
out.

It needs **no external install and no downloaded weights** — the ~30 published
parameters (from the open-access CC-BY paper) are built in — so it is a fast,
dependency-free, allele-independent baseline. It emits one `immunogenicity`
prediction per peptide (empty `allele`); `score > 0` leans immunogenic.

```python
from mhctools import Calis

predictor = Calis()
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score            # 0.30484 (higher = more immunogenic)
```

## PRIME

`PRIME` predicts CD8+ T-cell immunogenicity of class-I peptides by combining
MHC-I binding (via MixMHCpred, which it calls internally) with a
TCR-recognition propensity model. It emits one `immunogenicity` prediction per
(peptide, allele): `score` is the PRIME score (higher = more immunogenic) and
`percentile_rank` is the PRIME %Rank (lower = better).

PRIME is academic / non-commercial licensed, so mhctools shells out to an
install you provide rather than vendoring it.

```python
from mhctools import PRIME

predictor = PRIME(
    alleles=["HLA-A*02:01", "HLA-B*07:02"],
    program_name="PRIME",                    # or an absolute path
    mixmhcpred_path="/path/to/MixMHCpred",   # v3.0+, optional if on PATH
    timeout=300)
results = predictor.predict(["GILGFVFTL", "NLVPMVATV"])
results[0].immunogenicity.score
```

mhctools verifies MixMHCpred's reported version before PRIME inference and
rejects versions older than 3.0 or an unparseable version. The timeout covers
the PRIME process tree, including its nested MixMHCpred call.

## DeepImmuno

`DeepImmuno` predicts class-I CD8+ immunogenicity from the peptide and its
HLA-A/B/C allele with a small CNN (Li et al. 2021). It scores **9- and 10-mers
only** and supports a fixed set of ~62 alleles, snapping anything else to the
nearest it knows. It emits one `immunogenicity` prediction per (peptide,
allele); `score` is in 0–1 (higher = more immunogenic).

DeepImmuno ships its weights in-repo and is MIT-licensed, but its script loads
them with an old Keras 2 / TensorFlow stack, so mhctools shells out to
DeepImmuno's own CLI in a separate checkout. Run `mhctools fetch deepimmuno`,
or point at a manual clone with `DEEPIMMUNO_HOME`, and set `DEEPIMMUNO_PYTHON`
to an interpreter that has TensorFlow — with Keras 2, or newer TensorFlow plus
the `tf-keras` shim, since the wrapper sets `TF_USE_LEGACY_KERAS=1` for the
subprocess.

```python
from mhctools import DeepImmuno

DeepImmuno.fetch()
predictor = DeepImmuno(alleles=["HLA-A*02:01"])   # resolves DEEPIMMUNO_HOME / ~/DeepImmuno
results = predictor.predict(["NLVPMVATV", "GILGFVFTL"])
results[0].immunogenicity.score                   # 0.9568 (higher = more immunogenic)
```

## TLimmuno2

`TLimmuno2` is the odd one out: it predicts **class-II (CD4+)** immunogenicity
— the only class-II immunogenicity model here (`Calis`, `PRIME`, `BigMHC_IM` and
`DeepImmuno` are all class I).

It scores a peptide against a class-II allele (transfer-learned from class-II
binding) and emits one `immunogenicity` prediction per (peptide, allele):
`score` in 0–1 (higher = more immunogenic) and `percentile_rank` from its %Rank
against a background set, rescaled to 0–100 (lower = more immunogenic).

Native NetMHCIIpan-style keys (`DRB1_0803`, `HLA-DPA10103-DPB10101`) pass
through; common DR forms (`HLA-DRB1*08:03`) are converted; anything TLimmuno2
does not know raises. Its upstream license is ambiguous (an Apache-2.0 README
badge, no LICENSE file), which mhctools treats the same way as NetCleave: it
can fetch a pinned snapshot, but only when you confirm your own use is
authorized, so the first fetch requires `--accept-license`. `TLIMMUNO2_PYTHON`
names an interpreter that has TensorFlow (Keras 2, or newer TensorFlow plus
`tf-keras`).

```python
from mhctools import TLimmuno2

predictor = TLimmuno2(alleles=["DRB1_0803"])  # TLIMMUNO2_HOME, ~/TLimmuno2, then snapshot
results = predictor.predict(["FHTMWHVTRGAVLMY"])
results[0].immunogenicity.score                    # 0.9874 (higher = more immunogenic)
```

> ⚠️ TLimmuno2's %Rank is computed against ~90,000 background peptides **per
> distinct allele**, so a call costs about a minute per allele regardless of how
> many peptides you pass — batch peptides by allele. Class-II immunogenicity is
> noisier than class-I; a prioritization aid, not ground truth.
