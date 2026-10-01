# Peptidase cleavage evidence

Per-bond evidence about where peptidases can cut a peptide, from published
models and curated experimental observations. This is different from the
per-peptide cleavage *scores* in [antigen-processing predictors](../predictors/processing.md):
these APIs say which bond, by which enzyme, with what evidence, and they keep
categorical motif decisions separate from numerical scores.

## Which API do I use?

| I want… | Use | Read |
|---|---|---|
| Evidence for one or a few peptides, per bond and per enzyme | `predict_cleavage()` or `mhctools cleavage` | this page |
| Evidence for named epitopes with their flanks, or complete vaccine constructs, across tumor / APC / extracellular scenarios | `predict_cleavage_batch()` or `mhctools cleavage --input` | [batch assessments](batch.md) |
| A sequence-centered PDF for a vaccine construct | `mhctools vaccine-report` | [vaccine reports](../vaccine-reports.md) |
| To evaluate a model against measured data | `mhctools benchmark` | [benchmarks](../benchmarks.md) |
| One cleavage *score* per peptide | `Pepsickle`, `NetChop`, `NetCleave` | [processing predictors](../predictors/processing.md) |

Related pages: [models](models.md) (every model, what it assesses, how strict it
is), [reading the evidence](evidence.md) (what `matched`, `not_matched` and
`unsupported` mean), and [validation status](validation.md).

## Quickstart

```python
from mhctools import CleavageInput, DPP4qPISA

peptide = CleavageInput("HAEGTFTSD", source_id="GLP-1 fragment")
result = DPP4qPISA().predict(peptide)
print(result.sites[0].bond)   # 2: HA | EGTFTSD
print(result.sites[0].score)  # 2.1694, native qPISA score
print(result.to_dict())      # input, model, assay, limitations and source bonds
```

Or run the whole panel from the command line:

```sh
mhctools cleavage --list-models
mhctools cleavage --sequence RPPGFSPFR --model app2-xp --model cpn-basic
```

## Coordinates and chemistry


This API describes individual peptide bonds. Bond `b` splits
`sequence[:b] | sequence[b:]`; `source_bond` adds the zero-based
`source_start` offset. The per-peptide proteasome predictors remain available separately.

Inputs describe canonical, linear L-peptides. Strings assume free N and C
termini. Use `CleavageInput` to explicitly record `n_term="acetylated"` or
`c_term="amidated"`, or `unknown`. These forms are recorded but the initial
qPISA implementation abstains because they are outside its supported domain.
Other modifications, D-residues, cyclization and conjugates are unsupported;
do not strip their chemistry to obtain a score.

Exopeptidases assess exposed termini. An internal `HAE` is not an immediate
DPP4 site. To ask about a hypothesized product, use
`parent.fragment(start, end, n_term="free", c_term="free")` and predict that
fragment. The source offset is retained. This conditional analysis does not
predict formation of the fragment, cleavage order, rates or competition.
