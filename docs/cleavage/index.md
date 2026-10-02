# Peptidase activity

This API reports per-bond evidence about peptide hydrolysis by proteasomes and
other peptidases, using published models and curated experimental observations.
It covers enzymes in the cytosol, ER, endosomes, and extracellular settings.
See [antigen processing](../predictors/processing.md) for their relationship to
MHC presentation and [choosing models](choosing.md) for recommendations by
biological question.

Results identify the bond, enzyme, and evidence type. Some models produce a
numerical score; others report motif decisions or observations for an exact
substrate. These retain their own interpretation, as described in
[reading the evidence](evidence.md).

## Which API do I use?

| Task | Guide |
|---|---|
| Assess individual peptide bonds | [Quickstart](#quickstart) |
| Assess epitopes with flanks or complete vaccine constructs | [Batch assessments](batch.md) |
| Generate a sequence-centered PDF | [Vaccine reports](../vaccine-reports.md) |
| Evaluate a model against measurements | [Benchmarks](../benchmarks.md) |
| Get one cleavage score per peptide | [Processing predictors](../predictors/processing.md) |

Batch assessments cover tumor, APC, and extracellular scenarios. See
[models](models.md) for the supported enzymes and endpoints,
[reading the evidence](evidence.md) for result states, and
[validation status](validation.md) for the supporting data.

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
