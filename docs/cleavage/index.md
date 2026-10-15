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
from mhctools import DPP4qPISA, PeptideInput

peptide = PeptideInput("HAEGTFTSD", source_sequence_name="GLP-1 fragment")
result = DPP4qPISA().predict(peptide)
print(result.sites[0].bond)   # 2: HA | EGTFTSD
print(result.sites[0].score)  # 2.1694, native qPISA score
print(result.to_dict())      # input, model, assay, limitations and source bonds
```

For a motif rule, `site.match` is `1` when the documented pattern matches
and `0` when it does not. It is not the chance that the enzyme cuts the bond;
`site.score` remains reserved for a quantitative model's native output.
`result.pattern_mask` gives the same decisions aligned to every peptide bond:
entry `b - 1` describes bond `b`, and `None` means that bond was not assessed.
An internal rule fills the mask; a terminal rule assesses only the currently
exposed bond. Non-motif results have no pattern mask (`None`).

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
termini. Use the shared `PeptideInput` to record `n_term="acetylated"` or
`c_term="amidated"`, or `unknown`. `PeptideInput.fragment()` describes a
conditional product with explicit terminal chemistry and retained source
coordinates. Cleavage results preserve its context and source provenance;
context is descriptive and does not recalibrate a predictor. Attached peptides
are rejected because these models have not been validated for their chemistry.
Older cleavage result JSON with `source_id` is still readable. The initial
qPISA implementation abstains on modified termini because they are outside
its supported domain. Other
modifications, D-residues, cyclization and conjugates are unsupported; do not
strip their chemistry to obtain a score.

Exopeptidases assess exposed termini. An internal `HAE` is not an immediate
DPP4 site. To ask about a hypothesized product, use
`parent.fragment(start, end, n_term="free", c_term="free")` and predict that
fragment. The source offset is retained. This conditional analysis does not
predict formation of the fragment, cleavage order, rates or competition.

In mhctools 4, `PeptideInput` replaces `CleavageInput`. For Python callers,
rename `source_id` to `source_sequence_name`; the `fragment()` call stays the
same. Scalar cleavage JSON now uses schema version 4 and batch result JSON
uses version 3. Readers continue to accept older result files.
