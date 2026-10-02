# Known limits of these predictors

mhctools wraps published models. It does not improve them, and several of them
are weaker than their own papers suggest. This page collects every caveat
recorded elsewhere in these docs so you can find them before you rely on a
score, rather than after.

Each row links to the full explanation in context.

## At a glance

| Predictor | What to watch out for | Details |
|---|---|---|
| CD8 immunogenicity ([PRIME](predictors/immunogenicity.md#prime), [BigMHC](predictors/binding.md#bigmhc) IM, [DeepImmuno](predictors/immunogenicity.md#deepimmuno)) | Field-wide: ~AUC 0.5–0.65 on unseen tumor neoepitopes | [notes](predictors/immunogenicity.md#deepimmuno) |
| DeepImmuno | Scores 9- and 10-mers only, over a fixed set of about 62 alleles; other alleles are snapped to the nearest one it knows | [notes](predictors/immunogenicity.md#deepimmuno) |
| PRIME | Higher self-reported numbers are partly explained by documented train/test overlap; training positives are mostly viral | [notes](predictors/immunogenicity.md#read-this-before-trusting-a-score) |
| [Pepsickle](predictors/processing.md#pepsickle), [NetChop](predictors/processing.md#netchop) | The C-terminal score needs the residues after the peptide; with no `c_flanks` it is 0.0, which is not a prediction of no cleavage | [notes](predictors/processing.md#pepsickle) |
| [TLimmuno2](predictors/immunogenicity.md#tlimmuno2) | ~1 minute per distinct allele; class-II immunogenicity is noisier than class-I | [notes](predictors/immunogenicity.md#tlimmuno2) |
| [PlifePred2](predictors/peptide-pk.md#plifepred2) | Endpoint semantics are not established; units, transform, species and matrix are all inferred | [notes](predictors/peptide-pk.md#plifepred2) |
| [PeptiVerse](predictors/peptide-pk.md#peptiverse) | Fit on 130 examples, cross-validation only, no external test set; unsafe pickle serialization | [notes](predictors/peptide-pk.md#peptiverse) |
| [NetCleave](predictors/processing.md#netcleave) (class II) | Class-II C-terminal cleavage is a much weaker signal than class I (AUC ~0.66 vs ~0.91) | [notes](predictors/processing.md#netcleave) |
| [DeepTAP](predictors/processing.md#deeptap) | Self-reported evaluation; no independent TAP benchmark exists for any tool | [notes](predictors/processing.md#deeptap) |
| [ERAMER](predictors/processing.md#eramer) | Self-reported evaluation; ERAP1 trimming is intrinsically noisy | [notes](predictors/processing.md#eramer) |
| [CapHLA](predictors/binding.md#caphla) | Performance numbers are author-reported | [notes](predictors/binding.md#caphla) |
| [MixTCRpred](predictors/tcr.md#mixtcrpred) | Loading a PyTorch checkpoint can execute serialized code; use trusted sources | [notes](predictors/tcr.md#mixtcrpred) |
| [DPP4qPISA](cleavage/models.md#human-dpp4-qpisa) | Substrate-depletion estimates, not serum half-lives or probabilities | [peptidase activity guide](cleavage/index.md) |

## Three recurring themes

**Self-reported evaluation.** [DeepTAP](predictors/processing.md#deeptap), [ERAMER](predictors/processing.md#eramer), and [CapHLA](predictors/binding.md#caphla) are each
evaluated by their own authors, with no neutral benchmark to check them
against. For TAP and ERAP1 trimming this reflects the state of the field, not a
gap these particular tools left. Read those scores as pathway priors that help
prioritize, not as validated oracles.

**Generalization to novel neoepitopes.** CD8 immunogenicity predictors rank
well inside the regime they were trained on and fall toward chance outside it.
The [DeepImmuno notes](predictors/immunogenicity.md#deepimmuno) give the independent benchmark
numbers and the one neutral head-to-head comparison.

**Unestablished semantics.** [PlifePred2](predictors/peptide-pk.md#plifepred2) ships no
publication, training data, or target definition. mhctools reports its native
output and declines to claim a duration unless you explicitly opt in. This is
the clearest case of a general rule: where a transform to a physical unit is
unresolved, the wrapper leaves `value` empty rather than guessing.

## Conventions that are easy to over-read

A few properties of the output format look like stronger claims than they are.

**`score` is higher-is-better within one endpoint only.** It is a numerical
selection convention. A larger score is not universally better for a vaccine,
and scores from different predictors or endpoints are not interchangeable. See
[score, value, and percentile_rank](kinds.md#score-value-and-percentile_rank).

**`1-log50k` is not a probability.** For affinity predictions it is a monotone
rescaling useful for ordering. It is not calibrated and not inherently bounded
to 0–1.

**Sharing a unit is not sharing a measurement.** `pMHC_stability` (hours) is
the lifetime of a peptide-MHC complex; `peptide_half_life` (hours) is the
lifetime of the parent peptide. Different molecules, different assays.

**`mhctools ls` reports files, not working software.** `ready` means a path was
located. Use `mhctools predictors` for a capability report, where `LOCATED`,
`RUNNABLE`, and `REPRODUCED` are independent observations and `not checked` is
never promoted to success. See [inventory is not
capability](artifacts.md#inventory-is-not-capability).

**A motif non-match is not evidence of resistance.** In the cleavage API,
`not_matched` and `no_cleavage_detected` are not probabilities. See the
[peptidase activity guide](cleavage/evidence.md#interpreting-rule-based-evidence).

## Where the rest lives

- [Benchmark methodology and training overlap](benchmarks.md): how mhctools
  reports provenance, repeated measurements, and missing target-domain evidence
- [Optional backend conformance](optional-backends.md): what `verified` and
  `inference_reproduced` do and do not claim
- [Peptide PK, uptake, and tissue exposure](exposure-results.md): why these
  endpoints refuse a generic ordering
