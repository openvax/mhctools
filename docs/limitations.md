# Known limits of these predictors

mhctools wraps published models. It does not improve them, and several of them
are weaker than their own papers suggest. This page collects every caveat
recorded elsewhere in these docs so you can find them before you rely on a
score, rather than after.

Nothing here is new information — each row links to the full explanation in
context.

## At a glance

| Predictor | What to watch out for | Details |
|---|---|---|
| CD8 immunogenicity (`PRIME`, `BigMHC_IM`, `DeepImmuno`) | Field-wide: ~AUC 0.5–0.65 on unseen tumor neoepitopes | [notes](predictors.md#deepimmuno) |
| `TLimmuno2` | ~1 minute per distinct allele; class-II immunogenicity is noisier than class-I | [notes](predictors.md#tlimmuno2) |
| `PlifePred2` | Endpoint semantics are not established — units, transform, species and matrix are all inferred | [notes](predictors.md#plifepred2) |
| `PeptiVerse` | Fit on 130 examples, cross-validation only, no external test set; unsafe pickle serialization | [notes](predictors.md#peptiverse) |
| `NetCleave_II` | Class-II C-terminal cleavage is a much weaker signal than class I (AUC ~0.66 vs ~0.91) | [notes](predictors.md#netcleave) |
| `DeepTAP` | Self-reported evaluation; no independent TAP benchmark exists for any tool | [notes](predictors.md#deeptap) |
| `ERAMER` | Self-reported evaluation; ERAP1 trimming is intrinsically noisy | [notes](predictors.md#eramer) |
| `CapHLA` | Performance numbers are author-reported | [notes](predictors.md#caphla) |
| `MixTCRpred` | Loading a PyTorch checkpoint can execute serialized code — use trusted sources | [notes](predictors.md#mixtcrpred) |
| `DPP4qPISA` | Substrate-depletion estimates, not serum half-lives or probabilities | [cleavage guide](cleavage.md) |

## Three recurring themes

**Self-reported evaluation.** `DeepTAP`, `ERAMER`, and `CapHLA` are each
evaluated by their own authors, with no neutral benchmark to check them
against. For TAP and ERAP1 trimming this reflects the state of the field, not a
gap these particular tools left. Read those scores as pathway priors that help
prioritize, not as validated oracles.

**Generalization to novel neoepitopes.** CD8 immunogenicity predictors rank
well inside the regime they were trained on and fall toward chance outside it.
The [DeepImmuno notes](predictors.md#deepimmuno) give the independent benchmark
numbers and the one neutral head-to-head comparison.

**Unestablished semantics.** [`PlifePred2`](predictors.md#plifepred2) ships no
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
[cleavage guide](cleavage.md#interpreting-rule-based-evidence).

## Where the rest lives

- [Benchmark methodology and training overlap](benchmarks.md) — how mhctools
  reports provenance, repeated measurements, and missing target-domain evidence
- [Optional backend conformance](optional-backends.md) — what `verified` and
  `inference_reproduced` do and do not claim
- [Peptide PK, uptake, and tissue exposure](exposure-results.md) — why these
  endpoints refuse a generic ordering
