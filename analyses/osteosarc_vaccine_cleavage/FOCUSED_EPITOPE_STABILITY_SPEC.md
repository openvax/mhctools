# Focused epitope stability

## Problem and outcome

The cascade report gives arbitrary cut-weight scenarios and illustrative cuts
the visual prominence of predictor outputs. The user needs the fate of the
best qualifying mutant epitope, with parent disappearance as context.

## Implementation

1. Add reusable target-relative cleavage annotations that retain the native
   enzyme/model evidence and distinguish internal target loss, boundary release
   and flank trimming. Unavailable inputs remain explicit; a non-match never
   becomes protection. Keep extracellular mechanisms separate from processing.
2. Add a target-retention summary for simulated paths, with explicit loss causes,
   unknown outcomes and right censoring. A simulation median is conditional on
   supplied rates and cut weights, never a calibrated circulation half-life.
3. Build a reproducible focused report from frozen local inference artifacts.
   Lead with the best mapped-mutant candidate per MHC class and construct; ranks
   from different classes cannot establish one global winner. Show parent and
   released-target estimates side by side per estimator. Highlight the selected
   target and exact candidate bonds with large dashed red lines, label mechanisms
   and evidence type, and keep model availability visible. Put arbitrary 3x/10x
   trajectories and assumed unflagged cuts in optional supporting detail.
4. Preserve the original report and source manifests. Run focused regressions,
   lint, full tests, docs checks, self-review, CI, PR merge and PyPI release.
   Regenerate and render the new local PDF, inspect every page, and verify data
   against the frozen outputs. No patient-linked data goes into the PR.

## Scientific basis and limits

- Intact parent disappearance and surviving products are different endpoints:
  https://doi.org/10.1021/acs.jmedchem.1c00795
- Serum, anticoagulated plasma and blood can give different stability results:
  https://doi.org/10.1371/journal.pone.0178943
- PeptiVerse has 130 sequence half-life examples; long-vaccine accuracy is not
  established: https://doi.org/10.1038/s41467-026-74167-w
- Cavaco uses heterogeneous training assays and composition descriptors:
  https://doi.org/10.1111/cts.12985
- Chemistry, constrained structure and binding can prolong peptide lifetime;
  neither sequence-only estimator models these vaccine properties explicitly:
  https://doi.org/10.1021/ja405108p and the semaglutide FDA label.

No change to predictor scores or a new score-to-rate conversion is justified.
No absence of an enzyme flag is evidence that the target is protected. The
standalone target estimate assumes a free peptide present at time zero, not a
product's release time or lifetime inside its parent. Individual selected-target
loss does not establish loss of all possible alternative epitopes or presentation.
