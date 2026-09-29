# Contextual cleavage assessments

## Objective

Accept named epitope occurrences with native flanks and complete vaccine
sequences (including linkers/junctions), run explicitly selected cleavage
models, and overlay model-native evidence on intended epitopes. A batch can
compare tumor, APC and extracellular scenarios without treating those labels
as calibrated cell-specific predictors.

## Contract

- Preserve input occurrence IDs, source offsets, chemistry, sequence scope,
  epitope intervals, imported predictions and arbitrary source annotations.
  Native windows do not establish exposed molecular termini.
- Use zero-based, half-open epitope/fragment intervals and existing cleavage
  bond coordinates. Reject inconsistent peptide/sequence/interval combinations.
- Scenarios explicitly name models and compartments. Default panels report
  their limitations and missing biological coverage; optional models never
  silently disappear or substitute for requested models.
- Hypothetical fragments require explicit terminal chemistry and an assumption
  explaining their production. Their assessments remain conditional.
- Preserve categorical motif decisions, quantitative scores, substrate-only
  observations, unsupported inputs and runtime failures separately. Overlay
  internal bonds and the N/C boundaries without inventing terminal bonds.
- JSON round-trips retain the original inputs and evidence. Human reports
  consume these same records and do not assign aggregate protection scores.
- Expose upstream Pepsickle epitope and in-vitro model families with explicit
  constitutive/immunoproteasome selection, exact artifact identity and native
  score semantics. Reject settings the upstream model ignores.
- Review cathepsin/AEP evidence with pH, activation and assay scope. Include
  source-linked observations only where the original sequence and bond are
  verified. Missing transferable models remain explicit coverage gaps.

## Verification

Drive Python and CLI through the same mixed-context batch: duplicate peptides
with different flanks; multiple source occurrences; complete constructs with
junctions; conditional trimming; changed terminal chemistry; numeric and
categorical evidence; absent assets; save/reload and report rendering. Compare
proteasome variants against actual upstream inference. Separate published
model validation, adapter reproduction and any independently held-out data;
no zero labels inferred from unreported experimental bonds.

## Repository boundaries

mhctools owns predictor adapters and the cleavage evidence/report contract.
Generalized Exacto/LENS/pVACseq ingestion, evidence reconciliation and ranking
belong to Topiary (#365-370, #288). Vaxrank consumes those public APIs for
construct selection (#497) and exposes existing audit machinery in its CLI
(#445). This PR links those dependencies and does not duplicate their scoring
or antigen-source logic.
