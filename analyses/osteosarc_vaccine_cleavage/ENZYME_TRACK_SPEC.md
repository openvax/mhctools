# Clear enzyme tracks in the osteosarc atlas

## Problem

Tiny diamonds/triangles and match-only motif rows make enzyme evidence hard to read. An absent row or marker conflates an assessed non-match, an unsupported model, and an intact-terminal site outside a continuation segment. The FAP internal display alias also puts the cut between Gly and Pro, although the rule assesses Gly-Pro|non-Pro.

## Change

- Give the extracellular panel its own larger space and repeat the residue axis so sites can be read locally.
- Always show DPP4, MME, FAP internal and FAP N-terminal activity, plus the existing ANPEP/ENPEP rows. Split FAP modes; preserve their exact assessed bonds.
- Use large labeled motif-match blocks, open assessed-non-match markers, explicit unsupported/unassessed explanations, and a large DPP4 native-score callout connected only to intact SLP bond 2. Never project exposed-terminal scores onto continuation-page termini.
- State native DPP4 log2 depletion units, motif-only evidence, model applicability, and conditional enzyme exposure. Correct the FAP display alias without changing its rule or predictions.
- Re-render into a new timestamped run from the complete verified frozen source. Preserve all inference tables, source/weight/runtime provenance and the earlier artifacts. Update links/captions and bump the patch version.

## Verification

Regression checks cover fixed row visibility, DPP4 missing coefficients, MME length abstention, separate FAP modes, no-match versus unavailable, segment boundaries and unchanged native data. Render and inspect the atlas and manuscript pages. Run lint, full tests, docs checks and GitHub CI; merge and deploy through the repository release script.
