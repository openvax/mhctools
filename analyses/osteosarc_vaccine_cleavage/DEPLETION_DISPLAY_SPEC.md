# Simple loss labels and cleavage markers

## Change

- Replace DPP4 log2 callouts with whole-percent predicted loss in the source four-hour assay, calculated as `100 * (1 - 2**(-score))`. Keep native scores unchanged in CSV.
- Show positive estimates below one percent as `<1% predicted loss`, zero as `0% predicted loss`, and negative native estimates as `No predicted loss`; do not invent a zero-valued native score or fill unavailable results.
- Apply the same display wording to intact-SLP values on continuation pages without moving the terminal cut site.
- Replace the M blocks with large filled downward triangles in the enzyme row's colour. Keep assessed non-matches and unavailable states distinct. FAP/MME motifs do not acquire numerical loss scores.
- Keep one short legend for symbols and assay context; put the conversion and limitations in the report documentation.
- Apply the same loss labels and vector triangle/circle markers to the separate mechanism scorecards. Use one shared display helper and record its checksum in both reports' provenance.
- Re-render the complete atlas, manuscript set and mechanism scorecards from checksum-verified frozen predictions into new immutable runs, update links, and bump the patch release.

## Verification

Check percent conversion, tiny positive/negative/missing estimates and continuation scope. Verify source/native tables are unchanged; render and inspect all PDF pages. Run lint, focused/full tests, docs checks and GitHub CI, then merge and deploy from clean master.
