# Additional Sid stability-route predictions

## Objective

Run available, licensed enzyme-specific models on all 40 disclosed Sid vaccine peptide records, covering extracellular exposure and processing after uptake. Retain biological mechanisms, exact inputs, native scores, terminal topology, core coordinates, and unavailable-model reasons.

## Work

1. Verify the frozen source manifest and inventory. Audit additional pretrained enzyme-specific candidates for code/weight licensing, local inference, feature completeness, and applicability. Do not train replacement models or fabricate missing features.
2. Run all 18 installed CleaveNet MMP heads on every overlapping 10-residue window of every intact peptide. Preserve window coordinates and five-model spread; a window score does not identify its cleavage bond or imply a released fragment.
3. Run human Pepsickle constitutive and immunoproteasome digestion models, retaining their experimental status and separate native bond scores. Preserve the existing epitope-trained and NetChop tracks as different endpoints.
4. Expose the existing serum/plasma and extracellular peptidase assessments prominently with numerical versus motif evidence separated. List unresolved cathepsin/AEP or other requested enzyme coverage explicitly. Add further verified local pretrained models if the audit supports them.
   The audit supports bundled ITCell B/S internal and H initial-trimming profiles at 15/60/240 minutes, with explicit pH-6.5 assay scope and real author-code conformance. No structural features or replacement training are required.
5. Write new immutable native data/provenance and a readable PDF supplement organized by biological route. Use triangles/lines for site markers. Do not aggregate unlike outputs or convert substrate Z-scores/motif matches to percent peptide loss.
6. Verify real-model conformance, input/window/bond coordinates, source/core scope, counts, hashes and rendered PDF pages. Run lint, full tests and docs checks; ship a version-bumped PR only after CI is green, then deploy from clean master.

The primary question is where enzymatic susceptibility may occur on the route to antigen presentation. Enzyme exposure, concentration, repeated fragment trimming and clinical immune response are not inferred from sequence alone.

Patient-specific regenerated native outputs and PDF remain in ignored `local_results/`. The reusable model/data, verification tests, documentation and report generators ship via PR; the new patient-specific predictions are not published.
