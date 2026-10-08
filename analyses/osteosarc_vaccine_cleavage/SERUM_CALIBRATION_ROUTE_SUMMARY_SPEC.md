# Serum calibration and compact vaccine-route summary (#543)

## Goal

Audit reusable human serum/plasma perturbation measurements and expose the
smallest useful summary of exact-target preservation for three delivery routes.
Preserve measured kinetics, substrate discovery, recognition and assumptions as
distinct evidence. Intracellular antigen processing seeks productive epitope
generation and MHC loading; intact-parent persistence is not its objective.

## Deliverables

1. Source-linked, downloadable evidence inventory and reported measurement table
   with matrix, temperature, chemistry, concentration, perturbation, licensing,
   uncertainty and machine-readable reuse/identifiability decisions.
2. Assay-specific reference calculations only where a published measurement
   supports them. No arbitrary vaccine enzyme hazards, image-digitized fit, or
   transfer of hormone/protein kinetics to vaccine fragments. If no suitable
   human long-peptide inhibitor panel is found, record the search and minimum
   experiment needed to identify these rates.
3. Compact per-target route summary in the structured vaccine report: mRNA
   expression remains unassessed without RNA/construct/delivery evidence; APC
   processing distinguishes target-internal from boundary cleavage; serum
   reports calibrated target survival only when appropriate kinetics exist.
   Missing or below-threshold site evidence never becomes a protection verdict.
4. A one-page, readable route guide and reproducible JSON/CSV exports. The
   summary retains native model evidence and avoids a combined probability.

## Verification and release

Use primary sources before scientific edits. Test concentration/endpoint scope,
censoring, unresolved inhibitor attribution, boundary versus target-internal
cuts, and conditional cross-presentation. Inspect PDF renders. Run lint, full
tests, strict docs and links; self-review; bump patch version, PR, green CI,
merge, gated deployment from clean master, and verify the live PyPI artifacts.
