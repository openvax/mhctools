# Normalize vaccine window comparisons

Make the one-row-per-core CSV and Markdown window report lead with dimensionless
quantities: nodal APC arrivals as percent of administered free SLP copies (or
arrivals per input mRNA transcript), and time-averaged blood-core and nodal
surface-pMHC copies divided by input copies and observation duration. Define
these as mean copies per input molecule, not APC positivity or immune response.
Preserve raw copies and copy-hours in the machine-readable result and CSV for
audit. Leave zero-hour means undefined and missing-kinetics results unassessed.
Keep half-decline in hours because duration is an intrinsically timed endpoint.

Document the normalization and show an explicit numerical planning example from
the existing sensitivity script. Add focused checks for units, zero horizon,
missing data and mRNA ratios above one. Add a Research literature page with
primary dose-normalized node measurements, keeping route, formulation, analyte,
time, denominator and APC/core limitations beside each value. Bump patch
version, run lint/tests/docs, then PR, merge and release.
