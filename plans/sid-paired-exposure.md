# Sid peptide paired-route scenario report

Use the 40 publicly disclosed synthetic long peptide records in the frozen
2026-10-06 osteosarc inventory. Prefer the disclosed minimal epitope when it
occurs inside that specific SLP and has a matching class-I MHCflurry row;
otherwise select the lowest-percentile predicted class-I ligand from the
frozen table. Preserve its source rank and the selection basis; do not call
either a validated presented epitope or combine distinct cores.

Compare IM 1 µg of the disclosed free SLP with a counterfactual cytosolic mRNA
product containing a retained N-terminal initiator Met; shift the target-core
coordinates accordingly. Methionine excision is not modeled. Normalize RNA
outputs per input transcript and
per antigen copy produced; the full transcript, formulation, and RNA molecular
weight are unknown, so do not infer an equal-mass RNA comparison. Apply the
existing declared synthetic kinetics and three joint loss scenarios without
turning predictor scores into rate constants. Use 48-hour per-core APC arrival,
node pMHC, blood exposure, and time-varying core retention. Plot retention on a
logarithmic percent axis, excluding undefined zero-denominator RNA points.

Deliver a human-readable Markdown report, a machine-readable CSV, and a
reproducible plotting script. State that curves are conditional calculations,
that mRNA antigen can arise through direct APC expression and cross-dressing,
and that neither T-cell effect nor intrinsic circulation half-life is inferred.
