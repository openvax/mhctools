# Core-epitope delivery for candidate vaccine windows

Connect exact candidate sequences and all scored core occurrences to explicit
IM dose scenarios. Numerical outputs remain conditional on identifiable rates;
missing measurements must not acquire invented physiological defaults.

1. Extend the trajectory engine with first-entry accounting for core-bearing
   antigen in draining-node APCs. Keep uptake from node fluid, carriage by
   migrating APCs, endogenous production in nodal APCs and transferred pMHC
   distinct. Count each antigen lineage only once, while allowing RNA to produce
   multiple lineages. Report current blood core-bearing fragments, exact free
   cores, precursors and their integrated exposure. Describe the observed
   post-peak half-decline with its sampling bracket; do not call absorption- or
   production-limited IM kinetics an intrinsic blood half-life.
2. Support secreted and cytosolic mRNA antigen routing. Preserve direct APC
   translation and processing as well as secretion/reuptake. Add explicitly
   parameterized producer-cell MHC-I processing and pMHC transfer (cross-dressing).
   Source studies: Liang 2017 (PMID 28958578), Hassett 2024 (PMID 38161733),
   Rood 2026 (PMID 41545353), Jo 2026 (PMID 41986715). Mouse findings do not set
   universal human pathway fractions. Membrane topology, APC subsets and immune
   response are outside this effective compartment model.
3. Add a batch API/CLI/report for candidate windows and multiple exact cores.
   Default requested SLP amount is 1 microgram of the actual free peptide;
   calculate average molecular weight only for supported unmodified chemistry.
   RNA requires its actual transcript molecular weight. Preserve candidate
   identity, coordinates, core score contributions, source rates and missing
   inputs. Run each unique target occurrence/allele separately. Keep per-core
   delivery and score-weighted descriptive summaries separate from a molecule
   fraction; use explicit bounds for any-core arrival rather than double-count
   shared precursor molecules. Do not manufacture joint MHC competition.
4. Connect Vaxrank candidate windows to this batch contract through an opt-in
   assessment/report path, without silently changing established selection.
   Every retained scored core must be represented; record truncation or missing
   rates, and retain assay-specific half-life annotations separately from IM PK.

Validation: analytic competing hazards, irreversible first entry, migration,
RNA births and reuptake, producer processing, dose scaling, exact/free versus
precursor core accounting, multiple/overlapping cores, missing-input abstention,
CLI artifact preservation and a composed candidate-to-report demonstration.
For each repository release: version bump, lint/full tests, applicable report
smoke, green PR CI, merge and clean-main deployment, published-wheel verification.
Keep the shared Python environment and other active repository work untouched.
