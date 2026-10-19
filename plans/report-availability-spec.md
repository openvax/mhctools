# Fresh report availability (#533)

Generate fresh report availability from the prediction results and imported
assessment rows included in that run. Use concise `Included` / `Not included`
statuses for supported predictors; keep the existing predictor/status/reason/
sources output fields. Do not infer execution from model installation or from
an available model's presence in the catalog.

Split released PhageScout ELANE/CTSG sequence profiles from unavailable full
classifiers. Retain inferred-anchor, sparse-coverage and free-peptide transfer
limits, and the remaining CatL/AEP and blood/inflammation coverage gaps.
Archived reports remain immutable. Do not introduce new predictor runs.

Verify real bundled profile results and imported rows produce `Included`, absent
results produce `Not included`, and the full classifiers remain unavailable.
Exercise the fresh-run writer with synthetic data, inspect the coverage page,
and run lint/full required-backend tests and CI before release.
