# mRNA antigen presentation bibliography

Add a focused reference page under `docs/evidence/` covering direct APC
expression, antigen transfer and cross-presentation, and transfer of preformed
peptide–MHC complexes after mRNA vaccination. Summarize each primary paper in
one or two sentences and preserve species, formulation, route and endpoint.

Extract useful published measurements into a downloadable, source-located data
file. Keep observed expression, uptake, antigen concentration and immune
responses distinct; identify which model terms each observation could constrain
and which rate constants are not identifiable. Do not install literature values
as default human vaccine kinetics. Link original reusable source datasets where
available and retain the conditions needed to interpret numbers.

Link the page from the existing core-exposure and comparison guides and
MkDocs Research → Literature reviews navigation. Bump the patch version. Validate source claims,
data consistency, documentation links, strict MkDocs build, desktop and narrow
page layout, lint and the full test suite. Ship through green CI, merge and the
normal clean-main deployment; verify the live page and PyPI release.
