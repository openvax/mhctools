# Core-specific vaccine calibration audit (#549)

Audit published measurements against the observables required to estimate,
for a specified IM free-SLP window or mRNA construct, the first-arrival fraction
of intact core-bearing antigen in draining-node APCs, blood-core clearance, and
allele-specific surface pMHC persistence. Preserve source route, formulation,
species, analyte, site, dose and time. Include direct expression and acquired
antigen/pMHC pathways. Do not convert carrier label, bulk RNA, protease
activity, T-cell response or model-fitted values into a core-specific rate.

Publish a Research audit with the available assay-scoped observations, the
model's structural identifiability argument, explicit physical bounds, and
minimal measurements needed to distinguish uptake, continuing input and
clearance. Add the new mouse pMHC time course to the mRNA bibliography and its
measurement extract. Extend the machine-readable trajectory evidence inventory
with its endpoint and calibration limit. Cross-link the audit from the vaccine
guides and Research overview. Keep sensitivity scenarios labeled as assumptions;
do not add guessed physiological defaults.

Bump to 3.53.3 on a feature branch, run the repository lint and full test suite,
links and strict docs build, then PR, CI, merge and deploy from clean master.
