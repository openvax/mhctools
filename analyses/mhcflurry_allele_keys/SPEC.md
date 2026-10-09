# DLA allele-key resolution (#544)

Both MHCflurry adapters currently pass predictor keys through the legacy two-field normalizer. DLA-88*001:01 becomes DLA-88*01:01 and is rejected by exact membership.

Resolve against the loaded predictor supported keys: exact spelling first, then one unambiguous MHCgnomes Allele identity, with aliases disabled and all allele fields/annotations retained. Reject unknown keys, serotypes/groups and multiple equivalent dictionary keys unless the requested spelling is exact. Preserve input order, remove repeated resolved keys, and forward actual predictor keys into inference and calibration checks. Do not change shared normalization or substitute other alleles by sequence.

Keep resolution adapter-local and common to presentation/affinity. Add regression tests using injected predictors for canonical/legacy DLA, HLA, mouse, exact unconventional keys, duplicates/order, ambiguous/unknown/partial identities, and prediction key forwarding. Run focused tests, lint/full suite, real available DLA execution without implying biological calibration, docs checks and GitHub CI. Bump 3.47.16, merge then deploy from clean master and verify PyPI.

Primary nomenclature implementation: https://github.com/pirl-unc/mhcgnomes/blob/main/README.md (Allele fields, annotations, species and serotype distinctions).

Canine digit-width primary reference: https://www.georgehapp.com/Refs/Kennedy1999.pdf (three-digit major type, two-digit subtype); official database: https://www.ebi.ac.uk/ipd/mhc/group/DLA/.
