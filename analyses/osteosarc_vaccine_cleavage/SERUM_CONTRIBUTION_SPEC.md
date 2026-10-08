# Empirical serum contribution model

The requested endpoint is availability of an exact target sequence in either
the intact vaccine peptide or any successive fragment. Parent disappearance,
the enzyme making a cut, and the change after inhibiting that enzyme are three
different quantities.

## Implementation

1. Curate source-linked human serum/plasma substrate-specific kinetic evidence
   and broad plasma substrate-library observations. Preserve species, matrix,
   dilution, temperature, pH, chemistry, concentrations, and source locations.
   Library product counts are not molecular cleavage fractions or enzyme rates.
2. Add enzyme-labelled, absolute per-bond competing hazards to the existing
   target lineage simulator. Reassess each new fragment. Missing fragment rates
   remain unknown; explicitly assessed zero rates are distinct from missing data.
   Keep renal removal and uptake as separate optional competing hazards.
3. Calculate parent/target survival, destructive-cut attribution, and paired
   enzyme-removal sensitivity. Removing an enzyme subtracts its rates without
   reallocating them to the other enzymes. Sensitivities need not sum to 100%:
   successive enzymes can cooperate, or compete for the same precursor.
4. Include concentration-dependent reference calculations from measured human
   serum kinetics, plus a reproducible public-data example and an evidence
   coverage table. Do not fit universal human enzyme weights from animal
   product counts, purified-enzyme enrichment, or unrelated hormones.
5. Verify units, source values, mass conservation, sequential trimming,
   counterfactual signs, missing-state handling and chemistry. Run lint, full
   tests, docs checks, self-review, CI, merge and PyPI release.

## Evidence and boundaries

- Abid 2009 Table 1 reports apparent Km and serum-normalized Vmax for DPP4,
  aminopeptidase P and plasma kallikrein on mono-oxygenated NPY forms. Assay pH
  is 8 with added manganese, and the substrates have specific chemistry.
  https://doi.org/10.1074/jbc.M109.035253
- Kuoppala 2000 reports substrate-concentration-dependent ACE/CPN-like product
  shares in human plasma. These are product-pathway observations, not universal
  enzyme weights or an intact-vaccine half-life.
  https://doi.org/10.1152/ajpheart.2000.278.4.H1069
- Maffioli 2020 measures 228 synthetic 14-mers in diluted pig plasma. The
  observed product boundary spectrum can inform an explicitly exploratory
  terminal-trimming scenario, but cannot identify individual enzymes or supply
  calibrated human rates. Sequential cuts confound attribution of boundaries.
  https://doi.org/10.3390/molecules25184071

Measured reference calculations and transferred vaccine scenarios must remain
separate. No numerical FAP/MME/ACE rate is inferred from motif recognition alone.
This model addresses extracellular digestion, not systemic elimination or
productive antigen presentation. Patient-linked inputs and results stay local.
