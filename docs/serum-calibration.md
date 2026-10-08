# Human serum peptidase calibration audit

This audit addresses [#543](https://github.com/openvax/mhctools/issues/543).
It packages **12 study records and 27 reported hormone-stability measurements**,
with sources, assay conditions, perturbations, reuse terms and identifiability.
The research result is a bounded evidence inventory and minimum calibration
experiment. The audited inputs do not justify a new fit of transferable enzyme
hazards for human long vaccine peptides and their fragments.

## What was found

| Evidence | Available input | What it identifies |
|---|---|---|
| [Abid 2009](https://doi.org/10.1074/jbc.M109.035253), human serum oxidized NPY | Km/Vmax, inhibitors, products | Specific apparent DPP4/aminopeptidase-P/kallikrein reference rates; chemistry, pH 8 and manganese remain part of the assay |
| [Kuoppala 2000](https://doi.org/10.1152/ajpheart.2000.278.4.H1069), human plasma bradykinin | Concentration-dependent product attribution | ACE and CPN-like pathway dominance bounds, without transferable vaccine hazards |
| [Yi 2015](https://doi.org/10.1371/journal.pone.0134427), matched human collection tubes | Table 2: 27 half-life records; MS/antibody time-course figures; six TIF supplements | Assay-specific peptide-signal stability and cocktail protection; individual enzyme contributions remain unresolved |
| [Zhen 2016](https://doi.org/10.1042/BJ20151085), human plasma FGF21 | Selective FAP inhibitor, depletion and rescue; LC-MRM products and plotted trajectories | Causal support for FAP's C-terminal cleavage of FGF21; a 181-aa protein and plotted endpoints do not calibrate vaccine peptide hazards |
| [Bainbridge 2017](https://doi.org/10.1038/s41598-017-12900-8), FAP probes | Human-serum depletion, animal knockout and purified-enzyme comparisons | FAP reporter specificity; ordinary GP probes are also cleaved by PREP. Selective D-Ala probes have dyes and noncanonical chemistry |
| [Toräng 2016](https://doi.org/10.1152/ajpregu.00394.2015), human PYY | In-vitro products and in-vivo sitagliptin arm | Fragment and assay-endpoint differences; sampling inhibitors are quench reagents in vitro, and contaminated infusate prevents a reliable in-vivo conversion-rate estimate |
| [Nyimanu 2019](https://doi.org/10.1038/s41598-019-56157-9), human apelin | In-vivo identified metabolites | Products compatible with proposed MME activity; no matched MME perturbation to isolate its hazard |
| [Dufresne 2017](https://doi.org/10.1186/s12014-017-9176-7), broad human plasma peptidomics | Supplemental peptide/protein matrix, precursor intensities and identification frequency | Preanalytical proteolysis; ice/inhibitor versus warm samples confound temperature and inhibition, and products mix generation with loss |
| [Yi 2007](https://doi.org/10.1021/pr060550h), human serum/plasma | Broad inhibitor collection devices and endogenous fragments | Sequential processing and anticoagulant effects; numerical trajectories/licensing were not established in this audit |
| [FTMS 2007](https://doi.org/10.1016/j.ijms.2006.09.020) | Fractionated human peptide pool plus exogenous DPP4/APP2 | Substrate discovery and a purified DPP4 kinetic validation; intact-matrix kinetics remain unmeasured |
| [Maffioli 2020](https://doi.org/10.3390/molecules25184071) | Pig-plasma 228-peptide library and numerical products | Broad boundary evidence, with no selective named-enzyme calibration or validated human transfer |
| [Böttger 2017](https://doi.org/10.1371/journal.pone.0178943) | Mouse blood/plasma/serum comparison | Matrix/coagulation can alter mechanism and stability ranking; no human vaccine enzyme calibration |

The machine-readable inventory records inspected downloads and hashes where
available. API failures or unverified supplements stay explicit. A public
supplement's existence is distinguished from successful local inspection.
Searches covered human matrix degradation/inhibitors, DPP4/ACE/FAP/MME, matched
peptidomics and MSP-MS libraries; this audit is not proof that no other dataset
exists. Original figures and article prose are not redistributed.

## Reproduce reported reference clocks

```python
from mhctools import serum_calibration_evidence, serum_assay_parent_reference

records = serum_calibration_evidence()["measurements"]
reference = serum_assay_parent_reference("yi2015:t2:05", [0, 4, 8])
```

This record is GLP-1 G36A signal in human serum at room temperature, reported
as 4 ± 0.2 hours. Under the source's first-order assumption its reference signal
fractions are 1, 0.5 and 0.25. This is a reproduction of a reported parent clock,
without fitting raw trajectories or assigning enzyme hazards. Amidation,
matrix, assay and temperature remain explicit. Antibody-only and combined
MS/antibody rows retain their observable/concentration ambiguity.

`>96 h`, `<24 h` and donor ranges remain bounds with no midpoint. A replicate
SD is separate from those bounds and does not become a confidence interval.
The 96-hour horizon is the study's maximum; individual row horizons were not
reported in the table and stay unknown.
No bond allocation, target-survival trajectory or vaccine sequence lookup is
supplied by this reference. Original NPY enzyme-specific reference calculations
remain in [serum contributions](serum-contributions.md).

Source inconsistencies are recorded: Table 2, methods and body report different
room-temperature ranges; the glucagon EDTA range differs between table and
body; OXM-K33 in methods differs from generic OXM table labels. Table 2 values
are retained with these notes. Repeated assay/donor records are not independent
training examples.

```sh
python scripts/serum_calibration_review.py --output /tmp/serum-calibration \
  --source-xml /path/to/PMC4519045.xml
```

This audits every Table 2 row, including rowspans, censoring and uncertainty,
and exports source-linked JSON/CSV with checksums.

## Minimum vaccine-specific calibration experiment

Use the actual peptide chemistry, selected epitope and observed target-bearing
fragments. Compare vehicle with each selective inhibitor at the **same
temperature, substrate concentration, matrix fraction and donor preparation**.
Record serum/plasma processing, anticoagulant, pH, enzyme activation, inhibitor
identity/dose/selectivity and sample IDs. Include orthogonal depletion/rescue
and combination arms when needed to distinguish shared products and serial
pathways.

Measure intact parent and identified target-bearing products over sufficiently
dense early and late time points, using calibrated LC-MS, recovery controls,
stable-isotope standards and LOD/LOQ. Multiple donors and concentrations are
needed to assess donor variation and concentration dependence. Publish the
numerical trajectories and source product identities.

Fit only identifiable competing/serial substrate-scoped rates. Hold out donors
or sequences; separate parameter uncertainty from missing mechanisms. Report
changes in exact-target availability on enzyme removal without renormalizing
remaining hazards. MHC presentation and APC uptake need separate measurements.

For the shortest interpretation, use the [three-route summary](vaccine-target-summary.md).
