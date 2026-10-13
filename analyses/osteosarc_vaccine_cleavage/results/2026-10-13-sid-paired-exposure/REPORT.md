# Sid's disclosed SLPs: paired route scenarios

This is a **48-hour calculation under synthetic rates**, not a measured fraction
of Sid's dose, a clinically calibrated prediction, or a vaccine recommendation.
The 40 publicly disclosed SLP sequences are from [osteosarc.com/vaccines](https://osteosarc.com/vaccines/).
For each, the target is the disclosed minimal epitope if it occurs in
that SLP and has an exact frozen MHCflurry row; otherwise it is the
lowest-percentile class-I ligand in that table. These are candidate
ligands, **not observed presented epitopes**. The same amino-acid sequence
is then modeled as a free SLP or as a hypothetical cytosolic mRNA product.

[Log-scale core-retention time courses](core-retention-log.pdf) summarize all 40 peptides.
The SLP panel is the percent of dosed peptide molecules with this exact
core still present in any modeled fragment or pMHC state. The RNA panel
is the percent of antigen copies produced by that time still retaining
the core; at time zero its denominator is zero and the point is omitted.
These are per-core marginals, not an 'any of several cores' union.

## What the two routes mean

- **SLP:** 1 µg of each unmodified free peptide, using its sequence-derived
  molecular weight. Nodal APC arrival is a percent of input peptide copies.
- **mRNA:** one input transcript encoding the same peptide sequence. RNA
  arrival is antigen copies per input transcript, and also a percent of
  antigen copies produced by 48 hours. The real transcript, UTRs, signal
  peptide, LNP, and mass are unknown. Do not read it as an equal-µg comparison.
  This cytosolic route includes direct APC expression, antigen transfer
  and preformed pMHC transfer; secretion would require a different construct.
- The existing 19 disclosed mRNA context sequences are a different
  set of constructs and are outside this matched 40-SLP counterfactual.
- All transport, translation, processing and loading rates are declared
  synthetic assumptions from `examples/vaccine_sensitivity_scenarios.py`.
  The loss/clearance scenarios are jointly 0.1, 1, and 10 per hour.
  They are **not** physiological bounds. MHCflurry scores only select
  the displayed target; they do not set any rate or pMHC loading affinity.

## Per-peptide results

The three-number ranges below run from 0.1 to 10 /h joint-loss scenarios;
their endpoints are sensitivity results, not confidence intervals.
Node APC arrival counts first core-bearing entry. Mean pMHC is surface
copy-hours per input molecule divided by 48 hours; it is not a response
probability. RNA and SLP means have different input-molecule denominators.

| Gene / record | Chosen core / allele (rank %, basis) | SLP APC arrivals: central % (low–high loss) | RNA APC arrivals/transcript: central (low–high loss) | RNA APC arrivals / antigen made: central % | Mean node pMHC/input: SLP ; RNA |
| --- | --- | ---: | ---: | ---: | ---: |
| ADGRF5 / `ADGRF5-chr6-46856758:vaccine-peptide-1` | `EALLSKFSL` / HLA-B*08:01 (0.35, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| BMP1 / `BMP1-chr8-22192067:vaccine-peptide-3` | `KIILNFTTL` / HLA-C*01:02 (0.40, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| BTD / `BTD-chr3-15645183:vaccine-peptide-3` | `GVVDGLHTV` / HLA-C*01:02 (0.05, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| CD109 / `CD109-chr6-73818405:vaccine-peptide-1` | `FSGPGMSGM` / HLA-C*01:02 (0.30, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| CDC40 / `CDC40-chr6-110228838:vaccine-peptide-3` | `TLSPNGKSL` / HLA-C*01:02 (0.36, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| CENPL / `CENPL-chr1-173811223:vaccine-peptide-1` | `HLQKRLESV` / HLA-B*08:01 (0.04, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| CUL9 / `CUL9-chr6-43186384:vaccine-peptide-1` | `RMLDYYEEI` / HLA-C*01:02 (0.12, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| DYNC1H1 / `DYNC1H1-chr14-101980529:vaccine-peptide-3` | `KRFHATISF` / HLA-C*07:01 (0.00, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| DYNC1H1 / `DYNC1H1-chr14-101980529:vaccine-peptide-4` | `KRFHATISF` / HLA-C*07:01 (0.00, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| EPG5 / `EPG5-chr18-45925851:vaccine-peptide-1` | `SEIAVIRDWKK` / HLA-B*27:05 (0.98, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–18.4) | 27 | 2.77e-05 ; 0.00157 |
| EPG5 / `EPG5-chr18-45925851:vaccine-peptide-2` | `YLWQPSTSEI` / HLA-C*01:02 (1.49, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| EPG5 / `EPG5-chr18-45925851:vaccine-peptide-3` | `SVIRTLSTI` / HLA-B*27:05 (0.69, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| EXOC4 / `EXOC4-chr7-133274996:vaccine-peptide-1` | `SVIRTLSTI` / HLA-B*27:05 (0.69, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| EXOC4 / `EXOC4-chr7-133274996:vaccine-peptide-2` | `SVIRTLSTI` / HLA-B*27:05 (0.42, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| GLIS3 / `GLIS3-chr9-3856149:vaccine-peptide-1` | `ARNKLFFSTKW` / HLA-B*27:05 (0.19, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| H1-2 / `H1_2-chr6-26055824:vaccine-peptide-1` | `VKPKVVKPK` / HLA-B*27:05 (0.58, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| KDM3B / `KDM3B-chr5-138391529:vaccine-peptide-1` | `NLFLQPPKL` / HLA-C*01:02 (0.17, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| KDM3B / `KDM3B-chr5-138391529:vaccine-peptide-2` | `SSNPFLAFV` / HLA-C*01:02 (0.22, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–18.4) | 27 | 2.77e-05 ; 0.00157 |
| KIF1C / `KIF1C-chr17-5007293:vaccine-peptide-1` | `KTWEEKLRK` / HLA-B*27:05 (1.30, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| MAP2 / `MAP2-chr2-209694768:vaccine-peptide-3` | `RVVPFTKAL` / HLA-C*01:02 (0.04, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| MAP2 / `MAP2-chr2-209694768:vaccine-peptide-4` | `RVVPFTKAL` / HLA-C*01:02 (0.04, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| MT-ND5 / `MT_ND5-chrM-12994:vaccine-peptide-1` | `ATAGKSAQL` / HLA-C*01:02 (0.34, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| MYO15B / `MYO15B-chr17-75589953:vaccine-peptide-1` | `RAQAPTRVL` / HLA-C*01:02 (0.01, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| MYO9A / `MYO9A-chr15-71994488:vaccine-peptide-1` | `KKFDWIVFR` / HLA-B*27:05 (0.13, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| NME1 / `NME1-chr17-51154443:vaccine-peptide-1` | `ISSCDTRTM` / HLA-C*01:02 (1.07, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| NR2F2 / `NR2F2-chr15-96332299:vaccine-peptide-1` | `AGGQGGPGVPG` / HLA-B*27:05 (1.18, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| NR2F2 / `NR2F2-chr15-96332299:vaccine-peptide-2` | `AGGQGGPGVPG` / HLA-B*27:05 (1.18, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| PIP5K1A / `PIP5K1A-chr1-151242178:vaccine-peptide-3` | `AAPVATPAL` / HLA-C*01:02 (0.06, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| PTH1R / `PTH1R-chr3-46898732:vaccine-peptide-1` | `SFMLRAVSF` / HLA-B*08:01 (0.06, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| SLC25A12 / `SLC25A12-chr2-171813470:vaccine-peptide-3` | `ATAVYPTDL` / HLA-C*01:02 (0.93, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| SLC25A12 / `SLC25A12-chr2-171813470:vaccine-peptide-4` | `ATAVYPTDL` / HLA-C*01:02 (0.90, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| SMC5 / `SMC5-chr9-70298024:vaccine-peptide-3` | `RQKRIGNTR` / HLA-B*27:05 (0.20, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| TECPR1 / `TECPR1-chr7-98241127:vaccine-peptide-3` | `TIPAAPVGL` / HLA-C*01:02 (0.09, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| TECPR1 / `TECPR1-chr7-98241127:vaccine-peptide-4` | `SRWESAWSGL` / HLA-C*07:01 (0.25, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| TECPR1 / `TECPR1-chr7-98241127:vaccine-peptide-5` | `GRDRPWSGK` / HLA-B*27:05 (0.80, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| VPS13B / `VPS13B-chr8-99861923:vaccine-peptide-1` | `RQLPESLGEGL` / HLA-C*07:01 (0.05, predicted) | 4.81 (0.0598–37.7) | 1.52 (0.0237–20.2) | 27.1 | 2.92e-05 ; 0.00213 |
| VPS13B / `VPS13B-chr8-99861923:vaccine-peptide-2` | `RQGLFRLGI` / HLA-B*27:05 (0.17, predicted) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| VPS72 / `VPS72-chr1-151185547:vaccine-peptide-3` | `ERALLPLEL` / HLA-C*07:01 (0.25, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| WWC1 / `WWC1-chr5-168399513:vaccine-peptide-3` | `GHLQHELQF` / HLA-C*07:01 (0.08, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |
| ZNF436 / `ZNF436-chr1-23362762:vaccine-peptide-3` | `KSFGRSCHL` / HLA-C*01:02 (0.62, disclosed) | 4.81 (0.0598–37.7) | 1.51 (0.0237–19.2) | 27 | 2.52e-05 ; 0.00126 |

Across the 40 selected targets, the central-scenario median is
**4.81%** SLP dose arriving in nodal APCs, and **1.51** antigen arrivals
per input transcript. The shared assumed rates, rather than measured
sequence-specific kinetics, drive these values. SLP APC arrival is
identical to displayed precision across all rows by construction; small
RNA differences arise from the assumed route and target geometry.
This is **not a validated ranking of Sid's peptides**.

The CSV contains each scenario's core-retention percentages, blood
exposure and sampled post-peak half-decline bracket. This blood decline
is not an intrinsic circulation half-life while input continues.
The model cannot derive a T-cell response or vaccine efficacy from pMHC
alone. For empirical gaps and pathway evidence, see the
[calibration audit](../../../../docs/research/vaccine-calibration.md) and
[mRNA presentation research](../../../../docs/evidence/mrna-antigen-presentation.md).

## Provenance

Generated at 2026-10-13 15:51 UTC with mhctools 3.53.6.
- `analyses/osteosarc_vaccine_cleavage/results/2026-10-06T142010-273003-0400/tables/vaccine_sequence_inventory.csv`: SHA-256 `bd806a0d5fbb1fd7c3a6623e5cdf41ddaa9695c2fcf0a443e3aa2d518b869c85`
- `analyses/osteosarc_vaccine_cleavage/results/2026-10-06T142010-273003-0400/tables/slp_mhc_ligand_predictions.csv`: SHA-256 `bbe129bc8d5d6978e48db1138cf4ab16a87c3111c9d1df4e00fb053cd5b7bb6d`
- Script: `analyses/osteosarc_vaccine_cleavage/paired_exposure.py`.

From the repository root, install `.[vaccine-report]` and reproduce
into a **new** directory with:

```sh
PYTHONPATH=examples python analyses/osteosarc_vaccine_cleavage/paired_exposure.py \
  --output-dir analyses/osteosarc_vaccine_cleavage/results/NEW_RUN
```
