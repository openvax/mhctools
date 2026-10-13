# mRNA vaccine antigen presentation evidence and measurements

Intramuscular mRNA vaccination can put antigen on an APC through several routes:
the APC can translate vaccine RNA itself; it can take up antigen made by another
cell and process it; or it can acquire an already formed peptide–MHC-I complex
(cross-dressing). These routes have different kinetic inputs. Reporter-positive
cells, blood RNA, serum protein and T-cell priming are useful observations, but
none by itself measures the fraction of an encoded **epitope** reaching a
draining-node APC or its surface pMHC lifetime. The measurements below retain
their original species, formulation, dose and assay context.

## Direct expression and processing in APCs

- [Liang et al., 2017](https://doi.org/10.1016/j.ymthe.2017.08.006) gave
  modified mRNA–LNP encoding influenza H10 or a fluorescent reporter by IM or
  intradermal injection to rhesus macaques, then measured labeled LNP uptake and
  reporter translation separately in injection-site and draining-node immune
  cells. Monocytes and dendritic cells translated the reporter, while
  neutrophils took up LNP with little translation, so uptake alone is a poor
  proxy for productive expression ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5768558/)).
- [Hassett et al., 2024](https://doi.org/10.1016/j.omtn.2023.102083) followed
  reporter mRNA and protein after IM mRNA–LNP injection in rodents and
  cynomolgus macaques. RNA and expressed protein appeared at the injection site
  and in draining nodes, with macaque reporter expression especially frequent
  in node macrophages; these assays do not measure allele-specific peptide
  display ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10755037/)).
- [Takanashi et al., 2023](https://doi.org/10.1021/acs.molpharmaceut.2c01024)
  compared free and LNP-formulated reporter mRNA after IM injection in mice and
  measured expression in muscle and lymphoid organs. LNP formulation increased
  reporter signal in lymphoid tissues and improved responses to separately
  encoded vaccine antigens; stable reporter activity does not directly measure
  ongoing RNA translation ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10411422/)).
- [Rood et al., 2026](https://doi.org/10.1038/s41467-025-68080-x) compared
  endogenous and acquired influenza antigen presentation by murine APCs and
  selectively suppressed expression in hematopoietic cells or muscle cells after
  IM mRNA–LNP vaccination ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12873236/)).
  Direct APC expression mattered for the tested MHC-II and CD4 responses, while
  reducing muscle expression had a smaller effect; their [source data](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-68080-x/MediaObjects/41467_2025_68080_MOESM4_ESM.xlsx)
  contain assay-level comparisons, not a universal route fraction.
- [Muro et al., 2026](https://doi.org/10.1126/sciadv.aec7827) gave mice IM
  OVA mRNA–LNP together with reporter mRNA–LNP and stained draining-node
  cDC2 for surface SIINFEKL/H-2Kᵇ. The relative pMHC signal peaked at day 1,
  remained above control through day 5 and was not detected at day 7, while
  reporter expression followed a different time course; continuing production
  means this is not an intrinsic pMHC half-life
  ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13378541/)).

## Antigen made elsewhere and transferred to APCs

- [Lazzaro et al., 2015](https://doi.org/10.1111/imm.12505) tested a
  self-amplifying RNA vaccine in mice with in vivo and in vitro experiments that
  separated muscle expression from professional APC priming. Bone-marrow-derived
  APCs drove CD8 priming without direct transfection, supporting transfer and
  cross-presentation of muscle-derived antigen in this self-amplifying-RNA
  system.
- [Romani et al., 2017](https://doi.org/10.1038/s41598-017-11399-3) gave mice
  intravenous RNA lipoplexes encoding signal-tagged antigen and detected encoded
  p24 in serum after 24 hours. The experiment demonstrates secretion from an
  RNA-expressing cell, but its intravenous carrier and serum assay do not
  determine IM mRNA–LNP secretion yield or subsequent APC uptake
  ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5589809/)).
- [Jo et al., 2026](https://doi.org/10.1038/s41586-026-10353-6) dissected
  dendritic-cell pathways after mRNA–LNP vaccination in mice. Both cDC1 and cDC2
  could prime CD8 cells, and transfer of preformed peptide–MHC-I complexes from
  nonhematopoietic cells made a substantial, type-I-interferon-dependent
  contribution without the recipient APC processing the antigen again
  ([author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC13089314/)).

## Human blood kinetics and an existing quantitative model

- [Kent et al., 2024](https://doi.org/10.1021/acsnano.4c11652) serially
  measured vaccine RNA and ionizable lipid in plasma from 19 people after a
  bivalent mRNA-1273 booster. The similar **post-peak blood declines** of
  assay-linked RNA and lipid constrain circulating-carrier behavior for this
  formulation; the assay does not establish translation-competent RNA in APCs,
  encoded antigen concentration or intrinsic peptide elimination half-life
  ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11447892/)).
- [Dasti et al., 2025](https://doi.org/10.1002/psp4.70041) built a multiscale
  mRNA-vaccine model with APC uptake, intracellular translation and modeled MHC-II
  display, using macaque APC observations and human antibody data for
  calibration and evaluation. Its [published APC calibration data and code](https://github.com/cosbi-research/QSPmRNAVaccines/tree/main/Figure_reproduction/Figure%203a)
  provide a comparison for fitting strategies; antibody agreement does not
  independently identify every antigen-processing rate.

## Measurements worth carrying into a calibration dataset

The [machine-readable measurement extract](mrna-antigen-presentation-measurements.json)
records the values, units, assay, site, time and source figure for the
quantitative entries below. These are **observations**, not mhctools rate
defaults. The time courses and raw assay data, where available, are preferable
to fitting one rate to a reported peak or half-life.

| Experiment | Reported measurement | Potential use and boundary |
| --- | --- | --- |
| Hassett 2024, cynomolgus macaques, 300 µg reporter RNA per IM site, branched-DNA assay | First-observed RNA peaks: node 155 ng/g, injection site 98 ng/g, spleen 854 ng/g at 8 h; plasma 80 ng/mL at 2 h | Constrains combined RNA delivery/distribution at sampled sites. The first sample may miss the actual peak; tissue mass and plasma volume have different denominators. |
| Hassett 2024, same RNA assay | Exposure to last quantifiable sample: node 2,169 ng·h/g, injection site 1,352 ng·h/g, spleen 16,146 ng·h/g, plasma 905 ng·h/mL | Benchmark whole RNA trajectories after matching reporter, carrier and dose. It cannot isolate nodal APC entry or productive cytosolic RNA. |
| Hassett 2024, separate macaque EGFP cohort at 24 h | EGFP-positive node cells: macrophages 36.2%, monocytes 19.7%, conventional DCs 9.1%, plasmacytoid DCs 4.0% | Cell-specific reporter prevalence, not the fraction of the injected dose or the fraction of cells displaying a chosen epitope. |
| Kent 2024, 19 people after a bivalent mRNA-1273 booster | At 4 h, plasma RNA 6.5–112 copies/µL; peak up to 731 copies/µL around day 1–2 | Human blood RNA signal for a specific formulation and booster setting; short-amplicon PCR includes fragmented RNA. |
| Kent 2024, post-peak plasma | Assay-linked RNA slope 0.608/day and reported half-decline 1.14 days; ionizable lipid slope 0.607/day and reported half-decline 1.14 days | Check a modeled blood RNA/carrier tail. Continuing input and assay linkage prevent treating this as intracellular active-RNA decay or epitope half-life. |
| Muro 2026, mice, IM OVA and GFP mRNA–LNP, 1 µg RNA each | Relative SIINFEKL/H-2Kᵇ staining on GFP-positive node cDC2 peaked at day 1, remained above control through day 5 and was not detected by day 7 | Checks a mouse surface-pMHC *time pattern* for that construct and assay. The MFI ratio is not absolute pMHC copies or a measured complex-loss rate. |

Liang's flow-cytometry time points at 4 and 24 hours and day 9 distinguish
LNP-positive from translating cells, but do not establish absolute antigen
copies or surface pMHC per cell. Rood's measured presentation and perturbations
inform the direction of a **relative route contribution** for the tested mouse
constructs and assays. Jo's T-cell priming result identifies cross-dressing as
a route to include in relevant scenarios, not its transfer rate. The
[calibration audit](../research/vaccine-calibration.md) gives the
non-identifiability argument and the matched assays still needed.

For mhctools, match each observation to its modeled analyte before fitting:
extracellular carrier/RNA, active intracellular RNA, completed antigen,
core-containing fragments and surface peptide–MHC are distinct states. A
numerical APC-delivery fraction needs measurements of core-bearing fragments in
the relevant formulation and APC compartment; a pMHC lifetime additionally
needs allele-specific surface-complex measurements. See the [trajectory contract](../vaccine-trajectory.md#kinetic-contract),
[core-delivery guide](../vaccine-core-exposure.md) and
[calibration issue](https://github.com/openvax/mhctools/issues/549).
