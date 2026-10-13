# Dose-normalized vaccine delivery to draining lymph nodes

What fraction of a vaccine dose reaches a draining node? The answer depends on
**what was labeled, where it was injected, which node was collected, when it
was sampled, and what the denominator means**. This reference list collects
primary experiments near the [vaccine-window](../vaccine-windows.md) question.
It is a scale check for explicitly assumed scenarios, not a pooled estimate
for an IM SLP core or a default transport-rate distribution.

| Study | Route and material | Node measurement | Why it is not a core-to-APC fraction |
| --- | --- | --- | --- |
| [Wang et al., 2018](https://www.thno.org/v08p0223.htm) | Mouse tail-base SC; fluorescent iTEP polypeptide or albumin-binding ABD-iTEP | At 6 h, inguinal plus axillary nodes contained **0.36 ± 0.02%** and **1.31 ± 0.11% of injected label**, respectively (n = 3). | Whole-node carrier label is not intact vaccine core. |
| [Roces et al., 2019](https://strathprints.strath.ac.uk/69584/7/Roces_etal_MP_2019_Scale_independent_microfluidic_production_of_cationic_liposomal_adjuvants.pdf) | Mouse quadriceps IM; 5 µg radiolabeled H56 protein adsorbed to DDAB:TDB liposomes | Popliteal-node antigen label was **<0.3% of initial vaccine dose** across four conventional formulations, sampled at 6–192 h. | Radioiodine tracks antigen-associated label, not a sequenced epitope in APCs. The result is formulation- and node-specific. |
| [Zhu et al., 2017](https://www.nature.com/articles/s41467-017-01386-7) | Mouse SC; radiolabeled DNA–RNA nanocapsules carrying peptide neoantigen | At 24 h, **0.13% injected dose** was in inguinal and **0.11%** in axillary nodes. | The quantitative label is on the *carrier*. Peptide and carrier colocalization does not yield peptide copies per APC. |
| [Vicente et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24508040/) | Rabbit SC; indium-111-labeled polyglucosamine/squalene nanocapsules | Popliteal node had **3.8 ± 1.2% injected dose** at 48 h by gamma scintigraphy. | This is radiolabeled adjuvant carrier, not antigen or core peptide; some label dissociated in serum. |
| [Glycosylated OVA study, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8498032/) | Mouse hock SC; fluorescent OVA versus OVA–p(GluNAc), 5 µg antigen/hock | Calibrated whole-node imaging reports percent starting antigen dose/hock; conjugate gave **17-fold** more node signal at its 8 h maximum and **10-fold** greater node-signal AUC. | The fold changes compare formulations, not an IM dose fraction or intact-core/APC amount. |
| [Backlund et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9371699/) | Mouse tail-base SC; fluorescent gp100 peptide versus cell-penetrating-peptide conjugates | At 48 h, inguinal-node fluorescence was **3.1–5.8-fold** higher for most conjugates than free peptide. | Relative fluorescence has no injected-dose denominator; label-positive APCs are a cell fraction, not dose recovered in APCs. |
| [Hassett et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10755037/) | Macaque IM; reporter mRNA–LNP, 300 µg RNA/site | Reporter RNA in draining node measured **155 ng/g** at the first 8 h tissue sample; a separate cohort had **9.1% EGFP-positive cDCs** at 24 h. | Tissue concentration and positive-cell fraction have different denominators and cohorts; neither is peptide dose in APCs. |

Wang and Roces are the closest **antigen-associated, injected-dose
denominator** checks here. Wang used subcutaneous carrier polypeptides rather
than a free SLP; Roces used an intramuscular liposomal *protein* formulation.
Both report material present in a node at a sample time. Our
`node_apc_fraction_of_dosed_peptide` instead counts cumulative first arrival
of **core-bearing lineages inside nodal APCs** by the chosen horizon. A node
snapshot may include extracellular material, may miss material that passed
through earlier, and cannot tell whether a label still contains the target
sequence. No numerical conversion between these endpoints is justified.

## Related denominators that cannot be read as percent of dose

[A subcutaneous mouse mRNA–nanoparticle study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5835019/)
reported a node peak of about **170% injected dose per gram** at 8 h for
tritium-labeled RNA. This is a *concentration normalized by tissue mass*
(`%ID/g`); it is not 170% of the dose in one node. The paper also found
fewer than 5% nanoparticle-positive node cells, with up to 50% of those
positive cells classed as CD11c⁺/CD123⁺ DCs. Neither cell percentage supplies
RNA copies, translated peptide copies or a core-bearing APC fraction.

[Roces et al.](https://strathprints.strath.ac.uk/69584/7/Roces_etal_MP_2019_Scale_independent_microfluidic_production_of_cationic_liposomal_adjuvants.pdf)
also tested an engineered biotinylated liposome with an avidin predose: node
**lipid-label** recovery reached roughly 30% in inguinal and 40% in popliteal
nodes. The deliberate trapping intervention and lipid analyte make those
values a separate formulation experiment, not a generic IM protein baseline.
It illustrates how sharply formulation can change node recovery.

## What to measure next

To calibrate a window-specific probability, measure a known IM dose of the
actual SLP and timed, sequence-confirmed **intact core-bearing fragments**
in injection-site fluid, lymph, whole node, sorted APCs and blood. Report
both percent injected dose in each compartment and absolute copies per APC
subset; separately measure surface pMHC and cell counts. For mRNA, measure
active RNA, produced antigen, direct APC expression, extracellular antigen
transfer and preformed-pMHC transfer. Matched time courses are needed to
separate production from loss. The remaining gap is tracked in
[calibration issue #549](https://github.com/openvax/mhctools/issues/549)
and the [calibration audit](vaccine-calibration.md).
