# Antigen processing and T-cell responses

The causal chain is antigen availability → uptake or expression → generation
*or destruction* of a core → loading and surface peptide–MHC (pMHC) → T-cell
recognition and priming. Measuring a step does not measure all later steps:
binding affinity does not count presented complexes, and a responding T-cell
clone does not establish the amount of antigen that arrived at an APC.

## Primary experiments

| Study | Experiment and conclusion |
| --- | --- |
| [Manoury et al., 2002](https://www.nature.com/articles/ni754) | In human APCs processing myelin basic protein, inhibiting asparagine endopeptidase (AEP) increased presentation of the HLA-DR2 MBP(85–99) epitope; overexpressing AEP reduced it. AEP cleaves within the epitope, showing that more protease activity can *destroy* a target rather than generate it. |
| [Bijker et al., 2008](https://doi.org/10.1002/eji.200737995) | Compared minimal and extended OVA peptide vaccines with CpG in mice. Extended peptides yielded more sustained, dendritic-cell-focused presentation and stronger tumoricidal CTL activity; repeated minimal peptide dosing increased CTL counts without reproducing that effector quality. |
| [Rosalia et al., 2013](https://doi.org/10.1002/eji.201343324) | Compared whole protein with SLP uptake and processing in murine and human dendritic cells. SLPs reached a faster route to class-I presentation and increased CD8 activation in their in-vitro systems; protein and SLP therefore cannot share a generic processing multiplier. |
| [Rabu et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6422379/) | Varied cathepsin-sensitive linkers in artificial melanoma SLPs, then tested DC presentation, donor PBMC responses and a mouse tumor model. Some linkers improved class-I clone activation by up to about 100-fold and enhanced downstream responses, while class-II effects were smaller; sequence-specific processing can alter both endpoints. |
| [Influenza epitope quantification, 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6599079/) | Quantified 21 class-I influenza ligands after direct and cross-presentation and compared them with mouse CTL responses. The immunodominant NP and PA epitopes favored different routes; pMHC abundance contributed to response magnitude together with affinity and T-cell repertoire, rather than determining it alone. |

## Use in mhctools

The approximately 100-fold linker result is a *relative concentration-response
shift in one antigen, assay and APC system*. It is useful as a challenge case
for candidate-window ranking and to prioritize paired experiments, but it is
not a transferable per-hour cut rate or a dose fraction. Likewise, prolonged
in-vivo presentation after SLP vaccination includes continuing antigen input;
it does not identify an intrinsic surface-pMHC half-life.

Report the stage actually assessed: protease-site score, predicted binding,
measured surface pMHC, T-cell activation, or in-vivo response. The
[trajectory model](../vaccine-trajectory.md) and [window comparison](../vaccine-windows.md)
can propagate supplied rates and sensitivity scenarios, but their numerical
output cannot be calibrated from these response studies alone. The joint
[calibration issue](https://github.com/openvax/mhctools/issues/549) tracks the
missing matched, time-resolved core-bearing fragments and APC pMHC.
