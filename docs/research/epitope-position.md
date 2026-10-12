# Epitope position and peptide flanks

**Question:** does a core epitope near the N or C end of a protein or dosed
synthetic long peptide (SLP) reach MHC more often? The experiments below show
that *which residues flank the core* and *how the antigen enters the processing
pathway* matter. They do not establish a universal “put the core at the end”
rule. A position in a translated, folded protein is a different intervention
from moving the same core within an exogenous SLP.

## Controlled observations

| Primary experiment | What changed and what was measured | Result and scope |
| --- | --- | --- |
| [Gileadi et al., 1999](https://pubmed.ncbi.nlm.nih.gov/10427984/) | Moved two influenza nucleoprotein class-I epitopes to the N terminus of an expressed cytosolic protein; tested CTL recognition and vaccinia clearance. | N-terminal placement blocked presentation in this construct; adding at least two residues before the core or changing its N flank rescued it. This directly contradicts a general N-terminal advantage for *proteins*. |
| [Hearn et al., 2004](https://pubmed.ncbi.nlm.nih.gov/15342647/) | Compared EtxB-conjugated peptides ending at a class-I core with otherwise related C-terminally extended peptides; tested presentation kinetics and processing inhibition. | The C-terminal extension made presentation slower and proteasome-dependent in this delivery route. This supports a C-boundary constraint, not a rate for free IM SLP. |
| [ER trimming study, 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2855122/) | Synthesized N-extended SIINFEKL precursors and changed residues removed by ER aminopeptidases; measured trimming and presentation. | N-terminal precursor sequence changed generation of the final ligand; an N extension is not automatically fatal because trimming can occur. This is endogenous class-I processing, not extracellular dose survival. |
| [Rabu et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6422379/) | Connected melanoma class-I and class-II epitopes with alternative cathepsin-sensitive linkers in artificial SLPs; used human DC and T-cell clones, PBMC and transgenic mice. | Some linkers increased class-I cross-presentation by up to about 100-fold relative to a native long peptide; class-II presentation changed less. This is a within-system T-cell-assay potency comparison, not a 100-fold APC uptake fraction. |

An [ex-vivo long-peptide comparison](https://doi.org/10.1182/blood.V112.11.3904.3904)
reported about a ten-fold drop in T-cell activation for some C-terminally
extended CMV epitope variants versus N-terminally extended variants. It is a
conference abstract with antigen- and cell-dependent outcomes; it should
motivate matched experiments, not supply a universal correction factor.

## What this means for window selection

For each core and candidate window, retain the exact native or construct
flanks, core coordinates, class, formulation and processing route. Compare
candidate windows using the same model and assay context. For class I, a
correct C terminus must be generated or already exposed, while N-terminal
precursors may be trimmed; neither fact determines the full presentation
probability. For class II, a longer binding register and protection or
destruction during endosomal processing add different constraints. The
[processing-response review](processing-t-cell-responses.md) addresses why a
presentation advantage need not equal a T-cell advantage.

None of these studies measures the fraction of an IM free-SLP microgram
arriving in draining-node APCs as an intact, core-bearing lineage. Keep those
[candidate-window exposure](../vaccine-windows.md) outputs as explicit
sensitivity scenarios pending [calibration](https://github.com/openvax/mhctools/issues/549).
