# Core delivery and blood exposure after local vaccination

For a specified core, ask **what fraction of the starting peptide copies have
delivered that intact core into draining-node APCs by time T?** A delivered
fragment can be shorter than the input SLP. The reported fraction uses initial
peptide molecule counts, not the changing mass of trimmed products. It is
separate from the number still present in APCs at T and from surface pMHC.

`simulate_vaccine_trajectory` now tracks first presence independently of the
physical trajectory. `summarize_core_exposure(result)` exposes this accounting,
and `vaccine-trajectory --report` explains it alongside presentation exposure.
All rates and current-fragment assessments still require sources or explicit
scenario assumptions. Missing reachable kinetics give an unassessed result.

## Four distinct routes to nodal APC antigen

| Route | First-presence event |
|---|---|
| Uptake from node fluid | A core-bearing fragment enters a nodal APC |
| APC migration | A local APC carries core-bearing antigen or pMHC into the node |
| Endogenous expression | A nodal APC produces antigen from active vaccine RNA |
| Cross-dressing | A nodal APC acquires a preformed producer-cell pMHC-I complex |

The observer absorbs an antigen lineage at its first event. The physical
trajectory continues through processing, secretion or turnover; it cannot add
another arrival for the same lineage. RNA remains a repeatable source. Thus RNA
outputs report **antigen copies per input transcript**, which can exceed one,
and a separate fraction of antigen produced by the observation horizon. They
are never labeled a fraction of injected RNA delivered as peptide.

## Blood persistence of a core

Each time point reports blood antigen in two disjoint groups: the exact free
core and all longer fragments containing it. Their sum and integrated copy-hour
exposure remain available even when the original parent has disappeared.

After IM injection, blood antigen can continue arriving from the injection
site, lymph and, for secreted mRNA constructs, producing cells. The whole curve
therefore does not identify an intrinsic circulation elimination half-life.
The report gives a **post-peak half-decline bracket on the observation grid**:
the first downward crossing of half the largest sampled blood amount, measured
from that sampled peak. It neither extrapolates a censored tail nor claims that
the peak was resolved between observations. No blood signal or no observed
crossing is an explicit status. Refine the grid and inspect the full curve;
multiple input phases can cause later rebounds.

## mRNA biology and supported routes

Secretion followed by APC uptake is one pathway, not a universal serial model
of IM mRNA vaccines. [Liang 2017](https://pubmed.ncbi.nlm.nih.gov/28958578/) and
[Hassett 2024](https://pubmed.ncbi.nlm.nih.gov/38161733/) support direct mRNA
uptake/expression at the injection site and in draining nodes. In mouse studies,
[Rood 2026](https://pubmed.ncbi.nlm.nih.gov/41545353/) found an important role for
endogenous APC antigen in MHC-II responses, and
[Jo 2026](https://pubmed.ncbi.nlm.nih.gov/41986715/) identified a substantial
MHC-I cross-dressing contribution. These findings do not identify universal
human pathway proportions.

The [mRNA antigen-presentation bibliography](evidence/mrna-antigen-presentation.md)
summarizes these and related primary studies and provides an assay-scoped
measurement extract for future calibration.
The [core-specific calibration audit](research/vaccine-calibration.md)
separates measured pMHC time courses from intrinsic complex lifetimes and
gives the physical bounds on first-arrival fractions.

- `secreted_mrna` follows a supplied signal boundary, ER entry, secretory
  residence/cuts, secretion, uptake and subsequent APC processing. Failed ER
  entry and ER-associated retrotranslocation can also feed cytosolic processing.
- `cytosolic_mrna` places completed antigen in the transfected cell's cytosol.
  Direct APC processing requires no secretion. Producer-cell release is a
  separately supplied effective rate; no secretion signal is inferred.
- Both routes can include producer-cell MHC-I processing and transfer of its
  loaded complex to local or nodal APCs. The modeled donor allele is the
  specified allele. Transfer does not increment newly loaded APC complexes.

Producer-cell expression, release, processing and transfer require explicit
rates and fragment assessments. No MHC-II expression is inferred for non-APC
producer cells. Membrane topology, vesicles, debris, APC subsets, changing
inflammation and ligand competition are not separately resolved. These are
conditional kinetic scenarios, not a calibrated immune-response model.

Version 3.52.0 adds these routes. Older RNA input manifests may therefore need
additional explicit rates or declared exclusions. Never replace new unknown
rates with zero merely to obtain a numerical output.
