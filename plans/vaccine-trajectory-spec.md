# SLP and secreted-mRNA trajectory specification

Build an exploratory, single-target compartment model for local IM/SC delivery.
The endpoint is target loading and surface pMHC on local versus draining-node
APCs. Preserve the target in shorter fragments, and distinguish destruction,
clearance, RNA loss and pMHC turnover.

1. Free SLP starts in injection-site interstitium. mRNA starts as carrier-bound
   RNA, with competing local uptake, lymph-node drainage and loss. Separate
   producer-cell and APC endosomal escape, RNA decay and repeated translation.
2. Secreted mRNA explicitly supplies the full translated construct and signal
   peptide boundary. Model ER entry/signal removal, secretion, failed entry and
   ER-associated loss/retrotranslocation as competing effective steps.
3. Extracellular antigen branches through interstitium, afferent lymph,
   draining-node fluid and blood. Include local uptake and APC migration;
   neither modality must pass through blood before presentation.
4. MHC-I: endosomal escape, cytosolic cuts, TAP transport of eligible precursors,
   ER amino-terminal trimming, exact-ligand loading and surface trafficking.
   An optional exact-ligand vacuolar loading channel is separate. MHC-II:
   endolysosomal cuts and loading of an explicitly bounded ligand containing
   the exact selected span. Loading lumps groove availability/CLIP/DM editing.
5. Cleavage consumes a current fragment and retains its target-bearing product.
   Rates are supplied per compartment and current fragment; mechanisms and
   provenance survive the calculation. Bound pMHC leaves the free-peptide cut
   model, then undergoes separately supplied trafficking and turnover.
6. Solve a linear first-order expected-copy system without new heavy runtime
   dependencies. Translation creates antigen without consuming RNA. Report
   antigen conservation, cumulative loading and source-normalized units.
7. Missing reachable rates or cleavage assessments prevent numerical prediction.
   Explicit zero rates require a stated assumption. Native predictor outputs
   remain evidence, not loading probabilities or enzyme hazards. No biological
   rate defaults, cross-estimator averaging or serum-to-lymph transfer.
8. Deliver a public Python API, JSON CLI, complete mechanism/assumption audit,
   reproducible synthetic scenarios, route guide and meaningful analytical,
   conservation, fragment, loading-gate and missing-input tests.

Excluded: whole-body PBPK, sequence-derived transfection/secretion prediction,
active-enzyme calibration, saturation/MHC competition, cell-subset heterogeneity,
cross-dressing, antibody/T-cell response and clinical efficacy. These must stay
visible in the model audit. Synthetic numerical examples validate behavior and
illustrate sensitivity; they are not Sid or patient forecasts.
