# Known gaps

What mhctools does not do yet, and where the work is tracked. This is the one
page that carries issue numbers; the rest of the documentation describes what
exists. It is reviewed at release time.

## Cleavage validation and coverage

- **Cathepsin S/L/B and AEP (legumain) prediction for novel sequences.** The
  built-in panel has no transferable model; the gap is shown explicitly in batch
  coverage reports. Experimental observations can be imported now through
  `reference_panels`. [#470](https://github.com/openvax/mhctools/issues/470)
- **ProsperousPlus adapter**, blocked on the open-license requirement.
  [#281](https://github.com/openvax/mhctools/issues/281)
- **Serum and extracellular peptidases beyond the current candidates**: activated
  blood and inflammation-associated proteases (thrombin, plasmin, kallikrein,
  elastase, cathepsin G, proteinase 3), other exopeptidases and
  compartment-specific enzymes. Curate activation, inhibitors, tissue exposure
  and assay conditions before choosing a predictor; a generic Arg/Lys or
  hydrophobic-residue scan would not distinguish these enzymes.
  [#278](https://github.com/openvax/mhctools/issues/278)
- **MMP substrate predictions** (CleaveNet), which need their own native endpoint
  and assay validation rather than conversion to per-bond probabilities.
  [#476](https://github.com/openvax/mhctools/issues/476)
- **Remaining intracellular peptidases**, including LAP3 and BLMH, and the
  unsettled THOP1 epitope-substrate disagreement in the primary reports. TPP2
  and NPEPPS also need exact substrate observations so their permissive rules
  gain reproduction records.
  [#334](https://github.com/openvax/mhctools/issues/334)
- **Conflicting IRAP precursor sequence** in Georgiadou 2010 (`DIRSSVQNKL` in
  the text and Table I, `DIRSSQVNKL` in the Figure 3F caption) blocks curating
  that record. [#332](https://github.com/openvax/mhctools/issues/332)

## Benchmarks and validation

- **Held-out, assay-specific benchmarking of the shipped cleavage models**,
  including observed non-cleavages, terminal modifications, homologous-sequence
  leakage checks and abstention. Purified-enzyme turnover and disappearance of
  intact peptide in serum are separate endpoints. Reconstructing the Pepsickle
  paper's processed validation partition belongs here too.
  [#291](https://github.com/openvax/mhctools/issues/291)
- **Half-life endpoints**: auditing pepADMET endpoints, overlap with the Tan
  model, and missing local inference artifacts.
  [#294](https://github.com/openvax/mhctools/issues/294)

## Downstream integration

Cleavage overlays, source-table reconciliation and ranking belong to other
projects; mhctools provides the evidence contract.

- Topiary: [combining source tables using original evidence](https://github.com/openvax/topiary/issues/366)
  and [exposing extracellular cleavage and peptide half-life evidence](https://github.com/openvax/topiary/issues/288).
- Vaxrank: [adopting generalized Topiary source tables](https://github.com/openvax/vaxrank/issues/497)
  for vaccine construction.

## Maintenance

- Making the SMM subset audit independent of zlib compression differences.
  [#465](https://github.com/openvax/mhctools/issues/465)
