# Known gaps

What mhctools does not do yet, and where the work is tracked. This is the one
page that carries issue numbers; the rest of the documentation describes what
exists. It is reviewed at release time.

## Cleavage validation and coverage

- **Cathepsin L and AEP (legumain), and cellular validation of cathepsin profiles.**
  Bundled [ITCell B/S/H specificity profiles](cleavage/itcell.md) now score new
  canonical sequences under their published assay scope. CatL/AEP remain
  unavailable, and cellular/pH/exposure calibration is still a gap. Experimental observations can be imported through
  `reference_panels`. [#470](https://github.com/openvax/mhctools/issues/470)
- **ProsperousPlus adapter**, blocked on the open-license requirement.
  [#281](https://github.com/openvax/mhctools/issues/281)
- **Serum and extracellular peptidases beyond the current candidates**: activated
  blood and inflammation-associated proteases (thrombin, plasmin, kallikrein,
  proteinase 3), other exopeptidases and
  compartment-specific enzymes. Curate activation, inhibitors, tissue exposure
  and assay conditions before choosing a predictor; a generic Arg/Lys or
  hydrophobic-residue scan would not distinguish these enzymes.
  [#278](https://github.com/openvax/mhctools/issues/278). Native
  [PhageScout ELANE/CTSG sequence profiles](cleavage/phagescout.md) are available;
  independent inflammatory/whole-serum validation remains open.
- **Concrete additional protease candidates**, from the
  [primary-source and artifact audit](cleavage/candidates.md). The remaining
  tasks supplement the available PhageScout native sequence tracks:

  | Candidate | Next step | Tracking |
  |---|---|---|
  | PhageScout ELANE/CTSG classifiers | Obtain exact training medians and mature-protein normalization context; keep structural models separate | [#521](https://github.com/openvax/mhctools/issues/521) |
  | CatL | Add assay-scoped sequence-recognition and observed-site evidence | [#514](https://github.com/openvax/mhctools/issues/514) |
  | AEP/legumain | Add pH-scoped Asn/Asp recognition and observed sites | [#515](https://github.com/openvax/mhctools/issues/515) |
  | Thrombin | Curate extended cooperative recognition alternatives | [#516](https://github.com/openvax/mhctools/issues/516) |
  | Plasmin | Keep Lys/Arg-conditioned assay profiles and subsite interactions separate | [#517](https://github.com/openvax/mhctools/issues/517) |
  | Plasma kallikrein KLKB1 | Add human plasma-enzyme evidence separately from tissue KLKs | [#518](https://github.com/openvax/mhctools/issues/518) |
  | ELANE/PRTN3 | Add distinct inflammatory-enzyme recognition and independent benchmarks | [#519](https://github.com/openvax/mhctools/issues/519) |

  PhageScout native sequence artifacts are reproduced. Sequence-recognition
  rules can progress independently when their source scope is clear; external
  accuracy and whole-matrix validation depend on the benchmark work below.
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
  intact peptide in serum are separate endpoints. Reconstructing the [Pepsickle](predictors/processing.md#pepsickle)
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

- Making the [SMM](predictors/binding.md#smm-and-smm-pmbec) subset audit independent of zlib compression differences.
  [#465](https://github.com/openvax/mhctools/issues/465)
