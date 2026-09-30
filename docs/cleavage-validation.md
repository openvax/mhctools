# Cleavage coverage and validation status

This review accompanies the contextual batch API. Software conformance,
source-observation reproduction and independent prediction validation are
different claims. No new held-out biological performance number is claimed.

## Proteasome variants

The [Pepsickle author implementation](https://github.com/pdxgx/pepsickle)
provides epitope, gradient-boosted digestion and neural digestion families.
The [paper](https://doi.org/10.1093/bioinformatics/btab628) describes its
independent digestion evaluation. The batch adapter's live checks compare
human/all-mammal neural C/I profiles with the upstream functions on the same
sequence and assert that C/I can change the result. These are conformance
checks, not independent accuracy measurements.

The [paper repository's MASTER.sh](https://github.com/pdxgx/pepsickle-paper/blob/master/MASTER.sh)
expects held-out inputs under `data/validation_data/digestion_data/raw/` and
processed, overlap-filtered validation FASTAs. Those paths are absent from
the inspected repository tree at `c448c4db81925afad78477e74a7d25e0209d3bce`. Its available `data/raw/digestion_map_files`
instead feed the training-data extraction path. Reusing those maps and
calling them held-out would be incorrect. Reconstructing the paper's
held-out data and verifying study/sequence overlap remains [#469](https://github.com/openvax/mhctools/issues/469)
with the existing assay-aware evaluation machinery in [#291](https://github.com/openvax/mhctools/issues/291).

The gradient-boosted artifact records scikit-learn **0.23.2**. It fails under
current scikit-learn (`sklearn.ensemble._gb_losses` is missing). The isolated
Python 3.8.20/0.23.2 runtime added for [#471](https://github.com/openvax/mhctools/issues/471)
executes both C/I routes and matches direct upstream inference with networking
disabled. It records the actual subprocess's package, code and weight identity.
The MAGE-A3 sequence in this conformance test occurs in upstream training data;
it is explicitly **not** held-out validation. Neural and gradient-boosted
predictions are not silently substituted.

## APC endolysosomal enzymes

The current built-in panel has no transferable cathepsin S/L/B or AEP model.
IRAP is an exact-substrate source catalog, and NetCleave-II is a class-II
C-terminal processing proxy, not an enzyme-specific cathepsin predictor.
The absence is explicit in batch coverage reports and tracked in
[#470](https://github.com/openvax/mhctools/issues/470), a focused child of #334.

Primary sources for the next validation block:

| Source | Useful evidence | Applicability constraints |
|---|---|---|
| [CatS/AEP processing of MBP](https://pubmed.ncbi.nlm.nih.gov/11745393/) | Product mapping, epitope generation/destruction and protection by HLA-DR | Purified/experimental processing conditions; not a general peptide survival model |
| [DIPPS specificity profiling](https://pubmed.ncbi.nlm.nih.gov/28733325/) | Cathepsin profiles and pH-dependent legumain specificity | Denatured in-gel substrates, condition-specific preferences; cannot turn sequence logos into calibrated probabilities |
| [Legumain activation](https://pubmed.ncbi.nlm.nih.gov/22232165/) | Activation and substrate-dependent pH response | Both activation state and P1 residue affect the result; unconditional cut-after-Asn rules omit this context |
| [Cathepsin MSP-MS comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC10399199/) | Human B/L/S and other cathepsins on 228 14-mers at pH 4.6/7.2, 15/60 minutes, 25 C | Distinct CatB dipeptidyl-carboxypeptidase and endopeptidase activity; source-specific product detection thresholds |

The MSP-MS publication deposits mass-spectrometry data as **MSV000090043 /
PXD035641**. Its workflow uses quadruplicates and significance/fold-change
criteria for detected products. Missing products are not automatically
non-cleaved bonds. The available source supports assay-aware curation; it
does not supply an already validated predictor of long-vaccine processing.
No unverified cathepsin sequence/product rows have been transcribed into this PR.

The batch `reference_panels` input can preserve curated experiments from
these studies now, with each condition in its own named panel. Actual
novel-sequence inference needs a separately verified adapter/data model.
ProsperousPlus availability and runtime/licensing review remains #281.

## Source-backed long-peptide case study

[Wada et al. 2018](https://doi.org/10.1371/journal.pone.0199249) is explicitly
assigned to validation in Pepsickle's Table 1. This PR now curates all 47
detected products from Figure 2A and 2C: two 31-residue vaccine constructs
containing the same epitopes in different orders, joined by RR linkers.
The dataset retains first-detection times, figure row IDs, parent endpoints
and the purified murine-i20S assay conditions. The correction was checked.

The real gradient-boosted immunoproteasome model scores the complete constructs.
Its native scores are paired with 24 observed internal construct/bond pairs;
parent sequence ends never become cleavage labels. The pinned raw training-map
audit covers 79 files and 58 distinct source sequences. It finds no exact
construct or observed seven-residue cleavage-context overlap. Wada's DOI is
absent from the reported training DOI fields; unresolved source identifiers,
homology-family overlap and the original fitted partition remain explicit
limitations. This is a small study-held-out case study, not a new calibrated
performance estimate or a reconstruction of all 225 author validation windows.

The [fixture README](../tests/data/wada2018/README.md) gives the primary source,
license, curation scope and reproducible command. The report preserves observed
products through batch save/reload, with endpoint/boundary/internal overlays.
No missing product is relabeled as a negative, and no product detection time is
converted into a predicted cleavage rate or presentation outcome.

Additional extracellular coverage is tracked in [#476](https://github.com/openvax/mhctools/issues/476):
CleaveNet's MMP substrate scores need their own native endpoint and assay
validation rather than conversion to per-bond probabilities.

## Required validation before biological ranking

Use exact source substrates/products, enzymes, species, chemical forms,
activation, pH, temperature, time and detection limits. Separate positive
product observations from verified non-cleavage measurements. Audit training
overlap by study and sequence/family before claiming held-out performance.
Preserve unsupported inputs and failures in the denominator.

Evaluate purified-enzyme cleavage, whole-matrix degradation and antigen
presentation as distinct endpoints. In particular, [human-DC long-peptide
cross-presentation](https://pubmed.ncbi.nlm.nih.gov/24587108/) supports a
possible proteasome/TAP route; it does not calibrate every construct or imply
that cleavage scores alone predict presentation. [IFN-gamma tumor-organoid
experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC10655636/) likewise show
why a single sequence-only profile cannot encode the tumor's cellular state.
