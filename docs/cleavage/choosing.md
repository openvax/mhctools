# Choosing processing and peptidase models

Choose a biological endpoint and a plausible compartment before choosing an
enzyme. The recommendations below reflect the model's endpoint and available
coverage. They do not rank predictors by independently established accuracy.
See [antigen processing](../predictors/processing.md) for the pathway overview.

## Class I peptide generation

Start with [Pepsickle](../predictors/processing.md#pepsickle) for proteasomal
site prediction, with [NetChop](../predictors/processing.md#netchop) as another
model to compare. Retain the parent sequence and flanks when asking whether an
epitope boundary can be generated. The [Pepsickle paper](https://doi.org/10.1093/bioinformatics/btab628)
describes separate epitope-trained and digestion-trained families.

Use the epitope-trained family for an epitope-processing question. Use the
[constitutive and immunoproteasome digestion variants](models.md#proteasome-models)
when comparing those enzyme profiles. Select the profile from your biological
scenario; a tumor or APC label alone does not establish proteasome composition.
For full proteins and vaccine constructs, use [batch assessments](batch.md).
A cleavage-site score alone does not establish MHC presentation.

## ER trimming

Use [ERAMER](../predictors/processing.md#eramer) to assess a 9–16-residue
N-terminally extended precursor toward a target epitope. Choose its cascade
API for the aggregate trimming score, or [ERAMER's single-step model](models.md#erap1-and-existing-processing-models)
for the first exposed bond. These endpoints differ.

[ERAP2 motif evidence](models.md#the-models), available as `erap2-basic`,
can flag its Arg/Lys N-terminal preference. It is a limited rule, not a
full-context ERAP2 predictor. [ERAP1 experiments](https://pubmed.ncbi.nlm.nih.gov/12436110/)
show that trimming can generate or destroy epitopes, so more trimming is not
universally favorable. Pair this analysis with [TAP transport](../predictors/processing.md#deeptap)
and [MHC binding](../predictors/binding.md) as separate endpoints.

## Cytosolic trimming and degradation

For an exposed precursor terminus, select a
[cytosolic model](models.md#cytosolic-and-endosomal-candidates) matching the
question: aminopeptidase P for an N-terminal X-Pro bond, DPP8/9 for dipeptide
removal, or TPP2 for tripeptide removal. Their motif rules supply recognition
evidence. TPP2 and NPEPPS have permissive rules with weak selectivity information.

[THOP1 and neurolysin](models.md#cytosolic-and-endosomal-candidates) provide
exact-substrate experimental observations. Use these to inspect documented
substrates; they abstain on novel sequences. None of these models predicts
whole-cell peptide survival or competition with TAP and the proteasome.

## Class II endolysosomal processing

Use [NetCleave class II](../predictors/processing.md#netcleave) for a peptide's
C-terminal processing proxy, with at least three downstream residues. Its
[paper](https://pubmed.ncbi.nlm.nih.gov/34162981/) reports weaker performance for
class II than class I. It does not identify which cathepsin cuts the bond.

For enzyme-specific APC questions, cathepsins and AEP/legumain are relevant,
but mhctools currently has no transferable model for their activity on novel
sequences. [Cathepsin S experiments](https://pubmed.ncbi.nlm.nih.gov/9616206/)
establish a role in class II presentation; they do not validate a generic
sequence-only predictor. Use the [curated reference panels](validation.md#apc-endolysosomal-enzymes)
for their measured substrates and conditions. Keep these coverage gaps explicit
when combining results with [class II binding predictions](../predictors/binding.md#netmhciipan).

## Endosomal cross-presentation

[IRAP/LNPEP](models.md#cytosolic-and-endosomal-candidates) is an endosomal
aminopeptidase implicated in [class I cross-presentation](https://pubmed.ncbi.nlm.nih.gov/19498108/).
The `lnpep-observed` model is an exact-substrate lookup, not a general
cross-presentation predictor. Its source assay used purified enzyme at pH 8.0;
do not interpret it as calibration in an acidified endosome.

Use [batch scenarios](batch.md#scenarios-and-trimming) to record an explicit
routing hypothesis and any conditional fragments. Compartment labels filter
models by enzyme location; they do not predict where an antigen travels.

## Extracellular peptide degradation

For a human DPP4 substrate, use [DPP4 qPISA](models.md#human-dpp4-qpisa) to assess
removal of the first two residues from a free N-terminus. The
[qPISA study](https://doi.org/10.1038/s44320-024-00071-4) models substrate depletion
by purified enzyme under its assay conditions, not serum half-life.

For broader candidate screening, select
[extracellular models](models.md#serum-and-extracellular-candidates) by the bond:
CPN or activated CPB2 for a C-terminal Lys/Arg; ACE for ordinary C-terminal
dipeptide removal; aminopeptidase P for N-terminal X-Pro; or FAP for its
proline-associated activities. MME, ENPEP, and ANPEP provide narrower preference
flags. Preserve terminal chemistry and, for CPB2, the explicit active-enzyme
assumption.

Use `--compartment extracellular` for the broader annotated panel, or select
`serum` or `plasma` deliberately. These filters establish neither enzyme
concentration nor calibration in that fluid. For a duration endpoint, consult
[peptide half-life models](../predictors/peptide-pk.md); for delivery, consult
[uptake and exposure results](../exposure-results.md). Do not combine the
peptidase outputs into a serum-stability probability.

## Whole-substrate MMP scoring

[CleaveNet](cleavenet.md) supplies local whole-substrate MMP Z-scores and ensemble
spread, with optional ten-residue window scanning. Use its dedicated results
when you need this assay-scoped endpoint; it does not produce bond tracks.
