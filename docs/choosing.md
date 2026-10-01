# Choosing a predictor

Start from the question you are asking, then check the practical constraints
(license, installation, input) before the model details.

## By question

| I want to know… | Kind | Start with |
|---|---|---|
| Will this peptide bind my class I allele? | `pMHC_affinity` | `NetMHCpan` (DTU license) or `MHCflurry` (open); both report IC50 in nM, so the two are directly comparable |
| Will it be presented on the cell surface? | `pMHC_presentation` | `NetMHCpan41`/`42`, `MHCflurry`; add `MixMHCpred`, `BigMHC_EL` or `CapHLA` as independent opinions |
| The same for class II | `pMHC_affinity` / `pMHC_presentation` | `NetMHCIIpan` 4.3 and `MixMHC2pred`; `MixMHC2pred` and NetMHCIIpan were independently co-best in a 2024 class II benchmark |
| How long the pMHC complex lasts | `pMHC_stability` | `NetMHCstabpan` |
| Whether the peptide is cut out of the protein | `proteasome_cleavage` | `Pepsickle` (open) or `NetChop`; class II uses `NetCleave_II` |
| Whether it reaches the ER | `tap_transport` | `DeepTAP` |
| Whether ERAP1 trims a precursor | `erap_trimming` | `ERAMER` |
| Whether a T cell responds | `immunogenicity` | `Calis` as a baseline; `PRIME`, `BigMHC_IM`, `DeepImmuno`; `TLimmuno2` for class II |
| Whether a specific TCR recognises it | `pMHC_TCR_binding` | `NetTCR`, `Tulip`, `MixTCRpred` |
| How long a free peptide survives in serum | `peptide_half_life` | `PeptiVerse`, `PlifePred2` |
| Which peptidase cuts which bond | none | the [cleavage API](cleavage/index.md) |

The [predictor matrix](predictor-matrix.md) has the full list with inputs and
install routes.

## By constraint

**No license to request or sign.** `MHCflurry`, `CapHLA`, `SMM`/`SMMPMBEC`,
`Pepsickle`, `DeepTAP`, `DeepImmuno`, `Calis` and `RandomBindingPredictor` are
open or built in. The DTU tools need a license you request from DTU; the Gfeller
lab tools (`MixMHCpred`, `MixMHC2pred`, `PRIME`, `MixTCRpred`) and `BigMHC` and
`NetTCR` are academic, non-commercial. See [licensing](licensing.md).

**Nothing to download.** Only `Calis` and `RandomBindingPredictor`.

**Runs in your Python without a separate environment.** `MHCflurry`, `CapHLA`,
`SMM`, `Calis`, `Pepsickle` (its neural models) and `NetTCR`. The DTU tools and
the Gfeller tools are external executables. The torch and TensorFlow models
(`DeepTAP`, `DeepImmuno`, `TLimmuno2`, `MixTCRpred`, `Tulip`, `PeptiVerse`,
`PlifePred2`) run in a separate interpreter; see [environment
variables](env-vars.md).

**Peptides only, no allele.** `Calis`, `DeepTAP`, `ERAMER`, `PeptiVerse`,
`PlifePred2`, and the cleavage predictors (which also want flanks).

## Good habits

- **Use more than one predictor for the same kind** and compare. Everything in
  one table is the point of mhctools; see [recipes](recipes.md) and
  `mhctools predict-table`.
- **Compare like with like.** `value` is comparable across predictors of the
  same kind (IC50 in nM). `score` is predictor-specific and is not. See
  [prediction kinds](kinds.md#score-value-and-percentile_rank).
- **Do not rank neoepitopes by immunogenicity alone.** Every current CD8
  immunogenicity predictor falls toward chance on unseen tumor neoepitopes; see
  the [warning](predictors/immunogenicity.md#read-this-before-trusting-a-score).
- **Read the known limits** for the predictors you pick: [limitations](limitations.md).
