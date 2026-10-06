# ITCell cathepsin specificity profiles

ITCell provides sequence-specific predictions for **human cathepsins B and S
(internal cleavage)** and **cathepsin H (initial N-terminal aminopeptidase
trimming)**. The nine published 15/60/240-minute count matrices are bundled:
there is no external installation or structural feature requirement.

```python
from mhctools import ITCellCleavage, get_cleavage_model

internal = ITCellCleavage(enzyme="S", minutes=240).predict("ACDEFGHIKLMNPQRSTVWYV")
initial_trim = get_cleavage_model("itcell-cath-240").predict("MALWMRLLPLL")
for site in internal.sites:
    print(site.bond, site.score, site.reason)
```

The native score sums log2 ratios between the source P4-P4prime specificity
profile and published human background frequencies. The released implementation
adds one count and uses the maximum column count plus 20 as the shared
denominator. Missing terminal flanks contribute zero. We reproduce that
formula, including its unnormalized background frequencies, and retain all
scores. Author-code candidate thresholds are **strictly >3 for B/S** and
**strictly >2 for H**. The paper's prose describes the >3 threshold as
threefold; the released code uses log2, so that interpretation is inconsistent.
Use the native score and the stated source-code threshold.

These are specificity scores, **not probabilities, percentages lost, kinetic
rates or serum survival**. The profile's minutes select the source assay
observations; they do not predict degradation after that duration in a new
sample. The source used 228 tetradecapeptides with recombinant human enzymes
at 0.2 ug/mL, pH 6.5 and 1 mM TCEP. It substituted norleucine for methionine and
did not represent cysteine. Outputs explicitly flag those contexts and absent
flanks. Predictions require canonical sequence with assumed free termini.

H scores only the currently exposed first bond. B models internal specificity,
not its carboxydipeptidase activity. No fragment generation, repeated trimming,
uptake, cellular accessibility or enzyme activation is simulated.

The independent implementation is tested against native scores from the
released Perl primitive. The author's combined `cleave.pl` wrapper repeats
the 240-minute S profile in place of the 15/60-minute profiles; we run each
matrix independently, and [reported this upstream](https://github.com/salilab/itcell-lib/issues/1). This is not the
structure-dependent 2023 PCSS/SVM cathepsin predictor.

Sources: [ITCell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6219782/),
[author code](https://github.com/salilab/itcell-lib),
[released matrices](https://doi.org/10.5281/zenodo.3227044).
The bundled count matrices and background frequencies retain the author's
LGPL-2.1 license, copyright attribution, exact file hashes and archive identity
in `mhctools/data/itcell_profiles.json` and `ITCELL_LICENSE.txt`.
