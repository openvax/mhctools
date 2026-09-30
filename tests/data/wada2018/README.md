# Wada 2018 long-peptide digestion case study

This fixture transcribes **all 28 products in Figure 2A and all 19 products in
Figure 2C** from Wada H, Shimizu A, Osada T, Tanaka Y, Fukaya S, Sasaki E
(2018), *Development of a novel immunoproteasome digestion assay for synthetic
long peptide vaccine design*, PLoS ONE 13(7):e0199249.

- [Primary article and methods](https://doi.org/10.1371/journal.pone.0199249)
- [Original Figure 2](https://journals.plos.org/plosone/article/figure/image?size=large&id=10.1371/journal.pone.0199249.g002)
- [Correction](https://doi.org/10.1371/journal.pone.0205567): IFN-gamma typography;
  it does not replace the digestion maps.

Source material is copyright 2018 Wada et al., licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes here are
manual transcription into JSON, zero-based half-open product/epitope
coordinates, explicit source row IDs, and machine-readable assay annotations.
The downloaded original figure's SHA-256 is recorded in the JSON. Product
sequences and color-coded first-detection times were visually checked against
that figure; the parent sequences and epitope identities match Table 1.

The two constructs contain the same three epitopes and RR linkers in different
orders. The assay used purified **murine** 20S immunoproteasome at pH 7.8 and
37 C, with 20 micrograms of substrate and 2 micrograms of enzyme in 300
microliters. Purple/green/red products were first detected at 1/2/4 hours.
Product abundance and detection limits were not reported. Terminal chemical
states were not explicitly reported and remain unknown; Pepsickle inference
is explicitly conditional on its sequence-only representation.

Parent sequence ends are not cleavage events. Internal product boundaries
support observed cleavage, but a product's first detection does not date each
individual cut. Unreported products/bonds are not negative labels. Neither the
source's CTL assay results nor its claimed assay-to-CTL agreement are used as
cleavage labels here.

## Replay and training overlap

Run from a checkout with the legacy runtime provisioned:

```sh
python scripts/setup_test_backends.py pepsickle
source env/test-backends/activate.sh
python -m scripts.evaluate_wada_cleavage --out wada.json --html wada.html
```

To audit training overlap, add `--training-snapshot /path/to/pepsickle-paper`
pointing to the extracted author snapshot at
`c448c4db81925afad78477e74a7d25e0209d3bce`. The script requires an exact digest
of all 79 raw map files, including extensionless and multi-substrate files.
It audits study DOI, complete sequence and all possible seven-residue windows
centered on P1 with terminal padding. This raw-window inventory is a
conservative superset of selected training windows, not a reconstruction of
the original fitted partition.

`pepsickle-gb-evaluation.json` records one real inference run, all 47 products,
the 24 observed internal construct/bond pairs, exact runtime/code/weight
identity, and the completed raw-map audit. The audit finds no exact construct
or observed-window overlap among 58 distinct training source sequences. It
also retains map files with unresolved study identifiers; these do not become
evidence of complete study-level independence. CI replays the scores against
the pinned runtime, while source-data and report-integrity checks run offline.

[Pepsickle Table 1](https://doi.org/10.1093/bioinformatics/btab628) identifies
Wada 2018 among its six held-out studies. The study-held-out designation and
our explicit overlap audit are retained separately; family-level independence
and the absent processed validation partition remain unresolved. This is a
small, selected two-construct case study, not a reproduction of the paper's
full validation set. Native scores remain native scores. No AUC, specificity,
product-yield, serum-stability or presentation-accuracy estimate is claimed.
