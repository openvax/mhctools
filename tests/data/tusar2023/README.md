# Cathepsin L source-observation example

`catl-protected-peptides.json` preserves all labeled CatL products for **p3
(KDLLHPSP)** and **p4 (EIDLRNPKGN)** on page 1 of Supplementary Data 4 of
[Tusar et al. 2023](https://doi.org/10.1038/s42003-023-04772-8).
This is a selected two-substrate demonstration, not the complete study.

The chromatograms were visually inspected. Supplementary Table 17 independently
lists the same internal cuts: p3 bonds 3/4, p4 bonds 4/5/8. Supplementary Table 3
marks both substrates as protected; the article defines this as N-acetylation
and C-amidation. The peptide assay used 1 micromolar CatL at pH 5.5, 37 C for
2 hours. Its conditions must not be confused with the separate spike-protein
experiment at pH 6.5. Numeric substrate concentrations and detection limits
are unreported. No unreported site receives a negative label.

The article and supplements are CC BY 4.0. Changes are transcription to JSON,
zero-based half-open product coordinates, source IDs and condition annotations.
The source PDF hash is retained in the fixture. The intervals in `epitopes`
are product-overlay examples, **not assertions of MHC binding or presentation**.

Replay with no optional runtime:

```sh
mhctools cleavage --input tests/data/tusar2023/catl-protected-peptides.json \
  --out catl.json --html catl.html
```

The imported reference panel accepts only the observed sequence/chemical form.
Changing either terminal modification or the sequence causes abstention.
The original products and assay provenance survive save/reload in the inputs
and reference panel. This adds source evidence; it does not solve transferable
cathepsin S/L/B or AEP prediction (#470).
