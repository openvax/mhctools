# Enzyme-centered predictions

Use a named enzyme class when the question starts with an enzyme: **what activity
has been measured, and what can each available method say about my peptide?**
It discovers the enzyme's named methods and runs them together, while keeping
their endpoints separate. It does not turn a pattern match or native score
into a probability of cleavage.

```python
from mhctools import AssayConditions, DPP4, FAP

dpp4 = DPP4()
for method in dpp4.describe()["methods"]:
    print(method["name"], method["evidence_type"], method["description"])

peptide = "HAE" + "G" * 27  # synthetic example, not a measured substrate
report = dpp4.predict(peptide)
print(report.summary())

# Same dimensionless endpoint for each assay-backed enzyme:
for enzyme in (DPP4(), FAP()):
    estimate = enzyme.predict(peptide).assay_fractions[0]
    print(enzyme.enzyme, estimate["parent_remaining_fraction"])

# A sensitivity scenario: 60 minutes at 0.2 µM active enzyme.
scenario = AssayConditions(exposure_minutes=60, enzyme_concentration_um=0.2)
print(dpp4.predict(peptide, conditions=scenario).assay_fractions[0])
```

The DPP4 report for this example contains four **independent** answers:

| Method | What it says | Example result |
| --- | --- | --- |
| `dpp4-qpisa` | Predicted substrate depletion in [Gudipati et al.'s purified-human-DPP4 assay](https://doi.org/10.1038/s44320-024-00071-4); native log2 score | 2.1694 at the exposed second bond |
| `dpp4-p1-pro-ala` | [Kühn-Wache et al.'s](https://doi.org/10.1515/bc.2011.028) P1 Pro/Ala preference; a binary motif match | Matches `HA\|E`; other P1 residues can still be cleaved |
| `dpp4-p1prime-nonpro` | The [qPISA study's](https://doi.org/10.1038/s44320-024-00071-4) tested P1′ non-Pro preference; a separate binary motif match | Matches `HA\|E`; a tested P1′ Pro non-match is not a universal no-cut claim |
| `dpp4-qpisa-keane2011-parent-loss` | Estimated **full-length parent** half-life and fraction remaining in the [Keane et al. purified-human-DPP4 assay](https://doi.org/10.1111/j.1742-4658.2011.08051.x) | About 17-minute half-life and 8.8% remaining after the default 60 minutes, within the calibration domain |

`HA|E` marks the bond after the first two residues, not a new peptide
sequence. The two motifs come from distinct source observations. A peptide
can miss either preference and still receive a qPISA score; no motif is used
as a hard gate for the quantitative model or the assay estimate. The Keane
fraction uses a fitted half-life and first-order approximation, and may
return `unassessed` outside its measured length and score range.

The default 60-minute exposure is a **reporting reference**, not an
observation about vaccination. Both enzymes emit `parent_remaining_fraction`
between 0 and 1, `parent_lost_fraction`, and `estimated_half_life_minutes`.
By default each uses its own purified-assay conditions: DPP4 used 0.112 µM
enzyme at pH 7.6, FAP used 0.2 µM at pH 8.0, and both used 1 µg peptide in
15 µL at 37 °C. Thus the fractions share an endpoint and time but are not
intrinsic cross-enzyme rates. Setting one `AssayConditions` concentration for
both gives a sensitivity comparison at the same assumed active-enzyme
concentration; it still does not equalize pH, substrate molarity, or enzyme
accessibility. Concentration scaling assumes the first-order parent-loss rate
changes linearly with active enzyme. There are no kinetic measurements here
to vary substrate loading, pH, or temperature.

Inspect the full source, assay, compartment, activity and limitations for
each method with `dpp4.describe()`. Retrieve native results by model name:

```python
score = report.by_model["dpp4-qpisa"].sites[0].score
mask = report.by_model["dpp4-p1-pro-ala"].pattern_mask
fraction = report.assay_fractions[0]
rows = report.summary_rows()  # one endpoint and status per method
payload = report.to_dict()    # JSON-compatible, schema_version 1
```

`EnzymePredictor.available()` lists supported enzymes. `FAP()`
similarly reports its N-terminal dipeptidyl and internal endopeptidyl motifs
plus the separate [Keane assay estimate](models.md#the-corresponding-purified-fap-estimate).
For an enzyme without an assay calibration, the result has no invented
fraction; the summary says so. If an enzyme requires an explicitly active
form, pass `enzyme_state="active"` (for example, CPB2). Optional model
runtimes are included only with `include_optional=True`.

The older `predict(peptide, minutes=10)` spelling still selects a 10-minute
source-assay exposure. New code should use `AssayConditions` when changing
exposure time or enzyme concentration. These parameters do not change motif
matches or native cleavage scores. The `minutes` argument of
`ITCellCleavage(...)` has a different meaning: it selects a published
specificity profile and does not request a new exposure prediction.

The `enzymes` namespace gives every default-cataloged enzyme an object, including
those without a dedicated class. Its attribute names remove spaces:

```python
from mhctools import enzymes

ace = enzymes.ACE
cats = enzymes.CathepsinS
print(cats.describe()["methods"])
```

For configuration files and command-line integrations, the explicit string
registry resolves named classes. Other cataloged enzymes fall back to the
same generic predictor:

```python
from mhctools import DPP4, FAP, EnzymePredictor, enzyme_predictors

assert enzyme_predictors["dpp4"] is DPP4
assert isinstance(EnzymePredictor.from_name("fap"), FAP)
ace = EnzymePredictor.from_name("ace")
```

`EnzymePredictor(DPP4)` also accepts the class. The older
`EnzymePredictor("DPP4")` constructor remains valid for callers of 4.0.5;
new Python code should use `DPP4()` directly.

The fraction is **parent peptide lost in one purified-enzyme assay**. It is
not serum half-life, a fraction of intact core epitope surviving, an IM-vaccine
delivery fraction or an APC-presentation probability. For candidate windows
and those different endpoints, see [vaccine window exposure](../vaccine-windows.md)
and [vaccine core calibration](../research/vaccine-calibration.md).
