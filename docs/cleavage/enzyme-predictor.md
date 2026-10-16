# Enzyme-centered predictions

Use a named enzyme class when the question starts with an enzyme: **what activity
has been measured, and what can each available method say about my peptide?**
It discovers the enzyme's named methods and runs them together, while keeping
their endpoints separate. It does not turn a pattern match or native score
into a probability of cleavage.

```python
from mhctools import DPP4

dpp4 = DPP4()
for method in dpp4.describe()["methods"]:
    print(method["name"], method["evidence_type"], method["description"])

peptide = "HAE" + "G" * 27  # synthetic example, not a measured substrate
report = dpp4.predict(peptide, minutes=10)
print(report.summary())
```

The DPP4 report for this example contains four **independent** answers:

| Method | What it says | Example result |
| --- | --- | --- |
| `dpp4-qpisa` | Predicted substrate depletion in [Gudipati et al.'s purified-human-DPP4 assay](https://doi.org/10.1038/s44320-024-00071-4); native log2 score | 2.1694 at the exposed second bond |
| `dpp4-p1-pro-ala` | [Kühn-Wache et al.'s](https://doi.org/10.1515/bc.2011.028) P1 Pro/Ala preference; a binary motif match | Matches `HA\|E`; other P1 residues can still be cleaved |
| `dpp4-p1prime-nonpro` | The [qPISA study's](https://doi.org/10.1038/s44320-024-00071-4) tested P1′ non-Pro preference; a separate binary motif match | Matches `HA\|E`; a tested P1′ Pro non-match is not a universal no-cut claim |
| `dpp4-qpisa-keane2011-parent-loss` | Fraction of **full-length parent peptide** predicted lost after a specified time in the [Keane et al. purified-human-DPP4 assay](https://doi.org/10.1111/j.1742-4658.2011.08051.x) | About 33% lost after 10 minutes, within the calibration domain |

`HA|E` marks the bond after the first two residues, not a new peptide
sequence. The two motifs come from distinct source observations. A peptide
can miss either preference and still receive a qPISA score; no motif is used
as a hard gate for the quantitative model or the assay estimate. The Keane
fraction uses a fitted half-life and first-order approximation, and may
return `unassessed` outside its measured length and score range.

Inspect the full source, assay, compartment, activity and limitations for
each method with `dpp4.describe()`. Retrieve native results by model name:

```python
score = report.by_model["dpp4-qpisa"].sites[0].score
mask = report.by_model["dpp4-p1-pro-ala"].pattern_mask
fraction = report.by_model["dpp4-qpisa-keane2011-parent-loss"]
rows = report.summary_rows()  # one endpoint and status per method
payload = report.to_dict()    # JSON-compatible, schema_version 1
```

`EnzymePredictor.available()` lists supported enzymes. `FAP()`
similarly reports its N-terminal dipeptidyl and internal endopeptidyl motifs
plus the separate [Keane assay estimate](models.md#the-corresponding-purified-fap-estimate).
For an enzyme without an assay calibration, `minutes` produces no invented
fraction; the summary says so. If an enzyme requires an explicitly active
form, pass `enzyme_state="active"` (for example, CPB2). Optional model
runtimes are included only with `include_optional=True`.

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
