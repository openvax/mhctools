# Peptidase contributions to target availability

`simulate_enzyme_degradation` models successive cuts with **absolute,
enzyme-labelled rates**. `enzyme_removal_effects` asks how much more or less of
an exact target remains when one enzyme is removed. A surviving target counts
in either the full parent or any shorter fragment. Each new fragment is
reassessed, so newly exposed termini can change the next step.

This is a conditional extracellular-digestion model. It does not supply
universal human enzyme weights, patient circulation PK, uptake or presentation.
The bundled empirical evidence supports substrate-specific reference
calculations and a separate broad terminal-trimming scenario.

## What the available data support

| Evidence | Quantitative input | Supported use |
|---|---|---|
| [Abid 2009, human serum](https://doi.org/10.1074/jbc.M109.035253), Table 1 | Apparent Km and serum-normalized Vmax for DPP4, aminopeptidase P, and plasma kallikrein | Concentration-dependent local rates for the specific oxidized NPY forms; DPP4 and aminopeptidase P compete on the parent, kallikrein acts on a trimmed product |
| [Kuoppala 2000, human plasma](https://doi.org/10.1152/ajpheart.2000.278.4.H1069) | ACE pathway >90% at low nanomolar bradykinin; CPN-like pathway >90% at high micromolar bradykinin | Reference bounds demonstrating concentration-dependent dominance, without transferring these percentages to vaccine peptides |
| [Maffioli 2020, pig plasma](https://doi.org/10.3390/molecules25184071), Table S1 | 113 unique observed boundary contexts across a 228-peptide library; immediate N1/N2/C1/C2 boundary counts 46/30/7/7 | An explicitly exploratory terminal-only hazard allocation, with absolute kinetics supplied separately |
| [Semis 2019, mouse plasma peptides](https://pmc.ncbi.nlm.nih.gov/articles/PMC6693507/) | 244 possible ACE substrates/products; broad discovery used isolated plasma peptides plus purified ACE | Broad ACE substrate evidence, without a whole-plasma enzyme rate calibration |

The inventory includes assay species, matrix, chemistry, concentration range,
source locations and download hashes. Library M placeholders denote
**norleucine**, not methionine. The source's `Filtered` worksheet contains 114
qualifying product rows; its `Cleavage Sites` worksheet contains 113 unique
contexts. Context counts are not numbers of molecular cleavage events.

FAP, MME and activated CPB2 retain their separate specificity evidence. This
inventory supplies no applicable vaccine serum-rate calibration for them.
Enzyme abundance, recognition/enrichment scores, and absence of a motif are
insufficient to fill that gap. Intracellular processing stays outside this model.

## Measured reference calculations

```python
from mhctools import serum_reference_kinetics

reference = serum_reference_kinetics(
    "oxidized-NPY-1-36", concentration_um=5, serum_fraction=0.1)
for enzyme in reference["enzymes"]:
    print(enzyme["enzyme"], enzyme["local_rate_per_hour"],
          enzyme["share_of_reported_channels"])
```

For this reference, the local rates are **8.29/h for DPP4** and **0.446/h for
aminopeptidase P**: approximately 95% and 5% of these two reported channels.
Kallikrein is not included in that denominator because its measured substrate
was NPY(3-36), not NPY(1-36). Other observed pathways were not quantitatively
calibrated. These shares are calculated from the paper's Km/Vmax point
estimates, rather than being independent observed serum-loss percentages.

The conversion is `k = 60 * serum_fraction * Vmax / (Km + concentration)`.
Vmax is pmol/min/µL serum, equivalent to µM/min per unit serum fraction. Km and
concentration are in µM; k is in inverse hours. At finite concentration this
is a **local** rate: it changes as substrate concentration falls and is not a
constant half-life for the entire digestion trajectory. A mass-action
first-order approximation requires concentrations well below Km.

The study used mono-oxygenated peptide, 37 °C, pH 8 and added manganese.
Changing dilution assumes linear activity scaling. Calls outside the reported
5–70 µM range or 0.1 serum fraction are flagged as extrapolations. No chemical
sequence or motif lookup silently transfers these rates to a free vaccine
peptide. The NPY calculations and plasma boundary scenario remain separate.

## Compare enzyme removal

```python
from mhctools import DegradationTarget, EnzymeCutRate, enzyme_removal_effects

def fragment_rates(sequence):
    # Synthetic example; replace each state's rates with supported measurements
    # or explicitly documented assumptions. Never insert native model scores.
    if sequence == "ACDEFGHIK":
        return (
            EnzymeCutRate(2, "trimming route", 2.0, "synthetic assumption", "example"),
            EnzymeCutRate(4, "internal route", 0.5, "synthetic assumption", "example"),
        )
    if sequence == "DEFGHIK":
        return (EnzymeCutRate(3, "internal route", 0.2,
                              "synthetic assumption", "example"),)
    return None  # unavailable kinetics; do not claim survival

result = enzyme_removal_effects(
    "ACDEFGHIK", DegradationTarget("selected core", 2, 7), fragment_rates,
    enzymes=["trimming route", "internal route"], times_hours=[0.5, 1, 4],
    scenario="synthetic rate example", seed=42,
)
```

The baseline distinguishes full parent, target-bearing fragments, target
destruction, clearance, uptake and unknown outcomes. First-cut attribution
names the enzyme responsible for parent disappearance. Destructive-cut attribution
names the enzyme making the final cut; it does not identify every enzyme
required earlier in that path.

The removal effect is the change in target availability, in **percentage
points**, at each requested time. Removing an enzyme subtracts its hazards;
the other rates remain unchanged. Do not renormalize remaining rates to the
original total or add a second half-life-derived degradation hazard.

Effects can be negative: blocking a harmless trim can leave the precursor
exposed to a faster destructive route. Effects need not sum to 100% because
enzymes can compete or act successively. Compensation and inhibitor off-target
effects are outside this counterfactual. Common seeds make runs reproducible;
finite Monte Carlo sampling error remains separate from biological uncertainty.
`target_gain_mc_se_upper_bound_pp` reports `100 / sqrt(n_paths)` percentage
points: a worst-case bound obtained from each arm's maximum binomial variance,
valid under any coupling between arms. It stays positive even if a small sample
observes no events. Tiny effects comparable to that scale cannot establish a useful
enzyme ranking. This is not a biological uncertainty interval.

`None` from the rate provider means missing kinetics. An empty tuple means
explicitly assessed zero cleavage hazard in the supplied scenario. If unknown
states are reached, the result supplies lower/upper identification bounds and
withholds the point effect. These bounds cover missing states, not all
unmodeled biology, and are not confidence intervals.

## A broad terminal-trimming scenario

`empirical_terminal_cut_rates` uses the four immediate terminal boundary counts
as relative scenario weights. It does not assign individual enzyme identities.
More distant boundaries can arise through successive trimming, so they do not
become direct multi-residue jumps. Every fragment gets a new terminal assessment.

```python
import math
from mhctools import empirical_terminal_cut_rates

def fragment_rates(sequence):
    hours = half_lives_from_one_estimator.get(sequence)
    if hours is None or hours <= 0:
        return None
    return empirical_terminal_cut_rates(
        sequence, total_rate_per_hour=math.log(2) / hours,
        allow_transfer=True,  # explicit pig/14-mer-to-vaccine scenario
    )
```

This scenario is **terminal-only** and omits internal cleavage. Its unique
context counts are structural weights, not measured event frequencies.
Long-peptide sequence preferences, human transfer and fragment rates remain
unvalidated. It is unsuitable as a reassuring headline estimate of epitope
protection. The default requires 14-residue inputs; other lengths need explicit
transfer opt-in, and lengths below five remain unassessed. Run PeptiVerse and
Cavaco clocks separately, alongside other explicit cut-allocation assumptions.

## Reproduce the public example

```sh
python scripts/serum_contribution_example.py --output /tmp/serum-model
```

The example writes JSON and CSV for the NPY reference kinetics and synthetic
target trajectories. Its chart, when matplotlib is installed, shows parent
disappearance versus target survival in the terminal-only scenario. Synthetic
rates and transferred library weights are labelled as assumptions. To audit
the source counts, supply the downloaded Table S1:

```sh
python scripts/serum_contribution_example.py --output /tmp/serum-model \
  --library-table /path/to/Table_S1.xlsx
```

This requires pandas/openpyxl and verifies the published file hash and all
boundary counts. Analytic regression oracles check exponential survival,
serial cuts, nonadditive enzyme removal, protective trimming, competing removal,
missing fragment states and preservation of chemistry. These validate the
calculation, not vaccine biological accuracy. Further calibration is tracked
in [#541](https://github.com/openvax/mhctools/issues/541) and
[#291](https://github.com/openvax/mhctools/issues/291).
