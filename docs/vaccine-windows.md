# Assess candidate vaccine windows

`mhctools vaccine-windows` combines dose conversion, core-preserving degradation,
APC delivery and presentation for every supplied window and core/allele. A free,
unmodified SLP defaults to **1 µg of that window**, using its average molecular
weight. RNA needs a supplied amount: a mass dose requires the molecular weight
of the **actual full transcript**, including its chemistry, not the encoded protein.

This is a conditional calculation. It does **not** infer human IM transport,
uptake, enzyme exposure or clearance rates from sequence. Missing inputs produce
`unassessed`, with the missing keys listed; they never become zero exposure.
Serum half-life predictions and protease motif scores are separate evidence and
are not silently converted into IM rate constants.

## Run a complete example

From a source checkout:

```sh
python examples/vaccine_window_scenario.py > windows.json
mhctools vaccine-windows --input windows.json --output exposure.json \
  --table exposure.csv --report exposure.md
```

All rates in this example are explicitly synthetic. At 48 hours the two windows
have approximately 0.50 and 0.33 nodal APC arrival fractions per core because
we **assigned** different competing clearance rates. Those values are software
examples, not biological predictions for those sequences. Both windows contain
two scored core occurrences, which remain separate in the output.

The JSON retains each complete trajectory and its kinetic sources. The CSV has
one row per core occurrence/allele. The Markdown report explains dose, missing
inputs, arrival and presentation endpoints. Existing files are never overwritten.

## Input contract

Use `schema_version: 1`, a nonempty `scenario`, `times_hours`, and `candidates`.
Each candidate specifies:

- `id`: unique window identifier.
- `peptide`: sequence or an exact [PeptideInput](exposure-results.md) record.
- `delivery`: `slp`, `secreted_mrna` or `cytosolic_mrna`.
- `cores`: every contributing scored target, with `label`, `start`, `end`,
  `allele`, `mhc_class` and optional `score_contribution`/`score_source`.
  Coordinates are zero-based, half-open on this candidate. Repeated occurrences
  retain their distinct coordinates; the same core on different alleles also
  remains separate. Duplicate occurrence/allele records are rejected.
- `kinetics`: profiles keyed by `I` or `II`. A profile supplies `rates`,
  `tap_max_length` or `class_ii_max_length`, and `cleavage_assessments` in the
  [trajectory input format](vaccine-trajectory.md). Candidate profiles replace
  the corresponding request-level class profile; rates are not partially merged.
- Optional `dose`, `signal_end`, `administration`, `organism`, and `metadata`.
  Secreted mRNA uses the full encoded construct and requires its signal boundary.

The core interval is the target to preserve, not necessarily a minimum MHC
binding motif. If only a longer scored ligand is known, label it as that ligand;
do not invent a shorter binding core. All longer fragments retaining the exact
target occurrence contribute to delivery and blood exposure. Cuts strictly
inside that interval destroy that target; flank cuts can preserve it.

Cleavage tables are specific to each actual sequence and current fragment.
Unlisted reachable fragments remain unassessed. An empty channel list means an
explicit, sourced no-cut assumption; it is not the default. Modified termini or
attachments are preserved but remain unassessed by the current trajectory engine.

## Read the endpoints

- **Nodal APC fraction of dosed peptide:** cumulative first presence of the
  core-bearing lineage in a draining-node APC, divided by starting SLP molecules.
  This includes local APC migration. It is not the residual peptide mass fraction.
- **Arrivals per input RNA:** antigen copies first appearing in nodal APCs per
  input transcript; translation can make this exceed one. The JSON also provides
  arrivals divided by antigen copies produced by the observation horizon.
- **Nodal pMHC copy-hours:** integrated surface presentation in nodal APCs.
  Delivery does not imply successful loading, presentation or T-cell priming.
- **Blood core copy-hours:** integrated exact core plus longer core-bearing
  fragments. Exact free core and longer precursor counts are also separate.
- **Post-peak half-decline bracket:** sampled time interval for the blood
  core-bearing pool to fall to half its observed maximum. Continued absorption
  or production prevents interpreting this as an intrinsic circulation half-life.

For nonnegative supplied scores, the report computes a score-weighted mean of
per-core arrivals per input molecule, **only when every core is assessed**. This
is a design descriptor. It does not identify the fraction of unique molecules
carrying *any* core: multiple cores can share one precursor, or a cut can separate
them into different daughters. Separate marginal runs do not model their joint
fate or MHC competition. No automatic winner or ranking change is imposed.

For IM mRNA mechanisms, primary evidence, cross-dressing and direct APC
expression, see [core delivery and blood exposure](vaccine-core-exposure.md).
For calibration requirements, see [comparing exposure](vaccine-comparison.md).

## Explicit sensitivity scenarios

To run complete illustrative SLP, cytosolic mRNA and secreted mRNA cases:

```sh
python examples/vaccine_sensitivity_scenarios.py --output-dir sensitivity
```

Or apply the same declared assumptions to an exported candidate request:

```sh
python examples/vaccine_sensitivity_scenarios.py \
  --input windows.json --output-dir sensitivity
```

The script writes three sets of JSON, CSV and Markdown results. It changes all
unproductive loss/clearance rates together (0.1, 1 and 10 per hour), while keeping
transport, explicit boundary-cut cascades and RNA production assumptions fixed.
This is a joint scenario comparison, not an estimated biological range or a
one-parameter causal attribution. The chosen sequences are not validated binders;
the example signal boundary and transcript mass are hypothetical. Supplying
an input replaces its kinetic profiles with the stated synthetic profiles, while
preserving its candidates, scored target occurrences and dose/chemistry records.
Missing RNA doses and unsupported chemistry remain unassessed.
