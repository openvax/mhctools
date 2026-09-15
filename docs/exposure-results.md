# Peptide PK, uptake, and tissue-exposure results

`mhctools.pred` keeps pharmacokinetic and delivery endpoints distinct. It does
not expose a generic `delivery_score`, because systemic half-life, clearance,
exposure, CPP class confidence, quantitative cellular uptake, and tissue
concentration are not interchangeable measurements.

## Endpoint identity

The endpoint kinds are:

| Kind | Meaning |
|---|---|
| `systemic_elimination_half_life` | In-vivo terminal elimination half-life |
| `systemic_clearance` | Systemic or apparent clearance |
| `distribution_volume` | Systemic or apparent distribution volume |
| `systemic_exposure` | A study-defined systemic exposure quantity, such as AUC |
| `cpp_classification` | CPP class label and model confidence |
| `cellular_uptake` | Quantitative cellular uptake measurement or estimate |
| `tissue_concentration` | Concentration in a named tissue or compartment |

These remain separate from ex-vivo `serum_half_life`, `plasma_half_life`, and
`blood_half_life`, and from `pMHC_stability`, which describes dissociation of a
peptide-MHC complex. New plasma-assay results use `plasma_half_life`; they must
not be labeled as serum stability or in-vivo elimination.

## Measurement context

Every prediction of one of these endpoints requires an immutable
`MeasurementContext`. Its versioned fields preserve:

- whether the result was observed, fitted, simulated, or ML-predicted;
- whether it is available, unsupported, missing, out of domain, or failed;
- analyte, compartment, assay matrix, unit, and linear value transform;
- total versus unbound concentration and systemic versus apparent PK scope;
- class label and the meaning of a predictor-native score; and
- a time-series identifier, timepoint, unit, and time origin.

Unknown fields remain `None`; consumers must not replace them with biological
defaults. A stored `Prediction.value` is a linear physical value with an
explicit unit. A CPP class confidence belongs in `score`, alongside
`score_semantics` and `class_label`, and is not a duration, uptake amount, or
percentage delivered.

Unavailable results carry no stale numeric output. `status="unsupported"`,
`"missing"`, `"out_of_domain"`, or `"failed"` records the distinction and
`detail` may explain it.

Concentration-time points share a `series_id` and have explicit `timepoint`,
`time_unit`, and `time_origin`. These endpoints are MHC-independent, so their
predictions have an empty `allele`; they must not be copied once per HLA allele.

## Ordering and downstream projection

`best_direction()` and `PeptideResult.best_by*()` intentionally reject these
context-dependent kinds. Greater organ accumulation, longer circulation, or a
higher predictor-native output is not universally preferable. A downstream
Topiary or Vaxrank integration may project explicitly selected fields and
context, but it owns any use-specific ranking policy and clinical
interpretation. `mhctools` only transports endpoint data and provenance.

The contract follows the distinctions used in the systemic peptide PK review
by [Nordell et al.](https://doi.org/10.1007/s40262-025-01615-z) and the
cell-line-, cargo-, and assay-dependent quantitative uptake data in
[POSEIDON](https://doi.org/10.1186/s13321-024-00810-7).
