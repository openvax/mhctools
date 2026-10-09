# Vaccine target survival: three routes

**Track the exact target epitope and productive MHC loading.** The full parent
can disappear while a shorter fragment preserves the target. Processing is
also how the target is generated.

| Route | What needs to work | Most useful compact readout |
|---|---|---|
| **Expression from mRNA** | RNA delivery and translation; processing of the actual encoded construct, including junctions and routing | Expression yield + target presented on MHC. For cytosolic MHC-I: target-internal proteasome flags, C-terminal release and ER trimming evidence |
| **Peptide uptake by APCs** | Uptake into the relevant compartment; target release before destructive trimming; MHC loading | Target-internal versus boundary evidence for endolysosomal MHC-II processing, with MHC-I cross-presentation assessed separately |
| **Serum exposure** | Exact target persists in parent or fragments until uptake | Percentage retaining the target at the intended uptake time; parent half-life alongside it; measured, assumed or unknown kinetics |

A cut inside the selected target removes that exact sequence. A cut outside it
can retain the target and sometimes aid release. A boundary flag alone does not
establish successful liberation, and intact target retention does not establish
uptake or presentation. Alternative epitopes remain a separate question.
ERAP1 can generate or destroy epitopes, and human dendritic-cell experiments
demonstrate a cytosolic route for some SLPs. These routes depend on the construct,
cell and exposure. [ERAP1 experiments](https://doi.org/10.1038/ni860),
[human DC cross-presentation](https://doi.org/10.1371/journal.pone.0089897).

## One line per target

**Gene + mutation | exact target + HLA | route | internal-cut concern |
release evidence | target survival or unknown | evidence scope**

`compact_target_summary(construct)` and `target-summary.json` in
`mhctools vaccine-report` retain site flags, missing bonds and conditional
routes. Expression yield, uptake, MHC loading and serum target lifetime remain
unassessed when the input has no corresponding evidence. Predicted native MHC
ranks remain model-specific evidence.

For explicit local SLP/mRNA rate scenarios through loading and surface display,
see [APC pMHC trajectories](vaccine-trajectory.md). Missing trajectory inputs
remain unknown; the compact report does not infer rates from site flags.

The compact summary has no combined protection probability. "No flags" means
no flags in available assessments, with coverage shown; it does not mean
resistance. MHC-I ligand C-terminal predictors and named-enzyme digestion
models retain their different endpoints.

PeptiVerse and Cavaco parent estimates can be displayed side by side. Target
survival requires fragment-specific kinetics and a cut-path model; see
[target degradation](target-degradation.md). Serum digestion and patient
circulation PK remain distinct endpoints. Tissue access, binding and clearance
need additional evidence.

## Reproduce the compact guide and calibration audit

```sh
python scripts/serum_calibration_review.py --output /tmp/vaccine-route-review --pdf
```

The optional `--sid-input /path/to/focused_evidence.json` exports a compact
target table from a checksummed existing report without new predictor inference.
The [human perturbation audit](serum-calibration.md) documents which rates are
identifiable and the measurements still needed for vaccine-specific calibration.
