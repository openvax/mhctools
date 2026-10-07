# Target-preserving degradation: implementation spec

Build an explicitly exploratory, local peptide-fragment simulation. The primary
endpoint is retention of an exact, coordinate-defined target epitope, rather than
retention of the full vaccine peptide.

* Sample exponentially distributed waiting times using one half-life estimator
  at a time. Recompute half-life and extracellular recognition on every retained
  fragment; cache exact sequences without discarding occurrence coordinates.
* A flank/boundary cut retains the target-containing product. A cut strictly
  inside the target destroys that exact epitope. Stop when the target is lost,
  removed from circulation, taken up, or the observation horizon is reached.
* Cut-location weights, exponential kinetics, and enzyme exposure are explicit
  scenario assumptions. Native recognition scores are never treated as serum
  rates. Uniform/background cuts cover missing enzymes and coefficients.
* Keep PeptiVerse and Cavaco runs separate. Show their fragment estimates and
  trajectory curves side by side; scenario ranges are not confidence intervals.
* Allow separately supplied clearance and uptake rates, with competing hazards
  and distinct outcomes. Serum prediction alone does not specify either rate.
  Proteasome/cathepsin/CPP outputs do not determine extracellular event rates.
* Unsupported fragment estimates end in an explicit unknown outcome. Missing
  targets and ambiguous mutation mappings remain unknown. For class II, assess
  the mapped binding core separately from the longer predicted ligand. When a
  source-mapped mutation lies in its flank, track a core+mutation span so an
  intact wild-type binding core is not counted as mutant-target retention.
* Validate coordinate boundaries, sequential trimming/re-recognition, known
  kinetic limits, censoring, competing hazards, reproducibility, and independent
  estimator results. Use synthetic fixtures in the public repository. Keep all
  patient-linked outputs in ignored local results, preserving the frozen report.

Deliver the reusable simulation, documentation and tests through a versioned PR;
then run it locally on the Sid source targets and supported predicted mutant
ligands, with a readable comparison report and complete scenario provenance.
