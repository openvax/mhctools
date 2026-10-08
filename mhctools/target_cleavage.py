"""Relate native cleavage evidence to one exact target, without inventing rates."""

from .serum_degradation import DegradationTarget
from numbers import Integral


def annotate_target_cleavage(target, results):
    """Annotate native sites as target splitting, boundary release or trimming.

    Parameters
    ----------
    target : DegradationTarget
        Exact target interval in original-source coordinates. For reassessed
        fragments, ``CleavageInput.source_start`` locates their new termini.
    results : iterable of CleavageResult
        Native evidence for the same exact peptide, chemistry and occurrence.
        No predictor is run and no model score threshold is applied here.

    Returns
    -------
    tuple of dict
        One row per assessed bond, or a bond-less row for an unavailable or
        unresolved model. Every score keeps its endpoint, units and assay.
        ``target_effect`` describes what that cut WOULD do, not whether it
        occurs. Non-matches and unavailable results never establish protection.
        Compartments are retained, so callers can separate extracellular
        digestion from processing after uptake. This function never supplies
        cut probabilities or combines enzyme scores.
    """
    results = tuple(results)
    if not results:
        raise ValueError("At least one cleavage result is required")
    peptide = results[0].peptide
    offset = peptide.source_start
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in (target.start, target.end)):
        raise ValueError("Target requires integer coordinates")
    local_target = DegradationTarget(target.label, target.start - offset, target.end - offset)
    local_target.validate(peptide.sequence)
    rows = []
    for result in results:
        if result.peptide != peptide:
            raise ValueError("Do not mix peptide sequences, chemistry or source occurrences")
        model = result.model
        common = dict(
            model=model.name, model_version=model.version, enzyme=model.enzyme,
            evidence=model.evidence, motif_strictness=model.motif_strictness,
            scored_endpoint=model.scored_endpoint, score_name=model.score_name,
            score_units=model.score_units, compartments=model.compartments,
            assay=model.assay, limitations=model.limitations,
            references=model.references, conditions=result.conditions,
            substrate_observation=result.substrate_observation)
        if not result.sites:
            rows.append(dict(
                **common, bond=None, source_bond=None, bond_label=None,
                target_effect=None, score=None,
                status="unavailable" if result.unsupported_reason else "no_resolved_site",
                reason=result.unsupported_reason or "No site resolved; no target protection inferred"))
        for site in result.sites:
            bond = offset + site.bond
            effect = ("target_split" if target.start < bond < target.end else
                      "boundary_release" if bond in (target.start, target.end) else
                      "flank_trim")
            label = "%s%d|%s%d" % (
                peptide.sequence[site.bond - 1], bond,
                peptide.sequence[site.bond], bond + 1)
            rows.append(dict(**common, bond=site.bond, source_bond=bond,
                             bond_label=label, target_effect=effect,
                             status=site.status, score=site.score, reason=site.reason))
    return tuple(rows)
