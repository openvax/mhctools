"""Local target-centered evidence and explicit serum-degradation scenarios.

Source inputs remain local. No sequence is sent to a prediction service.
"""

from functools import lru_cache
import hashlib

from mhctools import DegradationTarget, predict_cleavage


# Exposure is assumed in these scenarios, rather than inferred from a motif.
# Activated CPB2 is deliberately omitted without an explicit activation state.
EXTRACELLULAR_MODELS = (
    "dpp4-qpisa", "ace-dipeptidyl", "cpn-basic", "app2-xp", "fap-endo-gp",
    "fap-dipeptidyl", "mme-hydrophobic", "anpep-ala", "enpep-acidic",
)


def positions(sequence, subsequence):
    return [i for i in range(len(sequence) - len(subsequence) + 1)
            if sequence.startswith(subsequence, i)]


def mutant_region(record, variant):
    """Map source-highlighted residues via exact unique sequence anchors.

    The source highlight is a source annotation, not an independent genomic
    validation. A source-core anchor can locate a point mutation despite altered
    vaccine flanks. No alignment, deletion-junction, or frameshift offsets are
    guessed when the exact source highlight is absent.
    """
    sequence = record["sequence"]
    context = variant.get("aa_sequence")
    start, end = variant.get("aa_hi_start"), variant.get("aa_hi_end")
    if not context or start is None or end is None:
        return None, "No source-highlighted mutant-region coordinates"
    if not 0 <= start < end <= len(context):
        raise ValueError("Invalid source mutant highlight")
    matches = positions(context, sequence)
    if len(matches) == 1:
        source_offset, vaccine_offset, length = matches[0], 0, len(sequence)
        anchor = "exact parent sequence"
    else:
        core = record.get("minimal_epitope")
        offset = record.get("minimal_epitope_offset")
        matches = positions(context, core) if core and offset not in (None, "") else []
        if len(matches) != 1:
            return None, "No unique exact parent/core anchor for source mutation"
        source_offset, vaccine_offset, length = matches[0], int(float(offset)), len(core)
        if sequence[vaccine_offset:vaccine_offset + length] != core:
            raise ValueError("Source-core sequence mismatch")
        anchor = "exact disclosed core; vaccine flanks not mapped"
    left, right = max(start, source_offset), min(end, source_offset + length)
    if left >= right:
        return None, "Source mutation does not overlap the exact mapped anchor"
    return (left - source_offset + vaccine_offset, right - source_offset + vaccine_offset), anchor


@lru_cache(maxsize=None)
def extracellular_evidence(sequence):
    """Reassess the current fragment; return native evidence without rates."""
    return predict_cleavage(sequence, models=EXTRACELLULAR_MODELS)


def recognition_bonds(sequence):
    """One flag per bond; overlapping rules never become independent votes."""
    bonds = set()
    for result in extracellular_evidence(sequence):
        for site in result.sites:
            if (site.status == "matched" or
                    result.model.name == "dpp4-qpisa" and site.status == "scored" and site.score > 0):
                bonds.add(site.bond)
    return bonds


def recognition_weights(sequence, multiplier):
    """Arbitrary sensitivity multiplier over a positive uniform background."""
    if multiplier < 1:
        raise ValueError("Recognition multiplier must be at least one")
    matched = recognition_bonds(sequence) if multiplier > 1 else set()
    return [float(multiplier if bond in matched else 1)
            for bond in range(1, len(sequence))]


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")


def target_inventory(records, variants, ligand_rows):
    """Source targets plus the best qualifying mutant ligand/core per class.

    Also retain ALL qualifying predicted candidates in an audit table. MHC-I
    candidates use the fresh presentation rank <= 2%; MHC-II candidates use
    rank <= 5% and require a uniquely mapped contiguous binding core. For a
    mutation in a class-II flank, retain the core AND the closest mapped mutant
    residue together: an intact wild-type core alone is not a mutant target.
    No claim about joint loss of all alternative candidates is inferred.
    """
    by_id = {record["sequence_record_id"]: record for record in records}
    mappings = {key: mutant_region(record, variants[record["variant_id"]])
                for key, record in by_id.items()}
    targets, candidates, coverage = [], [], []

    def add(record, label, start, end, kind, **details):
        target = DegradationTarget(label, start, end)
        target.validate(record["sequence"])
        region, reason = mappings[record["sequence_record_id"]]
        return dict(sequence_record_id=record["sequence_record_id"], gene=record["gene"],
                    target_label=label, start=start, end=end,
                    target_sequence=record["sequence"][start:end], kind=kind,
                    mutant_overlap=(max(start, region[0]) < min(end, region[1]) if region else None),
                    mutation_mapping=reason, **details)

    for record in records:
        key = record["sequence_record_id"]
        core, offset = record.get("minimal_epitope"), record.get("minimal_epitope_offset")
        has_core = bool(core and offset not in (None, ""))
        if has_core:
            start = int(float(offset))
            if record["sequence"][start:start + len(core)] != core:
                raise ValueError("Source-core mismatch")
            targets.append(add(record, "Source target", start, start + len(core), "source minimal epitope"))
        # A separate disclosed class-I target may be present even when the
        # variant's single minimal_epitope belongs to another vaccine construct.
        other = variants[record["variant_id"]].get("peptides", {}).get("CeGaT_Class_I")
        matches = positions(record["sequence"], other) if other else []
        if len(matches) == 1 and other != core:
            targets.append(add(record, "Source class-I candidate", matches[0], matches[0] + len(other),
                               "source class-I candidate"))
        region, reason = mappings[key]
        coverage.append(dict(sequence_record_id=key, gene=record["gene"],
                             source_core_available=has_core, mutant_region_start=region[0] if region else None,
                             mutant_region_end=region[1] if region else None, mutation_mapping=reason))

    for row in ligand_rows:
        if float(row["percentile_rank"]) > float(row["rank_threshold"]):
            continue
        record = by_id[row["sequence_record_id"]]
        start, end = int(row["start"]) - 1, int(row["end"])
        if record["sequence"][start:end] != row["peptide"]:
            raise ValueError("Predicted ligand coordinate mismatch")
        if row["mhc_class"] == "II":
            core = row["binding_core"]
            matches = positions(row["peptide"], core) if core else []
            if len(matches) != 1:
                candidates.append(dict(**row, target_assessment="Ambiguous or absent binding-core mapping"))
                continue
            start, end = start + matches[0], start + matches[0] + len(core)
        core_start, core_end = start, end
        region, _ = mappings[record["sequence_record_id"]]
        ligand_start, ligand_end = int(row["start"]) - 1, int(row["end"])
        mutant_ligand_overlap = (max(ligand_start, region[0]) < min(ligand_end, region[1])
                                 if region else None)
        mutant_core_overlap = (max(core_start, region[0]) < min(core_end, region[1])
                               if region else None)
        if row["mhc_class"] == "II" and mutant_ligand_overlap and not mutant_core_overlap:
            # Flanking residues can participate in class-II TCR recognition.
            # Preserve a contiguous core+mutation span, not just a WT core.
            if region[1] <= core_start:
                start = min(end, region[1]) - 1
            else:
                end = max(start, region[0]) + 1
        entry = add(record, "Predicted mutant " + row["mhc_class"], start, end,
                    "predicted class-I ligand" if row["mhc_class"] == "I" else
                    "predicted class-II binding core" if mutant_core_overlap else
                    "predicted class-II core + mutant flank" if mutant_ligand_overlap else
                    "predicted class-II binding core",
                    mhc_class=row["mhc_class"], allele=row["allele"], percentile_rank=float(row["percentile_rank"]),
                    full_ligand=row["peptide"], ligand_start=ligand_start, ligand_end=ligand_end,
                    binding_core_start=core_start if row["mhc_class"] == "II" else None,
                    binding_core_end=core_end if row["mhc_class"] == "II" else None,
                    binding_core=row.get("binding_core"), mutant_ligand_overlap=mutant_ligand_overlap,
                    mutant_binding_core_overlap=mutant_core_overlap if row["mhc_class"] == "II" else None,
                    predictor=row["predictor"], model_version=row["model_version"])
        entry["target_assessment"] = ("source-mapped mutant candidate" if entry["mutant_overlap"] else
                                      "mutation location unknown" if entry["mutant_overlap"] is None else
                                      "does not overlap mapped mutant region")
        candidates.append(entry)

    for key in by_id:
        for mhc_class in ("I", "II"):
            eligible = [r for r in candidates if r["sequence_record_id"] == key and
                        r.get("mhc_class") == mhc_class and r.get("mutant_overlap") is True]
            if eligible:
                best = min(eligible, key=lambda r: (r["percentile_rank"], r["start"], r["end"], r["allele"]))
                targets.append(best)
    return targets, candidates, coverage
