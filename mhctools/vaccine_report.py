"""Structured, route-aware vaccine peptide reports.

This module separates observed/model-native evidence from the biological
route policy used to display it. It never combines unlike model scores into a
single degradation, preservation, presentation, or placement probability.
"""

from dataclasses import asdict, dataclass
from datetime import datetime
import csv
import hashlib
import json
import math
from pathlib import Path

from .cleavage import AMINO_ACIDS


DELIVERY_MODALITIES = ("synthetic_long_peptide", "rna_encoded")
RNA_ROUTINGS = ("cytosolic", "secreted", "lysosomal_targeted")
TRACK_CONTEXTS = (
    "cytosolic_proteasome",
    "er_trimming",
    "endolysosomal",
    "extracellular_interstitial",
    "tumor_stroma",
    "circulation",
)
ROUTE_RELEVANCE = ("primary", "conditional", "not_applicable")


@dataclass(frozen=True)
class EpitopeWindow:
    """One-based inclusive intended target window."""

    label: str
    start: int
    end: int

    @classmethod
    def from_dict(cls, value):
        return cls(str(value["label"]), int(value["start"]), int(value["end"]))


@dataclass(frozen=True)
class MHCWindow:
    """One candidate ligand window with one predictor's native rank."""

    mhc_class: str
    allele: str
    start: int
    end: int
    percentile_rank: float
    predictor: str

    @classmethod
    def from_dict(cls, value):
        return cls(
            str(value["mhc_class"]),
            str(value["allele"]),
            int(value["start"]),
            int(value["end"]),
            float(value["percentile_rank"]),
            str(value["predictor"]),
        )


@dataclass(frozen=True)
class CleavageTrack:
    """One model-native per-bond track; ``None`` means unassessed."""

    name: str
    context: str
    evidence_type: str
    scores: tuple
    threshold: float = None
    units: str = "native model output"
    provenance: str = ""

    @classmethod
    def from_dict(cls, value):
        scores = tuple(
            None if item is None else float(item) for item in value["scores"]
        )
        threshold = value.get("threshold")
        return cls(
            name=str(value["name"]),
            context=str(value["context"]),
            evidence_type=str(value["evidence_type"]),
            scores=scores,
            threshold=None if threshold is None else float(threshold),
            units=str(value.get("units", "native model output")),
            provenance=str(value.get("provenance", "")),
        )


@dataclass(frozen=True)
class VaccineConstruct:
    """A delivered peptide or translated antigen construct."""

    identifier: str
    sequence: str
    delivery: str
    routing: str = None
    exposures: tuple = ()
    intended_epitopes: tuple = ()
    mhc_windows: tuple = ()
    cleavage_tracks: tuple = ()

    @classmethod
    def from_dict(cls, value):
        delivery = str(value["delivery"])
        routing = value.get("routing")
        if delivery == "rna_encoded" and routing is None:
            routing = "cytosolic"
        construct = cls(
            identifier=str(value["id"]),
            sequence=str(value["sequence"]),
            delivery=delivery,
            routing=None if routing is None else str(routing),
            exposures=tuple(str(item) for item in value.get("exposures", ())),
            intended_epitopes=tuple(
                EpitopeWindow.from_dict(item)
                for item in value.get("intended_epitopes", ())
            ),
            mhc_windows=tuple(
                MHCWindow.from_dict(item) for item in value.get("mhc_windows", ())
            ),
            cleavage_tracks=tuple(
                CleavageTrack.from_dict(item)
                for item in value.get("cleavage_tracks", ())
            ),
        )
        construct.validate()
        return construct

    def validate(self):
        if not self.identifier:
            raise ValueError("Construct id must not be empty")
        if (not self.sequence or self.sequence != self.sequence.upper() or
                not set(self.sequence) <= AMINO_ACIDS):
            raise ValueError(
                "%s must have a nonempty canonical uppercase peptide sequence"
                % self.identifier)
        if self.delivery not in DELIVERY_MODALITIES:
            raise ValueError(
                "Unknown delivery %r; choices: %s"
                % (self.delivery, ", ".join(DELIVERY_MODALITIES)))
        if self.delivery == "rna_encoded":
            if self.routing not in RNA_ROUTINGS:
                raise ValueError(
                    "RNA construct routing must be one of %s"
                    % ", ".join(RNA_ROUTINGS))
        elif self.routing is not None:
            raise ValueError("routing is only valid for rna_encoded constructs")
        length = len(self.sequence)
        for window in self.intended_epitopes:
            if not 1 <= window.start <= window.end <= length:
                raise ValueError("Invalid intended epitope interval in %s" % self.identifier)
        for window in self.mhc_windows:
            if window.mhc_class not in ("I", "II"):
                raise ValueError("MHC class must be I or II")
            if not 1 <= window.start <= window.end <= length:
                raise ValueError("Invalid MHC interval in %s" % self.identifier)
            if not math.isfinite(window.percentile_rank) or window.percentile_rank < 0:
                raise ValueError("MHC percentile rank must be finite and nonnegative")
        names = [track.name for track in self.cleavage_tracks]
        if len(names) != len(set(names)):
            raise ValueError("Cleavage track names must be unique within a construct")
        for track in self.cleavage_tracks:
            if track.context not in TRACK_CONTEXTS:
                raise ValueError("Unknown cleavage context %r" % track.context)
            if len(track.scores) != max(0, length - 1):
                raise ValueError(
                    "%s track %s needs exactly %d internal-bond values"
                    % (self.identifier, track.name, max(0, length - 1)))
            if track.threshold is not None and not math.isfinite(track.threshold):
                raise ValueError("Track thresholds must be finite")
            if any(value is not None and not math.isfinite(value)
                   for value in track.scores):
                raise ValueError("Track scores must be finite numbers or null")


@dataclass(frozen=True)
class VaccineReportInput:
    """Validated schema for the report CLI and Python API."""

    constructs: tuple
    schema_version: int = 1

    @classmethod
    def from_dict(cls, value):
        if value.get("schema_version") != 1:
            raise ValueError("vaccine report schema_version must be 1")
        constructs = tuple(
            VaccineConstruct.from_dict(item) for item in value.get("constructs", ())
        )
        if not constructs:
            raise ValueError("At least one vaccine construct is required")
        identifiers = [item.identifier for item in constructs]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("Vaccine construct ids must be unique")
        return cls(constructs=constructs)

    @classmethod
    def from_json(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def to_dict(self):
        return {
            "schema_version": self.schema_version,
            "constructs": [
                {
                    "id": construct.identifier,
                    "sequence": construct.sequence,
                    "delivery": construct.delivery,
                    "routing": construct.routing,
                    "exposures": list(construct.exposures),
                    "intended_epitopes": [asdict(item) for item in construct.intended_epitopes],
                    "mhc_windows": [asdict(item) for item in construct.mhc_windows],
                    "cleavage_tracks": [asdict(item) for item in construct.cleavage_tracks],
                }
                for construct in self.constructs
            ],
        }


@dataclass(frozen=True)
class RoutePolicy:
    context: str
    relevance: str
    rationale: str


def processing_route_policy(construct):
    """Return declared display relevance for every supported route context."""
    if construct.delivery == "synthetic_long_peptide":
        circulation = "circulation" in construct.exposures
        tumor_stroma = "tumor_stroma" in construct.exposures
        values = {
            "extracellular_interstitial": (
                "primary",
                "Injected free peptide is exposed before and during uptake.",
            ),
            "endolysosomal": (
                "primary",
                "Internalized SLP can be processed for class-II presentation.",
            ),
            "cytosolic_proteasome": (
                "conditional",
                "DC cytosolic export can support proteasome/TAP-dependent cross-presentation.",
            ),
            "er_trimming": (
                "conditional",
                "Relevant only after cross-presentation has produced an ER-accessible precursor.",
            ),
            "tumor_stroma": (
                "conditional" if tumor_stroma else "not_applicable",
                "Requires declared tumor/stromal exposure.",
            ),
            "circulation": (
                "conditional" if circulation else "not_applicable",
                "Requires declared blood/plasma exposure; injection alone does not establish it.",
            ),
        }
    elif construct.routing == "cytosolic":
        values = {
            "cytosolic_proteasome": (
                "primary", "Translated cytosolic antigen enters endogenous class-I processing."
            ),
            "er_trimming": (
                "primary", "ER trimming can act on transported class-I precursors."
            ),
            "endolysosomal": (
                "conditional", "Autophagy can route endogenous antigen toward class-II loading."
            ),
            "extracellular_interstitial": (
                "not_applicable", "No free extracellular peptide is declared."
            ),
            "tumor_stroma": (
                "not_applicable", "No free extracellular peptide is declared."
            ),
            "circulation": (
                "not_applicable", "No free circulating peptide is declared."
            ),
        }
    elif construct.routing == "secreted":
        values = {
            "cytosolic_proteasome": (
                "conditional", "Translation-associated or mislocalized products can enter class-I processing."
            ),
            "er_trimming": (
                "conditional", "Relevant only to ER-accessible class-I precursors."
            ),
            "endolysosomal": (
                "conditional", "Secreted antigen may be recaptured into endolysosomal processing."
            ),
            "extracellular_interstitial": (
                "primary", "The construct explicitly declares secretion."
            ),
            "tumor_stroma": (
                "conditional", "Depends on expression site and stromal exposure."
            ),
            "circulation": (
                "conditional", "Secretion does not by itself establish systemic exposure."
            ),
        }
    else:
        values = {
            "cytosolic_proteasome": (
                "conditional", "Translation still occurs before lysosomal targeting."
            ),
            "er_trimming": (
                "conditional", "Relevant only to class-I precursors escaping the targeted route."
            ),
            "endolysosomal": (
                "primary", "The encoded antigen is explicitly targeted to lysosomal processing."
            ),
            "extracellular_interstitial": (
                "not_applicable", "Lysosomal targeting does not declare peptide secretion."
            ),
            "tumor_stroma": (
                "not_applicable", "No free extracellular peptide is declared."
            ),
            "circulation": (
                "not_applicable", "No free circulating peptide is declared."
            ),
        }
    return tuple(RoutePolicy(context, *values[context]) for context in TRACK_CONTEXTS)


def _overlaps_intended(window, construct):
    return any(
        window.start <= epitope.end and epitope.start <= window.end
        for epitope in construct.intended_epitopes
    )


def select_mhc_windows(construct, mhc_class, maximum=10):
    """Pick at most ``maximum`` windows with rank, target, and allele diversity."""
    if maximum < 1:
        raise ValueError("maximum must be positive")
    ordered = sorted(
        (item for item in construct.mhc_windows if item.mhc_class == mhc_class),
        key=lambda item: (item.percentile_rank, item.allele, item.start, item.end),
    )
    selected = []
    seen = set()

    def add(item):
        span = (item.start, item.end)
        if span in seen or len(selected) >= maximum:
            return
        selected.append(item)
        seen.add(span)

    intended = [item for item in ordered if _overlaps_intended(item, construct)]
    if intended:
        add(intended[0])
    for allele in dict.fromkeys(item.allele for item in ordered):
        for item in ordered:
            if item.allele == allele:
                add(item)
                break
    for item in ordered:
        add(item)
    return tuple(selected)


def placement_assessments(construct):
    """Report boundary and internal evidence separately for each target/track."""
    policy = {item.context: item for item in processing_route_policy(construct)}
    rows = []
    for epitope in construct.intended_epitopes:
        for track in construct.cleavage_tracks:
            internal_bonds = range(epitope.start, epitope.end)
            boundary_bonds = tuple(
                bond for bond in (epitope.start - 1, epitope.end)
                if 1 <= bond < len(construct.sequence)
            )

            def observed(bonds):
                return [
                    (bond, track.scores[bond - 1])
                    for bond in bonds
                    if track.scores[bond - 1] is not None
                ]

            def supported(values):
                if track.threshold is None:
                    return []
                return [bond for bond, score in values if score >= track.threshold]

            internal = observed(internal_bonds)
            boundary = observed(boundary_bonds)
            rows.append({
                "construct_id": construct.identifier,
                "delivery": construct.delivery,
                "routing": construct.routing,
                "epitope": epitope.label,
                "epitope_start": epitope.start,
                "epitope_end": epitope.end,
                "track": track.name,
                "context": track.context,
                "route_relevance": policy[track.context].relevance,
                "threshold": track.threshold,
                "internal_assessed_bonds": [bond for bond, score in internal],
                "internal_supported_bonds": supported(internal),
                "boundary_assessed_bonds": [bond for bond, score in boundary],
                "boundary_supported_bonds": supported(boundary),
                "score_units": track.units,
            })
    return rows


def _timestamped_output_dir(base_dir, generated_at):
    if generated_at.tzinfo is None or generated_at.utcoffset() is None:
        raise ValueError("generated_at must include a timezone")
    return Path(base_dir) / generated_at.strftime("%Y-%m-%dT%H%M%S-%f%z")


def _lane_assign(windows):
    """Pack selected windows into non-overlapping lanes without dropping any.

    ``windows`` arrives in selection order, and the returned rank is that
    order's one-based index, so a figure label matches ``display_rank`` in
    ``selected-mhc-windows.csv`` for the same window. Lanes are added as
    needed: a window that was selected is always drawn, because a figure that
    quietly omitted one would disagree with its own audit table.
    """
    lane_ends = []
    assigned = []
    ordered = sorted(
        enumerate(windows, start=1),
        key=lambda item: (item[1].start, item[1].end, item[1].percentile_rank),
    )
    for rank, window in ordered:
        lane = next(
            (index for index, end in enumerate(lane_ends) if window.start > end),
            None,
        )
        if lane is None:
            lane_ends.append(window.end)
            lane = len(lane_ends) - 1
        else:
            lane_ends[lane] = window.end
        assigned.append((rank, window, lane))
    return assigned


def _render_report_pdf(report, path, maximum_mhc_windows):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.patches import FancyBboxPatch, Rectangle
    except ImportError as error:
        raise ImportError(
            "vaccine report PDFs require matplotlib; install mhctools[vaccine-report]"
        ) from error

    colors = {"I": "#1769aa", "II": "#6a3d9a"}
    context_colors = {
        "cytosolic_proteasome": "#009e73",
        "er_trimming": "#7b4ab5",
        "endolysosomal": "#5b3d91",
        "extracellular_interstitial": "#007b83",
        "tumor_stroma": "#b45f06",
        "circulation": "#8a5a00",
    }
    with PdfPages(path) as pages:
        for construct in report.constructs:
            policy = {item.context: item for item in processing_route_policy(construct)}
            displayed_tracks = [
                track for track in construct.cleavage_tracks
                if policy[track.context].relevance != "not_applicable"
            ]
            upper = [
                track for track in displayed_tracks
                if track.context in ("cytosolic_proteasome", "er_trimming")
            ]
            lower = [track for track in displayed_tracks if track not in upper]
            length = len(construct.sequence)
            # Lane count follows the selection, so the ligand band grows
            # instead of discarding windows the audit table still lists.
            lane_layout = {
                mhc_class: _lane_assign(select_mhc_windows(
                    construct, mhc_class, maximum=maximum_mhc_windows))
                for mhc_class in ("I", "II")
            }
            header_y = {
                mhc_class: 1.02 + 0.42 * (
                    1 + max((lane for _, _, lane in assigned), default=0)
                )
                for mhc_class, assigned in lane_layout.items()
            }
            upper_base = header_y["I"] + 0.78
            lower_base = -(header_y["II"] + 0.78)
            # Margin past the outermost track line, whose curve spans 0.48. A
            # side with no route-relevant track stops just past its heading
            # rather than reserving an empty band.
            ymax = (upper_base + (len(upper) - 1) * 0.72 + 0.75 if upper
                    else header_y["I"] + 0.55)
            ymin = (lower_base - (len(lower) - 1) * 0.72 - 0.75 if lower
                    else -(header_y["II"] + 0.55))
            # Page height follows the drawn extent at a fixed inches-per-unit
            # scale, so a sparse construct does not print a blank half page and
            # a dense one keeps the same band spacing.
            axes_fraction = 0.9
            figure_height = max(6.0, (ymax - ymin) * 1.15 / axes_fraction)
            fig = plt.figure(figsize=(16, figure_height))
            # A near-full-height axes keeps the band spacing that was chosen
            # in data units from being padded by default subplot margins.
            ax = fig.add_axes([0.12, 0.048, 0.86, axes_fraction])
            ax.set_xlim(0.1, length + 0.9)
            ax.set_ylim(ymin, ymax)
            ax.axis("off")

            ax.add_patch(Rectangle((0.53, -0.47), length - 0.06, 0.94,
                                   facecolor="#f3f5f7", edgecolor="none"))
            for epitope in construct.intended_epitopes:
                ax.add_patch(Rectangle(
                    (epitope.start - 0.46, -0.56),
                    epitope.end - epitope.start + 0.92,
                    1.12,
                    facecolor="#fff0b3", edgecolor="#c78c00", linewidth=1.7,
                ))
            for index, residue in enumerate(construct.sequence, start=1):
                ax.text(index, 0, residue, ha="center", va="center",
                        family="monospace", fontsize=21, fontweight="bold",
                        color="#18222d")
                if index in (1, length) or index % 5 == 0:
                    ax.text(index, -0.68, str(index), ha="center", va="top",
                            fontsize=8.5, color="#56616c")

            for mhc_class, sign in (("I", 1), ("II", -1)):
                assigned = lane_layout[mhc_class]
                for rank, window, lane in assigned:
                    y = sign * (1.02 + lane * 0.42)
                    patch = FancyBboxPatch(
                        (window.start - 0.43, y - 0.14),
                        window.end - window.start + 0.86,
                        0.28,
                        boxstyle="round,pad=0.02,rounding_size=0.08",
                        facecolor="white", edgecolor=colors[mhc_class], linewidth=1.2,
                    )
                    ax.add_patch(patch)
                    label = "%s%d %s %.2g%%" % (
                        mhc_class, rank, window.allele.replace("HLA-", ""),
                        window.percentile_rank,
                    )
                    ax.text((window.start + window.end) / 2, y, label,
                            ha="center", va="center", fontsize=7.5,
                            fontweight="bold", color=colors[mhc_class], clip_on=True)
                ax.text(0.35, sign * header_y[mhc_class],
                        "MHC-%s: %d drawn of %d eligible; bar number is the "
                        "selection rank in selected-mhc-windows.csv" % (
                            mhc_class, len(assigned),
                            sum(1 for item in construct.mhc_windows
                                if item.mhc_class == mhc_class)),
                        ha="left", va="center", fontsize=9,
                        fontweight="bold", color=colors[mhc_class])

            for direction, tracks, base in ((1, upper, upper_base), (-1, lower, lower_base)):
                for track_index, track in enumerate(tracks):
                    y = base + direction * track_index * 0.72
                    color = context_colors[track.context]
                    ax.plot([0.55, length + 0.45], [y, y], color="#c9d0d6", linewidth=0.7)
                    assessed = [
                        (bond + 0.5, score) for bond, score in enumerate(track.scores, start=1)
                        if score is not None
                    ]
                    if assessed:
                        values = [score for x, score in assessed]
                        low, high = min(values), max(values)
                        span = high - low or 1.0
                        xs = [x for x, score in assessed]
                        ys = [y + direction * 0.48 * ((score - low) / span)
                              for x, score in assessed]
                        ax.plot(xs, ys, color=color, linewidth=1.8)
                        ax.scatter(xs, ys, color=color, s=11)
                        if track.threshold is not None:
                            threshold_y = y + direction * 0.48 * (
                                (track.threshold - low) / span
                            )
                            ax.plot(
                                [0.55, length + 0.45],
                                [threshold_y, threshold_y],
                                color=color, linewidth=0.65, alpha=0.42,
                                linestyle="--",
                            )
                            hit_x = [
                                x for x, score in assessed if score >= track.threshold
                            ]
                            hit_y = [
                                y + direction * 0.48 * ((score - low) / span)
                                for x, score in assessed if score >= track.threshold
                            ]
                            ax.scatter(hit_x, hit_y, color=color, s=30,
                                       edgecolor="white", linewidth=0.5)
                    ax.text(0.35, y, "%s [%s]" % (
                        track.name, policy[track.context].relevance),
                        ha="right", va="center", fontsize=9, color=color)

            route = construct.delivery
            if construct.routing:
                route += " / " + construct.routing
            route = route.replace("_", " ")
            fig.suptitle("%s | %s | %d aa" % (
                construct.identifier, route, length),
                fontsize=16, fontweight="bold", y=0.965)
            fig.text(
                0.5, 0.025,
                "Each enzyme/model has its own native-scale track. Filled points meet that "
                "track's declared threshold; they are not comparable probabilities.",
                ha="center", va="bottom", fontsize=8.5, color="#4e5963",
            )
            pages.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def generate_vaccine_report(
        report, output_base, generated_at=None, maximum_mhc_windows=10):
    """Generate a timestamped PDF and machine-readable audit artifacts."""
    if isinstance(report, dict):
        report = VaccineReportInput.from_dict(report)
    if not isinstance(report, VaccineReportInput):
        raise TypeError("report must be VaccineReportInput or a schema dictionary")
    if generated_at is None:
        generated_at = datetime.now().astimezone()
    output_dir = _timestamped_output_dir(output_base, generated_at)
    output_dir.mkdir(parents=True, exist_ok=False)

    normalized_path = output_dir / "normalized-input.json"
    normalized_path.write_text(
        json.dumps(report.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    route_rows = []
    placement_rows = []
    selected_rows = []
    for construct in report.constructs:
        for route in processing_route_policy(construct):
            route_rows.append({
                "construct_id": construct.identifier,
                "delivery": construct.delivery,
                "routing": construct.routing or "",
                **asdict(route),
            })
        placement_rows.extend(placement_assessments(construct))
        for mhc_class in ("I", "II"):
            selected = select_mhc_windows(
                construct, mhc_class, maximum=maximum_mhc_windows
            )
            for display_rank, window in enumerate(selected, start=1):
                selected_rows.append({
                    "construct_id": construct.identifier,
                    "display_rank": display_rank,
                    **asdict(window),
                })

    def write_rows(path, rows, columns):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write_rows(
        output_dir / "processing-routes.csv", route_rows,
        ["construct_id", "delivery", "routing", "context", "relevance", "rationale"],
    )
    (output_dir / "placement-assessments.json").write_text(
        json.dumps(placement_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_rows(
        output_dir / "selected-mhc-windows.csv", selected_rows,
        ["construct_id", "display_rank", "mhc_class", "allele", "start", "end",
         "percentile_rank", "predictor"],
    )
    pdf_path = output_dir / "vaccine-processing-report.pdf"
    _render_report_pdf(report, pdf_path, maximum_mhc_windows)
    (output_dir / "README.md").write_text(
        "# Vaccine processing report\n\n"
        "This route-aware report preserves model-native scores. It does not calculate "
        "an ensemble probability or claim that a predicted cut occurs in vivo. RNA and "
        "synthetic-long-peptide delivery use distinct declared processing policies.\n",
        encoding="utf-8",
    )
    checksums = {}
    for artifact in sorted(output_dir.iterdir()):
        if artifact.is_file() and artifact.name != "SHA256SUMS.json":
            checksums[artifact.name] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    (output_dir / "SHA256SUMS.json").write_text(
        json.dumps(checksums, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return output_dir
