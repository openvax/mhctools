#!/usr/bin/env python3
"""Render mechanism-separated site scores from a verified frozen Sid run.

No prediction or network requests occur. CSV scores stay in native units;
DPP4 display uses predicted assay-loss percentages. Categorical evidence
stays categorical, and missing sites remain unassessed. Coordinates refer
to the intact disclosed vaccine input.
"""

import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
from xml.sax.saxutils import escape

from display_labels import dpp4_loss_label
from mhctools.eramer_cleavage import ERAMERCleavage
from mhctools.peptidases import cleavage_models, get_cleavage_model


PDF_NAME = "cleavage-mechanisms.pdf"
MAX_BONDS = 24
HUMAN = "pepsickle-in-vivo-human-only"
MAMMAL = "pepsickle-in-vivo-all-mammal"
ORGANISMS = ("human", "nonhuman", "mixed", "unknown")
INTERNAL_MOTIFS = ("prep-pro", "mme-hydrophobic", "fap-endo-gp")
PROCESSING = {
    HUMAN: (
        "Pepsickle human", "MHC-I processing proxy", "proteasome-associated; type agnostic",
        "Human epitope-trained neural ensemble; no constitutive/immunoproteasome assignment."),
    MAMMAL: (
        "Pepsickle all-mammal", "MHC-I processing proxy", "proteasome-associated; type agnostic",
        "Mammalian epitope-trained neural ensemble; no constitutive/immunoproteasome assignment."),
    "netchop-3.1-cterm-3.0": (
        "NetChop Cterm 3.0", "MHC-I processing proxy", "no individual enzyme assigned",
        "Trained on MHC-I ligand C-termini; processing/selection signal, not a purified-enzyme assay."),
    "netchop-3.1-20s-3.0": (
        "NetChop 20S 3.0", "Proteasome digestion", "20S proteasome; no catalytic subunit assigned",
        "Trained on in-vitro proteasome degradation experiments."),
    "netcleave-i-hla": (
        "NetCleave-I HLA", "MHC-I processing proxy", "no individual enzyme assigned",
        "MHC-I ligand-derived C-terminal processing; 8 upstream plus 3 downstream residues in this run."),
    "netcleave-ii-hla": (
        "NetCleave-II HLA", "MHC-II processing proxy", "no named cathepsin assigned",
        "MHC-II ligand-derived C-terminal processing; 13 upstream plus 3 downstream residues in this run."),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        raise ValueError("Cannot write an empty assessment table")
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verified_source(source):
    """Verify the complete immutable manifest, including its relative paths."""
    source = Path(source).resolve()
    manifest = json.loads((source / "SHA256SUMS.json").read_text())
    for relative, expected in manifest.items():
        path = (source / relative).resolve()
        if not path.is_relative_to(source) or not path.is_file() or sha256(path) != expected:
            raise ValueError("Source checksum verification failed: " + relative)
    required = {
        "provenance.json", "tables/vaccine_sequence_inventory.csv",
        "tables/model_catalog.csv", "tables/slp_quantitative_bond_scores.csv",
        "tables/slp_motif_assessments.csv",
    }
    if not required <= set(manifest):
        raise ValueError("Source manifest does not cover required inputs")
    return manifest


def preferred_pepsickle(organism):
    """Apply the requested display preference without asserting performance."""
    if organism not in ORGANISMS:
        raise ValueError("Unknown organism context: " + organism)
    return HUMAN if organism == "human" else MAMMAL


def mechanism_catalog(source_catalog, organism):
    """Annotate frozen model identities without loading predictor assets."""
    preferred = preferred_pepsickle(organism)
    builtin = {model.name: model for model in cleavage_models()}
    builtin[ERAMERCleavage.model.name] = ERAMERCleavage.model
    rows = []
    for original in source_catalog:
        name = original["model"]
        if name in PROCESSING:
            label, route, enzyme, assay = PROCESSING[name]
            topology, removed, strictness = "internal", "", ""
            exposure = (
                "Endolysosomal access after uptake; no named-enzyme exposure inferred"
                if name == "netcleave-ii-hla" else
                "Cytosolic access after SLP uptake/export; not direct extracellular cleavage"
            )
            limitations = "Native processing score; no kinetics, uptake, or enzyme/subunit attribution."
            references = original["reference"]
            if name == "netchop-3.1-cterm-3.0" or name == "netchop-3.1-20s-3.0":
                references += ";https://services.healthtech.dtu.dk/services/NetChop-3.1/"
            if name in (HUMAN, MAMMAL):
                references += ";https://github.com/pdxgx/pepsickle"
                limitations += " Human-only is experimental and trained on a smaller dataset; preference is a user policy."
        else:
            if name not in builtin:
                raise ValueError("Unannotated source model: " + name)
            metadata = builtin[name]
            enzyme, assay = metadata.enzyme, metadata.assay
            route = ";".join(metadata.compartments)
            label = enzyme + " / " + name
            strictness = metadata.motif_strictness or ""
            references = ";".join(metadata.references)
            limitations = metadata.limitations
            if metadata.evidence == "motif_rule":
                rule = get_cleavage_model(name)
                topology, removed = rule.topology, str(rule.removed)
            elif name == "dpp4-qpisa":
                topology, removed = "n_terminal", "2"
            elif name == "eramer-step":
                topology, removed = "n_terminal", "1"
            else:
                topology, removed = "whole_substrate_reference", ""
            exposure = "Requires enzyme access in " + route + "; location does not establish exposure"
            if topology in ("n_terminal", "c_terminal"):
                exposure += "; only the currently exposed terminus, no successive trimming simulated"
            if name == "cpb2-basic":
                exposure += "; frozen producer explicitly assumed activated CPB2"
        rows.append({
            "model": name, "label": label, "biological_route": route,
            "enzyme_attribution": enzyme, "topology": topology,
            "residues_removed": removed, "training_or_assay": assay,
            "evidence_type": original["evidence_type"],
            "score_units": original["score_units"], "motif_strictness": strictness,
            "exposure_requirement": exposure, "limitations": limitations,
            "pepsickle_preference": (
                "preferred" if name == preferred else "comparison" if name in (HUMAN, MAMMAL) else ""
            ),
            "source_included": original["included"],
            "original_biological_context": original["biological_context"],
            "primary_references": references,
        })
    return rows


def core_bounds(record):
    """Return verified zero-based half-open bounds of the frozen annotation."""
    core, offset = record["minimal_epitope"], record["minimal_epitope_offset"]
    if not core or not offset:
        return None
    start = int(float(offset))
    end = start + len(core)
    if record["sequence"][start:end] != core:
        raise ValueError("Core annotation does not match " + record["sequence_record_id"])
    return start, end


def enriched_assessments(records, quantitative, motifs, catalog):
    """Join coordinates and mechanisms while preserving raw score strings."""
    parents = {record["sequence_record_id"]: record for record in records}
    models = {row["model"]: row for row in catalog}
    rows, seen = [], set()
    for evidence, source_rows in (("quantitative_model", quantitative), ("motif_rule", motifs)):
        for original in source_rows:
            parent = parents[original["sequence_record_id"]]
            model = models[original["model"]]
            if model["evidence_type"] != evidence:
                raise ValueError("Source evidence type mismatch: " + original["model"])
            if original["sequence"] != parent["sequence"]:
                raise ValueError("Source sequence mismatch: " + parent["sequence_record_id"])
            value = original["bond"]
            bond = int(float(value)) if value else None
            if value and float(value) != bond:
                raise ValueError("Noninteger bond")
            key = parent["sequence_record_id"], original["model"], bond
            if key in seen:
                raise ValueError("Duplicate source assessment: " + str(key))
            seen.add(key)
            if bond is not None:
                sequence = parent["sequence"]
                if not 1 <= bond < len(sequence):
                    raise ValueError("Bond outside input")
                if (original["left_residue"], original["right_residue"]) != (sequence[bond - 1], sequence[bond]):
                    raise ValueError("Bond residues disagree with input")
                expected = (
                    int(model["residues_removed"]) if model["topology"] == "n_terminal" else
                    len(sequence) - int(model["residues_removed"]) if model["topology"] == "c_terminal" else bond
                )
                if bond != expected:
                    raise ValueError("Terminal model assessed an internal/unexposed bond")
            score = original.get("score", "")
            if evidence == "quantitative_model":
                state = "scored" if original["assessable"] == "True" else "unassessed"
                if state == "scored" and (bond is None or not math.isfinite(float(score))):
                    raise ValueError("Scored assessment needs a finite score and bond")
                if state == "unassessed" and score:
                    raise ValueError("Unassessed source has a numerical score")
            else:
                state = original["status"]
                if state not in ("matched", "not_matched", "unsupported"):
                    raise ValueError("Unknown motif state")
            bounds = core_bounds(parent)
            in_core = bool(bond is not None and bounds and bounds[0] < bond < bounds[1])
            rows.append({
                "sequence_record_id": parent["sequence_record_id"], "gene": parent["gene"],
                "protein_change": parent["protein_change"], "vaccines": parent["vaccines"],
                "sequence": parent["sequence"], "sequence_source_url": parent["source_url"],
                "bond": "" if bond is None else str(bond), "bond_label": original["bond_label"],
                "minimal_epitope": parent["minimal_epitope"],
                "minimal_epitope_offset": parent["minimal_epitope_offset"],
                "core_bond": str(bond - bounds[0]) if in_core else "",
                "inside_source_minimal_epitope": str(in_core),
                "model": original["model"], "assessment_state": state, "score": score,
                "score_units": original.get("score_units", "") if evidence == "quantitative_model" else "",
                "motif_strictness": original.get("motif_strictness", ""),
                "assessment_reason": original.get("reason", ""),
                "unsupported_reason": original.get("unsupported_reason", ""),
                "recorded_biological_context": original.get("biological_context", original.get("compartments", "")),
                "biological_route": model["biological_route"],
                "enzyme_attribution": model["enzyme_attribution"], "topology": model["topology"],
                "training_or_assay": model["training_or_assay"],
                "exposure_requirement": model["exposure_requirement"],
                "pepsickle_preference": model["pepsickle_preference"],
                "primary_references": model["primary_references"],
            })
    return rows


def bond_segments(length, maximum=MAX_BONDS):
    return [list(range(start, min(length, start + maximum))) for start in range(1, length, maximum)]


def site_lookup(rows):
    return {(row["sequence_record_id"], row["model"], row["bond"]): row for row in rows}


def assessment_at(lookup, record_id, model, bond):
    return lookup.get((record_id, model, str(bond))) or lookup.get((record_id, model, ""))


def cell_text(row):
    if row is None or row["assessment_state"] in ("unassessed", "unsupported"):
        return "NA"
    if row["assessment_state"] == "scored":
        if row.get("model") == "dpp4-qpisa":
            return dpp4_loss_label(float(row["score"]))
        return "%.3f" % float(row["score"])
    return "Match" if row["assessment_state"] == "matched" else "No match"


def render_pdf(path, records, rows, catalog, organism, source_name):
    """Draw exact site scorecards, keeping terminal models outside the grid."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.platypus import Paragraph, Table, TableStyle

    width, height = landscape(A4)
    canvas = Canvas(str(path), pagesize=(width, height))
    canvas.setTitle("Sid vaccine peptides: cleavage sites, scores and mechanisms")
    canvas.setAuthor("mhctools")
    ink, muted = colors.HexColor("#172c3a"), colors.HexColor("#54636c")
    style = ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12, textColor=ink)

    def motif_marker(row, x, y):
        if row is None or row["assessment_state"] in ("unassessed", "unsupported"):
            canvas.setFillColor(muted)
            canvas.setFont("Helvetica", 8)
            canvas.drawCentredString(x, y, "NA")
        elif row["assessment_state"] == "matched":
            colour = "#7b3e98" if row["model"].startswith("fap-") else "#007b83"
            canvas.setFillColor(colors.HexColor(colour))
            path = canvas.beginPath()
            path.moveTo(x - 6, y + 8)
            path.lineTo(x + 6, y + 8)
            path.lineTo(x, y - 3)
            path.close()
            canvas.drawPath(path, stroke=0, fill=1)
        else:
            canvas.setStrokeColor(muted)
            canvas.setLineWidth(.8)
            canvas.circle(x, y + 2, 3, stroke=1, fill=0)

    def paragraph(text, x, y, available=width - 72, size=9):
        para = Paragraph(text, ParagraphStyle("p", parent=style, fontSize=size, leading=size + 3))
        _, h = para.wrap(available, height)
        para.drawOn(canvas, x, y - h)
        return y - h

    def header(title, subtitle):
        canvas.setFillColor(ink)
        canvas.setFont("Helvetica-Bold", 18)
        canvas.drawString(36, height - 42, title)
        return paragraph(subtitle, 36, height - 54)

    def footer():
        canvas.setFillColor(muted)
        canvas.setFont("Helvetica", 7)
        canvas.drawString(36, 22, "Frozen inference: " + source_name + " | mechanism annotations and display policy are separate")
        canvas.drawRightString(width - 36, 22, str(canvas.getPageNumber()))
        canvas.showPage()

    def table(data, widths, y, font_size=8):
        formatted = [[Paragraph(escape(str(cell)), ParagraphStyle(
            "cell", parent=style, fontSize=font_size, leading=font_size + 2)) for cell in row] for row in data]
        item = Table(formatted, colWidths=widths, hAlign="LEFT")
        item.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eef1")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("LINEBELOW", (0, 0), (-1, -1), 0.3, colors.HexColor("#d3dce1")),
        ]))
        _, h = item.wrap(width - 72, height)
        if y - h < 38:
            raise ValueError("PDF table does not fit page")
        item.drawOn(canvas, 36, y - h)
        return y - h

    preferred = preferred_pepsickle(organism)
    models = {row["model"]: row for row in catalog}
    processing_order = [preferred, MAMMAL if preferred == HUMAN else HUMAN,
                        "netchop-3.1-cterm-3.0", "netchop-3.1-20s-3.0", "netcleave-i-hla", "netcleave-ii-hla"]
    processing_captions = {
        HUMAN: "MHC-I epitope proxy; proteasome type agnostic",
        MAMMAL: "MHC-I epitope proxy; proteasome type agnostic",
        "netchop-3.1-cterm-3.0": "MHC-I ligand C-terminal processing proxy",
        "netchop-3.1-20s-3.0": "In-vitro proteasome digestion model",
        "netcleave-i-hla": "MHC-I C-terminal processing proxy",
        "netcleave-ii-hla": "MHC-II processing proxy; no named cathepsin",
    }
    y = header("Cleavage sites, scores and mechanisms", "Sid's disclosed synthetic long vaccine peptides | human study | frozen local predictions")
    y = paragraph(
        f"<b>Coverage:</b> {len(records)} source records, {len({r['sequence'] for r in records})} unique SLP sequences. "
        "RNA-encoded sequences are outside this view. Gold marks the source-listed minimal epitope only where its recorded offset matches. "
        "That annotation is an mRNA minimal/candidate epitope and does not establish the experimentally used minimal target of every SLP.", 36, y - 18)
    y = paragraph(
        "<b>Coordinates:</b> bond b is between SLP residues b and b+1. An internal core cut lies strictly between the core boundaries. "
        "Numeric grid cells show native model scores, rounded to three decimals; DPP4 shows predicted percent loss in the source four-hour assay. "
        "CSV retains native scores. NA means unassessed/unsupported. Triangles mark motif matches; hollow circles mark non-matches.", 36, y - 14)
    y = paragraph(
        "<b>Mechanism:</b> ligand-trained models report processing proxies; they cannot name the enzyme responsible for a bond. "
        "20S reports an in-vitro proteasome model. Named peptidases have their own assay/motif evidence. "
        "For injected SLPs, cytosolic scores require uptake/export, ER trimming requires ER access, and extracellular activity requires enzyme exposure. "
        "No uptake, successive digestion, kinetics, MHC protection, or overall survival is predicted.", 36, y - 14)
    y = paragraph(
        f"<b>Pepsickle preference:</b> organism context = {organism}; {escape(models[preferred]['label'])} is preferred. "
        "All-mammal is preferred for nonhuman, mixed or unknown context. Both scores remain separate; this is a user-selected population policy, "
        "not evidence that human-only is more accurate. Upstream calls human-only experimental and notes its smaller training set. "
        "All-mammal training does not validate arbitrary non-mammalian organisms.", 36, y - 14)
    y = paragraph(
        "<b>NetChop performance:</b> DTU reports Cterm 3.0 performs best for CTL epitope boundaries. "
        "It is trained on MHC-I ligand C-termini; 20S uses in-vitro degradation data. "
        "Cterm is useful for ligand boundaries, while 20S has the more direct purified-proteasome interpretation. "
        "Neither score is a clinical cleavage probability. The 0.5 emphasis applies only to each 0-1 processing row; it is not a cross-model calibration.", 36, y - 14)
    paragraph("<b>Provenance:</b> source checksums verified; inference not repeated; no CleaveNet results included. "
              "Full model identities, assay scopes and primary references are in mechanism_catalog.csv; exact states and scores are in cleavage_sites.csv. "
              "The three exact-sequence reference models (THOP1, NLN, IRAP) were excluded in the source run and are not extrapolated here.", 36, y - 14)
    footer()

    y = header("Processing models: distinct endpoints", "The same bond can score differently because training endpoints differ; no ensemble score is computed.")
    data = [["Model", "Biological route / enzyme attribution", "Training endpoint and applicability"]]
    for name in processing_order:
        model = models[name]
        data.append([model["label"] + (" (preferred)" if name == preferred else ""),
                     model["biological_route"] + "; " + model["enzyme_attribution"], model["training_or_assay"]])
    y = table(data, [172, 245, width - 72 - 417], y - 18)
    y = paragraph("<b>DPP4 qPISA:</b> only the exposed N-terminal dipeptide bond (b=2) of the intact SLP. "
                  "Predicted percent loss in the source four-hour assay; native scores remain in CSV. Missing coefficients cause abstention. "
                  "<b>ERAMER:</b> ERAP1's initial N-terminal trimming step (b=1), only for 9-16-residue inputs; native PWM specificity.", 36, y - 16)
    paragraph("<b>Key primary sources:</b><br/>"
              "NetChop: https://services.healthtech.dtu.dk/services/NetChop-3.1/ ; DOI 10.1007/s00251-005-0781-7<br/>"
              "Pepsickle: DOI 10.1093/bioinformatics/btab628 ; https://github.com/pdxgx/pepsickle<br/>"
              "NetCleave: DOI 10.1038/s41598-021-92632-y<br/>"
              "DPP4: DOI 10.1038/s44320-024-00071-4 ; ERAMER: PMID 38925438", 36, y - 14, size=8)
    footer()

    y = header("Named peptidases: topology matters", "All assessments use the intact SLP with assumed free termini; repeated trimming and newly exposed fragments are not simulated.")
    data = [["Model / enzyme", "Topology", "Annotated locations", "Evidence"]]
    for model in catalog:
        if model["model"] in PROCESSING or model["source_included"] != "True":
            continue
        evidence = model["motif_strictness"] + " motif" if model["evidence_type"] == "motif_rule" else model["score_units"]
        data.append([model["model"] + " / " + model["enzyme_attribution"],
                     model["topology"].replace("_", " ") + ("; removes " + model["residues_removed"] if model["residues_removed"] not in ("", "0") else ""),
                     model["biological_route"], evidence])
    y = table(data, [205, 155, 165, width - 72 - 525], y - 10, font_size=7)
    paragraph("Required = clears this rule's recognition gate; preferred = incomplete preference; permissive = broad, weak flag. "
              "None establishes a cleavage rate. PREP: PMID 22750443; MME: PMID 6349683; FAP: PMID 16480718. "
              "CPB2 explicitly assumes activated enzyme in the frozen analysis. Every rule's complete primary references and limitations are in the catalog.",
              36, y - 10, size=8)
    footer()

    lookup, index = site_lookup(rows), []
    terminal_models = [model for model in catalog if model["topology"] in ("n_terminal", "c_terminal")]
    for record in records:
        bounds = core_bounds(record)
        segments = bond_segments(len(record["sequence"]))
        for part, bonds in enumerate(segments, start=1):
            page = canvas.getPageNumber()
            bookmark = record["sequence_record_id"] + "-" + str(part)
            canvas.bookmarkPage(bookmark)
            canvas.addOutlineEntry(record["gene"] + " / " + record["vaccines"] + f" / bonds {bonds[0]}-{bonds[-1]}", bookmark)
            index.append({"sequence_record_id": record["sequence_record_id"], "gene": record["gene"],
                          "vaccines": record["vaccines"], "first_bond": bonds[0], "last_bond": bonds[-1], "pdf_page": page})
            header(record["gene"] + " " + record["protein_change"],
                   escape(record["vaccines"]) + f" | {len(record['sequence'])} aa | bonds {bonds[0]}-{bonds[-1]} | segment {part}/{len(segments)}")
            segment = record["sequence"][bonds[0] - 1:bonds[-1] + 1]
            canvas.setFillColor(ink)
            canvas.setFont("Courier-Bold", 13)
            canvas.drawString(36, height - 95, segment)
            core_note = (
                "Source minimal epitope: " + record["minimal_epitope"] + f" (SLP {bounds[0]+1}-{bounds[1]})"
                if bounds else "No located source minimal epitope for this SLP; all internal bonds still shown"
            )
            paragraph(escape(core_note), 36, height - 106, size=8)
            paragraph("INTERNAL SITE SCORES - native 0-1 outputs; bold cells meet that model's 0.5 display threshold", 36, height - 122, size=8)
            left, col_width = 244, (width - 36 - 244) / MAX_BONDS
            y = height - 151
            for col, bond in enumerate(bonds):
                x = left + col * col_width
                in_core = bool(bounds and bounds[0] < bond < bounds[1])
                canvas.setFillColor(colors.HexColor("#f6dea0") if in_core else colors.HexColor("#edf1f3"))
                canvas.rect(x, y - 17, col_width - 1, 29, fill=1, stroke=0)
                canvas.setFillColor(ink)
                canvas.setFont("Helvetica", 7)
                canvas.drawCentredString(x + col_width / 2, y + 1, str(bond))
                canvas.setFont("Courier", 8)
                canvas.drawCentredString(x + col_width / 2, y - 12, record["sequence"][bond-1] + "|" + record["sequence"][bond])
            for row_number, name in enumerate(processing_order):
                model = models[name]
                row_y = y - 29 - row_number * 22
                label = model["label"] + (" [preferred]" if name == preferred else " [comparison]" if name in (HUMAN, MAMMAL) else "")
                canvas.setFillColor(ink)
                canvas.setFont("Helvetica-Bold" if name == preferred else "Helvetica", 8)
                canvas.drawString(36, row_y, label)
                canvas.setFillColor(muted)
                canvas.setFont("Helvetica", 6.5)
                canvas.drawString(36, row_y - 9, processing_captions[name])
                for col, bond in enumerate(bonds):
                    row = assessment_at(lookup, record["sequence_record_id"], name, bond)
                    high = bool(row and row["assessment_state"] == "scored" and float(row["score"]) >= .5)
                    x = left + col * col_width
                    fill = "#e4dff0" if name == "netcleave-ii-hla" else "#d6ece7" if name == "netchop-3.1-20s-3.0" else "#e2ebf4"
                    canvas.setFillColor(colors.HexColor(fill if high else "#f5f7f8"))
                    canvas.rect(x, row_y - 5, col_width - 1, 17, fill=1, stroke=0)
                    canvas.setFillColor(ink if row and row["assessment_state"] == "scored" else muted)
                    canvas.setFont("Helvetica-Bold" if high else "Helvetica", 7.5)
                    canvas.drawCentredString(x + col_width / 2, row_y, cell_text(row))
            y -= 169
            paragraph("INTERNAL ENZYME MOTIFS - triangle = match; hollow circle = no match; NA = outside scope", 36, y, size=8)
            for row_number, name in enumerate(INTERNAL_MOTIFS):
                row_y = y - 28 - row_number * 17
                canvas.setFillColor(ink)
                canvas.setFont("Helvetica", 8)
                canvas.drawString(36, row_y, models[name]["enzyme_attribution"] + " / " + name + " (" + models[name]["motif_strictness"] + ")")
                for col, bond in enumerate(bonds):
                    row = assessment_at(lookup, record["sequence_record_id"], name, bond)
                    x = left + col * col_width
                    canvas.setFillColor(colors.HexColor("#f5f7f8"))
                    canvas.rect(x, row_y - 5, col_width - 1, 15, fill=1, stroke=0)
                    motif_marker(row, x + col_width / 2, row_y)
            y -= 90
            paragraph("EXPOSED TERMINI ONLY - N1/N2/N3 remove 1/2/3 N-terminal residues; C1/C2 remove 1/2 C-terminal residues", 36, y, size=8)
            for i, model in enumerate(terminal_models):
                column, line = i // 8, i % 8
                x, row_y = 36 + column * 390, y - 28 - line * 14
                matches = [r for r in rows if r["sequence_record_id"] == record["sequence_record_id"] and r["model"] == model["model"]]
                if len(matches) != 1:
                    raise ValueError("Expected one terminal assessment per intact input")
                row = matches[0]
                suffix = " | " + row["bond_label"] if row["bond"] else " | no assessed bond"
                number = cell_text(row)
                motif = row["assessment_state"] in ("matched", "not_matched")
                if motif:
                    number = ""
                elif row["assessment_state"] == "scored" and model["model"] != "dpp4-qpisa":
                    number += " PWM specificity"
                topology_label = ("N" if model["topology"] == "n_terminal" else "C") + model["residues_removed"]
                enzyme = model["enzyme_attribution"] + (" (active assumption)" if model["model"] == "cpb2-basic" else "")
                text = enzyme + " [" + topology_label + "] / " + model["model"] + suffix + " : " + number
                canvas.setFillColor(ink)
                canvas.setFont("Helvetica", 7.5)
                canvas.drawString(x, row_y, text)
                if motif:
                    motif_marker(row, x + canvas.stringWidth(text, "Helvetica", 7.5) + 8, row_y)
            paragraph("Gold bond headers are strictly inside the located source epitope. Exposure and chemical assumptions are in the catalog. "
                      "CSV retains unassessed reasons and native precision. No score aggregation or resistance claim.", 36, 51, size=7)
            footer()
    canvas.save()
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; an existing directory is refused")
    parser.add_argument("--organism", choices=ORGANISMS, default="human", help="Pepsickle display preference; does not rerun inference")
    args = parser.parse_args()
    source = args.source_run.resolve()
    manifest = verified_source(source)
    records = [row for row in read_csv(source / "tables/vaccine_sequence_inventory.csv") if row["sequence_type"] == "synthetic_long_peptide"]
    catalog = mechanism_catalog(read_csv(source / "tables/model_catalog.csv"), args.organism)
    rows = enriched_assessments(records, read_csv(source / "tables/slp_quantitative_bond_scores.csv"),
                               read_csv(source / "tables/slp_motif_assessments.csv"), catalog)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "mechanism_catalog.csv", catalog)
    write_csv(output / "cleavage_sites.csv", rows)
    index = render_pdf(output / PDF_NAME, records, rows, catalog, args.organism, source.name)
    write_csv(output / "page_index.csv", index)
    source_provenance = json.loads((source / "provenance.json").read_text())
    provenance = {
        "schema_version": 1, "generated_at": datetime.now().astimezone().isoformat(),
        "prediction_execution": "not repeated; complete frozen manifest verified",
        "prediction_source_run": source.name, "source_manifest_sha256": sha256(source / "SHA256SUMS.json"),
        "source_files_verified": len(manifest), "source_provenance": source_provenance,
        "organism_context": args.organism, "preferred_pepsickle": preferred_pepsickle(args.organism),
        "score_display": "DPP4: whole-percent predicted loss in source 4 h assay; negative estimates: No predicted loss; small positive estimates: <1%; other numeric scores: three decimals; CSV retains native source strings",
        "display_labels_sha256": sha256(Path(__file__).with_name("display_labels.py")),
        "source_minimal_epitope_annotation": "source mRNA minimal/candidate epitope; only verified recorded offsets are highlighted",
        "renderer_sha256": sha256(Path(__file__)), "sequence_uploads": False,
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "git_worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], text=True)),
        "metadata_code_sha256": {name: sha256(Path(__file__).parents[2] / "mhctools" / name)
                                 for name in ("peptidases.py", "dpp4.py", "eramer_cleavage.py", "substrate_reference.py")},
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    (output / "SHA256SUMS.json").write_text(json.dumps({path.name: sha256(path) for path in sorted(output.iterdir()) if path.is_file()}, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
