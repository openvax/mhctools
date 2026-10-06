#!/usr/bin/env python3
"""Render native stability-route outputs; this command performs no inference."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
from xml.sax.saxutils import escape

from display_labels import dpp4_loss_label
from mechanism_scorecards import bond_segments, core_bounds, read_csv, sha256, write_csv
from stability_route import SERUM_MODELS


PDF_NAME = "stability-route-predictions.pdf"
RED = "#C5323B"
INK = "#182D3D"
GRAY = "#667A89"
GOLD = "#FFE6A0"
MEMBRANE = {"MMP14", "MMP15", "MMP16", "MMP24"}
GPI = {"MMP17", "MMP25"}
SERUM_LABELS = {
    "dpp4-qpisa": "DPP4 assay loss", "ace-dipeptidyl": "ACE / C-terminal 2 aa",
    "cpn-basic": "CPN / C-terminal 1 aa", "cpb2-basic": "CPB2 / activated",
    "app2-xp": "APP2 / N-terminal 1 aa", "fap-endo-gp": "FAP / internal",
    "fap-dipeptidyl": "FAP / N-terminal 2 aa", "mme-hydrophobic": "MME / internal",
    "anpep-ala": "ANPEP / N-terminal 1 aa", "enpep-acidic": "ENPEP / N-terminal 1 aa",
}


class Report:
    def __init__(self, output):
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen.canvas import Canvas
        self.output = output
        self.width, self.height = A4
        self.canvas = Canvas(str(output / PDF_NAME), pagesize=A4, pageCompression=1)
        self.canvas.setTitle("Sid vaccine peptides: stability-route predictions")
        self.canvas.setAuthor("mhctools")
        self.page = 0
        self.index = []

    def text(self, x, y, value, size=10, color=INK, font="Helvetica"):
        from reportlab.lib.colors import HexColor
        self.canvas.setFillColor(HexColor(color))
        self.canvas.setFont(font, size)
        self.canvas.drawString(x, y, str(value))

    def paragraph(self, x, top, value, width=None, size=10, color=INK):
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.platypus import Paragraph
        style = ParagraphStyle("body", fontName="Helvetica", fontSize=size,
                               leading=size * 1.35, textColor=color)
        paragraph = Paragraph(escape(value), style)
        width = width or self.width - x - 36
        _, height = paragraph.wrap(width, 1000)
        paragraph.drawOn(self.canvas, x, top - height)
        return top - height

    def site_marker(self, x, bottom, top):
        from reportlab.lib.colors import HexColor
        self.canvas.saveState()
        self.canvas.setStrokeColor(HexColor(RED))
        self.canvas.setLineWidth(1.6)
        self.canvas.setDash(3, 2)
        self.canvas.line(x, bottom, x, top)
        self.canvas.restoreState()

    def new_page(self, title, subtitle, record=None, kind="overview", segment=""):
        if self.page:
            self.canvas.showPage()
        self.page += 1
        self.index.append(dict(page=self.page, section=kind, segment=segment,
                               sequence_record_id=record["sequence_record_id"] if record else ""))
        self.text(36, 808, "STABILITY ROUTE / NATIVE MODEL EVIDENCE", 9, GRAY)
        self.text(36, 772, title, 24, font="Helvetica-Bold")
        self.paragraph(36, 752, subtitle, size=10)
        if record:
            self.text(36, 716, record["sequence_record_id"], 8, GRAY)
        self.text(36, 22, "mhctools / human context / conditional enzyme access", 8, GRAY)
        self.text(self.width - 65, 22, str(self.page), 9, GRAY)

    def sequence(self, record, y=690, start=0, end=None, x=36, width=None):
        from reportlab.lib.colors import HexColor
        sequence = record["sequence"]
        end = len(sequence) if end is None else end
        width = width or self.width - x - 36
        step = width / (end - start)
        bounds = core_bounds(record)
        if bounds:
            low, high = max(start, bounds[0]), min(end, bounds[1])
            if low < high:
                self.canvas.setFillColor(HexColor(GOLD))
                self.canvas.rect(x + (low - start) * step, y - 6, (high - low) * step, 23,
                                 fill=1, stroke=0)
        for index in range(start, end):
            self.text(x + (index - start + .5) * step - 4, y, sequence[index],
                      min(15, step * .72), font="Courier-Bold")
            if index == start or index == end - 1 or (index + 1) % 5 == 0:
                self.text(x + (index - start + .5) * step - 4, y + 25, index + 1, 8, GRAY)
        return step

    def mmp_page(self, record, rows):
        self.new_page(record["gene"] + " / MMP substrate windows",
                      "Strongest overlapping 10-mer for each of CleaveNet's 18 MMP heads. "
                      "The window locates susceptibility; the exact cut is unknown.", record, "mmp")
        sequence = record["sequence"]
        for start in range(0, len(sequence), 35):
            end = min(len(sequence), start + 35)
            self.sequence(record, y=683 - (start // 35) * 44, start=start, end=end)
        top = 626 - ((len(sequence) - 1) // 35) * 44
        columns = (36, 94, 165, 231, 371, 428, 483)
        for x, title in zip(columns, ("MMP", "LOCATION", "RESIDUES", "10-MER", "Z-SCORE", "5-MODEL SD", "CORE AA")):
            self.text(x, top, title, 7.5, GRAY, "Helvetica-Bold")
        grouped = defaultdict(list)
        for row in rows:
            grouped[row["enzyme"]].append(row)
        for index, enzyme in enumerate(sorted(grouped, key=lambda name: int(name[3:]))):
            strongest = max(grouped[enzyme], key=lambda row: float(row["z_score"]))
            y = top - 24 - index * 21
            values = (enzyme, "membrane" if enzyme in MEMBRANE else "GPI anchor" if enzyme in GPI else "soluble",
                      "%d-%d" % (int(strongest["window_start"]) + 1, int(strongest["window_end"])),
                      strongest["window_sequence"], "%.2f" % float(strongest["z_score"]),
                      "%.2f" % float(strongest["ensemble_sd"]), strongest["source_core_overlap_residues"] or "NA")
            for x, value in zip(columns, values):
                self.text(x, y, value, 10, font="Courier" if x == 231 else "Helvetica")
        bottom = top - 414
        self.paragraph(36, bottom, "Native Z-score is relative substrate signal; SD is ensemble spread. "
                       "Each row selects the maximum mean for this MMP and peptide; first window breaks ties. "
                       "All windows remain in the CSV. CORE AA counts source-core overlap, not cuts; NA means "
                       "unestablished core. No released fragment, exact cut, intact-SLP loss, or serum exposure is inferred.", size=9)

    def bond_page(self, record, bonds, scores, serum, cats, first):
        from reportlab.lib.colors import HexColor
        self.new_page(record["gene"] + " / sites and terminal enzymes",
                      "Proteasomes: cytosolic access after uptake/export. Cathepsins: endolysosomal access. "
                      "Dashed red lines mark candidate bonds; gold marks a source core.", record, "sites",
                      "%d-%d" % (bonds[0], bonds[-1]))
        start, end = bonds[0] - 1, bonds[-1] + 1
        self.text(36, 690, "RESIDUES %d-%d" % (start + 1, end), 9, GRAY, "Helvetica-Bold")
        step = self.sequence(record, y=662, start=start, end=end, x=143, width=416)
        tracks = (
            ("pepsickle-in-vitro-2-human-only-constitutive", "Constitutive", 0, 1, .5, "#246D8E"),
            ("pepsickle-in-vitro-2-human-only-immunoproteasome", "Immunoproteasome", 0, 1, .5, "#8A598B"),
            ("itcell-cats-240", "Cathepsin S / 240", -15, 5, 3, "#246D8E"),
            ("itcell-catb-240", "Cathepsin B / 240", -15, 5, 3, "#8A598B"),
        )
        for index, (model, label, low, high, threshold, color) in enumerate(tracks):
            base = 588 - index * 62
            self.paragraph(36, base + 37, label, width=100, size=10)
            self.text(36, base + 4, "%s to %s" % (low, high), 8, GRAY)
            # Repeat residues beside each model so a cut is visibly between its flanks.
            for residue_index in range(start, end):
                self.text(143 + (residue_index - start + .5) * step - 3.3,
                          base + 42, record["sequence"][residue_index],
                          min(11, step * .72), font="Courier-Bold")
            self.canvas.setStrokeColor(HexColor("#DDE5EA"))
            self.canvas.setLineWidth(.5)
            self.canvas.line(143, base, 559, base)
            threshold_y = base + 33 * (threshold - low) / (high - low)
            self.canvas.setStrokeColor(HexColor("#BFB4B5"))
            self.canvas.setDash(2, 3)
            self.canvas.line(143, threshold_y, 559, threshold_y)
            self.canvas.setDash()
            previous = None
            for bond in bonds:
                row = scores[(record["sequence_record_id"], model, str(bond))]
                value = float(row["score"])
                if not low <= value <= high:
                    raise ValueError("Native score outside declared chart scale")
                x = 143 + (bond - start) * step
                y = base + 33 * (value - low) / (high - low)
                self.canvas.setStrokeColor(HexColor(color))
                self.canvas.setLineWidth(1.4)
                if previous:
                    self.canvas.line(*previous, x, y)
                self.canvas.setFillColor(HexColor(color))
                self.canvas.circle(x, y, 1.6, fill=1, stroke=0)
                # Pepsickle display >=.5; ITCell author threshold is strictly >3.
                candidate = value >= threshold if model.startswith("pepsickle") else value > threshold
                if candidate:
                    self.site_marker(x, base, base + 52)
                    self.text(x - 7, base - 8, "%.2f" % value, 7.5, RED)
                previous = (x, y)
        self.text(36, 354, "Pepsickle: native 0-1 output, display >=0.5. Cat B/S: native log2 specificity, >3.", 8.5, GRAY)
        self.text(36, 340, "Cat 240 is the source count profile; all 15/60/240 profiles remain in the native tables.", 8.5, GRAY)
        if first:
            h = cats[(record["sequence_record_id"], "itcell-cath-240", "1")]
            value = float(h["score"])
            self.text(36, 318, "CATHEPSIN H / INITIAL N-TERMINAL TRIM", 9, GRAY, "Helvetica-Bold")
            self.text(306, 318, "bond 1: %.2f" % value, 10, RED if value > 2 else INK)
            for offset, residue in enumerate(record["sequence"][:4]):
                self.text(468 + offset * 24, 318, residue, 14, font="Courier-Bold")
            if value > 2:
                self.site_marker(484, 312, 334)
            self.text(36, 301, "Native log2 specificity; candidate >2. Assumed free N terminus; no repeated trimming.", 8.5, GRAY)
            self.text(36, 277, "SERUM / EXTRACELLULAR PANEL - CONDITIONAL EXPOSURE", 9, GRAY, "Helvetica-Bold")
            y = 256
            for model in SERUM_MODELS:
                rows = serum.get((record["sequence_record_id"], model), [])
                if not rows:
                    raise ValueError("Missing frozen extracellular model assessment")
                unsupported = next((r for r in rows if r["unsupported_reason"]), None)
                if unsupported:
                    value = "NA: " + unsupported["unsupported_reason"]
                elif model == "dpp4-qpisa":
                    value = "bond 2: " + dpp4_loss_label(float(rows[0]["score"])) + " (purified-enzyme assay)"
                else:
                    matched = [r["bond"] for r in rows if r["assessment_state"] == "matched"]
                    value = "Motif at bond " + ", ".join(matched) if matched else "No motif match"
                self.text(36, y, SERUM_LABELS[model], 9)
                bottom = self.paragraph(208, y + 9, value, width=351, size=8.5,
                                        color=GRAY if unsupported or value == "No motif match" else INK)
                if any(r["assessment_state"] == "matched" for r in rows):
                    self.site_marker(194, y - 3, y + 10)
                y = min(y - 18, bottom - 5)
            if y < 40:
                raise ValueError("Serum panel crosses footer")
        else:
            self.paragraph(36, 311, "Terminal enzyme assessments are on this record's first site page. "
                           "This continuation retains every internal bond without compressing the residue axis.", size=11)
            self.paragraph(36, 256, "These models assess conditional sequence susceptibility after uptake. "
                           "They do not predict cellular uptake, presentation, enzyme dose, whole-serum survival "
                           "or immune response. Full-precision native scores and context limitations are in the tables.", size=11)

    def finish(self):
        self.canvas.save()
        write_csv(self.output / "page_index.csv", self.index)


def render(output):
    if (output / "SHA256SUMS.json").exists():
        raise ValueError("Refusing to overwrite a frozen report")
    records = read_csv(output / "vaccine_records.csv")
    mmp = read_csv(output / "mmp_window_scores.csv")
    peps = read_csv(output / "proteasome_bond_scores.csv")
    cat = read_csv(output / "cathepsin_bond_scores.csv")
    serum = read_csv(output / "existing_enzyme_assessments.csv")
    audit = json.loads((output / "model_availability.json").read_text())
    report = Report(output)
    report.new_page("Stability on the route to presentation", "%d disclosed vaccine peptide records / "
                    "%d unique sequences / human organism preference" %
                    (len(records), len({r["sequence"] for r in records})))
    core_count = sum(core_bounds(r) is not None for r in records)
    y = report.paragraph(36, 697, "New predictions cover active matrix metalloproteinases, human constitutive "
                         "and immunoproteasome digestion, and human cathepsin B/S internal specificity and H "
                         "initial N-terminal trimming. Existing extracellular and serum/plasma enzyme "
                         "assessments are included beside each peptide. Verified source minimal/candidate "
                         "cores are available for %d/%d records; remaining core interiors are unassessed." %
                         (core_count, len(records)), size=13)
    for title, description in (
        ("BEFORE UPTAKE / MMPs", "%d overlapping ten-mers x 18 MMPs = %d native window scores. "
         "CleaveNet scores substrates; it does not locate a cut." % (len(mmp) // 18, len(mmp))),
        ("AFTER UPTAKE / PROTEASOMES", "%d bond scores from human Pepsickle digestion C/I models; "
         "cytosolic access after uptake/export is required. "
         "Separate from the existing epitope-trained Pepsickle and NetChop Cterm processing proxies." % len(peps)),
        ("AFTER UPTAKE / CATHEPSINS", "%d assessments from nine published ITCell B/S/H profiles. "
         "Main pages show 240-minute profiles; all three source profiles are retained." % len(cat)),
        ("EXPOSED TERMINI / BLOOD AND TISSUES", "DPP4 remains on its purified-enzyme assay-loss scale. "
         "ACE, CPN, activated CPB2, APP2, FAP, MME, ANPEP and ENPEP remain qualitative motif assessments."),
    ):
        report.text(36, y - 31, title, 10, GRAY, "Helvetica-Bold")
        y = report.paragraph(36, y - 43, description, size=12)
    report.paragraph(36, y - 28, "No combined serum-survival percentage or half-life is inferred. "
                     "Sequence susceptibility requires actual enzyme exposure, activation and accessible "
                     "substrate. Free termini are an assumption where terminal trimming is assessed. "
                     "Source cores are minimal/candidate annotations; not every vaccine's target is established.", size=11)
    report.new_page("How to read the models", "Native endpoints stay separate; dashed red lines mark candidate bonds.")
    y = 695
    for text in (
        "CleaveNet / active MMP substrate susceptibility: mean relative cleavage Z-score and five-model "
        "population SD. Higher score means higher predicted source-assay signal. No exact cut, survival "
        "probability or degradation time is supplied. Soluble and membrane MMPs require separate exposure assumptions.",
        "Human Pepsickle C/I / digestion output: every internal bond is scored. Each track repeats the sequence; "
        "a dashed red line between residues uses the 0.5 display threshold. Human-only models are experimental "
        "and use a smaller dataset; organism matching "
        "is the requested preference, not evidence of higher performance.",
        "ITCell / cathepsin specificity: source pH 6.5, recombinant human enzyme, 228 tetradecapeptides. "
        "Native log2 scores use author thresholds >3 for B/S and >2 for initial H trimming. Profile times are "
        "training observations, not predicted degradation times. Missing flanks contribute zero; source M "
        "represents norleucine and cysteine is unrepresented. Context flags are retained in the native data.",
        "Existing terminal enzyme evidence / intact input: positions mark only presently exposed N/C "
        "termini. DPP4 assay-loss percentage is relative signal in a purified-enzyme experiment, not loss in "
        "serum. Motif matches and non-matches remain qualitative. CPB2 assumes activated enzyme; no fragment "
        "cascade is simulated. NA gives the exact unavailable reason.",
        "Gold source core: a disclosed minimal/candidate epitope only where the sequence offset verifies. "
        "A cut inside this span may warrant investigation before binding; these predictions do not measure "
        "bound-MHC protection or prove epitope destruction. An overlapping MMP window is not a core cut.",
    ):
        y = report.paragraph(36, y, text, size=12) - 24
    report.new_page("Additional enzyme coverage", "What is executable now, and what remains unavailable.")
    y = 695
    for entry in audit:
        report.text(36, y, entry["predictor"], 12, font="Helvetica-Bold")
        y = report.paragraph(36, y - 12, entry["status"] + ": " + entry["reason"], size=10) - 18
    grouped_mmp, grouped_serum = defaultdict(list), defaultdict(list)
    for row in mmp:
        grouped_mmp[row["sequence_record_id"]].append(row)
    for row in serum:
        grouped_serum[(row["sequence_record_id"], row["model"])].append(row)
    score_lookup = {(r["sequence_record_id"], r["model"], r["bond"]): r for r in peps + cat}
    for record in records:
        report.mmp_page(record, grouped_mmp[record["sequence_record_id"]])
        for index, bonds in enumerate(bond_segments(len(record["sequence"]))):
            report.bond_page(record, bonds, score_lookup, grouped_serum, score_lookup, index == 0)
    report.new_page("Sources and reproducibility", "Full-precision evidence and exact source identities accompany this PDF.")
    y = 695
    for entry in audit:
        y = report.paragraph(36, y, entry["predictor"] + ": " + "; ".join(entry["sources"]), size=9) - 14
    for text in (
        "NetChop model distinction: https://services.healthtech.dtu.dk/services/NetChop-3.1/",
        "DPP4 assay: https://pmc.ncbi.nlm.nih.gov/articles/PMC11612144/",
        "Native data: native_predictions.json; mmp_window_scores.csv; proteasome_bond_scores.csv; "
        "cathepsin_bond_scores.csv; existing_enzyme_assessments.csv; existing_model_catalog.csv.",
        "Provenance: source and model/code/weight hashes, substrate identities, runtime versions, exposure "
        "assumptions and core scope. SHA256SUMS.json covers all final files; page_index.csv maps records to pages.",
        "Validation: source-code conformance and coordinate checks establish faithful execution. "
        "They do not supply an independently audited accuracy benchmark or whole-serum calibration.",
    ):
        y = report.paragraph(36, y, text, size=10) - 16
    report.finish()
    provenance = json.loads((output / "provenance.json").read_text())
    provenance["scripts"][Path(__file__).name] = sha256(__file__)
    provenance["pdf_pages"] = report.page
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(output / PDF_NAME)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    render(parser.parse_args().output.resolve())
