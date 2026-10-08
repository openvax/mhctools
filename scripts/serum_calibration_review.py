#!/usr/bin/env python3
"""Audit reported hormone measurements and export the compact route guide."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from mhctools import __version__
from mhctools.serum_calibration import (
    serum_assay_parent_reference, serum_calibration_evidence,
)


def reported_hormone_table(xml_path):
    """Expand the publisher's rowspans and preserve each Table 2 record."""
    root = ET.parse(xml_path).getroot()
    if root.findtext('.//article-id[@pub-id-type="doi"]') != "10.1371/journal.pone.0134427":
        raise ValueError("Expected Yi 2015 primary-source XML")
    table = root.find('.//table-wrap[@id="pone.0134427.t002"]')
    if table is None:
        raise ValueError("Source Table 2 is missing")
    active = {}
    rows = []
    for tr in table.findall('.//tbody/tr'):
        row = [None] * 5
        following = {}
        for column, (value, remaining) in active.items():
            row[column] = value
            if remaining > 1:
                following[column] = (value, remaining - 1)
        column = 0
        for cell in tr:
            while column < len(row) and row[column] is not None:
                column += 1
            value = " ".join("".join(cell.itertext()).split())
            row[column] = value
            rowspan = int(cell.attrib.get('rowspan', 1))
            if rowspan > 1:
                following[column] = (value, rowspan - 1)
            column += 1
        if any(value is None for value in row):
            raise ValueError("Incomplete source table row")
        active = following
        rows.append(row)
    if len(rows) != 27:
        raise ValueError("Expected all 27 Table 2 measurements")
    return rows


def measurement_records(rows):
    """Keep donor ranges, replicate SD and censoring distinct."""
    records = []
    for number, (substrate, sample, temperature, reported, method) in enumerate(rows, 1):
        if substrate.startswith('(') and substrate.endswith(')'):
            substrate = substrate[1:-1]
        values = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", reported)]
        point = sd = None
        if reported.startswith('>'):
            lower, upper, kind = values[0], None, "right_censored"
            uncertainty = "Lower half-life bound; no point estimate"
        elif reported.startswith('<'):
            lower, upper, kind = 0.0, values[0], "left_censored"
            uncertainty = "Upper half-life bound; no point estimate"
        elif '±' in reported:
            point, sd = values
            lower = upper = point
            kind = "mean_with_replicate_sd"
            uncertainty = "Mean +/- SD from one sample, three replicates; not a confidence interval"
        elif len(values) == 2:
            lower, upper = values
            kind = "donor_range"
            uncertainty = "Range across donor samples; not a confidence interval or paired inhibitor effect"
        else:
            point = lower = upper = values[0]
            kind = "reported_point"
            uncertainty = "No uncertainty supplied in Table 2"
        plasma = 'plasma' in sample.lower()
        serum = sample.lower() == 'serum'
        cocktail = sample.startswith(('P700', 'P800'))
        records.append(dict(
            measurement_id="yi2015:t2:%02d" % number, study_id="yi2015",
            source_location="Table 2, data row %d" % number,
            substrate=substrate, sample=sample,
            species="human", matrix="plasma" if plasma else "serum" if serum else "whole blood",
            temperature_label=temperature,
            temperature_celsius=24 if temperature == 'RT' else 4,
            temperature_source="Table 2 footnote: +/-1 C; methods say 25+/-1 C and text says 24+/-2 C",
            methods=method, reported_half_life_hours=reported,
            half_life_kind=kind, half_life_point_hours=point,
            half_life_lower_hours=lower, half_life_upper_hours=upper,
            replicate_sd_hours=sd, uncertainty=uncertainty,
            perturbation="proprietary inhibitor cocktail" if cocktail else "no added inhibitor cocktail",
            collection_additive=sample.split()[0] if not serum else "clotted serum",
            inhibitor_identity=sample.split()[0] if cocktail else None,
            inhibitor_dose=None, enzyme_attribution="unresolved multi-enzyme activity",
            observable="MS precursor signal" if method == 'MS' else
                       "Active GLP1 antibody signal" if method == 'ELISA' else
                       "MS precursor / active GLP1 antibody signals combined in source table",
            concentration_um=0.4 if method == 'MS' else 0.0004 if method == 'ELISA' else None,
            concentration_note="MS approximately 0.4 uM; antibody assay 400 pM; combined table rows do not resolve assay-specific half-lives",
            matrix_fraction=10 / 11 if method == 'MS' else None,
            ph=None, enzyme_activation="not independently measured",
            chemistry={
                'G36A': 'GLP1(7-36), C-terminal amide; isotope-labelled control is separate',
                'G37': 'GLP1(7-37); isotope-labelled control is separate',
                'GIP(1–42)': 'Synthetic GIP1-42; see source Table1; terminal chemistry not independently verified',
                'OXM(1–37)': 'OXM-K33 in methods; generic OXM1-37 in Table1/2; identity ambiguity retained',
                'Glucagon': 'Synthetic glucagon; isotope-labelled control is separate',
            }[substrate],
            lod=None, observed_horizon_hours=None, study_maximum_horizon_hours=96,
            raw_sample_id=None, product_ids=None,
            identifiable_enzyme_hazards=False, vaccine_transfer_validated=False,
        ))
    return records


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compact_route_guide(output):
    """Render one page centered on the selected target and productive loading."""
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('BodyCompact', fontName='Helvetica', fontSize=10,
                             leading=14, textColor=colors.HexColor('#233347')))
    styles.add(ParagraphStyle('CellCompact', fontName='Helvetica', fontSize=9,
                             leading=12, textColor=colors.HexColor('#233347')))
    def para(text, style='BodyCompact'):
        return Paragraph(text, styles[style])
    rows = [
        ['Route', 'What must happen', 'Main concern / useful readout'],
        ['mRNA expression',
         '<b>Translate the construct, then generate and load the target.</b><br/>Cytosolic antigen: proteasome, TAP and ER trimming for MHC-I. Secreted or lysosome-targeted constructs use different routes.',
         'Readout: target detected on MHC, plus expression yield.<br/><b>Watch:</b> cuts inside the target, wrong C terminus, over-trimming, linker/junction effects.<br/>Expression yield needs RNA and delivery evidence.'],
        ['Peptide uptake by APCs',
         '<b>Reach an APC and release a loadable target.</b><br/>Endolysosomal processing supports MHC-II; cytosolic export can enable MHC-I cross-presentation.',
         'Readout: target-specific MHC presentation after uptake.<br/><b>Watch:</b> destructive core cuts before loading or failure to access the required compartment.<br/>Uptake and processing need separate evidence.'],
        ['Serum exposure',
         '<b>Keep the target intact until uptake.</b><br/>The target can survive inside a shorter fragment even after its parent disappears.',
         'Readout: % retaining the exact target at the intended uptake time; parent half-life alongside it.<br/><b>Watch:</b> core cleavage or trimming into the target.<br/>Serum digestion alone does not determine circulation PK.'],
    ]
    table = Table([[para(value, 'CellCompact') for value in row] for row in rows],
                  colWidths=[94, 209, 225])
    table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#e5edf5')),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('LEFTPADDING', (0,0), (-1,-1), 10), ('RIGHTPADDING', (0,0), (-1,-1), 10),
        ('TOPPADDING', (0,0), (-1,-1), 10), ('BOTTOMPADDING', (0,0), (-1,-1), 10),
        ('LINEBELOW', (0,0), (-1,-1), .5, colors.HexColor('#d1dae5')),
    ]))
    story = [para('Will the vaccine target reach MHC?', 'Title'),
             para('<b>Track the target epitope, not just the full peptide.</b> '
                  'Cuts outside it can help processing; cuts through it remove that exact target. '
                  'Successful presentation also requires compatible HLA and MHC loading.'),
             Spacer(1,14), table, Spacer(1,16),
             para('<b>The smallest useful summary per target</b>'),
             para('Gene + mutation | exact target + HLA | route | internal-cut flags | '
                  'release/terminal evidence | target survival or unknown | evidence scope'),
             Spacer(1,12),
             para('<b>What the available evidence can tell us</b>'),
             para('Processing predictors provide site or ligand evidence. Native scores are not '
                  'a survival percentage. PeptiVerse and Cavaco parent half-life estimates stay '
                  'side by side. Target-survival trajectories and enzyme-removal effects require '
                  'fragment-specific rates or clearly named assumptions.'),
             Spacer(1,10),
             para('Human inhibitor data support specific hormone/protein mechanisms. FAP depletion '
                  'supports FGF21 cleavage; ordinary GP reporters can also be cut by PREP. '
                  'Inhibitor cocktails protect several hormones but cannot assign individual enzyme '
                  'rates. A transferable human long-vaccine-peptide calibration remains unavailable '
                  'in the audited evidence.'),
             Spacer(1,14),
             para('Sources: <a href="https://doi.org/10.1371/journal.pone.0089897">human DC cross-presentation</a>; '
                  '<a href="https://doi.org/10.1038/ni860">ERAP1 generation/destruction</a>; '
                  '<a href="https://doi.org/10.1371/journal.pone.0134427">human hormone stability</a>; '
                  '<a href="https://doi.org/10.1042/BJ20151085">FAP in human plasma</a>; '
                  '<a href="https://doi.org/10.1038/s41598-017-12900-8">FAP/PREP reporter selectivity</a>.', 'CellCompact')]
    def footer(canvas, doc):
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#65758a'))
        canvas.drawString(42, 25, 'mhctools %s | Evidence and assumptions stay separate | %s' % (__version__, doc.page))
    SimpleDocTemplate(str(output), pagesize=(612,792), leftMargin=42, rightMargin=42,
                      topMargin=35, bottomMargin=40, title='Vaccine target survival: three routes').build(
                          story, onFirstPage=footer, onLaterPages=footer)


def sid_target_rows(source):
    """Summarize a checksummed existing report without rerunning predictors."""
    manifest = json.loads(source.with_name('SHA256SUMS.json').read_text())
    entry = manifest['files'][source.name]
    if (hashlib.sha256(source.read_bytes()).hexdigest() != entry['sha256']
            or source.stat().st_size != entry['bytes']):
        raise ValueError('Focused report checksum mismatch')
    data = json.loads(source.read_text())
    rows = []
    for item in data['primary']:
        target = item['target']
        # Frozen target coordinates are zero-based, half-open. A physical
        # internal bond b is strictly between these boundaries.
        start, end = int(target['start']), int(target['end'])
        models = ('pepsickle-in-vivo-human-only', 'netchop-3.1-cterm-3.0') if target['mhc_class'] == 'I' else ('netcleave-ii-hla',)
        flags = []
        for model in models:
            sites = [s for s in item['processing_sites'] if s['model'] == model and start < s['bond'] < end]
            count = sum(bool(s['display_flag']) for s in sites)
            flags.append('%s: %d internal flags (%d assessed bonds)' % (model, count, len(sites))
                         if sites else '%s: unassessed inside target' % model)
        rows.append(dict(
            gene=target['gene'], mutation=item['protein_change'], target=target['target_sequence'],
            mhc_class=target['mhc_class'], allele=target['allele'],
            native_mhc_percentile_rank=target['percentile_rank'],
            mrna='Unassessed: actual encoded construct, junctions, RNA delivery and translation yield required',
            apc='Site evidence only; uptake/loading unassessed. ' + '; '.join(flags),
            serum='Exact-target survival uncalibrated; parent disappearance does not imply target destruction',
            peptiverse_parent_hours=item['estimates']['PeptiVerse']['parent'],
            cavaco_parent_hours=item['estimates']['Cavaco']['parent'],
            interpretation='Native site or ligand-terminus flags are model-specific evidence; zero flags do not establish protection',
            sequence_record_id=target['sequence_record_id'],
        ))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source-xml', type=Path, help='Audit Yi 2015 Table 2 against publisher XML')
    parser.add_argument('--pdf', action='store_true', help='Also create the one-page route guide')
    parser.add_argument('--sid-input', type=Path, help='Export compact rows from a checksummed focused report')
    args = parser.parse_args()
    evidence = serum_calibration_evidence()
    if args.source_xml:
        if measurement_records(reported_hormone_table(args.source_xml)) != evidence['measurements']:
            raise ValueError('Packaged measurements differ from primary Table 2')
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'evidence-inventory.json').write_text(json.dumps(evidence, indent=2) + '\n')
    write_csv(args.output / 'reported-measurements.csv', evidence['measurements'])
    references = [serum_assay_parent_reference(r['measurement_id'], [0, 1, 4, 24, 96])
                  for r in evidence['measurements']]
    (args.output / 'assay-parent-references.json').write_text(json.dumps(references, indent=2, allow_nan=False) + '\n')
    if args.sid_input:
        rows = sid_target_rows(args.sid_input)
        write_csv(args.output / 'sid-target-summary.csv', rows)
        (args.output / 'sid-target-summary.json').write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')
        (args.output / 'sid-source.json').write_text(json.dumps(dict(
            source=str(args.sid_input.resolve()), source_sha256=hashlib.sha256(args.sid_input.read_bytes()).hexdigest(),
            predictor_inference_repeated=False, new_vaccine_rate_fit=False,
            report_version=__version__), indent=2) + '\n')
    if args.pdf:
        compact_route_guide(args.output / 'vaccine-target-three-routes.pdf')
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(args.output.iterdir()) if p.is_file()}
    (args.output / 'SHA256SUMS.json').write_text(json.dumps(hashes, indent=2) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
