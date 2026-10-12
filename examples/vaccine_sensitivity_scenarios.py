"""Run explicit synthetic loss-rate scenarios for SLP and mRNA core delivery.

These are unit-scale numerical experiments, not estimated physiological ranges.
Use --input for an exported candidate request, including Vaxrank request.json.
All other profiles in that request are replaced and retained in input_request.
"""

import argparse
from copy import deepcopy
import csv
import json
from pathlib import Path

from mhctools import assess_vaccine_windows, format_vaccine_window_report, vaccine_route_steps, vaccine_window_rows

from vaccine_window_scenario import example_request


SOURCE = 'Synthetic sensitivity assumption; no human or animal kinetic calibration'


def synthetic_kinetics(candidate, loss):
    raw = candidate['peptide']
    sequence = raw['sequence'] if isinstance(raw, dict) else raw
    positions = {0, len(sequence), candidate.get('signal_end', 0)}
    for core in candidate['cores']:
        positions.update((core['start'], core['end']))
    positions = sorted(positions)
    compartments = ('interstitium', 'afferent_lymph', 'node_fluid', 'blood',
                    'producer_cytosol', 'producer_er', 'producer_secretory_er') + tuple(
        region + '_' + phase for region in ('local', 'node', 'systemic')
        for phase in ('endosome', 'cytosol', 'er', 'secretory_er'))
    assessments = []
    for start in positions:
        for end in positions:
            if end <= start:
                continue
            for compartment in compartments:
                # Explicit fragment cascade at target boundaries; ER trimming
                # excluded here. This does not predict protease recognition.
                bonds = [] if compartment.endswith('_er') else [p-start for p in positions if start < p < end]
                assessments.append({'compartment': compartment, 'start': start, 'end': end,
                                    'sequence': sequence[start:end], 'basis': 'assumed', 'source': SOURCE,
                                    'channels': [{'bond': bond, 'enzyme': 'synthetic boundary cleavage',
                                                  'rate_per_hour': 1, 'basis': 'assumed', 'source': SOURCE} for bond in bonds]})
    profiles = {}
    for mhc_class in sorted({core['mhc_class'] for core in candidate['cores']}):
        steps = vaccine_route_steps(candidate['delivery'], mhc_class)
        rates = {step.name: {'value_per_hour': 0, 'basis': 'disabled', 'source': SOURCE + '; explicit unused route'} for step in steps}
        selected = {'lymph_entry': .4, 'vascular_absorption': .2, 'lymph_transit': 1, 'node_outflow': .3,
                    'carrier_drainage': .2, 'producer_antigen_release': .1, 'cross_dressing_node': .1,
                    'cross_dressing_local': .1}
        for region in ('producer', 'local', 'node', 'systemic'):
            selected.update({region + '_' + key: value for key, value in {
                'antigen_uptake': .5, 'cross_escape': .2, 'tap': .7, 'er_loading': 1,
                'ii_loading': 1, 'surface_export': 1, 'autophagy': .2,
                'rna_uptake': .5, 'rna_escape': .3, 'rna_decay': .2, 'translation': 10,
                'signal_entry': 1, 'secretion': 1, 'failed_entry': .1, 'erad': .1,
                'secretory_to_endosome': .2}.items()})
        for step in steps:
            value = selected.get(step.name)
            if step.name.startswith('migration_'):
                value = .1
            if step.name.endswith(('_loss', '_clearance')):
                value = loss
            if value is not None:
                rates[step.name] = {'value_per_hour': value, 'basis': 'assumed', 'source': SOURCE}
        profiles[mhc_class] = {'rates': rates, 'cleavage_assessments': assessments,
                              'tap_max_length' if mhc_class == 'I' else 'class_ii_max_length': len(sequence)}
    return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', help='Optional candidate request; omit for built-in synthetic SLP and RNA examples')
    parser.add_argument('--output-dir', required=True, help='New directory')
    args = parser.parse_args()
    output = Path(args.output_dir)
    if output.exists():
        parser.error('Output directory already exists')
    request = json.loads(Path(args.input).read_text()) if args.input else example_request()
    if not args.input:
        for delivery in ('cytosolic_mrna', 'secreted_mrna'):
            candidate = deepcopy(request['candidates'][0])
            candidate.update(id=delivery, delivery=delivery, dose={
                'amount': 1, 'unit': 'ug', 'molecular_weight_g_per_mol': 1e6,
                'basis': 'assumed', 'source': SOURCE + '; hypothetical transcript mass'})
            if delivery == 'secreted_mrna':
                candidate['peptide'] = 'MAAAAA' + candidate['peptide']
                candidate['signal_end'] = 6  # Assumed boundary; not a predicted signal peptide.
                for core in candidate['cores']:
                    core['start'] += 6
                    core['end'] += 6
            request['candidates'].append(candidate)
    output.mkdir(parents=True)
    for loss in (.1, 1., 10.):
        scenario = deepcopy(request)
        scenario['scenario'] = 'Synthetic sensitivity: every unproductive loss/clearance rate = %g/h; NOT a physiological range' % loss
        scenario.pop('kinetics', None)
        for candidate in scenario['candidates']:
            candidate['kinetics'] = synthetic_kinetics(candidate, loss)
        result = assess_vaccine_windows(scenario)
        name = 'loss-%g' % loss
        (output / (name + '.json')).write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        (output / (name + '.md')).write_text(format_vaccine_window_report(result))
        rows = vaccine_window_rows(result)
        with (output / (name + '.csv')).open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else ['candidate_id', 'status'])
            writer.writeheader()
            writer.writerows(rows)
    print(output)


if __name__ == '__main__':
    main()
