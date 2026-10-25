#!/usr/bin/env python3
"""Check control coverage, workload and contrast algebra, not biological success."""

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime
from itertools import product
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / 'docs/协作/报告研究条件核验.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    RECORD.write_text(json.dumps({'status': 'running'}) + '\n')
    source_paths = [
        '设计/research-v3/condition-plan.tsv',
        '设计/research-v3/condition-plan.json',
        '报告/sections/research-extension.tex',
        'scripts/B-audit/check_research_condition_plan.py',
    ]
    inputs = {p: sha(ROOT / p) for p in source_paths}
    with (ROOT / source_paths[0]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    plan = json.loads((ROOT / source_paths[1]).read_text())
    assert len(rows) == len({r['condition_id'] for r in rows}) == 84
    assert {r['campaign'] for r in rows} == {'primary', 'confirmation'}
    checks = {'unique_condition_ids_and_campaigns': True}
    stimuli = ['vehicle', 'Torin1', 'MK8722']
    architectures = ['linker_EGFP', 'untagged']
    counts = {}
    for campaign, reagent, genotypes in [
        ('primary', 'TFEB_siRNA_1', ['WT', '3A', '3T']),
        ('confirmation', 'TFEB_siRNA_2', ['WT', '3A']),
    ]:
        expected = set()
        for g, a, s in product(genotypes, architectures, stimuli):
            expected.add(('research', reagent, g, a, s))
        for g, a, s in product(['WT', '3A'], architectures, stimuli):
            expected.add(('background_control', 'NT', g, a, s))
        for si, s in product(['NT', reagent], stimuli):
            expected.add(('no_rescue', si, 'none', 'nonexpressing_carrier', s))
            expected.add(('EGFP_control', si, 'none', 'EGFP_only', s))
        for s in stimuli:
            expected.add(('reagent_mock', 'none', 'none', 'none', s))
        current = [r for r in rows if r['campaign'] == campaign]
        actual = [tuple(r[k] for k in ('arm', 'tfeb_siRNA', 'genotype', 'architecture', 'stimulus')) for r in current]
        assert set(actual) == expected, campaign
        assert len(actual) == len(set(actual)), campaign
        counts[campaign] = len(actual)
    assert counts == {'primary': 45, 'confirmation': 39}
    checks['all_research_and_required_control_arms_complete_in_both_campaigns'] = True

    # Expand the written difference-of-differences algebra independently.
    weights = Counter()
    for a, s, sign in [('linker_EGFP', 'MK8722', 1), ('untagged', 'MK8722', -1),
                       ('linker_EGFP', 'Torin1', -1), ('untagged', 'Torin1', 1)]:
        for g, stimulus, coefficient in [('3A', s, 1), ('3A', 'vehicle', -1),
                                         ('WT', s, -1), ('WT', 'vehicle', 1)]:
            weights[(g, a, stimulus)] += sign * coefficient
    nonzero = {k: v for k, v in weights.items() if v}
    assert len(nonzero) == 8 and all(abs(v) == 1 for v in nonzero.values())
    assert sum(nonzero.values()) == 0
    assert all(k[2] != 'vehicle' for k in nonzero)
    # C must equal A_M - A_T with A_s = (3A-WT in E) - (3A-WT in U) at stimulus s.
    combined = Counter()
    for s, outer in (('MK8722', 1), ('Torin1', -1)):
        for g, gsign in (('3A', 1), ('WT', -1)):
            combined[(g, 'linker_EGFP', s)] += outer * gsign
            combined[(g, 'untagged', s)] -= outer * gsign
    assert {k: v for k, v in combined.items() if v} == nonzero
    for campaign in counts:
        available = {(r['genotype'], r['architecture'], r['stimulus']) for r in rows
                     if r['campaign'] == campaign and r['arm'] == 'research'}
        assert set(weights) <= available
    checks['primary_contrast_available_in_both_campaigns_and_vehicle_cancels'] = True

    n = plan['minimum_confirmatory_independent_days_per_campaign']
    wells = plan['culture_wells_per_condition_per_day']
    assert (n, wells) == (6, 2)
    for campaign in counts:
        assert plan[campaign + '_conditions'] == counts[campaign]
        assert plan[campaign + '_minimum_wells'] == counts[campaign] * n * wells
    assert sum(counts.values()) * n * wells == plan['combined_minimum_wells_per_assay_endpoint'] == 1008
    checks['culture_well_counts_include_controls_and_two_separate_campaigns'] = True
    r06 = plan['R06']
    assert r06['added_conditions_per_day'] == len(['WT', '3A']) * len(stimuli) == 6
    assert counts['primary'] + r06['added_conditions_per_day'] == r06['conditions_per_day'] == 51
    assert r06['conditions_per_day'] * n * wells == r06['minimum_wells_per_assay_endpoint'] == 612
    checks['R06_added_conditions_and_wells_and_C_equals_A_M_minus_A_T'] = True
    for p, h in inputs.items():
        assert sha(ROOT / p) == h, p
    result = {'status': 'passed', 'checked_at': datetime.now(ZoneInfo('America/New_York')).isoformat(),
              'checks': checks, 'inputs': inputs, 'conditions': counts,
              'minimum_wells_per_assay_endpoint': 1008,
              'primary_nonzero_coefficients': [{'genotype': k[0], 'architecture': k[1],
                                                'stimulus': k[2], 'weight': v} for k, v in sorted(nonzero.items())],
              'limitations': ['Counts and algebra only; no power calculation or experimental validation',
                              'Does not evaluate whether prose controls will work in cells']}
    RECORD.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'checks': len(checks), 'conditions': counts}))


if __name__ == '__main__':
    main()
