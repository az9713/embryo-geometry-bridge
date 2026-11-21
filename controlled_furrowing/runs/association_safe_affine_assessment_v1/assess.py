"""Replay a frozen rounded-affine LP with and without component channels."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SMALL = 'runs/material_cohort_descent_v1/results.json'
LARGE = 'runs/transfer_ambiguity_linear_v1/actual/results.json'
REVIEW = 'runs/transfer_ambiguity_linear_v1/actual_independent/receipt.json'
SOLVER = 'scripts/transfer_ambiguity_linear.py'
CHANNELS = ('maximum_shortening_um', 'region_averaged_shortening_um',
            'apical_max_shortening_um', 'basal_max_shortening_um')
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
SIGNS = {'max_contrast_change': 1, 'max_negative_contrast_change': -1}


def require(check, message):
    if not check:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, obj):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(obj, stream, indent=2, allow_nan=False)
        stream.write('\n')


def freeze_inputs():
    protocol = read(OUT / 'protocol.json')
    pins = protocol['study_relative_sha256']
    require(digest(Path(__file__)) == protocol['source_sha256'], 'assessment source mismatch')
    require(all(digest(ROOT / path) == value for path, value in pins.items()), 'study input hash mismatch')
    small, large, receipt = (read(ROOT / path) for path in (SMALL, LARGE, REVIEW))
    require(receipt['status'] == 'PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC'
            and all(receipt['checks'].values())
            and {p: pins[p] for p in (SMALL, LARGE, SOLVER)}
            == {p: receipt['bindings'][p] for p in (SMALL, LARGE, SOLVER)}
            and large['input_hashes'][SMALL] == pins[SMALL], 'original tangent authority mismatch')
    require(large['columns'] == protocol['parameter_columns'], 'parameter column order mismatch')
    require(protocol['allowance_um_exact'] == '1/1000', 'expected prior fixed allowance')
    baseline = [row for row in large['rows'] if row['mode'] == 'small_prediction_change'
                and Q(row['budget_exact']) == Q(protocol['allowance_um_exact'])]
    require(len(baseline) == 1 and Q(baseline[0]['allowance_um_exact']) == Q(1, 1000),
            'missing original small-prediction-change row')
    order = [{'record_id': 'optoGEF2_small_col_' + letter, 'observable': channel}
             for letter in 'BCDEFGHIJKLMNOPQR' for channel in CHANNELS]
    subset = {'source_joined34': [i for i, r in enumerate(order) if r['observable'] in CHANNELS[:2]],
              'ordered68': list(range(len(order)))}
    require(len(subset['source_joined34']) == 34 and len(order) == 68, 'invalid channel count')
    require(tuple(r['variant'] for r in small['candidates']) == VARIANTS, 'observer variant mismatch')
    matrices = {}
    for candidate in small['candidates']:
        affine = candidate['affine_LP']
        require(candidate['all68_constraint_order'] == order, 'same dose/channel order required')
        require(affine['trust_box_exact'] == protocol['trust_box_exact'], 'frozen trust box mismatch')
        matrix = affine['affine_derivative_coefficients_exact']
        require(len(matrix) == len(order) and all(len(row) == 5 for row in matrix), 'matrix dimensions')
        matrices[candidate['variant']] = matrix
    for first, second in zip(VARIANTS[::2], VARIANTS[1::2]):
        require(all(matrices[first][i] == matrices[second][i] for i in subset['source_joined34']),
                'component aggregation must not change the two total rows')
    return {'protocol_sha256': digest(OUT / 'protocol.json'), 'study_relative_sha256': pins,
            'allowance_um_exact': protocol['allowance_um_exact'], 'residual_offsets_exact': ['0'] * 68,
            'parameter_columns': large['columns'], 'trust_box_exact': protocol['trust_box_exact'],
            'observation_order': order, 'zero_based_row_indices': subset,
            'small_J_by_variant_exact': matrices, 'large_contrast_c_exact': large['rounded_affine_contrast_gradient_exact'],
            'original_disk_spatial_max_width_um_exact': baseline[0]['contrast_range_width_um_exact']}


def prove_endpoint(c, matrix, box, radius, witness):
    """Verify feasibility and the box-support upper bound from stored rationals."""
    c = list(map(Q, c))
    matrix = [tuple(map(Q, row)) for row in matrix]
    box = [tuple(map(Q, pair)) for pair in box]
    step = list(map(Q, witness['step_exact']))
    weights = list(map(Q, witness['multipliers_exact']))
    require(len(step) == len(c) == len(box) and len(weights) == 2 * len(matrix), 'witness dimensions')
    require(all(lo <= x <= hi for x, (lo, hi) in zip(step, box)), 'primal box infeasibility')
    require(all(abs(sum(x * y for x, y in zip(row, step))) <= radius for row in matrix),
            'primal observation infeasibility')
    require(all(w >= 0 for w in weights), 'negative dual multiplier')
    n = len(matrix)
    signed = [weights[i] - weights[i + n] for i in range(n)]
    reduced = [cj - sum(w * row[j] for w, row in zip(signed, matrix)) for j, cj in enumerate(c)]
    lower = sum(x * y for x, y in zip(c, step))
    upper = radius * sum(weights) + sum(max(v * lo, v * hi) for v, (lo, hi) in zip(reduced, box))
    require(0 <= lower <= upper, 'weak duality failure')
    require([lower, upper, upper - lower] == [Q(witness[key]) for key in
             ('primal_lower_exact', 'dual_upper_exact', 'duality_gap_exact')], 'witness scalar mismatch')
    return lower, upper


def evaluate(inputs, prior_results=None):
    producing = prior_results is None
    if producing:
        import numpy as np
        import scipy

        sys.path.insert(0, str(ROOT / 'scripts'))
        import transfer_ambiguity_linear as solver

        require(Path(solver.__file__).resolve() == ROOT / SOLVER, 'original LP implementation path')
    box = inputs['trust_box_exact']
    radius = Q(inputs['allowance_um_exact'])
    base_c = list(map(Q, inputs['large_contrast_c_exact']))
    interval = inputs['original_disk_spatial_max_width_um_exact']
    previous_low, previous_high = map(Q, interval)
    rows = []
    widths = {}
    for variant in VARIANTS:
        for group in ('source_joined34', 'ordered68'):
            indices = inputs['zero_based_row_indices'][group]
            matrix = [inputs['small_J_by_variant_exact'][variant][i] for i in indices]
            residual = [inputs['residual_offsets_exact'][i] for i in indices]
            require(set(residual) == {'0'}, 'diagnostic must use the baseline-change constraint')
            data = {'variant': variant, 'group': group, 'row_count': len(indices), 'endpoints': {}}
            prior = None if producing else prior_results['rows'][len(rows)]
            if prior is not None:
                require((prior['variant'], prior['group'], prior['row_count']) == (variant, group, len(indices)),
                        'stored result scope mismatch')
            bounds = []
            for name, sign in SIGNS.items():
                c = [sign * x for x in base_c]
                witness = (solver.max_affine(c, matrix, residual, radius, box) if producing
                           else prior['endpoints'][name])
                bounds.append(prove_endpoint(c, matrix, box, radius, witness))
                data['endpoints'][name] = witness
            low, high = (sum(end[j] for end in bounds) for j in range(2))
            data['width_um_exact'] = [str(low), str(high)]
            data['width_um_display'] = [float(low), float(high)]
            require(0 < low <= high and high - low < Q(2, 10**15), 'unresolved width enclosure')
            if variant == VARIANTS[0] and group == 'ordered68':
                require(max(previous_low, low) <= min(previous_high, high), 'archived LP disagrees')
            if producing:
                rounded = [[float(Q(x)) for x in row] for row in matrix]
                singular_values = np.linalg.svd(rounded, compute_uv=False)
                data['singular_values_floating_diagnostic'] = singular_values.tolist()
                data['numerical_rank_diagnostic'] = int(np.linalg.matrix_rank(rounded))
                require(data['numerical_rank_diagnostic'] == 5, 'unexpected diagnostic rank')
            else:
                require(prior['width_um_exact'] == data['width_um_exact']
                        and prior['width_um_display'] == data['width_um_display'], 'stored width mismatch')
            widths[variant, group] = (low, high)
            rows.append(data)
    growth = {}
    for variant in VARIANTS:
        fewer, more = (widths[variant, group] for group in ('source_joined34', 'ordered68'))
        ratio = [fewer[0] / more[1] - 1, fewer[1] / more[0] - 1]
        require(ratio[0] > 0, 'strict widening unresolved')
        growth[variant] = {'fractional_width_increase_exact': list(map(str, ratio)),
                           'percent_display': [float(100 * x) for x in ratio]}
    if not producing:
        require(len(prior_results['rows']) == len(rows) and prior_results['widening'] == growth,
                'stored comparison coverage mismatch')
        return None
    return {'status': 'COMPLETE_ROUNDED_AFFINE_CHANNEL_COMPARISON',
            'source_sha256': digest(Path(__file__)), 'input_snapshot_sha256': digest(OUT / 'inputs.json'),
            'protocol_sha256': inputs['protocol_sha256'], 'rows': rows, 'widening': growth,
            'numpy_version': np.__version__, 'scipy_version': scipy.__version__,
            'new_force_derivative_or_equilibrium_jobs': 0, 'larger_observation_data_read': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('run', 'verify'))
    mode = parser.parse_args().mode
    frozen = freeze_inputs()
    if mode == 'run':
        require(not any((OUT / name).exists() for name in ('inputs.json', 'results.json')), 'output exists')
        save(OUT / 'inputs.json', frozen)
        save(OUT / 'results.json', evaluate(frozen))
    results = read(OUT / 'results.json')
    require(read(OUT / 'inputs.json') == frozen and results['source_sha256'] == digest(Path(__file__))
            and results['input_snapshot_sha256'] == digest(OUT / 'inputs.json')
            and results['protocol_sha256'] == frozen['protocol_sha256'], 'input and result binding mismatch')
    evaluate(frozen, results)
    print(json.dumps({'verified_endpoint_certificates': 16, 'verified_widths': 8,
                      'results_sha256': digest(OUT / 'results.json')}))


if __name__ == '__main__':
    main()
