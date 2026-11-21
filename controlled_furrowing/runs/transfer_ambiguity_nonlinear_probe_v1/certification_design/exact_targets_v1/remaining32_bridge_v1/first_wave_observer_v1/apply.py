"""Scalar observer transfer for the four independently reviewed first-wave roots."""

import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path


OUT = Path(__file__).resolve().parent
BRIDGE = OUT.parent
EXACT = BRIDGE.parent
DESIGN = EXACT.parent
PROBE = DESIGN.parent
ROOT = PROBE.parents[1]
REPLAY = BRIDGE / 'first_wave_reduction_v1/independent_actual_preparation'
NATIVE = REPLAY / 'actual_replay'
PRIMARY = PROBE / 'independent_actual_preparation/actual_replay/receipt.json'
CACHED = DESIGN / 'full_outcome_coverage_v1/observer_application_v1'
BUDGET = DESIGN / 'budgets_v1'
SECONDARY = PROBE / 'observer_robustness'
GEOMETRY = EXACT / 'remaining32_geometry_v1/preparation/protocol.json'
PRODUCER = BRIDGE / 'first_wave_reduction_v1/actual/root_terminal_admission.json'
LENGTH = Q(67, 2)
RADIUS = Q(1, 2**34)
CASES = ('max_contrast_change_small_01', 'max_contrast_change_small_13',
         'max_negative_contrast_change_small_01', 'max_negative_contrast_change_small_16')
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
OBS = ('maximum_shortening_um', 'region_averaged_shortening_um',
       'apical_max_shortening_um', 'basal_max_shortening_um')
ROOT_STATUS = 'PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT'
SOURCE_REVIEW_STATUS = 'PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_OBSERVER_SOURCE'

PINS = {
    PRIMARY: 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    CACHED / 'results.json': '1df2bc372edd12687fad17edf63013509994e6bbe41c02150cefb25952eab9f6',
    CACHED / 'independent_readback/receipt.json': 'a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d',
    BUDGET / 'actual/results.json': '689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872',
    BUDGET / 'independent_actual/receipt.json': '978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c',
    SECONDARY / 'actual/results.json': '53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647',
    SECONDARY / 'independent_actual/receipt.json': '703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    GEOMETRY: '754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850',
    PRODUCER: '4c56643bdf0c4fdedf5f15f9d51ce89964e3f346f868eddff7ae6a41d3513452',
    REPLAY / 'replay_frozen.py': '4682165bf59bd633092c042a7df619c9bc1ad272dc80d1755404c2999a65287f',
    REPLAY / 'protocol.json': '9ca811a3e7f9e35346dd28d040ff3f370599e014ce3df942c48d347f6bd2f24a',
    REPLAY / 'preflight.json': '90242d2cde1a81c075cf342f716cf8a229f45b5ce3f73778fdd9db6ade249d91',
    REPLAY / 'root_source_review.json': '12e14fb46070d40694c09f3b2a3d76492b1dc85659c7a61c17df6b617f6915e0',
    ROOT / 'docs/NODAL-OBSERVER-ERROR-TRANSFER.md':
        'd0d5f711a6a8f2ebe46e8dcb76c70234753bda714836dbcf7d4be888ff7ba171',
    ROOT / 'evidence/nodal_observer_error_independent_review.json':
        '335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a',
}
OS_PIN = 'aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path, expected, hashes, parse=True):
    path = ROOT / path
    key = rel(path)
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == expected and hashes.setdefault(key, expected) == expected,
            'changed hash or conflicting authority: ' + key)
    if parse:
        return json.loads(raw)


def save(path, obj):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')
    return sha(path)


def exact(x):
    require(type(x) is str, 'expected exact rational string')
    return Q(x)


def pair(seq):
    lo, hi = map(exact, seq)
    require(lo <= hi, 'reversed interval')
    return lo, hi


def literal(x, expected):
    return type(x) is int and x == expected


def pick(rows, **match):
    [row] = [r for r in rows if all(r[k] == v for k, v in match.items())]
    return row


def native_summary(os_row, receipt, terminal, phase):
    require(os_row['status'] == 'PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT'
            and literal(os_row['pid'], 20640) and literal(os_row['unified_exec_session'], 78089)
            and literal(os_row['observed_exit_code'], 0) and os_row['direct_root_tool_observation'] is True
            and os_row['all_four_native_cases_and_672_component_hull_links_rehashed'] is True
            and os_row['full_five_gate_complete'] is False, 'observed successful independent replay required')
    require(terminal['status'] == 'FINISHED' and literal(terminal['exit_code'], 0)
            and terminal['receipt_sha256'] == os_row['receipt_sha256']
            and receipt['status'] == 'COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY'
            and literal(receipt['exact_spatial_roots'], 4) and receipt['archived_observers_replaced'] is False
            and receipt['full_five_gate_complete'] is False
            and all(literal(d['pid'], os_row['pid']) and all(d[k] == v for k, v in phase.items())
                    for d in (receipt, terminal)), 'phase and true successful terminal mismatch')
    rows = receipt['case_results']
    require(tuple(x['case_id'] for x in rows) == CASES, 'complete ordered four-case native batch')
    for row in rows:
        require(row['native_file'] == rel(NATIVE / (row['case_id'] + '.json'))
                and row['status'] == ROOT_STATUS and row['observer_state_budget_proved'] is True,
                'wrong native path or insufficient root budget')
    require(os_row['native_cases_sha256'] == {x['native_file']: x['native_sha256'] for x in rows},
            'observed and native case closures disagree')
    return rows


def native_index(hashes, os_sha):
    require(os_sha == OS_PIN, 'same root-observed native closure required')
    observed = read(NATIVE / 'os_exit_observation.json', os_sha, hashes)
    root_review = read(REPLAY / 'root_source_review.json', PINS[REPLAY / 'root_source_review.json'], hashes)
    require(root_review['status'] == 'PASS_ROOT_FIRST_WAVE_INDEPENDENT_GLOBAL_ENDPOINT_SOURCE'
            and root_review['checks'] and all(v is True for v in root_review['checks'].values()),
            'wrong native source approval')
    phase = {k: root_review[k] for k in ('source_sha256', 'protocol_sha256', 'preflight_sha256',
             'producer_root_terminal_admission_sha256', 'producer_reduction_input_admission_sha256')}
    phase['root_source_review_sha256'] = PINS[REPLAY / 'root_source_review.json']
    require(all(observed[k] == phase[k] for k in ('root_source_review_sha256',
            'producer_root_terminal_admission_sha256', 'producer_reduction_input_admission_sha256'))
            and phase['source_sha256'] == PINS[REPLAY / 'replay_frozen.py']
            and phase['protocol_sha256'] == PINS[REPLAY / 'protocol.json']
            and phase['preflight_sha256'] == PINS[REPLAY / 'preflight.json']
            and phase['producer_root_terminal_admission_sha256'] == PINS[PRODUCER], 'native phase pins')
    term = read(NATIVE / 'terminal.json', observed['terminal_sha256'], hashes)
    receipt = read(NATIVE / 'receipt.json', observed['receipt_sha256'], hashes)
    rows = native_summary(observed, receipt, term, phase)
    return rows, {**phase, 'pid': observed['pid']}


def plan(os_sha):
    hashes = {}
    for path, pin in PINS.items():
        read(path, pin, hashes, parse=False)
    rows, phase = native_index(hashes, os_sha)
    return {'status': 'FROZEN_FIRST_WAVE_EXACT_SPATIAL_OBSERVER_TRANSFER',
            'source_sha256': sha(Path(__file__)), 'case_order': list(CASES),
            'length_um_exact': str(LENGTH), 'radius_exact': str(RADIUS),
            'native_os_sha256': os_sha, 'native_phase': phase, 'native_cases': rows,
            'static_and_native_metadata_sha256': hashes, 'variants': list(VARIANTS),
            'observables': list(OBS), 'independent_source_gate_required': True,
            'scope': 'Four small fields only. No force, AD, spatial integration, refit or new observer.'}


def case_inputs(case, hashes, docs, phase, leaf):
    name = case['case_id']
    endpoint, number = name.rsplit('_small_', 1)
    checkpoint = 'small_' + number
    trial = pick(docs[GEOMETRY]['trials'], case_id=name)
    contract = read(trial['model_contract'], trial['model_contract_sha256'], hashes)
    target = leaf['target']
    keys = ('trial_box_sha256', 'model_contract_sha256', 'point_cache_sha256', 'R_binary_sha256',
            'center_binary_sha256', 'archived_center_shift_exact', 'exact_drive', 'energy_scale_binary_exact')
    require(leaf['status'] == ROOT_STATUS and leaf['case_id'] == target['case_id'] == contract['pilot_id'] == name
            and literal(leaf['pid'], phase['pid']) and all(leaf[k] == v for k, v in phase.items())
            and target['record_id'] == trial['record_id'] == contract['record_id']
            and all(target[k] == trial[k] for k in keys)
            and leaf['original_observer_field_replaced'] is False
            and leaf['closed_cube_invertible_R_from_contraction'] is True
            and exact(leaf['radius_exact']) == exact(trial['trial_radius_exact']) == RADIUS,
            'same native first-wave physical target, closed box, center and R required')
    box = read(trial['trial_box'], trial['trial_box_sha256'], hashes)
    require(box['record_id'] == target['record_id']
            and box['model_observation_contract_sha256'] == target['model_contract_sha256']
            and box['point_cache_sha256'] == target['point_cache_sha256']
            and exact(box['corrector_radius_exact']) == RADIUS, 'trial box target mismatch')

    prior = pick(docs[CACHED / 'results.json']['cases'], case_id=name)
    budget_doc = docs[BUDGET / 'actual/results.json']
    budget = pick(budget_doc['field_budgets'], endpoint_id=endpoint, checkpoint=checkpoint)
    budget_audit = pick(docs[BUDGET / 'independent_actual/receipt.json']['fields'],
                        endpoint_id=endpoint, checkpoint=checkpoint)
    secondary = docs[SECONDARY / 'independent_actual/receipt.json']
    old_review = docs[CACHED / 'independent_readback/receipt.json']
    s_path = SECONDARY / 'actual' / endpoint / (checkpoint + '.json')
    s_pin = pick(docs[SECONDARY / 'actual/results.json']['field_observer_receipts'], path=rel(s_path))['sha256']
    saved = read(s_path, s_pin, hashes)
    primary = docs[PRIMARY]
    model = budget['center_and_model']
    require(saved['role'] == budget['role'] == 'small_record' and saved['endpoint_id'] == endpoint
            and saved['checkpoint'] == checkpoint and saved['step'] == budget['step']
            and secondary['input_and_output_sha256'][rel(s_path)] == s_pin
            and all(saved[k] == budget[k] for k in ('field', 'field_sha256', 'verification', 'verification_sha256',
                    'record_id'))
            and prior['record_id'] == budget_audit['record_id'] == budget['record_id'] == contract['physical_source_record_id']
            and prior['field_sha256'] == budget_audit['field_sha256'] == saved['field_sha256']
            == target['original_observer_field_sha256'] == contract['archived_field_sha256']
            and saved['field'] == contract['archived_field'] == budget_audit['field']
            and saved['verification'] == contract['archived_verification']
            and saved['verification_sha256'] == contract['archived_verification_sha256']
            and target['observer_goal_exact'] == budget_audit['e_goal_dimensionless_exact']
            == budget['e_goal_dimensionless_exact'] == contract['observer_goal_exact']
            and target['original_addition_gap_exact'] == budget_audit['binary_gap_exact']
            == model['node_addition_gap_dimensionless_exact'] == contract['archived_addition_gap_exact']
            and target['archived_center_shift_exact'] == contract['archived_shift_exact']
            and prior['root']['sigma'] == contract['archived_shift_exact']
            and contract['archived_observer_variant'] == VARIANTS[0]
            and exact(contract['archived_observer_allowance_um_exact'])
            == exact(budget_doc['conditional_source_allowance_um_exact']),
            'same immutable field, verification, original budget, and recenter shift required')
    drive = contract['exact_requested_drive']
    require(exact(drive['optical_dose_exact']) == exact(saved['step']['dose_exact']) == exact(box['dose_exact'])
            and exact(drive['common_gain_exact']) == exact(model['gain_exact'])
            and exact(drive['activation_strength_exact']) == exact(model['exact_target_strength'])
            == exact(target['exact_drive']) == exact(drive['common_gain_exact']) * exact(drive['optical_dose_exact']),
            'original exact dose and material gain required')
    for path, pin in ((saved['field'], saved['field_sha256']), (saved['verification'], saved['verification_sha256'])):
        require(primary['all_hashes'][path] == old_review['current_input_sha256'][path]
                == secondary['input_and_output_sha256'][path] == pin, 'archived field or verification source review mismatch')
    state = read(saved['field'], saved['field_sha256'], hashes)
    verification = read(saved['verification'], saved['verification_sha256'], hashes)
    return budget, prior, saved, state, verification, box


def binary(x):
    require(type(x) in (int, float), 'binary coordinate expected')
    return Q(x)


def edge_constants(points):
    """Bounds on dR/dt and |dZ/dR| for one Lagrange P2 edge."""
    first = tuple(-3 * a + 4 * b - c for a, b, c in zip(*points, strict=True))
    last = tuple(a - 4 * b + 3 * c for a, b, c in zip(*points, strict=True))
    radial = min(first[0], last[0])
    require(radial > 0, 'archived trace must have strictly positive radial derivative')
    return radial, max(abs(first[1]), abs(last[1])) / radial


def graph_bound(radial, slope, epsilon):
    require(epsilon >= 0 and slope >= 0 and radial > 8 * epsilon, 'monotone trace guard must be strict')
    return radial - 8 * epsilon, Q(5, 4) * (1 + slope) * epsilon


def geometry(budget, prior, saved, state, box, sigma, xi, epsilon):
    model = budget['center_and_model']
    X = [binary(v) for node in state['reference_nodes'] for v in node]
    archived = [binary(v) for node in state['current_nodes'] for v in node]
    rim = binary(state['parameters']['outer_radius'])
    fixed = {i for i in range(650) if X[i - i % 2] == rim or (i % 2 == 0 and X[i] == 0)}
    free = [i for i in range(650) if i not in fixed]
    require(len(X) == len(archived) == 650 and len(free) == len(state['free_displacement']) == 611
            and free == model['free_map'] == box['free'] and all(type(i) is int for i in box['free'])
            and state['parameters'] == model['material_parameters']
            and state['triangles'] == box['triangles'] and len(state['triangles']) == 144,
            'archived mesh and exact-cube free map differ')
    y = X.copy()
    for index, u in zip(free, state['free_displacement'], strict=True):
        y[index] += binary(u)
    require(max(abs(ya - a) for ya, a in zip(y, archived, strict=True)) == xi
            and all(Q.from_float(float(ya)) == a for ya, a in zip(y, archived, strict=True)),
            'original X+u addition gap changed')
    cubecenter = [exact(v) for node in box['current_centers_exact'] for v in node]
    radii = [exact(v) for node in box['coordinate_radii_exact'] for v in node]
    require([exact(v) for node in box['reference_nodes_exact'] for v in node] == X
            and max(abs(a-b) for a, b in zip(y, cubecenter, strict=True)) == sigma
            and len(cubecenter) == len(radii) == 650
            and all(radii[i] == (0 if i in fixed else RADIUS) for i in range(650))
            and all(cubecenter[i] == y[i] == archived[i] for i in fixed), 'true center shift and fixed coordinates')
    roi = pair(budget['fixed_roi_um_exact'])
    packed = saved['observer_contract']
    packed_roi = tuple(Q(v['numerator'], v['denominator']) for v in packed['current_radius_ROI_um'])
    az, bz = (packed['reference_' + side + '_um'] for side in ('apical', 'basal'))
    require(roi == pair(prior['ROI_um_exact']) == packed_roi
            and roi[0] == 0 < roi[1] == LENGTH * binary(state['parameters']['patch_radius']) <= LENGTH * rim
            and Q(az['numerator'], az['denominator']) == exact(model['reference_height_um_exact']) == LENGTH
            and Q(bz['numerator'], bz['denominator']) == 0, 'same ROI and flat reference required')
    graph = {}
    old_epsilon = exact(prior['root']['epsilon_um'])
    for face, side in (('top', 'apical'), ('bottom', 'basal')):
        trace = model['fixed_trace_endpoints'][face]
        edges = trace['trace_nodes']
        dofs = [2 * edges[0][0], 2 * edges[-1][-1]]
        require(len(edges) == 12 and all(a[-1] == b[0] for a, b in zip(edges, edges[1:]))
                and dofs == trace['radial_dofs'] and set(dofs) <= fixed
                and [cubecenter[i] for i in dofs] == [0, rim]
                and pair(trace['proved_fixed_radii_um_exact']) == (0, LENGTH * rim), 'fixed conforming trace span')
        controls = [[tuple(LENGTH * archived[2 * node + i] for i in (0, 1)) for node in edge] for edge in edges]
        stored_controls = state['boundary_' + face + '_quadratic_controls']
        require(controls == [[tuple(LENGTH * binary(v) for v in node) for node in edge]
                             for edge in stored_controls],
                'independently archived trace differs from reviewed geometry')
        entries = [edge_constants(edge) for edge in controls]
        radial = min(v[0] for v in entries)
        slope = max(v[1] for v in entries)
        old = prior['apical_and_basal_graph_error_bounds'][side]
        base = budget['stored_face_derivatives'][face]
        require(radial == exact(base['radial_min_um_exact']) == exact(old['radial_derivative_margin_exact'])+8*old_epsilon
                and slope == exact(base['global_stored_slope_cap_exact']) == exact(old['base_graph_Lipschitz_upper_exact'])
                and exact(old['Eulerian_vertical_error_um_exact']) == Q(5, 4)*(1 + slope)*old_epsilon,
                'rederived face coefficients mismatch original independent geometry')
        margin, error = graph_bound(radial, slope, epsilon)
        graph[side] = {'stored_radial_min_um_exact': str(radial), 'slope_cap_exact': str(slope),
                       'strict_radial_margin_um_exact': str(margin), 'graph_error_um_exact': str(error)}
    return roi, graph


def posterior(arithmetic, transferred, target):
    eta, z = (exact(arithmetic[k]) for k in ('eta_upper_exact', 'z_upper_exact'))
    lhs = eta + z * RADIUS
    require(eta >= 0 and 0 <= z < 1 and lhs < RADIUS and arithmetic['strict_selfmap'] is True
            and arithmetic['strict_cube_exclusion'] is False
            and lhs == exact(arithmetic['selfmap_lhs_upper_exact']), 'strict exact spatial self-map required')
    e = eta / (1-z)
    sigma, xi, goal = (exact(target[k]) for k in ('archived_center_shift_exact', 'original_addition_gap_exact',
                        'observer_goal_exact'))
    epsilon = LENGTH * (sigma + e + xi)
    require(0 <= e < RADIUS and e == exact(arithmetic['posterior_error_upper_exact'])
            and sigma >= 0 and xi >= 0 and sigma + e < goal
            and transferred['archived_state_budget_proved'] is True
            and sigma + e == exact(transferred['archived_shift_plus_posterior'])
            and epsilon == exact(transferred['physical_nodal_error_um_upper_exact']), 'original observer budget mismatch')
    return e, sigma, xi, epsilon


def compare(stored, sample, error, allowance):
    lo, hi = stored
    require(lo <= hi and error >= 0 and allowance >= 0, 'invalid residual inputs')
    root = (lo - error, hi + error)
    left, right = root[0] - sample, root[1] - sample
    lower, upper = max(Q(0), left, -right), max(abs(left), abs(right))
    status = 'IN' if upper <= allowance else 'OUT' if lower > allowance else 'UNRESOLVED'
    return root, (lower, upper), status


def apply_case(case, hashes, docs, phase, allowance):
    leaf = read(case['native_file'], case['native_sha256'], hashes)
    budget, prior, saved, state, verification, box = case_inputs(case, hashes, docs, phase, leaf)
    e, sigma, xi, epsilon = posterior(leaf['arithmetic'], leaf['observer_and_optional_minimum'], leaf['target'])
    roi, graph = geometry(budget, prior, saved, state, box, sigma, xi, epsilon)
    ea, eb = (exact(graph[s]['graph_error_um_exact']) for s in ('apical', 'basal'))
    errors = dict(zip(OBS, (ea+eb, ea+eb, ea, eb), strict=True))
    require(verification['variant'] == VARIANTS[0] and set(saved['variant_observables_um_exact']) == set(VARIANTS),
            'full frozen observer catalogue required')
    rows = []
    for variant in VARIANTS:
        values = saved['variant_observables_um_exact'][variant]
        require(set(values) == set(OBS), 'full frozen channel collection required')
        for obs in OBS:
            stored = pair(values[obs])
            sample = exact(budget['step']['observed_exact'][obs])
            if variant == VARIANTS[0]:
                frozen = pick(prior['comparisons'], observable=obs)
                require(stored == pair(frozen['stored_interval_um_exact']) == pair(verification['observer_intervals_um_exact'][obs])
                        and sample == exact(frozen['observed_um_exact']), 'primary source interval changed')
            wide, error, standing = compare(stored, sample, errors[obs], allowance)
            rows.append({'variant': variant, 'observable': obs, 'stored_interval_um_exact': list(map(str, stored)),
                         'observed_um_exact': str(sample), 'extra_observer_error_um_exact': str(errors[obs]),
                         'exact_spatial_root_observer_interval_um_exact': list(map(str, wide)),
                         'absolute_residual_um_exact': list(map(str, error)), 'standing': standing})
    return {'case_id': case['case_id'], 'root_proof_sha256': case['native_sha256'],
            'record_id': budget['record_id'], 'original_field_sha256': saved['field_sha256'],
            'posterior_exact': str(e), 'archived_sigma_exact': str(sigma), 'original_xi_exact': str(xi),
            'physical_nodal_error_um_exact': str(epsilon), 'ROI_um_exact': list(map(str, roi)),
            'graph_premises': graph, 'comparisons': rows}


def gate(review, protocol, preflight, binding):
    require(review['status'] == SOURCE_REVIEW_STATUS and review['checks']
            and all(v is True for v in review['checks'].values())
            and all(review[k] == v for k, v in binding.items())
            and review['permitted_case_ids'] == list(CASES)
            and protocol['source_sha256'] == preflight['source_sha256'] == binding['source_sha256']
            and protocol['native_os_sha256'] == binding['native_os_sha256']
            and preflight['protocol_sha256'] == binding['protocol_sha256']
            and preflight['status'] == 'PASS' and literal(preflight['actual_observer_applications'], 0)
            and preflight['checks'] and all(v is True for v in preflight['checks'].values()),
            'same frozen source and independent prelaunch review required')


def fixtures():
    def reject(action):
        try:
            action()
        except (ValueError, KeyError):
            return True
        raise ValueError('invalid fixture admitted')

    checks = {}
    checks['IN_closed_OUT_strict'] = (compare((-2, -1), 0, 0, 2)[2] == 'IN'
        and compare((2, 3), 0, 0, 2)[2] == 'UNRESOLVED' and compare((Q(9, 4), 3), 0, 0, 2)[2] == 'OUT')
    checks['signed_zero_crossing'] = compare((-1, 2), 0, Q(1, 2), 2) == (
        (Q(-3, 2), Q(5, 2)), (Q(0), Q(5, 2)), 'UNRESOLVED')
    checks['P2_derivative_and_graph'] = (edge_constants(((0, 0), (Q(1, 2), 1), (1, 2))) == (1, 2)
        and graph_bound(1, 2, Q(1, 16)) == (Q(1, 2), Q(15, 64))
        and reject(lambda: graph_bound(1, 2, Q(1, 8))))
    z, eta = Q(1, 2), RADIUS / 8
    transfer = {'archived_state_budget_proved': True, 'archived_shift_plus_posterior': str(RADIUS / 2),
                'physical_nodal_error_um_upper_exact': str(LENGTH * (RADIUS / 2 + RADIUS / 16))}
    target = {'archived_center_shift_exact': str(RADIUS / 4), 'original_addition_gap_exact': str(RADIUS / 16),
              'observer_goal_exact': str(RADIUS)}
    arith = {'eta_upper_exact': str(eta), 'z_upper_exact': str(z), 'strict_selfmap': True,
             'strict_cube_exclusion': False, 'selfmap_lhs_upper_exact': str(eta + z * RADIUS),
             'posterior_error_upper_exact': str(RADIUS / 4)}
    checks['posterior_includes_sigma_and_xi'] = (
        posterior(arith, transfer, target)[-1] == Q(transfer['physical_nodal_error_um_upper_exact'])
        and reject(lambda: posterior(arith, transfer, {**target, 'archived_center_shift_exact': '0'})))
    checks['typed_process_identifiers'] = not literal(False, 0) and not literal(20640.0, 20640)
    binding = {'source_sha256': 's', 'native_os_sha256': 'o', 'protocol_sha256': 'p', 'preflight_sha256': 'f'}
    proto = {'source_sha256': 's', 'native_os_sha256': 'o'}
    pre = {'source_sha256': 's', 'protocol_sha256': 'p', 'status': 'PASS',
           'actual_observer_applications': 0, 'checks': {'fixture': True}}
    rev = {**binding, 'status': SOURCE_REVIEW_STATUS, 'checks': {'fixture': True}, 'permitted_case_ids': list(CASES)}
    gate(rev, proto, pre, binding)
    checks['cold_phase_roundtrip_and_reject_mixed_hash'] = reject(lambda: gate(
        {**rev, 'native_os_sha256': 'wrong'}, proto, pre, binding))
    require(all(v is True for v in checks.values()), 'cold fixture failed')
    return checks


def prepare(os_sha):
    protocol = plan(os_sha)
    hashes = {}
    rows, phase = native_index(hashes, os_sha)
    docs, _ = archive_authority(hashes)
    for row in rows:
        leaf = read(row['native_file'], row['native_sha256'], hashes)
        _, _, saved, state, verification, box = case_inputs(row, hashes, docs, phase, leaf)
        require(all(len(box[k]) == 325 for k in ('reference_nodes_exact', 'current_centers_exact',
                'coordinate_radii_exact')) and len(state['current_nodes']) == 325
                and all(len(state['boundary_' + face + '_quadratic_controls']) == 12
                        for face in ('top', 'bottom'))
                and set(saved['variant_observables_um_exact']) == set(VARIANTS)
                and verification['variant'] == VARIANTS[0], 'complete archived metadata and scope')
    checks = {**fixtures(), 'all_four_real_native_and_archive_case_bindings': True}
    ph = save(OUT / 'protocol.json', protocol)
    pf = save(OUT / 'preflight.json', {'status': 'PASS', 'source_sha256': protocol['source_sha256'],
              'protocol_sha256': ph, 'checks': checks, 'cold_native_and_archive_reads_sha256': hashes,
              'actual_observer_applications': 0})
    print(json.dumps({'source_sha256': protocol['source_sha256'], 'protocol_sha256': ph, 'preflight_sha256': pf}))


def archive_authority(hashes):
    docs = {path: read(path, pin, hashes) for path, pin in PINS.items() if path.suffix == '.json'}
    ca, ci = (docs[CACHED / name] for name in ('results.json', 'independent_readback/receipt.json'))
    bd, bi = (docs[BUDGET / name] for name in ('actual/results.json', 'independent_actual/receipt.json'))
    sa, si = (docs[SECONDARY / name] for name in ('actual/results.json', 'independent_actual/receipt.json'))
    allowance = exact(ca['frozen_allowance_um_exact'])
    require(ci['status'] == 'PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS'
            and ci['current_input_sha256'][rel(CACHED / 'results.json')] == PINS[CACHED / 'results.json']
            and bi['status'] == 'PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS'
            and bi['source_and_input_hashes'][rel(BUDGET / 'actual/results.json')] == PINS[BUDGET / 'actual/results.json']
            and si['status'] == 'PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP'
            and si['input_and_output_sha256'][rel(SECONDARY / 'actual/results.json')] == PINS[SECONDARY / 'actual/results.json']
            and allowance > 0 and allowance == exact(bd['conditional_source_allowance_um_exact'])
            == exact(ci['frozen_model_discrepancy_um_exact']) == exact(si['source_configured_delta_um_exact'])
            == exact(sa['frozen_primary_source_allowance_um_exact']), 'original reviewed observer authorities required')
    return docs, allowance


def run(args):
    hashes = {}
    p = read(OUT / 'protocol.json', args.protocol_sha256, hashes)
    pf = read(OUT / 'preflight.json', args.preflight_sha256, hashes)
    review = read(OUT / 'independent_source_review/receipt.json', args.review_sha256, hashes)
    binding = {'source_sha256': sha(Path(__file__)), 'native_os_sha256': args.native_os_sha256,
               'protocol_sha256': args.protocol_sha256, 'preflight_sha256': args.preflight_sha256}
    gate(review, p, pf, binding)
    require(plan(args.native_os_sha256) == p, 'frozen plan changed')
    rows, phase = native_index(hashes, args.native_os_sha256)
    docs, allowance = archive_authority(hashes)
    dest = OUT / 'actual'
    dest.mkdir(exist_ok=False)
    mark = {**binding, 'source_review_sha256': args.review_sha256, 'pid': os.getpid()}
    save(dest / 'started.json', {'status': 'RUNNING', **mark})
    terminal = {'status': 'FAILED', 'exit_code': 1, **mark}
    try:
        result = [apply_case(case, hashes, docs, phase, allowance) for case in rows]
        require(all(hashes[key] == pin for key, pin in pf['cold_native_and_archive_reads_sha256'].items()),
                'changed cold native or archived read')
        counts = {v: dict(Counter(c['standing'] for r in result for c in r['comparisons'] if c['variant'] == v))
                  for v in VARIANTS}
        require(len(result) == 4 and all(sum(c.values()) == 16 for c in counts.values()), 'complete 64-comparison coverage')
        receipt_sha = save(dest / 'results.json', {'status': 'FIRST_WAVE_EXACT_SPATIAL_OBSERVERS_AWAIT_INDEPENDENT_READBACK',
                   **mark, 'hashes': hashes, 'allowance_um_exact': str(allowance), 'cases': result, 'counts': counts,
                   'total_primary_comparisons': 16, 'total_alternative_comparisons': 48,
                   'new_force_AD_or_integral_jobs': 0, 'full_calibration_or_source_fit_claimed': False})
        terminal.update(status='FINISHED', exit_code=0, results_sha256=receipt_sha)
        print(json.dumps({'results_sha256': receipt_sha, 'counts': counts}))
    except Exception as error:
        save(dest / 'failure.json', {'error': str(error), **mark})
        raise
    finally:
        save(dest / 'terminal.json', terminal)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    prep, actual = sub.add_parser('prepare'), sub.add_parser('run')
    for command in (prep, actual):
        command.add_argument('--native-os-sha256', required=True)
    for k in ('protocol', 'preflight', 'review'):
        actual.add_argument('--' + k + '-sha256', required=True)
    args = parser.parse_args()
    prepare(args.native_os_sha256) if args.mode == 'prepare' else run(args)
