"""Extract conditional observer error budgets for two frozen nonlinear branches."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import os
import traceback
from unittest.mock import patch

import nodal_observer_error_transfer as transfer
import transfer_ambiguity_observer_sweep as sweep


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
PRIMARY = sweep.PRIMARY
DESIGN = f'{PRIMARY}/certification_design'
OUT = ROOT / DESIGN / 'budgets_v1'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
PRIMARY_PROTOCOL = sweep.PRIMARY_PROT
NATIVE_REVIEW = f'{PRIMARY}/independent_actual_preparation/actual_replay/receipt.json'
MANIFEST = f'{PRIMARY}/field_manifest.json'
THEOREM_REVIEW = 'evidence/nodal_observer_error_independent_review.json'
T_SOURCE = 'scripts/nodal_observer_error_transfer.py'
S_SOURCE = 'scripts/transfer_ambiguity_observer_sweep.py'
PINS = {
    PRIMARY_PROTOCOL: '9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94',
    THEOREM_REVIEW: '335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a',
    T_SOURCE: '001769f23ba386fe6bde4530c12e011f72f8733bb1a4a2923182796e54da35de',
    S_SOURCE: '9eacb1de680cefb84d90d603cc1aca4dd6ef4892372a91cffa3480651c5874a8',
    f'{DESIGN}/README.md': '4f9ee03c3723f74b5cddecbb237424c75daed6962b331a77862866d24d3a3793',
    f'{DESIGN}/receipt.json': 'db730a0406337d2b779a6ca47b0668917bb8ecc01b177dd4272b0f154b6f6421',
}
ACTUAL_PINS = {
    NATIVE_REVIEW: 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    MANIFEST: '308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0',
    sweep.PRIMARY_RESULT: '2ddc2fb3f0da426551bf74b6e2aac7bb489dfe848f8d9e73cfc42fddecbe99d6',
    sweep.PRIMARY_TERMINAL: '84aafea79a60878e139ed034c199e4875e91bac21fa650354f12f944e87fd0b9',
}
LENGTH = Q(67, 2)
SIGMA = Q(1, 2)
RADIAL_FRACTION = Q(1, 2)
FACES = ('top', 'bottom')
OBS = sweep.contract.OBSERVABLES
ENDPOINTS = sweep.ENDPOINTS
SCOPE = ('Conditional dimensionless displacement-radius budgets from the exact stored-field margins. '
         'These thresholds specify sufficient bounds that future independent exact FE root certificates would need to establish. '
         'A certified root, minimum, static thickness of a future root, source uncertainty, loading branch and PDE accuracy remain separate obligations.')
require = sweep.require
sha = sweep.sha
pinned = sweep.pinned
save = sweep.save


def protocol():
    primary = pinned(ROOT / PRIMARY_PROTOCOL, PINS[PRIMARY_PROTOCOL])
    for path, expected in PINS.items():
        require(sha(ROOT / path) == expected, f'fixed design input hash changed: {path}')
    require(Path(transfer.__file__).resolve() == ROOT / T_SOURCE
            and Path(sweep.__file__).resolve() == ROOT / S_SOURCE,
            'resolved import paths differ from the frozen theorem or admission source')
    review = pinned(ROOT / THEOREM_REVIEW, PINS[THEOREM_REVIEW])
    require(review['status'] == 'PASS_CONDITIONAL_EXACT_NODAL_TO_OBSERVER_TRANSFER'
            and sweep.all_true(review['checks']) and review['reviewed_hashes'][T_SOURCE] == PINS[T_SOURCE],
            'conditional nodal observer theorem review failed')
    prereqs = {**sweep.PINS, **review['reviewed_hashes'], **primary['source_hashes']}
    for path, expected in prereqs.items():
        require(sha(ROOT / path) == expected, f'bound model, theorem or admission source changed: {path}')
    require(primary['observable_variant'] == sweep.VARIANTS[0]
            and Q(primary['large_observable']['length_unit_um_exact']) == LENGTH,
            'primary observable convention or length changed')
    plan = sweep.plan_checkpoints(primary)
    selected = [r for r in plan if r['role'] in ('small_record', 'large_final')]
    require(len(selected) == 36 and sum(r['role'] == 'small_record' for r in selected) == 34,
            'expected 34 small fields and two final large fields')
    return {
        'status': 'FROZEN_EXACT_OBSERVER_BUDGET_DESIGN', 'source_sha256': sha(SOURCE),
        'fixed_design_hashes': PINS, 'dependencies': prereqs, 'terminal_authority_hashes': ACTUAL_PINS,
        'scope': SCOPE, 'frozen_44_checkpoint_plan': plan, 'selected_36_fields': selected,
        'length_unit_um_exact': str(LENGTH),
        'observer': primary['observable_variant'],
        'frozen_conditional_allowance_um_exact': primary['selection']['conditional_source_allowance_um_exact'],
        'sigma_exact': str(SIGMA), 'radial_margin_spending_fraction_exact': str(RADIAL_FRACTION),
        'center': 'Unchanged exact rational X + injection(u) from the stored binary model coordinates and displacements.',
        'small_budget_formula': 'min over four observables of sigma*(delta - worst stored residual)/A',
        'large_budget_formula': 'sigma*distance of stored signed gap from zero / (2*(L0+C_apical))',
        'strict_radius_formula': 'e_goal = min(observation budget, min Rprime/(16 L0)) - worst binary addition gap',
        'pilot': 'Smallest positive e_goal among all 17 in each endpoint, tie by record ID, plus both final large fields. '
                 'A missing positive target blocks a complete four-pilot recommendation.',
        'launch': 'External SHA of independent own prelaunch, native independent actual primary receipt and whole 44-field manifest required.',
        'scientific_threshold': 'The prespecified conditional discrepancy allowance is an analysis parameter.',
        'root_error_bounds_proved': 0, 'new_force_evaluations': 0, 'new_equilibrium_solves': 0,
    }


def admit(prepared, review_sha, manifest_sha):
    require(review_sha == ACTUAL_PINS[NATIVE_REVIEW] and manifest_sha == ACTUAL_PINS[MANIFEST],
            'supplied actual admission digests differ from the frozen reviewed pair')
    review = pinned(ROOT / NATIVE_REVIEW, review_sha)
    manifest = pinned(ROOT / MANIFEST, manifest_sha)
    terminal = pinned(ROOT / sweep.PRIMARY_TERMINAL, ACTUAL_PINS[sweep.PRIMARY_TERMINAL])
    require(terminal['status'] == 'FINISHED', 'primary nonlinear terminal is incomplete')
    results = pinned(ROOT / sweep.PRIMARY_RESULT, ACTUAL_PINS[sweep.PRIMARY_RESULT])
    rows = sweep.planned_manifest(results, prepared['frozen_44_checkpoint_plan'])
    sweep.require_primary_review(review, manifest, ACTUAL_PINS[sweep.PRIMARY_RESULT],
                                 ACTUAL_PINS[sweep.PRIMARY_TERMINAL], rows)
    require(manifest['independent_primary_review_sha256'] == review_sha
            and type(manifest['independent_replay_observed_os_exit_code']) is int
            and manifest['independent_replay_observed_os_exit_code'] == 0,
            'whole manifest lacks a successful matching independent replay process')
    replay_terminal = sweep.read((ROOT / NATIVE_REVIEW).with_name('terminal.json'))
    require(replay_terminal['status'] == 'FINISHED' and replay_terminal['receipt_sha256'] == review_sha
            and type(replay_terminal['exit_code']) is int and replay_terminal['exit_code'] == 0,
            'independent replay terminal is incomplete or mixed')
    sweep.verify_primary_supervisor(terminal, results, ACTUAL_PINS[sweep.PRIMARY_RESULT],
                                    review['all_hashes'])
    # All 44 metadata entries and 88 independently reviewed digests are now admitted.
    # Only 36 selected final fields and their verification payloads are opened.
    selected = [r for r in rows if r['role'] in ('small_record', 'large_final')]
    require(len(selected) == 36, 'result collection does not match frozen final fields')
    return results, review, selected


def binary(value):
    require(type(value) in (int, float), 'numeric stored binary coordinate required')
    return Q(value)


def addition_gaps(reference, current, free, vector):
    require(len(reference) == len(current) and len(free) == len(vector)
            and len(set(free)) == len(free)
            and all(type(i) is int and 0 <= i < len(reference) for i in free),
            'invalid nodal dimensions or injection map')
    exact = list(map(binary, reference))
    for i, value in zip(free, vector, strict=True):
        exact[i] += binary(value)
    stored = list(map(binary, current))
    require(all(binary(float(x)) == y for x, y in zip(exact, stored, strict=True)),
            'stored coordinate does not equal binary64 rounding of X + injection(u)')
    gaps = [abs(x - y) for x, y in zip(exact, stored)]
    return exact, gaps


def model_receipt(model, state, step, materials, numerical_policy):
    require(state['reached_target'] is True and state['converged'] is True
            and state['parameters'] == materials['parameters']
            and state['integration_policy'] == numerical_policy['solve_quadrature']
            and state['reference_nodes'] == model.reference.tolist()
            and state['triangles'] == model.triangles.tolist()
            and state['activation_elements'] == model.activation_elements.tolist()
            and sum(state['activation_elements']) == 24,
            'field does not reproduce the frozen converged model and activation mask')
    ref, current = state['reference_nodes'], state['current_nodes']
    require(len(ref) == len(current) == 325 and all(len(x) == 2 for x in ref + current)
            and len(state['triangles']) == 144 and all(len(t) == 6 for t in state['triangles']),
            'unexpected P2 target geometry or 325/144 dimensions')
    free = [int(i) for i in model.free]
    rim = binary(state['parameters']['outer_radius'])
    prescribed = {2 * n + k for n, (r, _) in enumerate(ref) for k in (0, 1)
                  if binary(r) == rim or (binary(r) == 0 and k == 0)}
    require(len(free) == 611 and free == [i for i in range(650) if i not in prescribed],
            'reviewed sorted 611 free map disagrees with exact axis/rim rule')
    exact, gaps = addition_gaps([v for node in ref for v in node], [v for node in current for v in node],
                                free, state['free_displacement'])
    require(all(gaps[i] == 0 for i in prescribed), 'fixed node has a rounding or injection discrepancy')
    dose = Q(step['dose_exact'])
    gain = Q(materials['gain_exact'])
    requested = gain * dose
    actual = binary(state['activation_strength'])
    require(requested > 0 and requested == Q(step['requested_strength_exact'])
            and actual == binary(float(requested)) == Q(step['solver_strength_binary64_exact'])
            and actual - requested == Q(step['binary64_rounding_difference_exact']),
            'physical target gain/dose and stored rounded solver strength differ from frozen plan')
    controls, height = sweep.observer._state_controls(state, LENGTH)
    require(height == binary(state['parameters']['height']) * LENGTH, 'mismatched stored reference height')
    fixed = {}
    for name, ids in (('top', model.top_edges.tolist()), ('bottom', model.bottom_edges.tolist())):
        require(len(ids) == len(controls[name]) and all(len(edge) == 3 for edge in ids)
                and all(ids[i][-1] == ids[i + 1][0] for i in range(len(ids) - 1)),
                'observer model trace fails ordered P2 connectivity')
        observed = [[tuple(binary(x) * LENGTH for x in current[i]) for i in edge] for edge in ids]
        require(observed == controls[name], 'observer trace uses a different ordered node map')
        radial_dofs = [2 * ids[0][0], 2 * ids[-1][-1]]
        require(all(i in prescribed and gaps[i] == 0 for i in radial_dofs)
                and [exact[i] for i in radial_dofs] == [Q(0), rim],
                'actual axis and rim radii are not both exactly fixed')
        fixed[name] = {'trace_nodes': ids, 'radial_dofs': radial_dofs,
                       'proved_fixed_radii_um_exact': ['0', str(rim * LENGTH)]}
    return controls, {
        'free_map': free, 'fixed_trace_endpoints': fixed, 'node_addition_gap_dimensionless_exact': str(max(gaps)),
        'largest_gap_flat_indices': [i for i, gap in enumerate(gaps) if gap == max(gaps)],
        'reference_height_um_exact': str(height), 'dose_exact': str(dose), 'gain_exact': str(gain),
        'exact_target_strength': str(requested), 'solver_binary64_strength_exact': str(actual),
        'binary64_strength_rounding_error_exact': str(actual - requested),
        'material_parameters': materials['parameters'], 'material_parameter_storage_exact': materials['parameter_storage_exact'],
        'active_elements': [i for i, selected in enumerate(state['activation_elements']) if selected]}


def edge_facts(controls):
    require(bool(controls) and all(len(edge) == 3 for edge in controls)
            and all(left[-1] == right[0] for left, right in zip(controls, controls[1:])),
            'connected conforming P2 trace required')
    entries = []
    for edge, (p0, pm, p1) in enumerate(controls):
        radial, vertical = zip(p0, pm, p1)
        def endpoints(a):
            x, m, y = map(Q, a)
            return -3*x + 4*m - y, x - 4*m + 3*y
        d_r, d_z = endpoints(radial), endpoints(vertical)
        radial_min = min(d_r)
        require(radial_min > 0, 'nonpositive stored radial derivative')
        slope = max(map(abs, d_z)) / radial_min
        entries.append({'edge': edge, 'radial_start_end_um_exact': [str(radial[0]), str(radial[-1])],
                        'Rprime_local_endpoints_um_exact': list(map(str, d_r)),
                        'Zprime_local_endpoints_um_exact': list(map(str, d_z)),
                        'radial_derivative_min_um_exact': str(radial_min), 'stored_slope_cap_exact': str(slope)})
    return entries


def geometry_targets(controls, roi):
    require(0 <= roi[0] < roi[1], 'nonempty covered observer ROI required')
    faces = {}
    for name in FACES:
        edges = edge_facts(controls[name])
        radial = min(Q(e['radial_derivative_min_um_exact']) for e in edges)
        slope = max(Q(e['stored_slope_cap_exact']) for e in edges)
        proof = transfer.curve_bound(controls[name], Q(0), fixed_radial_endpoints=True, roi_um=roi)
        require(proof['status'] == 'CONDITIONAL_EXACT_ERROR_TRANSFER'
                and Q(proof['radial_derivative_margin_exact']) == radial
                and Q(proof['base_graph_Lipschitz_upper_exact']) == slope,
                'edgewise derivative and slope proof does not match the reviewed implementation')
        c = Q(5, 4) * LENGTH * (1 + slope)
        faces[name] = {'edges': edges, 'radial_min_um_exact': str(radial), 'global_stored_slope_cap_exact': str(slope),
                       'C_um_per_dimensionless_coordinate_exact': str(c), 'zero_radius_transfer_check': proof}
    t_radial = RADIAL_FRACTION * min(Q(f['radial_min_um_exact']) for f in faces.values()) / (8 * LENGTH)
    return faces, t_radial


def stored_small_residuals(record, rows, intervals, delta):
    require(set(intervals) == set(OBS) == {r['observable'] for r in rows} and len(rows) == 4,
            'need all four independent stored observer intervals for this record')
    found = {r['observable']: r for r in rows}
    values = {}
    for name in OBS:
        r = found[name]
        observed = Q(record['observed_exact'][name])
        pred = list(map(Q, intervals[name]))
        lo, hi = sweep.contract.interval_residual(observed, pred)
        require(r['record_id'] == record['record_id'] and observed == Q(r['observed_exact'])
                and pred == list(map(Q, r['prediction_interval_exact']))
                and (lo, hi) == (Q(r['minimum_interval_residual_exact']), Q(r['maximum_interval_residual_exact']))
                and delta == Q(r['discrepancy_allowance_exact']),
                'frozen row, reviewed stored interval or exact discrepancy differs')
        values[name] = {'observed_um_exact': str(observed), 'stored_interval_um_exact': list(map(str, pred)),
                        'stored_worst_residual_um_exact': str(hi), 'margin_um_exact': str(delta - hi)}
    return values


def strict_goal(allocations, coefficients, xi, radial_t):
    require(0 <= xi and radial_t > 0 and set(allocations) == set(coefficients)
            and bool(allocations) and all(v > 0 for v in coefficients.values()), 'invalid design coefficients')
    if min(allocations.values()) <= 0:
        return {'status': 'UNRESOLVED_NONPOSITIVE_OBSERVER_MARGIN', 'e_goal_dimensionless_exact': None}
    per_output = {k: allocations[k] / a for k, a in coefficients.items()}
    limiting = min(per_output, key=per_output.get)
    combined = min(per_output[limiting], radial_t)
    goal = combined - xi
    basis = {'observation_total_dimensionless_cap_by_observable_exact': {k: str(v) for k, v in per_output.items()},
             'limiting_observable': limiting, 'monotonicity_total_dimensionless_cap_exact': str(radial_t),
             'limiting_guard': 'radial_monotonicity' if radial_t < per_output[limiting] else 'observer_error',
             'total_dimensionless_coordinate_cap_exact': str(combined)}
    if goal <= 0:
        return {**basis, 'status': 'UNRESOLVED_BINARY_ADDITION_GAP_EXHAUSTS_BUDGET',
                'e_goal_dimensionless_exact': None}
    return {**basis, 'status': 'POSITIVE_CONDITIONAL_ERROR_TARGET', 'e_goal_dimensionless_exact': str(goal),
            'e_goal_dimensionless_display': float(goal), 'physical_nodal_error_goal_um_exact': str(LENGTH * combined),
            'conditional_observer_error_caps_um_exact': {k: str(a * combined) for k, a in coefficients.items()}}


def small_budget(rows, faces, xi, radial_t, delta):
    apical, basal = (Q(faces[name]['C_um_per_dimensionless_coordinate_exact']) for name in FACES)
    coefficients = dict(zip(OBS, (apical + basal, apical + basal, apical, basal)))
    allocations = {k: SIGMA * Q(rows[k]['margin_um_exact']) for k in OBS}
    goal = strict_goal(allocations, coefficients, xi, radial_t)
    content = {'stored_residuals': rows, 'small_margin_fraction_exact': str(SIGMA),
               'per_observer_allocations_um_exact': {k: str(v) for k, v in allocations.items()},
               'conditional_observer_coefficient_um_exact': {k: str(v) for k, v in coefficients.items()}}
    if goal['e_goal_dimensionless_exact'] is not None:
        retained = {}
        for k, cap in goal['conditional_observer_error_caps_um_exact'].items():
            worst = Q(rows[k]['stored_worst_residual_um_exact'])
            margin = Q(rows[k]['margin_um_exact'])
            remaining = delta - worst - Q(cap)
            require(remaining >= (1 - SIGMA) * margin > 0, 'exact small budget does not preserve strict allowance')
            retained[k] = str(remaining)
        content['conditional_remaining_allowance_um_exact'] = retained
    return {**content, **goal}


def separation_bounds(result, independent):
    def ordered(pair):
        a, b = map(Q, pair)
        require(a <= b, 'inverted stored large observable interval')
        return a, b
    intervals = [ordered(r['large_Q_um_exact']) for r in result['endpoints']]
    audited = [ordered(r['large_Q_um_exact']) for r in independent['endpoints']]
    pair = ordered(result['large_Q_positive_direction_minus_negative_direction_um_exact'])
    require(intervals == audited and pair == ordered(independent['stored_large_Q_signed_difference_um_exact'])
            and pair == (intervals[0][0] - intervals[1][1], intervals[0][1] - intervals[1][0]),
            'stored large signed gap differs from its independent exact interval')
    distance = max(Q(0), pair[0], -pair[1])
    return pair, distance


def large_budget(stored, faces, xi, radial_t, distance):
    apical = Q(faces['top']['C_um_per_dimensionless_coordinate_exact'])
    allocated = SIGMA * distance / 2
    goal = strict_goal({'Q': allocated}, {'Q': LENGTH + apical}, xi, radial_t)
    return {'stored_Q_um_exact': list(map(str, stored)), 'stored_gap_distance_from_zero_lower_um_exact': str(distance),
            'balanced_error_allocation_um_exact': str(allocated),
            'conditional_Q_error_coefficient_um_exact': str(LENGTH + apical), **goal}


def validate_design_transfer(controls, roi, faces, goal):
    if goal['e_goal_dimensionless_exact'] is None:
        return {}
    epsilon = Q(goal['physical_nodal_error_goal_um_exact'])
    total = Q(goal['total_dimensionless_coordinate_cap_exact'])
    proof = transfer.endpoint_bounds(controls['top'], controls['bottom'], epsilon,
                                     fixed_radial_endpoints=True, roi_um=roi)
    require(proof['status'] == 'CONDITIONAL_EXACT_ERROR_TRANSFER', 'design target does not satisfy reviewed transfer')
    for name, side in (('top', 'apical'), ('bottom', 'basal')):
        radial = Q(faces[name]['radial_min_um_exact'])
        require(Q(proof[side]['radial_derivative_margin_exact']) == radial - 8 * epsilon >= radial/2 > 0
                and Q(proof[side]['Eulerian_vertical_error_um_exact'])
                == total * Q(faces[name]['C_um_per_dimensionless_coordinate_exact']),
                'conditional curve error or strict derivative guarantee does not agree')
    if 'Q' in goal['conditional_observer_error_caps_um_exact']:
        require(Q(goal['conditional_observer_error_caps_um_exact']['Q'])
                == epsilon + Q(proof['apical']['Eulerian_vertical_error_um_exact']),
                'large Q budget must include fixed-axis value and apical same-radius height')
    else:
        require(proof['endpoint_error_um_exact'] == goal['conditional_observer_error_caps_um_exact'],
                'small observable error coefficients disagree with reviewed theorem')
    return proof


def choose_pilots(rows):
    selections = []
    for eid in ENDPOINTS:
        endpoint = [r for r in rows if r['endpoint_id'] == eid]
        small = [r for r in endpoint if r['role'] == 'small_record']
        large = [r for r in endpoint if r['role'] == 'large_final']
        require(len(small) == 17 and len(large) == 1, 'incomplete pilot selection set')
        if any(r['e_goal_dimensionless_exact'] is None for r in endpoint):
            return []
        limiting = min(small, key=lambda r: (Q(r['e_goal_dimensionless_exact']), r['record_id']))
        raw = min(small, key=lambda r: (min(Q(v['margin_um_exact']) for v in r['stored_residuals'].values()),
                                       r['record_id']))
        for row in (limiting, large[0]):
            out = {key: row[key] for key in ('endpoint_id', 'checkpoint', 'role', 'record_id',
                                            'field', 'field_sha256', 'e_goal_dimensionless_exact')}
            if row['role'] == 'small_record':
                out['record_with_smallest_raw_stored_observer_margin'] = raw['record_id']
            selections.append(out)
    return selections


def calculate(prepared, result, independent, selected):
    from transfer_pair_pilot import AxisymmetricModel, MAIN_RULE, SOLVER, SPARSE, Parameters
    primary = pinned(ROOT / PRIMARY_PROTOCOL, PINS[PRIMARY_PROTOCOL])
    params = {r['endpoint_id']: r for r in primary['endpoints']}
    output = {r['endpoint_id']: r for r in result['endpoints']}
    delta = Q(prepared['frozen_conditional_allowance_um_exact'])
    plan = {(r['endpoint_id'], r['checkpoint']): r['step'] for r in prepared['selected_36_fields']}
    gap, distance = separation_bounds(result, independent)
    models = {}
    rows = []
    for row in selected:
        eid, label = row['endpoint_id'], row['checkpoint']
        small = row['role'] == 'small_record'
        level = 'small' if small else 'large'
        material = params[eid]['materials'][level]
        key = (eid, level)
        if key not in models:
            models[key] = AxisymmetricModel(Parameters(**material['parameters']), SOLVER, MAIN_RULE, SPARSE)
        step = plan[eid, label]
        state = pinned(ROOT / row['field'], row['field_sha256'])
        verified = pinned(ROOT / row['verification'], row['verification_sha256'])
        controls, provenance = model_receipt(models[key], state, step, material, primary['numerical_policy'])
        xi = Q(provenance['node_addition_gap_dimensionless_exact'])
        right = binary(material['parameters']['patch_radius']) * LENGTH if small else (
            Q(primary['large_observable']['shoulder_radius_um_exact']))
        roi = (Q(0), right)
        faces, radial_t = geometry_targets(controls, roi)
        if small:
            rid = step['record_id']
            require(verified['variant'] == prepared['observer'], 'primary observer variant drift')
            residuals = [r for r in output[eid]['small_residual']['all68_residuals']['constraint_rows']
                         if r['record_id'] == rid]
            small_rows = stored_small_residuals(step, residuals, verified['observer_intervals_um_exact'], delta)
            goal = small_budget(small_rows, faces, xi, radial_t, delta)
        else:
            rid = primary['matched_input']['large_record_id']
            stored = list(map(Q, verified['surface']['Q_um_exact']))
            require(stored == list(map(Q, output[eid]['large_Q_um_exact'])),
                    'selected large field verification differs from result interval')
            goal = large_budget(stored, faces, xi, radial_t, distance)
        proof = validate_design_transfer(controls, roi, faces, goal)
        rows.append({**row, 'record_id': rid, 'center_and_model': provenance,
                     'fixed_roi_um_exact': list(map(str, roi)), 'stored_face_derivatives': faces,
                     **goal, 'conditional_transfer_at_design_cap': proof})
    require(len(rows) == 36 and sum(len(r.get('stored_residuals', {})) for r in rows) == 136,
            'complete 36-field and 136-observable budget coverage required')
    pilots = choose_pilots(rows)
    all_small = [value for row in rows for value in row.get('stored_residuals', {}).values()]
    minimum_stored_margin = min(Q(v['margin_um_exact']) for v in all_small)
    require(delta - minimum_stored_margin == Q(result['both_endpoints_required_allowance_um_exact'][1]),
            'maximum stored residual across 136 rows differs from the signed-off primary bound')
    large = [r for r in rows if r['role'] == 'large_final']
    gap_inflation = None
    if all(r['e_goal_dimensionless_exact'] is not None for r in large):
        extra = sum(Q(r['conditional_observer_error_caps_um_exact']['Q']) for r in large)
        post = (gap[0] - extra, gap[1] + extra)
        retained = max(Q(0), post[0], -post[1])
        require(retained >= (1-SIGMA) * distance > 0, 'conditional large Q gap reserve failed')
        gap_inflation = {'gap_if_each_future_error_meets_its_cap_um_exact': list(map(str, post)),
                         'remaining_distance_from_zero_lower_um_exact': str(retained)}
    return {'status': 'COMPLETE_POSITIVE_CONDITIONAL_36_FIELD_BUDGETS' if pilots else 'COMPLETE_WITH_UNRESOLVED_BUDGETS',
            'scope': SCOPE, 'field_budgets': rows, 'four_pilot_selection': pilots,
            'stored_primary_constraint_count': 136, 'minimum_stored_margin_um_exact': str(minimum_stored_margin),
            'conditional_source_allowance_um_exact': str(delta), 'stored_large_signed_gap_um_exact': list(map(str, gap)),
            'conditional_large_gap_reserve': gap_inflation,
            'exact_root_error_bounds_proved': 0, 'new_force_evaluations': 0, 'new_equilibrium_solves': 0}


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def fixtures(prepared):
    exact, gaps = addition_gaps([1., 0.], [1., 0.], [0], [2.**-54])
    graph = [[(Q(0), Q(0)), (Q(1, 2), Q(0)), (Q(1), Q(0))],
             [(Q(1), Q(0)), (Q(3, 2), Q(2)), (Q(2), Q(4))]]
    controls = {'top': graph, 'bottom': [[(r, z - 3) for r, z in edge] for edge in graph]}
    roi = (Q(0), Q(2))
    facts, t = geometry_targets(controls, roi)
    epsilon = Q(1, 64)
    q = transfer.curve_bound(graph, epsilon, fixed_radial_endpoints=True, roi_um=roi)
    # A right-shift of the shared knot transfers r=1+eps/2 from the baseline
    # second edge to the perturbed first edge, where the height stays zero.
    crossing_error = 2 * epsilon
    xi = Q(1, 2**54)
    goal = strict_goal({'Q': Q(10)}, {'Q': Q(2)}, xi, t)
    assumed = Q(goal['e_goal_dimensionless_exact'])
    delta = Q(2)
    margin = Q(1, 50)
    four_rows = {k: {'stored_worst_residual_um_exact': str(delta - margin), 'margin_um_exact': str(margin)}
                 for k in OBS}
    small = small_budget(four_rows, facts, xi, t, delta)
    large = large_budget((Q(1), Q(2)), facts, xi, t, Q(4))
    small_proof = validate_design_transfer(controls, roi, facts, small)
    large_proof = validate_design_transfer(controls, roi, facts, large)
    return {
        'planned_44_checkpoints_and_exactly_36_final_fields': (
            len(prepared['frozen_44_checkpoint_plan']) == 44 and len(prepared['selected_36_fields']) == 36),
        'addition_gap_uses_binary_before_rational_sum': exact[0] == 1 + xi and gaps == [xi, Q(0)],
        'reject_mismatched_stored_binary_rounding': rejects(lambda: addition_gaps([1.], [2.], [0], [2.**-54])),
        'negative_and_zero_margin_block_budget': all(strict_goal({'Q': a}, {'Q': Q(1)}, Q(0), t)
            ['e_goal_dimensionless_exact'] is None for a in (Q(-1), Q(0))),
        'zero_budget_after_rounding_gap_is_unresolved': strict_goal({'Q': t}, {'Q': Q(1)}, t, t)
            ['status'] == 'UNRESOLVED_BINARY_ADDITION_GAP_EXHAUSTS_BUDGET',
        'moving_knot_uses_global_edgewise_slope': (
            facts['top']['global_stored_slope_cap_exact'] == '4'
            and Q(5, 4)*epsilon < crossing_error <= Q(q['Eulerian_vertical_error_um_exact'])),
        'strict_eight_epsilon_boundary_rejected': transfer.curve_bound(
            graph, Q(1, 8), fixed_radial_endpoints=True, roi_um=roi)['status'] == 'UNRESOLVED_RADIAL_MONOTONICITY',
        'actual_goal_includes_gap_and_half_radial_reserve': (
            assumed + xi == t and 8 * LENGTH * (assumed + xi) == Q(1, 2) and assumed > 0),
        'small_retains_half_limiting_margin_and_correct_surface_composition': (
            min(map(Q, small['conditional_remaining_allowance_um_exact'].values())) == SIGMA * margin
            and small_proof['endpoint_error_um_exact'] == small['conditional_observer_error_caps_um_exact']),
        'large_equal_split_before_radial_guard_covers_axis_plus_same_radius': (
            Q(large['balanced_error_allocation_um_exact']) == 1
            and Q(large['conditional_observer_error_caps_um_exact']['Q'])
            == Q(large_proof['apical']['Eulerian_vertical_error_um_exact'])
            + Q(large['physical_nodal_error_goal_um_exact'])),
        'uncovered_roi_rejected': rejects(lambda: geometry_targets(controls, (Q(0), Q(3)))),
    }


def prepare():
    PREP.mkdir(parents=True, exist_ok=False)
    with (PREP / 'transfer_ambiguity_error_budgets_frozen.py').open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        original = Path.open
        def guard(path, *args, **kwargs):
            resolved = path.resolve()
            require(not resolved.is_relative_to(ROOT / sweep.PRIMARY_ACTUAL)
                    and not resolved.is_relative_to((ROOT / NATIVE_REVIEW).parent)
                    and resolved != ROOT / MANIFEST, 'preparation opened actual result or field bytes')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', guard):
            prepared = protocol()
            h = save(PREP / 'protocol.json', prepared)
            checks = fixtures(prepared)
        checks['no_actual_budget_payloads'] = not ACTUAL.exists()
        status = 'PASS' if sweep.all_true(checks) else 'FAIL'
        pf = save(PREP / 'preflight.json', {'status': status, 'checks': checks, 'source_sha256': sha(SOURCE),
                  'protocol_sha256': h, 'primary_protocol_sha256': PINS[PRIMARY_PROTOCOL],
                  'primary_actual_fields_opened': 0, 'actual_budgets_computed': 0,
                  'new_force_evaluations': 0, 'new_equilibrium_solves': 0})
        require(status == 'PASS', 'conditional budget preflight failed')
        print({'status': status, 'source_sha256': sha(SOURCE), 'protocol_sha256': h, 'preflight_sha256': pf}, flush=True)
    except BaseException as error:
        save(PREP / 'failure.json', {'type': type(error).__name__, 'error': str(error), 'traceback': traceback.format_exc()})
        raise


def launch(args):
    binding = {'source_sha256': sha(SOURCE), 'protocol_sha256': sha(PREP / 'protocol.json'),
               'preflight_sha256': sha(PREP / 'preflight.json'), 'primary_protocol_sha256': PINS[PRIMARY_PROTOCOL],
               'primary_actual_review_sha256': ACTUAL_PINS[NATIVE_REVIEW],
               'whole_44_field_manifest_sha256': ACTUAL_PINS[MANIFEST]}
    reviewer = pinned(args.prelaunch_review, args.prelaunch_review_sha256)
    require(sweep.good_review(reviewer) and all(reviewer.get(k) == v for k, v in binding.items()),
            'independent budget prelaunch review not admitted')
    preflight = pinned(PREP / 'preflight.json', binding['preflight_sha256'])
    require(preflight['status'] == 'PASS' and sweep.all_true(preflight['checks'])
            and preflight['protocol_sha256'] == binding['protocol_sha256']
            and preflight['source_sha256'] == binding['source_sha256']
            and preflight['primary_actual_fields_opened'] == 0
            and sha(PREP / 'transfer_ambiguity_error_budgets_frozen.py') == binding['source_sha256'],
            'preflight or immutable source snapshot mismatch')
    prepared = pinned(PREP / 'protocol.json', binding['protocol_sha256'])
    require(prepared == protocol(), 'current budget protocol differs from frozen preparation')
    result, independent, selected = admit(prepared, args.actual_review_sha256, args.field_manifest_sha256)
    binding['prelaunch_review_sha256'] = args.prelaunch_review_sha256
    ACTUAL.mkdir(exist_ok=False)
    terminal = {'status': 'FAILED', 'pid': os.getpid(), 'bindings': binding}
    save(ACTUAL / 'started.json', {**terminal, 'status': 'RUNNING', 'utc': sweep.stamp()})
    try:
        output = calculate(prepared, result, independent, selected)
        terminal['results_sha256'] = save(ACTUAL / 'results.json', {**output, 'bindings': binding})
        terminal['status'] = 'FINISHED'
    except BaseException as error:
        terminal['error'] = {'type': type(error).__name__, 'message': str(error)}
        save(ACTUAL / 'failure.json', {'error': terminal['error'], 'traceback': traceback.format_exc()})
        raise
    finally:
        terminal['utc'] = sweep.stamp()
        save(ACTUAL / 'terminal.json', terminal)
        print(terminal, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('prepare')
    stage = commands.add_parser('run')
    stage.add_argument('--prelaunch-review', required=True, type=Path)
    for name in ('prelaunch-review', 'actual-review', 'field-manifest'):
        stage.add_argument('--' + name + '-sha256', required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare()
    else:
        launch(args)


if __name__ == '__main__':
    main()

