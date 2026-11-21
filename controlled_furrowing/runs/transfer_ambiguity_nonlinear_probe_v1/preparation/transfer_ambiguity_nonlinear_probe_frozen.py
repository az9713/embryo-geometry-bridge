"""Prepare and run one frozen nonlinear diagnostic of local transfer ambiguity."""
import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from dataclasses import asdict
from unittest.mock import patch

import numpy as np
import torch

import material_cohort_direct as direct
import transfer_pair_pilot as pilot


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
OUT = ROOT / 'runs/transfer_ambiguity_nonlinear_probe_v1'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
LINEAR = 'runs/transfer_ambiguity_linear_v1'
LINEAR_RESULT = f'{LINEAR}/actual/results.json'
LINEAR_REVIEW = f'{LINEAR}/actual_independent/receipt.json'
LINEAR_TERMINAL = f'{LINEAR}/actual/terminal.json'
LINEAR_STARTED = f'{LINEAR}/actual/started.json'
LINEAR_OS_EXIT = f'{LINEAR}/actual/os_exit_observation.json'
PILOT = 'runs/transfer_pair_pilot_v1/preparation/protocol.json'
DIRECT = 'runs/material_cohort_direct_v1/protocol.json'
SMALL = 'data/derived/calibration_small_18um.csv'
LARGE_INPUTS = 'data/derived/spatial_inputs_only.json'
INPUT_HASHES = {
    LINEAR_RESULT: '68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f',
    LINEAR_REVIEW: '8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072',
    LINEAR_TERMINAL: 'd49142d9a24f4da88d610acf9979ca942639ba67a513a97245bd35d46af67a9d',
    LINEAR_STARTED: '5c71c1728570c6537213ba314b19da2caededec0ad5d1ca671bd047c98537df1',
    LINEAR_OS_EXIT: 'eadc64673469079511f034afc17406b839d7959ab818f82cebead8e0a916a826',
    PILOT: 'a2745a7e50e9a64ee6e5e597f04b30981c7910ad5e3ca2d3a16e6d4d0f3de927',
    DIRECT: '1ed458c8f0f7b9c9bacfdf0c04b907d474d06e5607ced5956ecc2d2f1adac953',
    SMALL: '1aa53b6e6ce0973faf0b038ae8d4059059b05c02bb87ed9698e24d1dea7b1518',
    LARGE_INPUTS: '799afd087541f1c6b18da34f83dce7c865ed5413ca7b34813d4987d1866b4326',
}
DIRECT_RUNTIME = ('material_cohort_direct.py', 'joint_gain_extension.py',
                  'joint_gain_calibration.py', 'joint_calibration_contract.py',
                  'piecewise_gain_minimax.py')
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
COLUMNS = ('log_apical_area', 'log_basal_area', 'apical_tension_additive',
           'log_basal_tension', 'log_common_gain')
BUDGET = Q(1, 100)
PRIMARY = 'current_disk__spatial_max'
SCOPE = ('Numerical nonlinear replay at two fixed materials on stored P2 fields, '
         'with exact stored-geometry and observer enclosures and two quadrature residual checks. '
         'Compatibility is conditional on the declared observation allowance. '
         'Source validity, certified FE roots, PDE error and biological prediction are outside this scope.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, payload):
    encoded = json.dumps(payload, allow_nan=False, indent=2) + '\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(encoded)
    return sha(path)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def binary(value):
    return Q.from_float(float(value))


def all_true(checks):
    return isinstance(checks, dict) and bool(checks) and all(v is True for v in checks.values())


def verify_hashes(bindings):
    for name, expected in bindings.items():
        require(sha(ROOT / name) == expected, f'changed input or source bytes: {name}')


def runtime_sources(pilot_protocol, direct_protocol):
    sources = dict(pilot_protocol['source_hashes'])
    require(set(sources) == {f'scripts/{name}' for name in pilot.SOURCE_NAMES},
            'unexpected pilot source inventory')
    for name in DIRECT_RUNTIME:
        sources[f'scripts/{name}'] = direct_protocol['source_hashes'][name]
    verify_hashes(sources)
    for name in sources:
        loaded = sys.modules[Path(name).stem]
        require(Path(loaded.__file__).resolve() == ROOT / name, f'import path drift: {name}')
    return sources


def admit_linear(result, review, terminal, started, observed_exit):
    require(review['status'] == 'PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC'
            and all_true(review['checks']), 'independent linear actual checks did not pass')
    for name in (LINEAR_RESULT, LINEAR_TERMINAL, LINEAR_STARTED, SMALL, PILOT):
        require(review['bindings'][name] == INPUT_HASHES[name], f'linear review binding drift: {name}')
    require(result['status'] == 'COMPLETE_LOCAL_AFFINE_DIAGNOSTIC_NOT_NONLINEAR_WITNESSES'
            and result['heldout_response_tables_read'] is False
            and result['new_equilibrium_solves'] is False, 'linear result scope or status drift')
    require(terminal['status'] == 'FINISHED' and type(terminal['exit_code']) is int
            and terminal['exit_code'] == 0 and terminal['pid'] == started['pid']
            and terminal['results_sha256'] == INPUT_HASHES[LINEAR_RESULT]
            and observed_exit['status'] == 'OBSERVED_UNIFIED_EXEC_PROCESS_EXIT_ZERO'
            and type(observed_exit['observed_exit_code']) is int
            and observed_exit['observed_exit_code'] == 0
            and observed_exit['worker_terminal_pid'] == terminal['pid']
            and observed_exit['worker_terminal_sha256'] == INPUT_HASHES[LINEAR_TERMINAL]
            and observed_exit['results_sha256'] == INPUT_HASHES[LINEAR_RESULT],
            'linear result lacks matching terminal and observed OS exit')


def sorted_small_rows(rows, archived):
    ids = {r['record_id'] for r in rows}
    require(len(rows) == len(ids) == len(archived) == 17
            and ids == {r['record_id'] for r in archived}, 'need all 17 archived small IDs once')
    matches = {r['record_id']: r for r in archived}
    prepared = []
    for row in sorted(rows, key=lambda r: (Q(r['activation_magnitude']), r['record_id'])):
        rid = row['record_id']
        expected = matches[rid]
        dose = Q(row['activation_magnitude'])
        observed = {key: str(Q(row[key])) for key in direct.OBSERVABLES}
        require(dose > 0 and dose == Q(expected['dose_exact'])
                and observed == {key: str(Q(expected['observed_exact'][key]))
                                 for key in direct.OBSERVABLES}, f'CSV/archive mismatch: {rid}')
        require(row['nominal_diameter_um'] == '18', f'unexpected small diameter label: {rid}')
        prepared.append({'record_id': rid, 'dose_exact': str(dose),
                         'observed_exact': observed, 'source_csv_dose_literal': row['activation_magnitude'],
                         'reported_diameter_um_exact': str(Q(row['reported_diameter_um']))})
    return prepared


def strength_receipt(dose, gain):
    requested = dose * gain
    rounded = binary(requested)
    require(dose >= 0 and gain > 0 and (dose == 0 or rounded > 0), 'invalid loading strength')
    return {'dose_exact': str(dose), 'requested_strength_exact': str(requested),
            'solver_strength_binary64_exact': str(rounded),
            'binary64_rounding_difference_exact': str(rounded - requested)}


def transformed(case, step):
    base = case['parameters']
    base_gain = Q(case['gain_exact'])
    params, gain = direct.proposal(base, base_gain, step)
    binary_step = list(map(binary, step))
    exp_factor = {i: binary(math.exp(float(binary_step[i]))) for i in (0, 1, 3, 4)}
    keys = pilot.MATERIAL_KEYS
    require(all(params[k] == v for k, v in base.items() if k not in keys)
            and gain == base_gain * exp_factor[4], 'transformation changed model, mesh or gain rule')
    unrounded = {keys[0]: binary(base[keys[0]]) * exp_factor[0],
                 keys[1]: binary(base[keys[1]]) * exp_factor[1],
                 keys[2]: binary(base[keys[2]]) + binary_step[2],
                 keys[3]: binary(base[keys[3]]) * exp_factor[3]}
    require(all(binary(params[name]) == binary(exact) for name, exact in unrounded.items()),
            'material multiplication/addition did not round as declared')
    return {'parameters': params, 'gain_exact': str(gain),
            'parameter_storage_exact': {key: str(binary(value)) if type(value) is float else str(value)
                                        for key, value in params.items()},
            'coordinate_binary64_exact': list(map(str, binary_step)),
            'coordinate_rounding_error_exact': [str(a - b) for a, b in zip(binary_step, step)],
            'exp_of_binary_coordinate_binary64_exact': {COLUMNS[i]: str(v)
                                                       for i, v in exp_factor.items()},
            'material_before_final_rounding_exact': {k: str(v) for k, v in unrounded.items()},
            'material_final_rounding_error_exact': {k: str(binary(params[k]) - v)
                                                    for k, v in unrounded.items()}}


def select_pairs(result):
    require(tuple(result['columns']) == COLUMNS, 'unexpected affine coordinate convention')
    rows = [row for row in result['rows'] if row['mode'] == 'source_residual_allowance'
            and row['budget_exact'] == str(BUDGET)]
    require(len(rows) == 1, 'exactly one source allowance row at the prescribed budget required')
    row = rows[0]
    delta = Q(result['source_baseline_residual_um_exact']) + BUDGET
    require(delta == Q(row['allowance_um_exact']) and delta > 0, 'affine source allowance drift')
    gradient = list(map(Q, result['rounded_affine_contrast_gradient_exact']))
    require(len(gradient) == len(COLUMNS), 'unexpected affine gradient width')
    pairs = []
    for name in ENDPOINTS:
        selected = row[name]
        require(selected['exact_primal_feasible'] is True and selected['exact_weak_duality'] is True,
                f'unverified linear endpoint: {name}')
        full_step = list(map(Q, selected['step_exact']))
        require(len(full_step) == len(COLUMNS), f'bad endpoint coordinate coverage: {name}')
        step = [value / 2 for value in full_step]
        shift = sum(g * value for g, value in zip(gradient, step))
        sign = 1 if name == ENDPOINTS[0] else -1
        require(2 * sign * shift == Q(selected['primal_lower_exact']),
                f'endpoint does not reproduce its affine objective: {name}')
        pairs.append((name, full_step, step, shift))
    width = pairs[0][3] - pairs[1][3]
    affine_lower, affine_upper = map(Q, row['contrast_range_width_um_exact'])
    require(width > 0 and width * 2 == affine_lower <= affine_upper,
            'selected endpoint contrast orientation or affine width drift')
    return pairs, delta, width


def prepare_protocol():
    verify_hashes(INPUT_HASHES)
    result, review, terminal, started, observed_exit = (
        read(ROOT / path) for path in (LINEAR_RESULT, LINEAR_REVIEW, LINEAR_TERMINAL,
                                      LINEAR_STARTED, LINEAR_OS_EXIT))
    admit_linear(result, review, terminal, started, observed_exit)
    pilot_protocol, direct_protocol = read(ROOT / PILOT), read(ROOT / DIRECT)
    sources = runtime_sources(pilot_protocol, direct_protocol)
    with (ROOT / SMALL).open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    small = sorted_small_rows(rows, direct_protocol['records'])
    with (ROOT / LARGE_INPUTS).open(encoding='utf-8') as handle:
        large = json.load(handle, parse_float=str)
    selected = pilot.select_input(rows, large)
    require(all(pilot_protocol['input_selection'][k] == v for k, v in selected.items()),
            'matched outcome-free candidate input drift')
    require({row['reported_diameter_um_exact'] for row in small}
            == {selected['matched_small_reported_diameter_um_exact']}, 'small intervention geometry drift')
    shoulder = Q(pilot_protocol['observable']['shoulder_radius_um_exact'])
    dose = Q(selected['matched_dose_exact'])
    require(2 * shoulder == Q(selected['large_reported_diameter_um_exact'])
            and shoulder == Q('24.94698765'), 'fixed current radius drift')
    baselines = {case['diameter_level']: case for case in pilot_protocol['cases']
                 if case['candidate'] == 'baseline'}
    require(set(baselines) == {'small', 'large'}, 'expected exactly two baseline geometries')
    for name, policy in (('solve', pilot.SOLVER), ('sparse', pilot.SPARSE),
                         ('solve_quadrature', pilot.MAIN_RULE), ('cross_force_quadrature', pilot.CROSS_RULE)):
        require(asdict(policy) == pilot_protocol['numerical_policy'][name], f'numerical policy drift: {name}')
    baseline = baselines['small']
    require(baseline['parameters'] == direct_protocol['base_parameters']
            and Q(baseline['gain_exact']) == Q(direct_protocol['base_gain_exact'])
            and baseline['gain_exact'] == baselines['large']['gain_exact']
            and {k for k, v in baseline['parameters'].items()
                 if baselines['large']['parameters'][k] != v} == {'patch_radius'},
            'pilot and small-cohort baseline differ outside intervention radius')
    for name, case in baselines.items():
        require(binary(Q(case['diameter_um_exact']) / (2 * pilot.LENGTH))
                == binary(case['parameters']['patch_radius']), f'bad patch geometry: {name}')
    pairs, delta, width = select_pairs(result)
    cases = []
    for name, full_step, step, shift in pairs:
        materials = {level: transformed(case, step) for level, case in baselines.items()}
        gain = Q(materials['small']['gain_exact'])
        require(gain == Q(materials['large']['gain_exact']), 'small/large gain mismatch')
        checkpoints = pilot.dose_schedule(dose, gain)
        require(checkpoints == [strength_receipt(Q(i, 3) * dose, gain) for i in range(4)],
                'large checkpoint exact strength mismatch')
        cases.append({'endpoint_id': name, 'full_linear_step_exact': list(map(str, full_step)),
                      'half_step_exact': list(map(str, step)),
                      'affine_large_Q_change_um_exact': str(shift), 'materials': materials,
                      'small_zero': strength_receipt(Q(0), gain),
                      'small_records': [{**r, **strength_receipt(Q(r['dose_exact']), gain)} for r in small],
                      'large_checkpoints': checkpoints})
    return {
        'version': 1, 'status': 'FROZEN_OUTCOME_FREE_NONLINEAR_DIAGNOSTIC',
        'scope': SCOPE, 'source_sha256': sha(SOURCE), 'source_hashes': sources,
        'input_hashes': INPUT_HASHES, 'runtime': pilot.runtime(), 'coordinates': list(COLUMNS),
        'selection': {'mode': 'source_residual_allowance', 'increment_um_exact': str(BUDGET),
                      'scale_each_frozen_step_exact': '1/2',
                      'affine_baseline_source_residual_um_exact': result['source_baseline_residual_um_exact'],
                      'conditional_source_allowance_um_exact': str(delta),
                      'half_pair_affine_large_Q_range_um_exact': str(width),
                      'half_pair_affine_large_Q_range_um_display': float(width),
                      'allowance_origin': 'Threshold for this frozen affine diagnostic, without an experimental noise model.'},
        'baseline': baselines, 'matched_input': selected, 'observable_variant': PRIMARY,
        'large_observable': pilot_protocol['observable'], 'numerical_policy': pilot_protocol['numerical_policy'],
        'endpoints': cases,
        'loading': {'serial_torch_threads': 1, 'small': 'Fresh zero equilibrium, then exact ascending small optical dose.',
                    'large': 'Separate fresh large zero, then 1/3, 2/3 and 1 times the matched dose.',
                    'zero_solve': 'Same fixed Q8 zero attempt used by inherited continuation, saved before replay checks.',
                    'retained_final_outcome_fields_per_endpoint': 18, 'retained_checkpoints_per_endpoint': 22,
                    'rounding': 'Each binary64 material, dose product and implied reached dose is recorded exactly.',
                    'stop_policy': 'A numerical or physical verification failure stops the serial queue and preserves attempted stages. '
                                   'An allowance outcome leaves both complete endpoint measurements in scope.'},
        'authority': {'independent_linear_actual_review_sha256': INPUT_HASHES[LINEAR_REVIEW],
                      'requires_independent_own_prelaunch_review': True,
                      'requires_external_own_prelaunch_review_sha256': True},
        'heldout_response_tables_read': False, 'retuning_or_fitting': False,
        'source_witness_certified': False, 'exact_integral_FE_root_certified': False,
        'continuum_accuracy_certified': False,
    }


def residual_summary(records, predictions, delta):
    observations = {row['record_id']: row['observed_exact'] for row in records}
    require(len(records) == len(observations) == len(predictions) == 17, '17-record predictions required')
    result = direct.evaluate(observations, predictions, {k: str(delta) for k in direct.OBSERVABLES})
    require(result['constraints'] == 68, 'incomplete all-observer residual coverage')
    bounds = (max(map(Q, result['discrepancy_lower_requirements_exact'].values())),
              max(map(Q, result['discrepancy_sufficient_requirements_exact'].values())))
    status = 'ALLOWANCE_SATISFIED' if bounds[1] <= delta else (
        'ALLOWANCE_EXCEEDED' if bounds[0] > delta else 'UNRESOLVED')
    require(result['status'] == {'ALLOWANCE_SATISFIED': 'IN', 'ALLOWANCE_EXCEEDED': 'OUT',
                                 'UNRESOLVED': 'UNRESOLVED'}[status], 'allowance predicates disagree')
    return {'status': status, 'required_allowance_um_exact': list(map(str, bounds)),
            'required_allowance_um_display': list(map(float, bounds)),
            'declared_allowance_um_exact': str(delta), 'all68_residuals': result}


def interval_difference(after, before):
    a, b = map(Q, after)
    c, d = map(Q, before)
    require(a <= b and c <= d, 'ordered interval bounds required')
    return [str(a - d), str(b - c)]


def checkpoint_plan(endpoint):
    yield 'small_zero', endpoint['small_zero'], 'small_zero'
    for index, row in enumerate(endpoint['small_records']):
        yield f'small_{index + 1:02d}', row, 'small_record'
    for index, step in enumerate(endpoint['large_checkpoints']):
        yield f'large_{index:02d}', step, 'large_final' if index == 3 else (
            'large_zero' if index == 0 else 'large_checkpoint')


def review_bindings(review, expected):
    return review.get('status') == 'PASS' and all_true(review.get('checks')) and all(
        review.get(key) == value for key, value in expected.items())


def load_launch(review_path, review_sha):
    protocol_path = PREP / 'protocol.json'
    preflight_path = PREP / 'preflight.json'
    require(sha(review_path) == review_sha, 'independent prelaunch receipt SHA mismatch')
    bindings = {'source_sha256': sha(SOURCE), 'protocol_sha256': sha(protocol_path),
                'preflight_sha256': sha(preflight_path),
                'linear_actual_review_sha256': INPUT_HASHES[LINEAR_REVIEW]}
    require(review_bindings(read(review_path), bindings), 'independent prelaunch review gate did not pass')
    preflight = read(preflight_path)
    require(preflight['status'] == 'PASS' and all_true(preflight['checks'])
            and preflight['source_sha256'] == bindings['source_sha256']
            and preflight['protocol_sha256'] == bindings['protocol_sha256']
            and preflight['nonlinear_jobs_launched'] == 0, 'invalid preparation preflight')
    require(sha(PREP / 'transfer_ambiguity_nonlinear_probe_frozen.py') == sha(SOURCE),
            'frozen driver snapshot drift')
    protocol = read(protocol_path)
    require(protocol == prepare_protocol(), 'runtime input or frozen protocol semantics drift')
    return protocol, {**bindings, 'prelaunch_review_sha256': review_sha}


def verify_field(model, state, step):
    expected = Q(step['solver_strength_binary64_exact'])
    require(state['reached_target'] is True and state['converged'] is True
            and binary(state['activation_strength']) == expected, 'requested strength did not verify')
    require(math.isfinite(state['scaled_free_force_linf'])
            and 0 <= state['scaled_free_force_linf'] < pilot.SOLVER.force_tolerance,
            'reported Q8 force failed')
    vector = torch.tensor(state['free_displacement'], dtype=torch.float64)
    nodes = model.coordinates(vector).detach().numpy()
    require(np.isfinite(nodes).all()
            and np.array_equal(model.reference, state['reference_nodes'])
            and np.array_equal(model.triangles, state['triangles'])
            and np.array_equal(nodes, state['current_nodes'])
            and np.array_equal(nodes[model.top_edges], state['boundary_top_quadratic_controls'])
            and np.array_equal(nodes[model.bottom_edges], state['boundary_bottom_quadratic_controls']),
            'stored field nodes, mesh or controls differ from the solver vector')


def new_zero(model):
    state = model._attempt(0., np.zeros(len(model.free)))
    state['requested_activation_strength'] = 0.
    state['reached_target'] = state['converged'] is True and binary(state['activation_strength']) == 0
    state['zero_initialization'] = 'Independent zero vector in this endpoint and geometry.'
    return state


def save_field(path, state):
    try:
        return save(path, state)
    except ValueError:
        # Preserve solver diagnostics even when nonfinite entries disqualify a field.
        raw = path.with_name('nonfinite_solver_field.txt')
        with raw.open('x', encoding='utf-8') as handle:
            handle.write(json.dumps(state, allow_nan=True, indent=2) + '\n')
        raise


def certify_small(model, state):
    verification, variants = direct.certify_state(model, state)
    cross = verification['Q12_cross_evaluation']['scaled_force_linf']
    geometry = verification['geometry']
    require(math.isfinite(cross) and 0 <= cross < pilot.SOLVER.force_tolerance
            and geometry['all_exact_J_positive'] is True and Q(geometry['J_lower_exact']) > 0
            and geometry['global_embedding_certified'] is True, 'Q12 or exact stored geometry gate failed')
    values = variants[PRIMARY]
    require(set(values) == set(direct.OBSERVABLES)
            and all(len(pair) == 2 and Q(pair[0]) <= Q(pair[1]) for pair in values.values()),
            'incomplete or invalid exact observation intervals')
    return {'verification': verification, 'variant': PRIMARY, 'observer_intervals_um_exact': values}


def solve_endpoint(protocol, endpoint, folder):
    materials = endpoint['materials']
    models = {level: pilot.build_models(pilot.Parameters(**material['parameters']))
              for level, material in materials.items()}
    meshes = {}
    for level, (model, cross) in models.items():
        pilot.verify_mesh_pair(model, cross)
        meshes[level] = pilot.mesh_receipt(model)
    save(folder / 'meshes.json', meshes)
    previous = {}
    predictions = {}
    checkpoints = []
    gain = Q(materials['small']['gain_exact'])
    shoulder = Q(protocol['large_observable']['shoulder_radius_um_exact'])
    for label, step, role in checkpoint_plan(endpoint):
        verify_hashes(protocol['input_hashes'])
        verify_hashes(protocol['source_hashes'])
        require(sha(SOURCE) == protocol['source_sha256'], 'runtime driver byte drift')
        level = 'small' if role.startswith('small') else 'large'
        model, cross = models[level]
        stage = folder / label
        stage.mkdir(exist_ok=False)
        save(stage / 'requested.json', {'level': level, 'role': role, 'step': step,
                                        'previous': previous[level][1] if level in previous else None})
        strength = float(Q(step['solver_strength_binary64_exact']))
        if level in previous:
            state = model.advance(strength, previous[level][0], output_directory=stage / 'solver')
        else:
            require(role in ('small_zero', 'large_zero') and strength == 0, 'fresh zero was skipped')
            state = new_zero(model)
        field_path = stage / 'field.json'
        field_sha = save_field(field_path, state)
        verify_field(model, state, step)
        if level == 'small':
            verified = certify_small(model, state)
        else:
            verified = pilot.checkpoint(model, cross, state, step, shoulder, role == 'large_final')
        verification_path = stage / 'verification.json'
        verification_sha = save(verification_path, verified)
        if level == 'large':
            require(verified['status'] == 'PASS' and all_true(verified['checks']),
                    f'large checkpoint failed: {label}')
            require(all(math.isfinite(verified[k]['scaled_force_linf'])
                        and 0 <= verified[k]['scaled_force_linf'] < pilot.SOLVER.force_tolerance
                        for k in ('Q8_mechanical_evaluation', 'Q12_cross_evaluation')),
                    f'nonfinite Q8/Q12 force at {label}')
        if role == 'small_record':
            predictions[step['record_id']] = verified['observer_intervals_um_exact']
        reached = binary(state['activation_strength'])
        continuation = state.get('robust_continuation', {})
        row = {'checkpoint': label, 'role': role, 'step': step,
               'actual_strength_binary64_exact': str(reached),
               'actual_dose_implied_by_binary_strength_exact': str(reached / gain),
               'field': relative(field_path), 'field_sha256': field_sha,
               'verification': relative(verification_path), 'verification_sha256': verification_sha,
               'accepted_path_implied_optical_doses': pilot.implied_dose(continuation.get('accepted_increments', []), gain),
               'attempts': len(continuation.get('attempts', [])),
               'rejected_attempts': sum(not item['accepted'] for item in continuation.get('attempts', []))}
        checkpoints.append(row)
        previous[level] = (state, {'field': row['field'], 'field_sha256': field_sha,
                                   'verification': row['verification'], 'verification_sha256': verification_sha})
        save(stage / 'passed.json', row)
        print(json.dumps({'endpoint': endpoint['endpoint_id'], 'checkpoint': label, 'status': 'VERIFIED'}), flush=True)
    delta = Q(protocol['selection']['conditional_source_allowance_um_exact'])
    residual = residual_summary(endpoint['small_records'], predictions, delta)
    require(len(checkpoints) == 22, 'incomplete stored checkpoint inventory')
    return {'status': 'COMPLETE_CONDITIONAL_STORED_FIELDS', 'endpoint_id': endpoint['endpoint_id'],
            'checkpoints': checkpoints, 'verified_final_outcome_fields': 18, 'small_residual': residual,
            'large_Q_um_exact': verified['surface']['Q_um_exact'], 'scope': SCOPE}


def worker(review_path, review_sha, endpoint_id, supervisor_pid):
    protocol, bindings = load_launch(review_path, review_sha)
    launch = read(ACTUAL / 'started.json')
    require(launch['pid'] == supervisor_pid == os.getppid() and launch['bindings'] == bindings
            and endpoint_id in launch['serial_endpoint_order'], 'worker requires matching live supervisor')
    endpoint, = [item for item in protocol['endpoints'] if item['endpoint_id'] == endpoint_id]
    folder = ACTUAL / endpoint_id
    require(folder.is_dir() and not (folder / 'worker_started.json').exists(), 'fresh supervisor stage required')
    pid = os.getpid()
    save(folder / 'worker_started.json', {'pid': pid, 'supervisor_pid': supervisor_pid,
                                         'endpoint_id': endpoint_id, 'bindings': bindings, 'utc': stamp()})
    terminal = {'status': 'FAILED', 'pid': pid, 'supervisor_pid': supervisor_pid, 'endpoint_id': endpoint_id,
                'binds': bindings}
    try:
        result = solve_endpoint(protocol, endpoint, folder)
        terminal['results_sha256'] = save(folder / 'results.json', result)
        terminal['status'] = 'FINISHED'
    except BaseException as error:
        terminal['error'] = {'type': type(error).__name__, 'message': str(error)}
        save(folder / 'failure.json', {'error': terminal['error'], 'traceback': traceback.format_exc()})
        raise
    finally:
        terminal['utc'] = stamp()
        save(folder / 'terminal.json', terminal)


def admit_worker(terminal, observed, bindings, expected_pid, endpoint_id):
    require(type(observed) is int and observed == 0 and terminal['status'] == 'FINISHED'
            and terminal['pid'] == expected_pid and terminal['endpoint_id'] == endpoint_id
            and terminal['binds'] == bindings and len(terminal['results_sha256']) == 64,
            f'worker failed or terminal binding mismatch: {endpoint_id}')


def summarize_pair(results):
    require([row['endpoint_id'] for row in results] == list(ENDPOINTS), 'exact ordered pair required')
    require(all(row['verified_final_outcome_fields'] == 18
                and row['small_residual']['all68_residuals']['constraints'] == 68 for row in results),
            'incomplete endpoint measurements')
    # Positive affine direction minus negative affine direction.
    separation = interval_difference(results[0]['large_Q_um_exact'], results[1]['large_Q_um_exact'])
    spans = [list(map(Q, row['small_residual']['required_allowance_um_exact'])) for row in results]
    required = [str(max(endpoints)) for endpoints in zip(*spans)]
    return {'large_Q_positive_direction_minus_negative_direction_um_exact': separation,
            'large_Q_signed_difference_um_display': list(map(lambda v: float(Q(v)), separation)),
            'large_Q_stored_intervals_disjoint': Q(separation[0]) > 0 or Q(separation[1]) < 0,
            'both_endpoints_required_allowance_um_exact': required,
            'both_endpoints_satisfy_declared_allowance': all(
                row['small_residual']['status'] == 'ALLOWANCE_SATISFIED' for row in results),
            'verified_outcome_fields': 36, 'small_observer_intervals': 136}


def launch(review_path, review_sha):
    protocol, bindings = load_launch(review_path, review_sha)
    ACTUAL.mkdir(exist_ok=False)
    pid = os.getpid()
    save(ACTUAL / 'started.json', {'pid': pid, 'bindings': bindings,
                                   'serial_endpoint_order': list(ENDPOINTS), 'utc': stamp()})
    terminal = {'status': 'FAILED', 'pid': pid, 'bindings': bindings, 'workers': []}
    results = []
    try:
        for endpoint in protocol['endpoints']:
            endpoint_id = endpoint['endpoint_id']
            folder = ACTUAL / endpoint_id
            folder.mkdir(exist_ok=False)
            command = [sys.executable, '-B', str(SOURCE), 'worker', '--review', str(review_path),
                       '--review-sha256', review_sha, '--endpoint', endpoint_id,
                       '--supervisor-pid', str(pid)]
            with (folder / 'stdout.log').open('x', encoding='utf-8') as log:
                child = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                save(folder / 'supervisor_child_started.json', {'pid': child.pid, 'supervisor_pid': pid,
                                                                'command': command, 'utc': stamp()})
                code = child.wait()
            child_terminal = folder / 'terminal.json'
            exit_receipt = {'pid': child.pid, 'supervisor_pid': pid, 'observed_exit_code': code,
                            'terminal_sha256': sha(child_terminal) if child_terminal.is_file() else None,
                            'endpoint_id': endpoint_id, 'utc': stamp()}
            exit_receipt_sha = save(folder / 'supervisor_os_exit.json', exit_receipt)
            terminal['workers'].append({**exit_receipt, 'os_exit_receipt_sha256': exit_receipt_sha})
            require(exit_receipt['terminal_sha256'] is not None, f'worker terminal absent: {endpoint_id}')
            receipt = read(child_terminal)
            admit_worker(receipt, code, bindings, child.pid, endpoint_id)
            child_started = read(folder / 'worker_started.json')
            require(child_started['pid'] == child.pid and child_started['endpoint_id'] == endpoint_id
                    and child_started['supervisor_pid'] == pid and child_started['bindings'] == bindings
                    and receipt['supervisor_pid'] == pid
                    and sha(folder / 'results.json') == receipt['results_sha256'],
                    'worker supervisor or result binding mismatch')
            results.append(read(folder / 'results.json'))
        summary = summarize_pair(results)
        terminal['results_sha256'] = save(ACTUAL / 'results.json', {
            'status': 'COMPLETE_CONDITIONAL_STORED_FIELD_COMPARISON', 'bindings': bindings,
            'endpoints': results, **summary, 'scope': SCOPE,
            'source_witness_certified': False, 'exact_integral_FE_root_certified': False,
            'continuum_accuracy_certified': False, 'heldout_response_tables_read': False})
        terminal['status'] = 'FINISHED'
    except BaseException as error:
        terminal['error'] = {'type': type(error).__name__, 'message': str(error)}
        raise
    finally:
        terminal['utc'] = stamp()
        save(ACTUAL / 'terminal.json', terminal)
        print(json.dumps(terminal), flush=True)


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def synthetic_checks(protocol):
    cases = protocol['endpoints']
    plans = [list(checkpoint_plan(case)) for case in cases]
    records = cases[0]['small_records']
    observations = {row['record_id']: row['observed_exact'] for row in records}
    delta = Q(1, 10)
    prediction = {rid: {key: (str(Q(value) - delta), str(Q(value) + delta))
                        for key, value in observed.items()} for rid, observed in observations.items()}
    at_limit = residual_summary(records, prediction, delta)
    rid, key = next(iter(prediction)), direct.OBSERVABLES[0]
    value = Q(observations[rid][key])
    prediction[rid][key] = (str(value + delta), str(value + 2 * delta))
    unresolved = residual_summary(records, prediction, delta)
    prediction[rid][key] = (str(value + 2 * delta), str(value + 3 * delta))
    exceeds = residual_summary(records, prediction, delta)
    compatibility_is_outcome = summarize_pair([
        {'endpoint_id': e, 'verified_final_outcome_fields': 18, 'small_residual': verdict, 'large_Q_um_exact': q}
        for e, verdict, q in zip(ENDPOINTS, (exceeds, at_limit), (['3', '5'], ['-2', '-1']))
    ])
    partial = prediction.copy()
    partial.pop(rid)
    exit_binding = {'source_sha256': 'b' * 64}
    terminal = {'status': 'FINISHED', 'pid': 10, 'endpoint_id': ENDPOINTS[0],
                'binds': exit_binding, 'results_sha256': 'c' * 64}
    admit_worker(terminal, 0, exit_binding, 10, ENDPOINTS[0])
    return {
        'exact_17_reconciled_and_ordered_doses_per_endpoint': all(
            [Q(row['dose_exact']) for row in case['small_records']] ==
            sorted(Q(row['dose_exact']) for row in case['small_records'])
            and len({row['record_id'] for row in case['small_records']}) == 17 for case in cases),
        'exact_half_steps_and_Ta_stays_zero': all(
            [Q(x) / 2 for x in case['full_linear_step_exact']] == list(map(Q, case['half_step_exact']))
            and Q(case['half_step_exact'][2]) == 0 for case in cases),
        'separate_zeros_and_complete_36_final_fields': all(
            len(plan) == 22 and plan[0][2] == 'small_zero' and plan[18][2] == 'large_zero'
            and sum(role in ('small_record', 'large_final') for _, _, role in plan) == 18 for plan in plans),
        'exact_checkpoint_strength_and_rounding': all(
            step == {**step, **strength_receipt(Q(step['dose_exact']),
                                               Q(case['materials']['small']['gain_exact']))}
            for case, plan in zip(cases, plans) for _, step, _ in plan),
        'exact_allowance_equality_overlap_and_strict_failure': (
            at_limit['status'] == 'ALLOWANCE_SATISFIED' and unresolved['status'] == 'UNRESOLVED'
            and exceeds['status'] == 'ALLOWANCE_EXCEEDED'
            and all(r['all68_residuals']['constraints'] == 68 for r in (at_limit, unresolved, exceeds))),
        'allowance_failure_retains_complete_pair': (
            compatibility_is_outcome['verified_outcome_fields'] == 36
            and compatibility_is_outcome['small_observer_intervals'] == 136
            and compatibility_is_outcome['both_endpoints_satisfy_declared_allowance'] is False),
        'reject_missing_small_observation': rejects(lambda: residual_summary(records, partial, delta)),
        'signed_interval_difference_and_bad_bounds': (
            interval_difference(['3', '5'], ['-2', '-1']) == ['4', '7']
            and rejects(lambda: interval_difference(['5', '3'], ['-2', '-1']))),
        'reject_bad_review_bytes_empty_checks_and_truthy_sentinels': (
            review_bindings({'status': 'PASS', 'checks': {'fixture': True}, **exit_binding}, exit_binding)
            and not review_bindings({'status': 'PASS', 'checks': {'fixture': True}, **exit_binding},
                                    {'source_sha256': 'd' * 64})
            and not review_bindings({'status': 'PASS', 'checks': {}, **exit_binding}, exit_binding)
            and not review_bindings({'status': 'PASS', 'checks': {'fixture': 1}, **exit_binding}, exit_binding)),
        'reject_worker_failure_cross_endpoint_and_pid_drift': (
            rejects(lambda: admit_worker(terminal, 1, exit_binding, 10, ENDPOINTS[0]))
            and rejects(lambda: admit_worker(terminal, 0, exit_binding, 10, ENDPOINTS[1]))
            and rejects(lambda: admit_worker(terminal, 0, exit_binding, 11, ENDPOINTS[0]))),
    }


def reference_checks(protocol):
    checks = {}
    receipts = {}
    shoulder = Q(protocol['large_observable']['shoulder_radius_um_exact'])
    for case in protocol['endpoints']:
        for level, material in case['materials'].items():
            label = f"{case['endpoint_id']}_{level}"
            model, cross = pilot.build_models(pilot.Parameters(**material['parameters']))
            pilot.verify_mesh_pair(model, cross)
            mesh = pilot.mesh_receipt(model)
            state = pilot.identity_state(model)
            geom = pilot.exact_geometry(state)
            surface = pilot.surface_receipt(state, shoulder)
            main_force = pilot.evaluate(model, state)['scaled_force_linf']
            cross_force = pilot.evaluate(cross, state)['scaled_force_linf']
            verification_fixture = {**state, 'reached_target': True, 'converged': True,
                                    'scaled_free_force_linf': main_force}
            zero_step = strength_receipt(Q(0), Q(material['gain_exact']))
            verify_field(model, verification_fixture, zero_step)
            changed = [row.copy() for row in state['current_nodes']]
            changed[0][1] += .125
            checks[f'{label}_reject_corrupt_serialized_node'] = rejects(
                lambda: verify_field(model, {**verification_fixture, 'current_nodes': changed}, zero_step))
            checks[label] = (geom['all_exact_J_positive'] is True and Q(geom['J_lower_exact']) > 0
                             and geom['global_embedding_certified'] is True and surface['usable_stored_graph'] is True
                             and surface['Q_um_exact'] == ['0', '0'] and mesh['patch_radius_is_vertex_aligned']
                             and mesh['axial_activation_faces_aligned'] and mesh['has_axis_and_outer_rim_fixed_components']
                             and 0 < mesh['activation_element_count'] < mesh['triangle_count']
                             and 0 <= main_force < pilot.SOLVER.force_tolerance
                             and 0 <= cross_force < pilot.SOLVER.force_tolerance)
            receipts[label] = {**mesh, 'exact_J_lower': geom['J_lower_exact'],
                               'exact_J_upper': geom['J_upper_exact'], 'Q_um_exact': surface['Q_um_exact'],
                               'unsolved_identity_Q8_force': main_force, 'unsolved_identity_Q12_force': cross_force}
    return checks, receipts


def prepare():
    require(not ACTUAL.exists(), 'actual output already exists')
    PREP.mkdir(parents=True, exist_ok=False)
    snapshot = PREP / 'transfer_ambiguity_nonlinear_probe_frozen.py'
    with snapshot.open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        protocol = prepare_protocol()
        protocol_sha = save(PREP / 'protocol.json', protocol)
        with (patch.object(pilot.AxisymmetricModel, '_attempt', side_effect=AssertionError('solve forbidden')),
              patch.object(pilot.AxisymmetricModel, 'advance', side_effect=AssertionError('advance forbidden')),
              patch.object(pilot.CrossModel, '_attempt', side_effect=AssertionError('cross solve forbidden')),
              patch.object(pilot.CrossModel, 'advance', side_effect=AssertionError('cross advance forbidden')),
              patch.object(subprocess, 'Popen', side_effect=AssertionError('subprocess forbidden'))):
            checks = synthetic_checks(protocol)
            reference, receipts = reference_checks(protocol)
        checks.update(reference)
        checks['no_actual_job_and_single_torch_thread'] = not ACTUAL.exists() and torch.get_num_threads() == 1
        status = 'PASS' if all_true(checks) else 'FAIL'
        digest = save(PREP / 'preflight.json', {
            'status': status, 'source_sha256': sha(SOURCE), 'protocol_sha256': protocol_sha,
            'linear_actual_review_sha256': INPUT_HASHES[LINEAR_REVIEW], 'checks': checks,
            'unsolved_zero_reference_receipts': receipts, 'nonlinear_jobs_launched': 0,
            'heldout_response_tables_read': False,
            'scope': 'Exact input and synthetic gate checks, plus zero identity field geometry and residual evaluation.'})
        require(status == 'PASS', 'nonlinear diagnostic preparation preflight failed')
        print(json.dumps({'status': status, 'source_sha256': sha(SOURCE),
                          'protocol_sha256': protocol_sha, 'preflight_sha256': digest}), flush=True)
    except BaseException as error:
        save(PREP / 'failure.json', {'error': type(error).__name__, 'message': str(error),
                                     'traceback': traceback.format_exc(), 'utc': stamp()})
        raise


def main():
    torch.set_num_threads(1)
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest='stage', required=True)
    stages.add_parser('prepare')
    for mode in ('launch', 'worker'):
        sub = stages.add_parser(mode)
        sub.add_argument('--review', required=True, type=Path)
        sub.add_argument('--review-sha256', required=True)
        if mode == 'worker':
            sub.add_argument('--endpoint', required=True, choices=ENDPOINTS)
            sub.add_argument('--supervisor-pid', required=True, type=int)
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    elif args.stage == 'launch':
        launch(args.review.resolve(), args.review_sha256)
    else:
        worker(args.review.resolve(), args.review_sha256, args.endpoint, args.supervisor_pid)


if __name__ == '__main__':
    main()

