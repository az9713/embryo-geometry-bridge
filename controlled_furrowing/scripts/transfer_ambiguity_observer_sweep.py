"""Secondary observer sweep of the frozen nonlinear pair, using only stored fields."""

import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import re
import traceback
from unittest.mock import patch

import annulus_observer as intervals
import joint_calibration_contract as contract
import quadratic_endpoint_observer as observer


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
PRIMARY = 'runs/transfer_ambiguity_nonlinear_probe_v1'
OUT = ROOT / PRIMARY / 'observer_robustness'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
PRIMARY_PROT = f'{PRIMARY}/preparation/protocol.json'
PRIMARY_PREFLIGHT = f'{PRIMARY}/preparation/preflight.json'
PRIMARY_PRELAUNCH = f'{PRIMARY}/preparation/prelaunch_independent_review.json'
PRIMARY_DRIVER = 'scripts/transfer_ambiguity_nonlinear_probe.py'
PRIMARY_ACTUAL = f'{PRIMARY}/actual'
PRIMARY_RESULT = f'{PRIMARY_ACTUAL}/results.json'
PRIMARY_TERMINAL = f'{PRIMARY_ACTUAL}/terminal.json'
REPLAY = f'{PRIMARY}/independent_actual_preparation'
DIRECT_PROT = 'runs/material_cohort_direct_v1/protocol.json'
SMALL_CSV = 'data/derived/calibration_small_18um.csv'
PINS = {
    PRIMARY_PROT: '9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94',
    PRIMARY_PREFLIGHT: 'f225616c4274a856576f587d40a93c89c64f88a9e4b3f1522d9f026434477211',
    PRIMARY_PRELAUNCH: 'df6daf32ca582e03ce7d4ab2ecc97e4ce44ca6420fe85134f352b5a64aa21b26',
    PRIMARY_DRIVER: '87ca09b46f9569903276923d937bef63d00bba9ad1b99d78b7e1708168f7e631',
    f'{REPLAY}/replay_frozen.py': 'e1ff63fddd6981f9ef8f05183c884537546cb66d59b9ec063911b22a1c50d27d',
    f'{REPLAY}/protocol.json': 'e1abffa662a3a64f024e0c6f54504af256cd04dcac01b44144b3757e62ffde99',
    f'{REPLAY}/preflight.json': '8884d9b6aeb78aa1f423a6a05094bd452b1684722ae6b7c23cb8a976d4a5d199',
    DIRECT_PROT: '1ed458c8f0f7b9c9bacfdf0c04b907d474d06e5607ced5956ecc2d2f1adac953',
    SMALL_CSV: '1aa53b6e6ce0973faf0b038ae8d4059059b05c02bb87ed9698e24d1dea7b1518',
}
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')
VARIANTS = ('current_disk__spatial_max', 'current_disk__regional_mean',
            'current_radial_line__spatial_max', 'current_radial_line__regional_mean')
LENGTH = Q('33.5')
TOLERANCE = Q('0.001')
ROOT_STEPS = 80
CELLS = 8192
MANIFEST_STATUS = 'COMPLETE_44_PRIMARY_CHECKPOINTS'
REVIEW_STATUS = 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
SCOPE = ('Secondary sensitivity to four fixed static observer conventions on the same 34 saved small fields. '
         'The minimum required common allowance refers to those stored fields and the given recorded values. '
         'The primary decision remains the value in its own receipt. '
         'Other observation operators, source registration and timing, geometry uncertainty, FE root error and PDE error lie outside this analysis.')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return digest(path.read_bytes())


def read(path):
    return json.loads(path.read_bytes())


def pinned(path, expected):
    require(isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected) is not None,
            f'malformed SHA-256 for {path.name}')
    payload = path.read_bytes()
    require(digest(payload) == expected, f'SHA-256 mismatch for {path.name}')
    return json.loads(payload)


def save(path, document):
    encoded = json.dumps(document, indent=2, allow_nan=False) + '\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(encoded)
    return sha(path)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def all_true(checks):
    return isinstance(checks, dict) and bool(checks) and all(value is True for value in checks.values())


def good_review(review):
    status = review.get('status')
    return (isinstance(status, str) and (status == 'PASS' or status.startswith('PASS_INDEPENDENT_'))
            and all_true(review.get('checks')))


def binary(value):
    require(type(value) in (int, float), 'stored scalar must be numeric')
    return Q(value)


def plan_checkpoints(primary):
    require(tuple(case['endpoint_id'] for case in primary['endpoints']) == ENDPOINTS,
            'changed primary endpoint order')
    plan = []
    for case in primary['endpoints']:
        endpoint = case['endpoint_id']
        require(len(case['small_records']) == 17 and len({r['record_id'] for r in case['small_records']}) == 17
                and len(case['large_checkpoints']) == 4, 'incomplete primary frozen plan')
        entries = [('small_zero', 'small_zero', case['small_zero'])]
        entries.extend((f'small_{i:02d}', 'small_record', row)
                       for i, row in enumerate(case['small_records'], 1))
        large_roles = ('large_zero', 'large_checkpoint', 'large_checkpoint', 'large_final')
        entries.extend((f'large_{i:02d}', role, step)
                       for i, (role, step) in enumerate(zip(large_roles, case['large_checkpoints'])))
        for label, role, step in entries:
            stage = f'{PRIMARY_ACTUAL}/{endpoint}/{label}'
            plan.append({'endpoint_id': endpoint, 'checkpoint': label, 'role': role, 'step': step,
                         'field': f'{stage}/field.json', 'verification': f'{stage}/verification.json'})
    require(len(plan) == 44 and sum(row['role'] == 'small_record' for row in plan) == 34,
            'checkpoint plan failed completeness check')
    return plan


def prepare_protocol():
    for path, expected in PINS.items():
        require(sha(ROOT / path) == expected, f'changed fixed input: {path}')
    primary = pinned(ROOT / PRIMARY_PROT, PINS[PRIMARY_PROT])
    legacy = pinned(ROOT / DIRECT_PROT, PINS[DIRECT_PROT])
    require(primary['source_sha256'] == PINS[PRIMARY_DRIVER]
            and primary['observable_variant'] == VARIANTS[0]
            and tuple(legacy['observation_variants']) == VARIANTS
            and Q(primary['large_observable']['length_unit_um_exact']) == LENGTH,
            'observer authority or common length drift')
    sources = {'scripts/quadratic_endpoint_observer.py': observer,
               'scripts/annulus_observer.py': intervals, 'scripts/joint_calibration_contract.py': contract}
    source_hashes = {name: primary['source_hashes'][name] for name in sources}
    for name, module in sources.items():
        require(Path(module.__file__).resolve() == ROOT / name
                and sha(ROOT / name) == source_hashes[name], f'changed observer dependency: {name}')
    extractor_authority = 'scripts/joint_gain_calibration.py'
    source_hashes[extractor_authority] = primary['source_hashes'][extractor_authority]
    require(sha(ROOT / extractor_authority) == source_hashes[extractor_authority],
            'historical four-variant extraction source drift')
    observations = {rid: {key: str(value) for key, value in observed.items()}
                    for rid, observed in contract.load_observations().items()}
    for case in primary['endpoints']:
        require(observations == {r['record_id']: r['observed_exact'] for r in case['small_records']},
                'source rows differ from primary frozen observations')
    return {
        'version': 1, 'status': 'PREDECLARED_SECONDARY_OBSERVER_SWEEP', 'scope': SCOPE,
        'source_sha256': sha(SOURCE), 'source_hashes': source_hashes, 'input_hashes': PINS,
        'endpoint_order': list(ENDPOINTS), 'variant_order': list(VARIANTS),
        'primary_kept': VARIANTS[0],
        'observer': {'length_unit_um_exact': str(LENGTH), 'tolerance_um_exact': str(TOLERANCE),
                     'maximum_cells': CELLS, 'root_steps': ROOT_STEPS,
                     'reference': 'Flat reference height from each unchanged field, and basal reference zero.',
                     'ROI': 'Closed current-radius interval from axis to binary64 initial patch radius times 33.5 µm.',
                     'signed_components': True, 'one_static_state_per_record': True,
                     'shared_maximum': 'Same signed Eulerian height supremum in all four variants.',
                     'mean_alternative': 'Line or disk integration for means; shared maxima for spatial_max.'},
        'conditional_source_allowance_um_exact': primary['selection']['conditional_source_allowance_um_exact'],
        'threshold_rule': 'Enclose max absolute residual over all 68 entries, then take max across endpoints. '
                          'Also give the max over the four predeclared conventions for this fixed pair.',
        'frozen_checkpoint_plan': plan_checkpoints(primary),
        'admission': {'independent_actual_primary_review_sha256': 'required external CLI pin',
                      'whole_44_field_manifest_sha256': 'required external CLI pin',
                      'independent_secondary_prelaunch_review_sha256': 'required external CLI pin',
                      'manifest_status': MANIFEST_STATUS,
                      'analysed_roles': ['small_record']},
        'additional_force_evaluations': False, 'additional_equilibrium_solves': False,
        'retuning': False, 'heldout_large_response_read': False,
    }


def thin_checkpoint(row):
    return {name: row[name] for name in ('endpoint_id', 'checkpoint', 'role', 'step', 'field', 'field_sha256',
                                        'verification', 'verification_sha256')}


def planned_manifest(results, plan):
    require(results['status'] == 'COMPLETE_CONDITIONAL_STORED_FIELD_COMPARISON'
            and results['verified_outcome_fields'] == 36 and results['small_observer_intervals'] == 136
            and tuple(case['endpoint_id'] for case in results['endpoints']) == ENDPOINTS,
            'primary results incomplete or from another pair')
    reviewed = []
    for case in results['endpoints']:
        require(case['status'] == 'COMPLETE_CONDITIONAL_STORED_FIELDS'
                and case['verified_final_outcome_fields'] == 18 and len(case['checkpoints']) == 22
                and case['small_residual']['all68_residuals']['constraints'] == 68,
                'primary endpoint lacks complete successful terminal data')
        reviewed.extend({'endpoint_id': case['endpoint_id'], **row} for row in case['checkpoints'])
    for expected, row in zip(plan, reviewed, strict=True):
        require(all(row[k] == v for k, v in expected.items()),
                f'mixed record/checkpoint in primary result: {expected["checkpoint"]}')
        require(all(isinstance(row[key], str) and re.fullmatch('[0-9a-f]{64}', row[key])
                    for key in ('field_sha256', 'verification_sha256')), 'malformed field or verification digest')
    return list(map(thin_checkpoint, reviewed))


def require_primary_review(review, manifest, results_sha, terminal_sha, expected_rows):
    require(review.get('status') == REVIEW_STATUS and review['checkpoints'] == 44
            and review['independently_recomputed_small_residuals'] == 136
            and review['new_equilibrium_solves'] is False and review['heldout_response_tables_read'] is False,
            'independent primary actual review did not pass')
    hashes = review['all_hashes']
    require(review['actual_results_sha256'] == results_sha and hashes[PRIMARY_RESULT] == results_sha
            and hashes[PRIMARY_TERMINAL] == terminal_sha
            and all(hashes.get(path) == h for path, h in PINS.items()),
            'independent primary review does not bind the primary results, terminal and frozen inputs')
    require(manifest.get('status') == MANIFEST_STATUS
            and manifest.get('primary_results_sha256') == results_sha
            and manifest.get('primary_terminal_sha256') == terminal_sha,
            'field manifest belongs to another or incomplete primary run')
    require(len(manifest['fields']) == len(expected_rows) == 44
            and list(map(thin_checkpoint, manifest['fields'])) == expected_rows,
            '44-field manifest does not exactly match the frozen primary checkpoint order and paths')
    reviewed_fields = []
    require(tuple(endpoint['endpoint_id'] for endpoint in review['endpoints']) == ENDPOINTS,
            'independent replay endpoint order drift')
    for endpoint in review['endpoints']:
        require(len(endpoint['checkpoints']) == 22, 'independent primary replay coverage incomplete')
        reviewed_fields.extend((endpoint['endpoint_id'], row['checkpoint'], row['field_sha256'])
                               for row in endpoint['checkpoints'])
    require(reviewed_fields == [(row['endpoint_id'], row['checkpoint'], row['field_sha256']) for row in expected_rows]
            and all(hashes.get(row[kind]) == row[kind + '_sha256']
                    for row in expected_rows for kind in ('field', 'verification')),
            'manifest field or verification lacks matching independent replay digest')


def verify_primary_supervisor(terminal, results, actual_sha, reviewed_hashes):
    require(terminal['status'] == 'FINISHED' and terminal['results_sha256'] == actual_sha
            and terminal['bindings'] == results['bindings'], 'primary supervisor did not finish')
    links = results['bindings']
    expected = {'source_sha256': PINS[PRIMARY_DRIVER], 'protocol_sha256': PINS[PRIMARY_PROT],
                'preflight_sha256': PINS[PRIMARY_PREFLIGHT], 'prelaunch_review_sha256': PINS[PRIMARY_PRELAUNCH]}
    require(all(links.get(key) == value for key, value in expected.items()), 'mixed primary run bindings')
    start_path = f'{PRIMARY_ACTUAL}/started.json'
    started = pinned(ROOT / start_path, reviewed_hashes[start_path])
    require(started['pid'] == terminal['pid'] and started['bindings'] == links
            and started['serial_endpoint_order'] == list(ENDPOINTS), 'primary supervisor start drift')
    for endpoint_id, worker in zip(ENDPOINTS, terminal['workers'], strict=True):
        folder = f'{PRIMARY_ACTUAL}/{endpoint_id}'
        require(worker['endpoint_id'] == endpoint_id and worker['supervisor_pid'] == terminal['pid']
                and type(worker['observed_exit_code']) is int and worker['observed_exit_code'] == 0,
                'primary worker did not exit zero')
        exit_path = f'{folder}/supervisor_os_exit.json'
        require(worker['os_exit_receipt_sha256'] == reviewed_hashes[exit_path], 'worker OS exit digest drift')
        exit_observed = pinned(ROOT / exit_path, reviewed_hashes[exit_path])
        require(exit_observed == {k: v for k, v in worker.items() if k != 'os_exit_receipt_sha256'},
                'mixed primary worker OS exit')
        terminal_path = f'{folder}/terminal.json'
        require(exit_observed['terminal_sha256'] == reviewed_hashes[terminal_path], 'worker terminal digest drift')
        child = pinned(ROOT / terminal_path, reviewed_hashes[terminal_path])
        require(child['status'] == 'FINISHED' and child['pid'] == worker['pid']
                and child['supervisor_pid'] == started['pid'] and child['binds'] == links
                and child['endpoint_id'] == endpoint_id, 'mixed or failed primary worker terminal')
        result_path = f'{folder}/results.json'
        require(child['results_sha256'] == reviewed_hashes[result_path], 'worker result digest drift')
        reply = pinned(ROOT / result_path, reviewed_hashes[result_path])
        require(reply == next(case for case in results['endpoints'] if case['endpoint_id'] == endpoint_id),
                'primary worker and supervisor results disagree')


def admit_actual(protocol, review_path, review_sha, manifest_path, manifest_sha):
    review = pinned(review_path, review_sha)
    require(review.get('status') == REVIEW_STATUS, 'independent primary actual review did not pass')
    independent_terminal = read(review_path.with_name('terminal.json'))
    require(independent_terminal['status'] == 'FINISHED' and type(independent_terminal['exit_code']) is int
            and independent_terminal['exit_code'] == 0 and independent_terminal['receipt_sha256'] == review_sha,
            'independent primary replay terminal does not authenticate the reviewed receipt')
    manifest = pinned(manifest_path, manifest_sha)
    results_sha = review['actual_results_sha256']
    terminal_sha = review['all_hashes'][PRIMARY_TERMINAL]
    terminal = pinned(ROOT / PRIMARY_TERMINAL, terminal_sha)
    require(terminal['status'] == 'FINISHED', 'primary terminal did not finish')
    results = pinned(ROOT / PRIMARY_RESULT, results_sha)
    expected = planned_manifest(results, protocol['frozen_checkpoint_plan'])
    require_primary_review(review, manifest, results_sha, terminal_sha, expected)
    verify_primary_supervisor(terminal, results, results_sha, review['all_hashes'])
    # No field payload is read until the whole successful primary run and its independent manifest match.
    for row in expected:
        for kind in ('field', 'verification'):
            require(sha(ROOT / row[kind]) == row[kind + '_sha256'], f'changed whole-manifest bytes: {row[kind]}')
    return results, expected, {'primary_actual_review_sha256': review_sha,
                               'whole_44_field_manifest_sha256': manifest_sha,
                               'primary_results_sha256': results_sha, 'primary_terminal_sha256': terminal_sha}


def extract_variants(observed):
    require(observed['status'] == 'CERTIFIED_STATIC_INTERVALS'
            and observed['maximum_shortening_certificate']['requested_width_achieved'] is True,
            'static observer did not certify requested maximum width')
    max_height = observed['four_endpoints']['maximum_Eulerian_height_shortening']
    require(max_height.hi - max_height.lo <= TOLERANCE, 'maximum observer tolerance exceeded')
    variants = {}
    for variant in VARIANTS:
        weighting, aggregation = variant.split('__')
        means = observed['mean_alternatives'][weighting]
        if aggregation == 'spatial_max':
            a = observed['spatial_max_components']['apical_downward']
            b = observed['spatial_max_components']['basal_upward']
        else:
            a = means['apical_regional_mean_downward']
            b = means['basal_regional_mean_upward']
        selected = (max_height, means['height_regional_mean_shortening'], a, b)
        require(all(i.lo <= i.hi and i.lc and i.rc for i in selected), 'nonclosed or reversed observer interval')
        variants[variant] = {key: [str(v.lo), str(v.hi)] for key, v in zip(contract.OBSERVABLES, selected)}
    return variants


def measure_field(state, frozen_parameters, step):
    require(state['converged'] is True and state['reached_target'] is True
            and state['parameters'] == frozen_parameters
            and binary(state['activation_strength']) == Q(step['solver_strength_binary64_exact']),
            'stored small field does not reproduce its frozen material and loading')
    roi = (Q(0), binary(frozen_parameters['patch_radius']) * LENGTH)
    ref = binary(frozen_parameters['height']) * LENGTH
    result = observer.observe_state(state, roi_um=roi, length_unit_um=LENGTH,
                                    reference_apical_um=ref, reference_basal_um=Q(0),
                                    weighting='current_disk', component_aggregation='spatial_max',
                                    tolerance_um=TOLERANCE, maximum_cells=CELLS, root_steps=ROOT_STEPS)
    return {'variant_observables_um_exact': extract_variants(result),
            'observer_contract': intervals.serial(result['contract']),
            'maximum_shortening_certificate': intervals.serial(result['maximum_shortening_certificate']),
            'positive_thickness_lower_um_exact': str(result['positive_thickness_lower_um'])}


def threshold(records, predictions, delta):
    require(len(records) == len(predictions) == 17, '17 predictions required per endpoint and variant')
    obs = {r['record_id']: r['observed_exact'] for r in records}
    result = contract.evaluate(obs, predictions, {key: str(delta) for key in contract.OBSERVABLES})
    require(len(obs) == 17 and result['constraints'] == 68, 'missing exact observer constraint')
    bounds = [max(map(Q, result[key].values())) for key in (
        'discrepancy_lower_requirements_exact', 'discrepancy_sufficient_requirements_exact')]
    status = {'IN': 'ALLOWANCE_SATISFIED', 'OUT': 'ALLOWANCE_EXCEEDED', 'UNRESOLVED': 'UNRESOLVED'}[result['status']]
    return {'status': status, 'required_allowance_um_exact': list(map(str, bounds)),
            'required_allowance_um_display': list(map(float, bounds)),
            'declared_allowance_um_exact': str(delta), 'all68_residuals': result}


def combined_requirement(rows, delta):
    pairs = [list(map(Q, row['required_allowance_um_exact'])) for row in rows]
    require(bool(pairs) and all(0 <= a <= b for a, b in pairs), 'bad common allowance intervals')
    lower, upper = map(max, zip(*pairs))
    status = 'ALLOWANCE_SATISFIED' if upper <= delta else (
        'ALLOWANCE_EXCEEDED' if lower > delta else 'UNRESOLVED')
    return {'status_at_frozen_allowance': status,
            'required_allowance_um_exact': [str(lower), str(upper)],
            'required_allowance_um_display': [float(lower), float(upper)]}


def run(protocol, primary_results, manifest_rows, bindings):
    primary = pinned(ROOT / PRIMARY_PROT, PINS[PRIMARY_PROT])
    cases = {case['endpoint_id']: case for case in primary['endpoints']}
    primary_receipts = {case['endpoint_id']: case['small_residual'] for case in primary_results['endpoints']}
    fixed_plan = {(r['endpoint_id'], r['checkpoint']): r for r in protocol['frozen_checkpoint_plan']}
    delta = Q(protocol['conditional_source_allowance_um_exact'])
    per_endpoint = {}
    receipts = []
    for eid in ENDPOINTS:
        case = cases[eid]
        predictions = {v: {} for v in VARIANTS}
        for row in manifest_rows:
            if row['endpoint_id'] != eid or row['role'] != 'small_record':
                continue
            step = fixed_plan[eid, row['checkpoint']]['step']
            state = pinned(ROOT / row['field'], row['field_sha256'])
            measurement = measure_field(state, case['materials']['small']['parameters'], step)
            rid = step['record_id']
            for variant, values in measurement['variant_observables_um_exact'].items():
                predictions[variant][rid] = values
            path = ACTUAL / eid / f"{row['checkpoint']}.json"
            h = save(path, {**row, 'record_id': rid, **measurement})
            receipts.append({'path': path.relative_to(ROOT).as_posix(), 'sha256': h})
        summaries = {variant: threshold(case['small_records'], predictions[variant], delta)
                     for variant in VARIANTS}
        require(summaries[VARIANTS[0]] == primary_receipts[eid],
                f'primary observer replay does not reproduce its exact archived result: {eid}')
        per_endpoint[eid] = summaries
    require(len(receipts) == 34, 'expected 34 measured stored small fields')
    per_variant = {v: combined_requirement([per_endpoint[e][v] for e in ENDPOINTS], delta) for v in VARIANTS}
    return {'status': 'COMPLETE_SECONDARY_STATIC_OBSERVER_SWEEP', 'scope': SCOPE, 'bindings': bindings,
            'by_endpoint_and_variant': per_endpoint, 'common_pair_requirement_by_variant': per_variant,
            'common_pair_requirement_across_all_four_conventions': combined_requirement(list(per_variant.values()), delta),
            'frozen_primary_source_allowance_um_exact': str(delta),
            'primary_threshold_outcome_by_endpoint': {eid: primary_receipts[eid]['status'] for eid in ENDPOINTS},
            'primary_replay_exactly_matches': True, 'measured_small_fields': 34,
            'small_constraint_evaluations': 2 * len(VARIANTS) * 68, 'field_observer_receipts': receipts,
            'additional_equilibria': 0, 'additional_force_evaluations': 0, 'retuning': False}


def prelaunch_bindings(path, review_sha):
    expected = {'source_sha256': sha(SOURCE), 'protocol_sha256': sha(PREP / 'protocol.json'),
                'preflight_sha256': sha(PREP / 'preflight.json'), 'primary_protocol_sha256': PINS[PRIMARY_PROT]}
    approved = pinned(path, review_sha)
    require(good_review(approved) and all(approved.get(key) == v for key, v in expected.items()),
            'invalid independent secondary prelaunch bindings')
    receipt = pinned(PREP / 'preflight.json', expected['preflight_sha256'])
    require(receipt['status'] == 'PASS' and all_true(receipt['checks'])
            and receipt['source_sha256'] == expected['source_sha256']
            and receipt['protocol_sha256'] == expected['protocol_sha256']
            and receipt['primary_actual_fields_opened'] == 0, 'secondary preparation preflight drift')
    require(sha(PREP / 'transfer_ambiguity_observer_sweep_frozen.py') == expected['source_sha256'],
            'secondary frozen driver bytes changed')
    return {**expected, 'secondary_prelaunch_review_sha256': review_sha}


def launch(args):
    prelaunch = prelaunch_bindings(args.prelaunch_review, args.prelaunch_review_sha256)
    protocol = pinned(PREP / 'protocol.json', prelaunch['protocol_sha256'])
    require(protocol == prepare_protocol(), 'frozen secondary protocol drift')
    primary_results, manifest, primary_bindings = admit_actual(
        protocol, args.actual_review, args.actual_review_sha256, args.field_manifest, args.field_manifest_sha256)
    bindings = {**prelaunch, **primary_bindings}
    ACTUAL.mkdir(exist_ok=False)
    pid = os.getpid()
    start = save(ACTUAL / 'started.json', {'pid': pid, 'bindings': bindings, 'utc': stamp()})
    terminal = {'status': 'FAILED', 'pid': pid, 'start_sha256': start, 'bindings': bindings}
    try:
        terminal['results_sha256'] = save(ACTUAL / 'results.json', run(protocol, primary_results, manifest, bindings))
        terminal['status'] = 'FINISHED'
    except BaseException as error:
        terminal['error'] = {'type': type(error).__name__, 'message': str(error)}
        save(ACTUAL / 'failure.json', {'error': terminal['error'], 'traceback': traceback.format_exc()})
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


def fixture(protocol):
    # On [0,1], these controls give apical shortening r² and basal 2-r/2.
    top = [[(0, 10), (Q(1, 2), Q(39, 4)), (1, 9)]]
    bottom = [[(0, 2), (Q(1, 2), Q(7, 4)), (1, Q(3, 2))]]
    measured = extract_variants(observer.observe_graphs(
        top, bottom, roi_um=(0, 1), reference_apical_um=10, tolerance_um=TOLERANCE))
    keys = contract.OBSERVABLES
    targets = (('5/2', '13/6', '1', '2'), ('5/2', '13/6', '1/2', '5/3'),
               ('5/2', '25/12', '1', '2'), ('5/2', '25/12', '1/3', '7/4'))
    plan = protocol['frozen_checkpoint_plan']
    rows = [thin_checkpoint({**row, 'field_sha256': digest(row['field'].encode()),
                            'verification_sha256': digest(row['verification'].encode())}) for row in plan]
    rs, ts = ('a' * 64, 'b' * 64)
    small_replay = [
        {'endpoint_id': e, 'checkpoints': [{'checkpoint': row['checkpoint'],
                                          'field_sha256': row['field_sha256']}
                                         for row in rows if row['endpoint_id'] == e]}
        for e in ENDPOINTS]
    hashes = {**PINS, PRIMARY_RESULT: rs, PRIMARY_TERMINAL: ts,
              **{r[kind]: r[kind + '_sha256'] for r in rows for kind in ('field', 'verification')}}
    review = {'status': REVIEW_STATUS, 'actual_results_sha256': rs,
              'checkpoints': 44, 'independently_recomputed_small_residuals': 136,
              'new_equilibrium_solves': False, 'heldout_response_tables_read': False,
              'all_hashes': hashes, 'endpoints': small_replay}
    manifest = {'status': MANIFEST_STATUS, 'primary_results_sha256': rs,
                'primary_terminal_sha256': ts, 'fields': rows}
    require_primary_review(review, manifest, rs, ts, rows)
    actual = {'status': 'COMPLETE_CONDITIONAL_STORED_FIELD_COMPARISON',
              'verified_outcome_fields': 36, 'small_observer_intervals': 136,
              'endpoints': [
                  {'endpoint_id': e, 'status': 'COMPLETE_CONDITIONAL_STORED_FIELDS',
                   'verified_final_outcome_fields': 18, 'small_residual': {'all68_residuals': {'constraints': 68}},
                   'checkpoints': [{k: r[k] for k in r if k != 'endpoint_id'}
                                   for r in rows if r['endpoint_id'] == e]} for e in ENDPOINTS]}
    swapped = [{**plan[0], 'step': plan[1]['step']}, *plan[1:]]
    primary = read(ROOT / PRIMARY_PROT)
    source = primary['endpoints'][0]['small_records']
    preds = {row['record_id']: {k: [str(Q(v) + 1)] * 2 for k, v in row['observed_exact'].items()}
             for row in source}
    exact1 = threshold(source, preds, Q(1))
    strict = threshold(source, preds, Q(99, 100))
    record_ids = {row['record_id'] for row in source}
    return {
        'analytic_four_variants_distinguish_weights_and_aggregations': all(
            measured[variant] == {k: [v, v] for k, v in zip(keys, expected)}
            for variant, expected in zip(VARIANTS, targets)),
        'combined_maximum_is_not_sum_of_component_maxima': (
            measured[VARIANTS[0]][keys[0]] == ['5/2', '5/2']
            and sum(Q(measured[VARIANTS[0]][k][0]) for k in keys[2:]) == 3),
        'whole_44_manifest_and_34_small_fields_fixed': len(rows) == 44
            and sum(row['role'] == 'small_record' for row in rows) == 34 and len(record_ids) == 17
            and planned_manifest(actual, plan) == rows and rejects(lambda: planned_manifest(actual, swapped)),
        'exact_signed_residual_and_threshold_cases': exact1['required_allowance_um_exact'] == ['1', '1']
            and exact1['status'] == 'ALLOWANCE_SATISFIED' and strict['status'] == 'ALLOWANCE_EXCEEDED'
            and strict['all68_residuals']['constraints'] == 68,
        'max_across_endpoints_or_variants_uses_both_bounds': combined_requirement([
            {'required_allowance_um_exact': ['1', '4']}, {'required_allowance_um_exact': ['3', '7/2']}], Q(3))
            == {'status_at_frozen_allowance': 'UNRESOLVED', 'required_allowance_um_exact': ['3', '4'],
                'required_allowance_um_display': [3., 4.]},
        'reject_partial_and_mixed_field_manifest': rejects(lambda: require_primary_review(
            review, {**manifest, 'fields': rows[:-1]}, rs, ts, rows)) and rejects(
            lambda: require_primary_review(review, {**manifest, 'fields': rows[:-1] + [rows[0]]}, rs, ts, rows)),
        'reject_primary_review_with_wrong_terminal_or_manifest': rejects(lambda: require_primary_review(
            {**review, 'all_hashes': {**hashes, PRIMARY_TERMINAL: rs}}, manifest, rs, ts, rows))
            and rejects(lambda: require_primary_review(review, {**manifest, 'primary_results_sha256': ts}, rs, ts, rows)),
        'reject_failed_review_or_wrong_reviewed_field_digest': rejects(lambda: require_primary_review(
            {**review, 'status': 'FAIL'}, manifest, rs, ts, rows)) and rejects(lambda: require_primary_review(
            {**review, 'all_hashes': {**hashes, rows[0]['verification']: rs}}, manifest, rs, ts, rows)),
    }


def prepare():
    PREP.mkdir(parents=True, exist_ok=False)
    snapshot = PREP / 'transfer_ambiguity_observer_sweep_frozen.py'
    with snapshot.open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        # Ensure preparation does not open even a primary started/terminal file.
        original_open = Path.open

        def permitted_open(path, *args, **kwargs):
            resolved = path.resolve()
            require(not resolved.is_relative_to(ROOT / PRIMARY_ACTUAL), 'preflight opened primary actual payload')
            return original_open(path, *args, **kwargs)

        with (patch.object(Path, 'open', permitted_open),
              patch.object(observer, 'observe_state', side_effect=AssertionError('stored-field observer forbidden'))):
            protocol = prepare_protocol()
            protocol_sha = save(PREP / 'protocol.json', protocol)
            checks = fixture(protocol)
        checks['actual_output_absent'] = not ACTUAL.exists()
        status = 'PASS' if all_true(checks) else 'FAIL'
        preflight = save(PREP / 'preflight.json', {
            'status': status, 'source_sha256': sha(SOURCE), 'protocol_sha256': protocol_sha,
            'primary_protocol_sha256': PINS[PRIMARY_PROT], 'checks': checks, 'primary_actual_fields_opened': 0,
            'actual_observer_loops_run': 0, 'additional_equilibria': 0, 'additional_force_evaluations': 0})
        require(status == 'PASS', 'secondary observer sweep preflight failed')
        print(json.dumps({'status': status, 'source_sha256': sha(SOURCE),
                          'protocol_sha256': protocol_sha, 'preflight_sha256': preflight}), flush=True)
    except BaseException as error:
        save(PREP / 'failure.json', {'error': type(error).__name__, 'message': str(error),
                                     'traceback': traceback.format_exc(), 'utc': stamp()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest='stage', required=True)
    stages.add_parser('prepare')
    actual = stages.add_parser('run')
    for name in ('prelaunch-review', 'actual-review', 'field-manifest'):
        actual.add_argument('--' + name, required=True, type=Path)
        actual.add_argument('--' + name + '-sha256', required=True)
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        launch(args)


if __name__ == '__main__':
    main()

