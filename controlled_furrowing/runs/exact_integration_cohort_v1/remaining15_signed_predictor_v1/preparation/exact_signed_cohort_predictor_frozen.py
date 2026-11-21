"""Prepare/replay the signed exact-integral predictor residual for record D."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
from time import monotonic
import traceback
from unittest.mock import patch

from flint import ctx as precision
import numpy as np

import exact_cohort_context as C
import exact_signed_predictor_L as signed
import parametric_root_cells as P
import review_cohort_exact_record as D


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
BASE = 'runs/exact_integration_cohort_v1'
OUT = ROOT / BASE / 'remaining15_signed_predictor_v1'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
RECORD = 'optoGEF2_small_col_D'
ORIGINAL = D._record_dir(RECORD)
REDUCTION = f'{BASE}/remaining15_centered_reduction_preparation'
GLOBAL = f'{REDUCTION}/independent_global_actual/{RECORD}'
GLOBAL_PREP = f'{REDUCTION}/independent_global_preparation'
PRODUCER = f'{REDUCTION}/{RECORD}/results.json'
PRODUCER_TERMINAL = f'{REDUCTION}/{RECORD}/terminal.json'
REVIEW = f'{GLOBAL}/receipt.json'
OS_EXIT = f'{GLOBAL}/os_exit_observation.json'
INPUTS = f'{GLOBAL}/signed_inputs.json'
ADMISSION = f'{REDUCTION}/{RECORD}_input_admission.json'
PIN = {
    'scripts/exact_signed_predictor_L.py': '00cebb37ed7cd23152ee46a659ff9878a32f62733faa2964d3b051fa7678dd1b',
    'scripts/parametric_root_cells.py': '362a19548aa60c4234017316e22deef4a80e93b73c725f47132573a1a9103b60',
    'scripts/algebraic_root_feasibility.py': '7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
    'scripts/exact_cohort_context.py': 'c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f',
    'scripts/review_cohort_exact_record.py': '67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285',
    f'{ORIGINAL}/protocol.json': 'e3165e9989d6e43b6f8b6819ea22496dc1f8598c091ffa02a02183ee6476c522',
    PRODUCER: 'a30067ae17c8d85418bf6002822a665a1737bfd6103aa92fa7ceef2c785e450e',
    PRODUCER_TERMINAL: '39f2f72cc770f13a2754d5652aeb752e534baff937ff8b2bff6790dbdf3f542d',
    REVIEW: 'ac2e091abaf8175d3d50b952b30bda8ebfcfed1cdb73ae3297bb8813df62ae4c',
    OS_EXIT: '517ba38a95d16e0714ce998f03e7cdc930943b3b0655e14ca23b738202df17aa',
    INPUTS: 'f6414069a462400cf32b55ddad59fce4625a06a1231b45027d6af95c643fb176',
    ADMISSION: '6355d543a739aeaf83d6592acd0a2e46a671ebbe16e6183a2f96dac01d9001e3',
    f'{GLOBAL_PREP}/review_centered_cohort_global_frozen.py':
        'a0090d6c28e798a220de09ecc91b255ae1d00ffa096a218a13885b664a34d5e9',
    f'{GLOBAL_PREP}/protocol.json': '47b00fd7aa4c744e36b4f859b144351eec98fca4ea259f7260d184beea88c028',
    f'{GLOBAL_PREP}/preflight.json': 'bdb446c4f9cf661455fd00017e7e77df75bb6f557ab8b2d32374e5a2e43e11d3',
}
RANK = 611
GAIN = 'independent_affine_gain_RF_intervals_exact'
SHARED = ('point_bulk_RF_intervals_exact', 'retained_surface_RF_intervals_exact',
          'RH_bulk_S_signed_intervals_exact', 'closed_five_offsets_exact', 'parameter_order')
SCOPE = ('Record D only. One uniform signed RF_I enclosure at the parameter-matched predictors on the original '
         'closed cell C. The whole-cube derivative bound comes from the independently reviewed full 168-component '
         'assembly. A strict exclusion applies only to the original radius cube around each predictor on this C.')


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def bootstrap():
    precision.prec = 640
    P.core.A.cache_clear()
    D.require(signed.P is P and signed.D is D and signed.C is C and precision.prec == 640,
              'shared pure signed-kernel aliases and Arb640 required')


def _bound(name, links, parse=True):
    return D._bound(name, PIN[name], links, parse=parse)


def check_authority(review, os_exit, terminal, producer, payload, record):
    D._terminal(terminal, PIN[PRODUCER])
    D.require(type(terminal['exit_code']) is int and terminal['exit_code'] == 0,
              'successful producer terminal at original D required')
    expected_os = {'status': 'OBSERVED_ZERO_EXIT_INDEPENDENT_GLOBAL_CENTERED_REPLAY',
                   'record_id': record, 'review_file': REVIEW, 'review_sha256': PIN[REVIEW],
                   'signed_inputs_file': INPUTS, 'signed_inputs_sha256': PIN[INPUTS],
                   'producer_results_sha256': PIN[PRODUCER], 'producer_terminal_sha256': PIN[PRODUCER_TERMINAL],
                   'producer_pid': terminal['pid']}
    D.require(all(os_exit[k] == v for k, v in expected_os.items())
              and type(os_exit['observed_exit_code']) is int and os_exit['observed_exit_code'] == 0,
              'same-record independently observed OS exit required')
    expected_review = {'record_id': record, 'producer_results_sha256': PIN[PRODUCER],
                       'producer_terminal_sha256': PIN[PRODUCER_TERMINAL],
                       'future_signed_inputs_sha256': PIN[INPUTS], 'actual_admission_sha256': PIN[ADMISSION],
                       'review_source_sha256': PIN[f'{GLOBAL_PREP}/review_centered_cohort_global_frozen.py'],
                       'review_protocol_sha256': PIN[f'{GLOBAL_PREP}/protocol.json'],
                       'review_preflight_sha256': PIN[f'{GLOBAL_PREP}/preflight.json'],
                       'global_F_entries': RANK, 'H_upper_entries': 6676}
    D.require(all(review[k] == v for k, v in expected_review.items())
              and review['status'] == 'PASS_INDEPENDENT_CENTERED_REDUCTION_UNRESOLVED'
              and review['root_certified'] is False and review['no_new_actual_mechanical_jobs'] is True
              and review['first_144_bulk_and_retained_24_surfaces'] is True
              and producer['status'] == 'UNREVIEWED_CENTERED_SUFFICIENT_TEST_UNRESOLVED'
              and producer['record_id'] == record and producer['bulk_components'] == 144
              and producer['retained_surface_components'] == 24
              and producer['admission_sha256'] == PIN[ADMISSION]
              and producer['actual_packet_review_sha256'] == review['actual_packet_review_sha256'],
              'same reviewed full centered reduction expected')
    D.require(payload['status'] == 'INDEPENDENT_EXACT_CENTERED_INPUT_ARRAYS' and payload['record_id'] == record
              and payload['reducer_exact_array_digest'] == review['exact_array_sha256']
              and all(D._digest(producer[k]) == h for k, h in review['exact_array_sha256'].items())
              and all(D._digest(producer[k]) == D._digest(payload[k]) for k in SHARED),
              'signed array links or complete reviewed producer array bytes differ')
    D.require(review['gain_slots'] == review['producer_gain_contains_exact_rational_nodal_gain_slots'] == 1656
              and review['producer_gain_overlapping_only_slots'] == 0
              and payload['raw_gain_operand_proven_to_cover_exact_box'] is True
              and payload['noncontaining_raw_gain_diagnostics'] == []
              and D._digest(payload['verified_producer_gain_RF_intervals_exact'])
                    == D._digest(producer['affine_gain_RF_intervals_exact'])
              and payload['declared_outward_dyadic_bits'] == 800
              and 3 < signed._interval(payload['gain_rational_pi_enclosure'])[0]
                    <= signed._interval(payload['gain_rational_pi_enclosure'])[1] < 4,
              'same independently rebuilt gain-site coverage on entire exact nodal box required')


def parse_operand(payload, context):
    keys = tuple(P.NAMES)
    D.require(payload['record_id'] == context.record_id
              and payload['same_closed_C_sha256'] == context.row['parameter_cell_sha256']
              and D._fraction(payload['same_radius_exact']) == context.old_radius
              and tuple(payload['parameter_order']) == keys == tuple(C.ORDER),
              'same D cell identity, original radius and five ordered parameters required')
    delta = [signed._interval(v) for v in payload['closed_five_offsets_exact']]
    cell = P.check_cell(context.row['parameter_cell'])
    D.require(cell == context.cell and len(delta) == len(keys) == 5
              and delta == [(cell[k][0] - context.base[k], cell[k][1] - context.base[k]) for k in keys]
              and all(lo <= 0 <= hi for lo, hi in delta)
              and context.drive == context.base['common_gain'] * context.dose,
              'exact original closed C, base, dose or offsets differ')
    bulk, surface, gain = [[signed._interval(v, dyadic=True) for v in payload[k]]
                            for k in (SHARED[0], SHARED[1], GAIN)]
    matrix = [[signed._interval(v, dyadic=True) for v in row] for row in payload[SHARED[2]]]
    D.require(len(bulk) == len(surface) == len(gain) == len(matrix) == RANK
              and all(len(row) == 5 for row in matrix), 'all 611 signed five-column rows required')
    return bulk, surface, gain, matrix, delta


def input_bundle():
    bootstrap()
    links = {}
    for name, module in ((rel(signed.__file__), signed), (rel(C.__file__), C), (rel(P.__file__), P),
                         (rel(P.core.__file__), P.core), (rel(D.__file__), D)):
        _bound(name, links, parse=False)
        D.require(sha(module.__file__) == PIN[name], 'imported review primitives drifted')
    original = _bound(f'{ORIGINAL}/protocol.json', links)
    context = C.load_context(RECORD)
    D.require(context.record_id == RECORD and context.row == original['row']
              and context.R.shape == (RANK, RANK) and context.u.shape == (RANK,)
              and context.S.shape == (RANK, 5) and all(np.isfinite(v).all() for v in (context.u, context.S, context.R)),
              'matched frozen D context, computed R and selected affine predictor required')
    for name, digest in context.input_hashes.items():
        D._bound(name, digest, links, parse=False)
    box = D._bound(original['coordinate_box_file'], original['coordinate_box_sha256'], links)
    theory = C.exact_node_box(context, context.old_radius)
    D.require(all(D._digest(box[k]) == D._digest(v) for k, v in theory.items()),
              'old D predictor and its explicit exact nodal box differ')
    old_eta, old_z, radius, provider = D._old_premises(context, original, links)
    review, os_exit, terminal, producer, payload = [_bound(name, links) for name in (REVIEW, OS_EXIT, PRODUCER_TERMINAL,
                                                                                     PRODUCER, INPUTS)]
    for name in (ADMISSION, f'{GLOBAL_PREP}/review_centered_cohort_global_frozen.py',
                 f'{GLOBAL_PREP}/protocol.json', f'{GLOBAL_PREP}/preflight.json'):
        _bound(name, links, parse=False)
    check_authority(review, os_exit, terminal, producer, payload, RECORD)
    operands = parse_operand(payload, context)
    totals = review['scalar_exact_readback']
    b = Q(totals['same_fullbox_b_exact'])
    z = old_z + b
    norm = Q(totals['exact_R_infinity_norm'])
    D.require(b >= 0 and norm > 0 and 0 <= z < 1
              and Q(payload['same_radius_exact']) == radius == context.old_radius
              and Q(payload['same_full_box_b_exact']) == b and Q(payload['same_total_Z_upper_exact']) == z
              and Q(payload['normalization_norm_exact']) == norm
              and Q(totals['old_cached_Z_exact']) == old_z and Q(totals['old_cached_eta_exact']) == old_eta
              and Q(totals['same_radius_exact']) == radius and Q(totals['total_Z_upper_exact']) == z
              and all(Q(producer[k]) == Q(v) for k, v in totals.items())
              and producer['same_closed_C_sha256'] == context.row['parameter_cell_sha256']
              and Q(producer['same_exact_source_dose']) == context.dose,
              'same independently reviewed original cached bounds and whole-box full-integral Z required')
    args = {'point_certificate_file': context.row['point_cache_file'],
            'predictor_file': str((ROOT / context.row['predictor_file']).resolve()),
            'dose_exact': str(context.dose), 'base_parameters_exact': dict(context.row['base_parameters_exact'])}
    D.require(context.row['point_cache_file'] in links and context.row['predictor_file'] in links,
              'exact selected point-cache and predictor bytes were not bound')
    return links, context, operands, args, old_z, radius, b, z, norm, provider


def classify(intervals, z, radius):
    D.require(bool(intervals) and all(a <= b for a, b in intervals), 'ordered signed row intervals required')
    upper = max(map(D._abs_upper, intervals))
    lower_by_row = [max(Q(0), a, -b) for a, b in intervals]
    row = max(range(len(intervals)), key=lower_by_row.__getitem__)
    lower = lower_by_row[row]
    lhs, exists, error = signed.strict_bound(upper, z, radius)
    threshold = (1 + z) * radius
    excluded = lower > threshold
    D.require(not (excluded and exists), 'inconsistent contraction and exclusion intervals')
    return {'eta_I_upper_exact': upper, 'uniform_RF_I_lower_exact': lower, 'lower_witness_row_zero_based': row,
            'strict_eta_budget_exact': (1 - z) * radius, 'selfmap_lhs_exact': lhs,
            'exclusion_threshold_exact': threshold, 'exclusion_margin_exact': lower - threshold,
            'sufficient_existence_test_passes': exists, 'sufficient_old_cube_exclusion_passes': excluded,
            'conditional_root_error_exact': error,
            'status': 'CANDIDATE_UNIFORM_EXACT_FE_ROOT_UNREVIEWED' if exists else (
                'CANDIDATE_ORIGINAL_MATCHED_CUBES_EXCLUDED_UNREVIEWED' if excluded else 'UNRESOLVED_SIGNED_PREDICTOR')}


def prep_doc(bundle):
    links, context, (_, _, _, _, offsets), args, old_z, radius, b, z, norm, provider = bundle
    return {'status': 'FROZEN_D_SIGNED_PREDICTOR_PREPARATION', 'record_id': RECORD, 'scope': SCOPE,
            'source_sha256': sha(SOURCE), 'input_hashes': links, 'precision_bits': 640,
            'parameter_order': list(P.NAMES), 'same_closed_C_sha256': context.row['parameter_cell_sha256'],
            'original_parameter_cell': context.row['parameter_cell'],
            'closed_parameter_offsets_exact': offsets, 'point_kernel_arguments': args,
            'unchanged_original_radius_exact': radius, 'unchanged_cached_Z_exact': old_z,
            'unchanged_fullbox_b_exact': b, 'unchanged_fullbox_z_exact': z,
            'exact_R_infinity_norm': norm, 'independent_cached_bound_provider': provider,
            'signed_gain_operand': GAIN, 'expected_finite_cache_site_count': 9504,
            'rigorous_existence_bound': 'max row upper |RF_I(p)| + z*r < r and z < 1',
            'rigorous_nonexistence_bound': 'one uniform row |RF_I(p)| lower > (1+z)*r',
            'exclusion_argument': 'A root in the same convex cube would require ||RF_I(p)||∞ <= (1+z)r by the whole-cube '
                                  'Lipschitz bound on I-RH_I, since RF_I=I-T for T=I-RF_I.',
            'conditional_cached_SPD_bindings': {key: D.PREPARATION_PINS[key] for key in (D.ENERGY_RESULT, D.ENERGY_REVIEW)},
            'bulk_and_retained_surface_semantics': 'Point bulk and signed bulk H*S on fixed C, independent exact-gain '
                                                  'force derivative, and whole-box retained surface defect. Cached anchor and '
                                                  'all cached and integral-defect slopes sum before applying each common offset.',
            'cached_kernel': 'Explicit D cache and S passed as arguments to prior reviewed pure point and remainder helper.',
            'actual_cached_site_assembly_executed': False}


def _raises(call):
    try:
        call()
    except ValueError:
        return True
    return False


def fixture(bundle):
    from copy import deepcopy
    _, c, _, _, _, r, _, z, _, _ = bundle
    payload = json.loads((ROOT / INPUTS).read_bytes())
    changed_cell = {**payload, 'same_closed_C_sha256': '0' * 64}
    changed_record = {**payload, 'record_id': 'optoGEF2_small_col_C'}
    changed_radius = dict(payload, same_radius_exact=str(2*r))
    changed_offsets = deepcopy(payload['closed_five_offsets_exact'])
    changed_offsets[2][0] = str(Q(changed_offsets[2][0]) - 1)
    changed_endpoint = [[*x] for x in payload[GAIN]]
    changed_endpoint[0] = ['1/3', '1']
    pass_case = classify([(Q(-1, 8), Q(1, 16))], Q(1, 2), Q(1))
    tie = classify([(Q(3, 2), Q(2))], Q(1, 2), Q(1))
    exclude = classify([(Q(-2), Q(-13, 8)), (Q(-1), Q(1))], Q(1, 2), Q(1))
    terminal = json.loads((ROOT / PRODUCER_TERMINAL).read_bytes())
    review = json.loads((ROOT / REVIEW).read_bytes())
    os_exit = json.loads((ROOT / OS_EXIT).read_bytes())
    producer = json.loads((ROOT / PRODUCER).read_bytes())
    changed_shared = {**payload, SHARED[0]: payload[SHARED[0]][:-1]}
    checks = signed.tests()
    checks.update({
        'authoritatively_selected_radius_and_positive_wholebox_bound': (
            Q(payload['same_radius_exact']) == r and Q(payload['same_total_Z_upper_exact']) == z and z < 1),
        'reject_cross_record_or_cell': _raises(lambda: parse_operand(changed_record, c))
            and _raises(lambda: parse_operand(changed_cell, c)),
        'reject_changed_closed_offset': _raises(lambda: parse_operand(
            {**payload, 'closed_five_offsets_exact': changed_offsets}, c)),
        'reject_non_dyadic_or_missing_signed_force_row': (
            _raises(lambda: parse_operand({**payload, GAIN: changed_endpoint}, c))
            and _raises(lambda: check_authority(review, os_exit, terminal, producer, changed_shared, RECORD))),
        'reject_changed_radius': _raises(lambda: parse_operand(changed_radius, c)),
        'strict_root_exclusion_and_boundary_classification': (
            pass_case['sufficient_existence_test_passes'] and pass_case['conditional_root_error_exact'] == Q(1, 4)
            and not tie['sufficient_old_cube_exclusion_passes'] and tie['exclusion_margin_exact'] == 0
            and exclude['sufficient_old_cube_exclusion_passes'] and exclude['exclusion_margin_exact'] == Q(1, 8)),
        'reject_mixed_review_terminal_or_merely_overlapping_gain': (
            _raises(lambda: check_authority(review, os_exit, {**terminal, 'pid': terminal['pid'] + 1}, producer, payload, RECORD))
            and _raises(lambda: check_authority({**review, 'producer_gain_overlapping_only_slots': 1},
                                                os_exit, terminal, producer, payload, RECORD))),
    })
    return checks


def prepare():
    PREP.mkdir(parents=True, exist_ok=False)
    with (PREP / 'exact_signed_cohort_predictor_frozen.py').open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        forbidden = AssertionError('preparation entered the actual cached point or variation')
        with (patch.object(P, 'build_point', side_effect=forbidden),
              patch.object(P, 'variation', side_effect=forbidden),
              patch.object(P.core, 'assemble_point', side_effect=forbidden),
              patch.object(P.core, 'local_terms', side_effect=forbidden)):
            bundle = input_bundle()
            checks = fixture(bundle)
        checks['actual_output_absent'] = not ACTUAL.exists()
        status = 'PASS' if D._yes(checks) else 'FAIL'
        digest = D._write(PREP / 'protocol.json', prep_doc(bundle))
        pf = D._write(PREP / 'preflight.json', {'status': status, 'record_id': RECORD, 'checks': checks,
                      'source_sha256': sha(SOURCE), 'protocol_sha256': digest,
                      'actual_cached_sites_assembled': 0, 'added_exact_integral_component_jobs': 0})
        D.require(status == 'PASS', 'D signed predictor preparation failed')
        print(json.dumps({'status': status, 'source_sha256': sha(SOURCE),
                          'protocol_sha256': digest, 'preflight_sha256': pf}), flush=True)
    except BaseException as error:
        D._write(PREP / 'failure.json', {'error': str(error), 'traceback': traceback.format_exc()})
        raise


def run(review_sha):
    current = {'source_sha256': sha(SOURCE), 'protocol_sha256': sha(PREP / 'protocol.json'),
               'preflight_sha256': sha(PREP / 'preflight.json'), 'global_reduction_review_sha256': PIN[REVIEW],
               'independent_signed_inputs_sha256': PIN[INPUTS], 'global_reduction_results_sha256': PIN[PRODUCER]}
    gate = D._bound(rel(PREP / 'prelaunch_independent_review.json'), review_sha, {})
    D.require(gate['status'] == 'PASS' and D._yes(gate['checks'])
              and all(gate[k] == v for k, v in current.items()), 'same independent signed-D prelaunch approval required')
    pf = D._bound(rel(PREP / 'preflight.json'), current['preflight_sha256'], {})
    D.require(pf['status'] == 'PASS' and D._yes(pf['checks']) and pf['actual_cached_sites_assembled'] == 0
              and pf['protocol_sha256'] == current['protocol_sha256'] and pf['source_sha256'] == current['source_sha256']
              and sha(PREP / 'exact_signed_cohort_predictor_frozen.py') == current['source_sha256'],
              'current frozen preflight drifted')
    bundle = input_bundle()
    proto = D._bound(rel(PREP / 'protocol.json'), current['protocol_sha256'], {})
    D.require(proto == json.loads(json.dumps(prep_doc(bundle), default=str)), 'bound D signed proof inputs changed')
    links, context, (bulk, surface, gain, bulk_s, offsets), args, old_z, radius, b, z, norm, provider = bundle
    import torch
    torch.set_num_threads(1)
    ACTUAL.mkdir(exist_ok=False)
    start = monotonic()
    link = {**current, 'prelaunch_review_sha256': review_sha}
    D._write(ACTUAL / 'started.json', {**link, 'pid': os.getpid(), 'record_id': RECORD, 'utc': now()})
    try:
        g, cache_m, remainder, min_j, center_z = signed._point_and_remainder(context, args, offsets)
        origins = [signed._add(x, y, s) for x, y, s in zip(g, bulk, surface, strict=True)]
        slopes = signed.combine_slopes(cache_m, bulk_s, gain, tuple(P.NAMES).index('common_gain'))
        intervals = signed.affine_rows(origins, slopes, offsets, remainder)
        classification = classify(intervals, z, radius)
        why = (D._minimum(context, links, old_z, z, radius, norm, provider)
               if classification['sufficient_existence_test_passes'] else None)
        result = {**link, **classification, 'record_id': RECORD, 'scope': SCOPE, 'input_hashes': links,
                  'same_closed_C_sha256': context.row['parameter_cell_sha256'],
                  'closed_parameter_offsets_exact': offsets, 'r_exact': radius, 'b_upper_exact': b,
                  'old_cached_Z_upper_exact': old_z, 'z_upper_exact': z, 'exact_R_infinity_norm': norm,
                  'cached_predictor_path_min_J_exact': min_j[0], 'cached_center_defect_Z0_exact': center_z[1],
                  'cached_signed_anchor_gR': g, 'cached_signed_affine_M': cache_m,
                  'nonnegative_cached_predictor_remainder_upper_exact': remainder,
                  'signed_combined_origins': origins, 'signed_combined_five_slopes': slopes,
                  'signed_full_integral_RF_predictor_intervals': intervals,
                  'conditional_curvature': why, 'actual_independent_scope_review_required': True,
                  'global_model_absence_claimed': False, 'independently_certified_root_or_exclusion': False,
                  'seconds': monotonic() - start}
        h = D._write(ACTUAL / 'results.json', result)
        D._write(ACTUAL / 'terminal.json', {**link, 'status': 'FINISHED', 'record_id': RECORD,
                  'results_sha256': h, 'pid': os.getpid(), 'exit_code': 0, 'utc': now()})
        print(json.dumps({'status': result['status'], 'results_sha256': h}), flush=True)
    except BaseException as error:
        D._write(ACTUAL / 'failure.json', {**link, 'error': str(error), 'traceback': traceback.format_exc()})
        D._write(ACTUAL / 'terminal.json', {**link, 'status': 'FAILED', 'record_id': RECORD,
                  'pid': os.getpid(), 'exit_code': 1, 'utc': now()})
        raise


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    stages = cli.add_subparsers(dest='stage', required=True)
    stages.add_parser('prepare')
    actual = stages.add_parser('run')
    actual.add_argument('--review-sha256', required=True)
    options = cli.parse_args()
    if options.stage == 'prepare':
        prepare()
    else:
        run(options.review_sha256)

