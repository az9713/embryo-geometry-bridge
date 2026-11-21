"""Later-wave adapter for the two reviewed signed global reduction kernels."""

import argparse
from collections import Counter
from fractions import Fraction as Q
import importlib.util
from itertools import combinations_with_replacement
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace


SOURCE = Path(__file__).resolve()
HERE = SOURCE.parent
LATER = HERE.parent
BRIDGE = LATER.parent
ROOT = next(parent for parent in HERE.parents if (parent / 'PLAN.md').is_file())
sys.path.insert(0, str(LATER))
import wave_execution as X
import hull_replay as J
import transfer_ambiguity_remaining32_first_wave_reduce as D


S, W, V, K = D.S, D.W, D.V, D.K
read, save, sha, rel, require = D.read, D.save, D.sha, D.rel, D.require
OLD_REPLAY = BRIDGE / 'first_wave_reduction_v1/independent_actual_preparation/replay_frozen.py'
AUTH = HERE / 'independent_source_review/receipt.json'
N = 611
RADIUS = Q(1, 2**34)
LEAF_COUNT = 168
RELEASE = '002b471dceed8666aa0f70e1e9673744f49d0bf3a70da65a194c3faf34b0adbe'
PINS = {
    J.SOURCE: 'd5cb7ccaba8526e9ca3f9d413ba0f6ea2ae24328f43cc41825af5314044adbfc',
    J.PREP / 'protocol.json': 'd7c4d611814350dff83aba7559ebfaaaff852bb232c2a0199b903a4740d8d2b3',
    J.PREP / 'preflight.json': 'f51a9789357238b6e39f6fe1806ffcafd5a66f693f6f3bc9919bba5516757dce',
    J.REVIEW: '1c6d2060393903f0fbebd14e17cf7d7cf588d473ab2368df19a48e450dd8b20c',
    J.PREP / 'hull_replay_frozen.py': 'd5cb7ccaba8526e9ca3f9d413ba0f6ea2ae24328f43cc41825af5314044adbfc',
    D.SOURCE: '58de9f5e96b43005b0ba83dd0f4d8cd1b35f0690eef58290a75ef5abce692d92',
    Path(K.__file__).resolve(): 'd31ca58afe2e1b935d723513fb9fa18f2f978885f8f6ec32e48766f8b15bf0e8',
    OLD_REPLAY: '4682165bf59bd633092c042a7df619c9bc1ad272dc80d1755404c2999a65287f',
    OLD_REPLAY.with_name('root_source_review.json'):
        '12e14fb46070d40694c09f3b2a3d76492b1dc85659c7a61c17df6b617f6915e0',
    D.SPD / 'actual_replay/os_exit_observation.json':
        'fd8ac1cb3b920b57096649544bba105442404fe19576aae99c53f849fe15a8d5',
}
SPD = {'observed_os_exit_sha256': PINS[D.SPD / 'actual_replay/os_exit_observation.json'],
       'replay_phase': {'source_sha256': '28c511fc6087d130f353f67c3cb655aa729031b557ae0bd84c9714b632435a13',
                        'protocol_sha256': 'a740f6a3305421de51a8fa7ce1e5da94bf308ae35ab7226c0ef3336e1beed8c7',
                        'preflight_sha256': 'c8150950f934719f6212b579bcac066573c94b09ac5943cd7f4c9b5d489339d7'}}
ZERO = {'actual_hull_leaves_opened': 0, 'actual_global_reductions': 0, 'new_force_AD_or_spatial_integrals': 0}
SOURCE_PASS = 'PASS_INDEPENDENT_LATER_WAVE_GLOBAL_REDUCTION_SOURCE'
HULL_OS = 'PASS_ROOT_OBSERVED_LATER_WAVE_HULL_REPLAY_EXITS'
PROD_OS = 'PASS_ROOT_OBSERVED_LATER_WAVE_GLOBAL_PRODUCER_EXIT'
PROD_STATUS = 'COMPLETE_LATER_WAVE_GLOBAL_PRODUCER_AWAITS_INDEPENDENT_ENDPOINT_REPLAY'
INDEP_STATUS = 'COMPLETE_INDEPENDENT_LATER_WAVE_SIGNED_ENDPOINT_REPLAY'


def endpoint_kernel():
    require(sha(OLD_REPLAY) == PINS[OLD_REPLAY], 'reviewed endpoint implementation')
    spec = importlib.util.spec_from_file_location('_prior_signed_endpoint_functions', OLD_REPLAY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.D is D and module.S is S and module.BITS == 768 and module.N == N,
            'shared pinned endpoint and physical authority aliases')
    return module


def static():
    require(all(sha(p) == h for p, h in PINS.items()), 'frozen local proof or reviewed arithmetic drift')
    require(J.X is X and X.S is J.S is S and J.W is W and S.RADIUS == RADIUS, 'same reviewed stages')
    profiles, release = J.static()
    require(release['later_waves_release_sha256'] == RELEASE and len(profiles) == 7, 'matching later-stage release')
    rows, _, _ = D.optional_spd(SPD)
    require(tuple(rows) == V.A.IDS and all(row['status'] == D.SPD_PASS for row in rows.values()),
            'already observed all32 SPD proof with unchanged data')
    p = {'status': 'FROZEN_LATER_WAVE_GLOBAL_REDUCTION_ADAPTER', 'source_sha256': sha(SOURCE),
         'source_pins': {rel(p): h for p, h in {**J.PINS, **PINS}.items()},
         'later_waves_release_sha256': RELEASE, 'profiles_sha256': X.PROFILE_SHA,
         'wave_case_order': {wave: profiles[wave]['case_order'] for wave in X.WAVES},
         'radius_exact': str(RADIUS), 'same_R_spd': SPD, 'producer_precision_bits': 640,
         'independent_endpoint_precision_bits': 768, **ZERO}
    return p, profiles


def approved(p, f, review, expected):
    require(review['status'] == SOURCE_PASS and V.yes(review['checks'])
            and all(review[k] == v for k, v in expected.items())
            and review['permitted_waves'] == list(X.WAVES)
            and review['approved_modes'] == ['producer', 'independent']
            and review['new_spatial_integrations_authorized'] is False
            and p['source_sha256'] == f['source_sha256'] == expected['source_sha256']
            and f['protocol_sha256'] == expected['protocol_sha256'] and f['status'] == 'PASS'
            and V.yes(f['checks']) and all(V.zero(x[k]) for x in (p, f) for k in ZERO),
            'separate frozen later signed reduction source approval')


def known_exit(row):
    handle = row.get('unified_exec_session')
    chunk = row.get('unified_exec_chunk_id')
    return (V.positive_int(row['pid']) and V.zero(row['observed_exit_code'])
            and (handle is None or V.positive_int(handle)) and
            (chunk is None or type(chunk) is str and bool(chunk.strip()))
            and (handle is not None or chunk is not None))


def hull_observation(doc, wave, cases, pins):
    require(doc['status'] == HULL_OS and doc['wave_id'] == wave
            and all(doc[k] == v for k, v in pins.items())
            and tuple(r['case_id'] for r in doc['cases']) == tuple(cases)
            and all(known_exit(r) for r in doc['cases'])
            and len({r['pid'] for r in doc['cases']}) == 4,
            'complete same-wave externally observed native 0 closure')
    return doc['cases']


def check_native(start, terminal, receipt, observed, phase, target, plan, raw):
    require(all(V.positive_int(r['pid']) and all(r[k] == v for k, v in phase.items())
                for r in (start, terminal, receipt))
            and start['status'] == 'RUNNING' and terminal['status'] == 'FINISHED'
            and V.zero(terminal['exit_code']) and terminal['receipt_sha256'] == observed['receipt_sha256']
            and receipt['status'] == J.RAW_PASS and receipt['target'] == target
            and receipt['record_id'] == target['record_id'] and Q(receipt['radius_exact']) == RADIUS
            and receipt['full_independent_hulls_computed'] is True
            and receipt['root_asserted'] is False and receipt['overlap_or_producer_containment_used_as_proof'] is False
            and receipt['admit_hull_only_global_reduction_after_separate_review'] is True
            and set(receipt['disjoint_counts']) == set(D.CHANNELS)
            and all(V.zero(v) for v in receipt['disjoint_counts'].values())
            and receipt['selector_plan_sha256'] == plan['plan_sha256']
            and receipt['plan_review_sha256'] == plan['independent_plan_review_sha256']
            and V.positive_int(receipt['producer_observed_pid']) and receipt['producer_observed_pid'] == raw['pid']
            and V.positive_int(receipt['producer_observed_unified_exec_session'])
            and receipt['producer_observed_unified_exec_session'] == raw['unified_exec_session']
            and receipt['raw_root_observed_mapping_sha256'] == S.E._digest(raw['sha256'])
            and type(receipt['raw_root_observed_component_files']) is int
            and receipt['raw_root_observed_component_files'] == len(raw['sha256']) == 715,
            'complete native hull outcome and actual raw/selector identity required')
    require(tuple((m['kind'], m['slot']) for m in receipt['case_components']) == S.OPS
            and all(type(m['slot']) is int and m['status'] == D.LEAF_PASS for m in receipt['case_components']),
            'all 168 local native proof rows, with no HOLD')


def open_wave(args, static_profile):
    pins = {'raw_wave_admission_sha256': args.raw_wave_admission_sha256,
            'raw_root_observed_exit_sha256': args.raw_root_observed_exit_sha256}
    wave = args.wave
    root = X.wave_folder(wave)
    observed = read(root / 'admission/global_hull_root_admission.json', args.hull_admission_sha256)
    cases = static_profile['case_order']
    observed_rows = hull_observation(observed, wave, cases, pins)
    xa = J.context_args(wave, args.raw_wave_admission_sha256)
    profile, admission, parts, raw = X.open_raw_phase(xa, args.raw_root_observed_exit_sha256)
    require(profile == static_profile, 'same frozen phase and full native raw observation')
    native = {}
    for obs in observed_rows:
        case = obs['case_id']
        ja = SimpleNamespace(wave=wave, case_id=case, protocol_sha256=PINS[J.PREP / 'protocol.json'],
                             preflight_sha256=PINS[J.PREP / 'preflight.json'], review=J.REVIEW,
                             review_sha256=PINS[J.REVIEW])
        hprofile, binding = J.gate(ja)
        require(hprofile == profile, 'same native hull source and later-wave scope')
        folder = root / 'independent_hulls' / case
        target = profile['raw_worker_targets'][case]
        plan, = [row for row in admission['plans'] if row['case_id'] == case]
        start, end, receipt = (read(folder / (name + '.json'), obs[name + '_sha256'])
                               for name in ('started', 'terminal', 'receipt'))
        phase = {**binding, **pins, 'pid': obs['pid'], 'case_id': case, 'wave_id': wave}
        check_native(start, end, receipt, obs, phase, target, plan, raw)
        native[case] = folder, receipt, phase, obs['receipt_sha256']
    # Verify the full sealed 672-file hull catalogue before parsing its signed operands.
    for case, (folder, receipt, _, _) in native.items():
        for meta in receipt['case_components']:
            kind, slot = meta['kind'], meta['slot']
            raw_dir, part = parts[case, kind]
            item = part['components'][slot]
            leaf_path, raw_path = folder / item['file'], raw_dir / item['file']
            require(meta['native_file'] == rel(leaf_path) and meta['packet_file'] == rel(raw_path)
                    and meta['packet_sha256'] == item['sha256'] == raw['sha256'][rel(raw_path)]
                    and sha(leaf_path) == meta['native_sha256'], 'exact 672 observed leaf digest closure')
    return admission, parts, native


def channel_operands(mode, I, meta, leaf, packet, part, row, proof, phase, tails):
    kind, slot = meta['kind'], meta['slot']
    target = row['target']
    ids = proof['global_indices']
    J.exact_worker_fields(phase['wave_id'], phase['case_id'], kind, packet, part, row, proof)
    W.validate_packet(packet, proof, target, slot, kind, tails)
    require((leaf['kind'], leaf['slot'], leaf['global_indices']) == (kind, slot, ids)
            and type(leaf['slot']) is int and ids == meta['global_indices']
            and len(set(ids)) == len(ids) and all(type(i) is int and 0 <= i < N for i in ids)
            and leaf['status'] == D.LEAF_PASS and leaf['record_id'] == target['record_id']
            and leaf['phase'] == phase and leaf['packet_sha256'] == meta['packet_sha256']
            and leaf['selector_leaf_content_sha256'] == W.B.digest(proof)
            and leaf['only_independently_valid_enclosures_used_for_hull_truth'] is True
            and leaf['overlap_not_used_as_proof'] is True and leaf['root_asserted'] is False
            and leaf['checked_cache_sites'] == (64 if kind == 'bulk' else 12)
            and (leaf['point_force_scope'], leaf['whole_H_scope']) == (S.POINT, S.BOX)
            and Q(leaf['radius_exact']) == RADIUS and set(leaf['signed_channels']) == set(D.CHANNELS),
            'same independently valid typed union proof, no unresolved signed scope')
    F, H = packet['point_force'], packet['whole_box_hessian']
    local = list(range(len(ids)))
    upper = list(combinations_with_replacement(local, 2))
    art = leaf['independent_point_and_whole_tails_exact']
    correction = leaf['independent_same_Y_whole_H_correction_exact']
    require(len(art['F']) == len(ids) and all(0 <= Q(v) <= tails[0] for v in art['F'])
            and leaf['producer_force_tails_exact'] == F['force_tail_upper_exact']
            and all(type(v) is list and len(v) == 3 and type(v[0]) is type(v[1]) is int for v in art['H'])
            and all(type(v) is list and len(v) == 4 and type(v[0]) is type(v[1]) is int for v in correction)
            and [tuple(v[:2]) for v in art['H']] == [tuple(v[:2]) for v in correction] == upper
            and all(0 <= Q(v[2]) <= tails[1] for v in art['H']), 'local separate F/H tails and finite correction')
    for value in correction:
        D.checked_interval(value[2:])
    raw = (F['F_I_minus_F_cache_intervals_exact'], [v[2:] for v in H['H_I_minus_H_cache_upper_exact']],
           F['Qstar_minus_cache_correction_exact'])
    values, misses = [], []
    for name, indices, original in zip(D.CHANNELS, (local, list(map(list, upper)), local), raw, strict=True):
        if mode == 'producer':
            endpoints, count = D.hull_channel(leaf['signed_channels'][name], indices, original)
        else:
            endpoints, count = I.validated_channel(leaf, name, indices, original)
        values.append(endpoints)
        misses.append(count)
    return {'kind': kind, 'slot': slot, 'global_indices': ids}, upper, values[:2], misses


def assembled(mode, I, case, raw_admission, parts, admitted_hulls):
    folder, receipt, phase, rh = admitted_hulls[case]
    row, d, tails, proofs = X.exact_context(raw_admission, case)
    target = row['target']
    plan = read(S.PLAN / case / 'plan.json', row['plan_sha256'])
    ctx = d['ctx']
    R = [[V.A.T.rational(v) for v in rr] for rr in ctx.R]
    norm = max(sum(abs(v) for v in rr) for rr in R)
    zq = Q(plan['sizing']['ZQ_exact'])
    require(len(R) == N and all(len(rr) == N for rr in R)
            and norm == Q(plan['sizing']['R_inf_exact']) > 0 and 0 <= zq < 1
            and V.A.a_sha(ctx.R) == target['R_binary_sha256'] and V.A.a_sha(ctx.u) == target['center_binary_sha256']
            and ctx.row['point_cache_sha256'] == target['point_cache_sha256']
            and ctx.drive == Q(target['exact_drive'])
            and Q(plan['whole_box_radius_exact']) == RADIUS and plan['trial_box_sha256'] == target['trial_box_sha256']
            and tuple(map(Q, plan['sizing']['caps_exact'])) == tails
            and len(plan['signed_cached_RF_center_rows_exact']) == N, 'same exact R, scope and signed cached anchor')
    anchor = [D.checked_interval(v) for v in plan['signed_cached_RF_center_rows_exact']]
    routes, packets, pins, counts = [], [], {}, Counter()
    force, hullH = [(Q(0), Q(0)) for _ in range(N)], {}
    for meta in receipt['case_components']:
        key = meta['kind'], meta['slot']
        raw_folder, part = parts[case, key[0]]
        filename = part['components'][key[1]]['file']
        path = folder / filename
        leaf = read(path, meta['native_sha256'])
        raw = read(raw_folder / filename, meta['packet_sha256'])
        route, upper, (f, h), misses = channel_operands(mode, I, meta, leaf, raw, part, row, proofs[key], phase, tails)
        routes.append(route)
        pins[rel(path)] = meta['native_sha256']
        for name, n in zip(D.CHANNELS, misses, strict=True):
            counts[name] += n
        if mode == 'producer':
            packets.append({**route, 'point_force': {'scope': S.POINT, 'F_I_minus_F_cache_intervals_exact': f},
                            'whole_box_hessian': {'scope': S.BOX, 'radius_exact': str(RADIUS),
                            'H_I_minus_H_cache_upper_exact': [[i, j, *v] for (i, j), v in zip(upper, h, strict=True)]}})
        else:
            I.scatter_add(force, hullH, route['global_indices'], f, h)
    require(dict(counts) == receipt['producer_containment_miss_counts']
            and len(routes) == len(pins) == LEAF_COUNT
            and routes == [{k: entry[k] for k in ('kind', 'slot', 'global_indices')}
                           for entry in d['entry']['component_catalogue']], 'complete proved signed free-map catalogue')
    if mode == 'producer':
        force, hullH = K.assemble(packets, routes, N)
    fields = {'target': target, 'radius_exact': str(RADIUS), 'cache_anchor_already_preconditioned': True,
              'R_infinity_norm_exact': str(norm), 'whole_cube_ZQ_exact': str(zq), 'maps': routes,
              'signed_cached_RF_exact': [list(map(str, v)) for v in anchor],
              'signed_global_point_force_defect_exact': [list(map(str, v)) for v in force],
              'signed_global_whole_cube_H_defect_upper_exact': [[i, j, *map(str, v)] for (i, j), v in sorted(hullH.items())],
              'same_selector_and_cached_geometry_bindings': plan['bindings'],
              'independent_local_hull_receipt_sha256': rh, 'independent_local_hull_leaf_sha256': pins,
              'producer_containment_miss_counts': dict(counts)}
    return row, d, plan, R, norm, zq, anchor, force, hullH, fields


def producer_observed(args, output, cases):
    root = output / 'producer'
    obs = read(root / 'root_observed_exit.json', args.producer_os_sha256)
    require(obs['status'] == PROD_OS and known_exit(obs) and obs['wave_id'] == args.wave
            and obs['hull_admission_sha256'] == args.hull_admission_sha256,
            'observed exact matching successful producer process required')
    start, end, result = (read(root / (k + '.json'), obs[k + '_sha256'])
                          for k in ('started', 'terminal', 'results'))
    stamp = {'mode': 'producer', 'wave_id': args.wave, 'pid': obs['pid'],
             'hull_admission_sha256': args.hull_admission_sha256,
             'raw_wave_admission_sha256': args.raw_wave_admission_sha256,
             'raw_root_observed_exit_sha256': args.raw_root_observed_exit_sha256,
             'source_sha256': sha(SOURCE), 'protocol_sha256': args.protocol_sha256,
             'preflight_sha256': args.preflight_sha256, 'source_review_sha256': args.review_sha256}
    require(all(V.positive_int(d['pid']) and all(d[k] == v for k, v in stamp.items()) for d in (start, end, result))
            and start['status'] == 'RUNNING' and end['status'] == 'FINISHED'
            and V.zero(end['exit_code']) and end['results_sha256'] == obs['results_sha256']
            and result['status'] == PROD_STATUS and tuple(r['case_id'] for r in result['cases']) == tuple(cases),
            'same complete producer wave, phase, true OS and all assembled outputs')
    for item in result['cases']:
        folder = root / item['case_id']
        require(all(sha(folder / name) == item[k] for name, k in (('operands.json', 'operands_sha256'),
                        ('results.json', 'case_results_sha256'))), 'seal all eight case outputs')
    return {r['case_id']: r for r in result['cases']}


def apply_case(mode, I, case, raw, parts, native, spd, dest, stamp, old):
    row, d, plan, R, norm, zq, anchor, f, H, fields = assembled(mode, I, case, raw, parts, native)
    target = row['target']
    if mode == 'producer':
        premise = D.min_premise(spd, {'opts': {**row, 'case_id': case}, 'd': d, 'plan': plan}, norm)
        K.ctx.prec = 640
        arithmetic = K.contract(R, f, H, anchor, zq, radius=RADIUS)
        transfer = D.transfer(arithmetic, target, premise, norm)
        status = 'LATER_WAVE_GLOBAL_PRODUCER_CANDIDATE_AWAITS_INDEPENDENT_REPLAY'
        comparison = None
    else:
        pdir = dest.parent / 'producer' / case
        metadata = old[case]
        require(read(pdir / 'operands.json', metadata['operands_sha256']) == fields,
                'independent exact scatter must reproduce all assembled signed operands')
        prior = read(pdir / 'results.json', metadata['case_results_sha256'])
        require(prior['status'] == 'LATER_WAVE_GLOBAL_PRODUCER_CANDIDATE_AWAITS_INDEPENDENT_REPLAY'
                and prior['root_asserted'] is False and prior['target'] == target
                and prior['operands_sha256'] == metadata['operands_sha256'], 'same observed producer arithmetic scope')
        premise = I.own_spd(spd, target, d['ctx'], plan, norm)
        arithmetic = I.classify(f, H, anchor, R, zq, RADIUS)
        transfer = I.observer_minimum(arithmetic, target, premise is not None, norm)
        status = ('PASS_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_ROOT' if arithmetic['strict_selfmap'] else
                  'SUFFICIENT_INDEPENDENT_EXCLUSION_OF_ONLY_THIS_CUBE' if arithmetic['strict_cube_exclusion'] else
                  'UNRESOLVED_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_TEST')
        comparison = I.compare_producer(arithmetic, transfer, prior)
    require(premise is not None, 'same exact R with independent positive symmetric part required')
    folder = dest / case
    folder.mkdir(exist_ok=False)
    oh = save(folder / 'operands.json', fields)
    root = mode == 'independent' and arithmetic['strict_selfmap']
    result = {'status': status, **stamp, 'case_id': case, 'target': target, 'radius_exact': str(RADIUS),
              'operands_sha256': oh, 'arithmetic': arithmetic, 'same_R_SPD_premise': premise,
              'observer_and_optional_minimum': transfer, 'producer_comparison': comparison,
              'root_asserted': root, 'archived_observer_field_replaced': False, 'full_five_gate_complete': False}
    rh = save(folder / 'results.json', result)
    return {'case_id': case, 'status': status, 'operands_sha256': oh, 'case_results_sha256': rh,
            'exact_spatial_root_certified': root,
            'exact_spatial_minimum_certified': mode == 'independent' and transfer['strict_minimum_proved'],
            'shift_inclusive_observer_budget_proved': mode == 'independent' and transfer['archived_state_budget_proved']}


def run(args):
    p, profiles = static()
    want = {'source_sha256': sha(SOURCE), 'protocol_sha256': args.protocol_sha256,
            'preflight_sha256': args.preflight_sha256, 'later_waves_release_sha256': RELEASE,
            'profiles_sha256': X.PROFILE_SHA}
    approved(read(HERE / 'protocol.json', args.protocol_sha256), read(HERE / 'preflight.json', args.preflight_sha256),
             read(AUTH, args.review_sha256), want)
    require(read(HERE / 'protocol.json', args.protocol_sha256) == p, 'runtime cold plan drift')
    output = X.wave_folder(args.wave) / 'global_reduction_v1'
    profile = profiles[args.wave]
    cases = profile['case_order']
    old = producer_observed(args, output, cases) if args.mode == 'independent' else None
    raw, parts, native = open_wave(args, profile)
    I, spd = endpoint_kernel(), D.optional_spd(SPD)
    dest = output / args.mode
    dest.mkdir(parents=True, exist_ok=False)
    stamp = {**want, 'pid': os.getpid(), 'mode': args.mode, 'wave_id': args.wave,
             'source_review_sha256': args.review_sha256, 'hull_admission_sha256': args.hull_admission_sha256,
             'raw_wave_admission_sha256': args.raw_wave_admission_sha256,
             'raw_root_observed_exit_sha256': args.raw_root_observed_exit_sha256}
    if args.mode == 'independent':
        stamp['producer_observed_os_sha256'] = args.producer_os_sha256
    save(dest / 'started.json', {'status': 'RUNNING', **stamp})
    end = {'status': 'FAILED', 'exit_code': 1, **stamp}
    completed = []
    try:
        for case in cases:
            completed.append(apply_case(args.mode, I, case, raw, parts, native, spd, dest, stamp, old))
            print(json.dumps(completed[-1]), flush=True)
        rh = save(dest / 'results.json', {'status': PROD_STATUS if args.mode == 'producer' else INDEP_STATUS,
                   **stamp, 'cases': completed, 'signed_hull_components': 672, 'new_force_AD_or_spatial_integrals': 0,
                   'full_five_gate_complete': False})
        end.update(status='FINISHED', exit_code=0, results_sha256=rh)
        print(json.dumps({'results_sha256': rh}), flush=True)
    except Exception as error:
        save(dest / 'failure.json', {'error': str(error), 'completed': completed, **stamp})
        raise
    finally:
        save(dest / 'terminal.json', end)


def fixtures(protocol):
    I = endpoint_kernel()
    scopes = {'point_force': {'scope': S.POINT, 'F_I_minus_F_cache_intervals_exact': [['3/64', '3/64'], ['1/32', '1/32']]},
              'whole_box_hessian': {'scope': S.BOX, 'radius_exact': str(RADIUS),
                  'H_I_minus_H_cache_upper_exact': [[0, 0, '1/16', '1/16'], [0, 1, '-1/32', '-1/32'],
                                                  [1, 1, '1/64', '1/64']]}}
    proj = {'kind': 'bulk', 'slot': 0, 'global_indices': [1, 0]}
    refF, refH = K.assemble([{**proj, **scopes}], [proj], 2)
    f, h = [(Q(0), Q(0)) for _ in range(2)], {}
    I.scatter_add(f, h, [1, 0], [(Q(3, 64),)*2, (Q(1, 32),)*2],
                 [(Q(1, 16),)*2, (Q(-1, 32),)*2, (Q(1, 64),)*2])
    require((f, h) == (refF, refH), 'exact signed symmetric scatter including reversed map')
    R = [[Q(1), Q(-1)], [Q(0), Q(1)]]
    anchor = [(Q(1, 64),)*2, (Q(-3, 64),)*2]
    K.ctx.prec = 640
    one = K.contract(R, f, h, anchor, Q(1, 4), radius=RADIUS)
    two = I.classify(f, h, anchor, R, Q(1, 4), RADIUS)
    require(one['signed_RF_I_rows_exact'] == two['signed_R_F_I_rows_exact'] == [['0', '0'], ['0', '0']]
            and Q(one['b_upper_exact']) == Q(two['b_H_upper_exact']) == Q(9, 64),
            'independent signed R and mirror exactly once, then cached RF anchor exactly once')
    target = {'archived_center_shift_exact': '1/8', 'original_addition_gap_exact': '1/16', 'observer_goal_exact': '1/4'}
    require(I.observer_minimum(two, target, True, Q(2))['physical_nodal_error_um_upper_exact'] == str(Q(201, 32)),
            'original archive receives sigma, posterior and xi')
    case_ids = next(iter(protocol['wave_case_order'].values()))
    record = {'status': HULL_OS, 'wave_id': 'wave_02',
              'cases': [{'case_id': c, 'pid': 100 + i, 'observed_exit_code': 0,
                         'unified_exec_session': 200 + i} for i, c in enumerate(case_ids)]}
    require(len(hull_observation(record, 'wave_02', case_ids, {})) == 4
            and known_exit({'pid': 101, 'unified_exec_chunk_id': '0d1a50', 'observed_exit_code': 0})
            and V.A.reject(lambda: hull_observation({**record, 'cases': record['cases'][:-1]}, 'wave_02', case_ids, {}))
            and V.A.reject(lambda: hull_observation({**record, 'wave_id': 'wave_03'}, 'wave_02', case_ids, {}))
            and not known_exit({'pid': 101.0, 'unified_exec_session': 4, 'observed_exit_code': 0})
            and not known_exit({'pid': 101, 'unified_exec_session': 4, 'observed_exit_code': False}),
            'typed observed sync/async process and complete matching four targets')
    return {'frozen_all32_same_R_SPD_with_all_seven_profiles': True, 'typed_complete_wave_observed_admission': True,
            'tiny_signed_oracle_reversed_projection_symmetric_H': True, 'preconditioned_cached_anchor_added_once': True,
            'shift_posterior_addition_gap_and_strict_scopes': True}


def prepare():
    p, _ = static()
    checks = fixtures(p)
    h = save(HERE / 'protocol.json', p)
    f = save(HERE / 'preflight.json', {'status': 'PASS', 'source_sha256': sha(SOURCE), 'protocol_sha256': h,
                                       'checks': checks, **ZERO})
    print(json.dumps({'source_sha256': sha(SOURCE), 'protocol_sha256': h, 'preflight_sha256': f}))


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    commands = cli.add_subparsers(dest='mode', required=True)
    commands.add_parser('prepare')
    for name in ('producer', 'independent'):
        stage = commands.add_parser(name)
        stage.add_argument('--wave', required=True, choices=X.WAVES)
        for key in ('protocol', 'preflight', 'review', 'hull-admission', 'raw-wave-admission', 'raw-root-observed-exit'):
            stage.add_argument('--' + key + '-sha256', required=True)
        if name == 'independent':
            stage.add_argument('--producer-os-sha256', required=True)
    args = cli.parse_args()
    prepare() if args.mode == 'prepare' else run(args)
