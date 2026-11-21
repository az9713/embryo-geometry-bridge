"""Review frozen later-wave selector admission and arithmetic without real plans."""

from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
PREP = HERE.parent
sys.path.insert(0, str(PREP))
SOURCE = 'b8cf7a1ec27f4a342473b40bdd8ec894f1183ab26c4a4bf2976115bbf52b6c63'
PROTOCOL = 'bf9f215e8adbfbd4cb01444177d38e9ff6e877e77bbcf0f4268782b743adba04'
PREFLIGHT = 'e0c01b71bb7355e50b940d84d3babe32ef943ee227eee0fa3fb82096b6aba433'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def reject(call, phrase=None):
    try:
        call()
    except (ValueError, KeyError) as error:
        require(phrase is None or phrase in str(error), f'wrong rejection: {error}')
        return True
    raise ValueError('expected rejected fixture')


def check_input_arithmetic(R, M, catalogue, joined, profiles, protocol):
    waves = [f'wave_{i:02d}' for i in range(2, 9)]
    require(list(profiles) == waves, 'all seven later waves')
    cases = [case for profile in profiles.values() for case in profile['case_order']]
    require(len(cases) == len(set(cases)) == 28 and set(cases).isdisjoint(M.FIRST)
            and set(cases) | set(M.FIRST) == set(catalogue) == set(joined), 'exact partition of original 32')
    require(protocol['ordered_case_ids_by_wave'] == {w: profiles[w]['case_order'] for w in waves},
            'protocol profile case order')
    modes = Counter()
    total_routes = 0
    for w in waves:
        is_wave = set(profiles[w]['case_order'])
        require(len(is_wave) == 4, 'four exact cases per wave')
        for case in is_wave:
            context = catalogue[case]
            cached = joined[case]
            trial = context['trial']
            require(R.route(profiles, w, case) == R.E.W.case_payload(context)
                    and cached['trial_box_sha256'] == trial['trial_box_sha256']
                    and Q(trial['trial_radius_exact']) == M.R == Q(1, 2**34)
                    and context['cached_input_kind'] == cached['cached_input_kind'],
                    'same trial box, record and cached authority')
            modes[cached['cached_input_kind']] += 1
            routes = context['component_catalogue']
            require(len(routes) == 168 and all(type(r['slot']) is int for r in routes), 'complete typed routes')
            require((context['point_force_scope'], context['whole_H_scope']) == (M.POINT, M.BOX),
                    'point and whole cube roles')
            f = [0] * 611
            h = [0] * 611
            for route, slot in zip(routes, M.ORDER, strict=True):
                ids = route['global_indices']
                require((route['kind'], route['slot']) == slot
                        and len(ids) == len(set(ids)) and all(type(i) is int and 0 <= i < 611 for i in ids),
                        'validated exact-degree-of-freedom incidence')
                for i in ids:
                    f[i] += 1
                    h[i] += len(ids)
            require(min(f) > 0 and min(h) > 0 and max(f) == 6 and max(h) == 72, 'whole projected coverage')
            z, norm = map(Q, (cached['reviewed_ZQ_upper_exact'], cached['R_inf_exact']))
            force, Hessian = (1-z)*M.R/4, (1-z)/2
            explicit = [force/(10000*norm*max(f)), Hessian/(10000*norm*max(h))]
            replay = M.independently_size(routes, z, norm)
            require(0 <= z < 1 and norm > 0 and min(explicit) > 0
                    and list(replay['local']) == list(map(Q, cached['localized_force_and_H_tail_caps_exact'])) == explicit
                    and force + (z + Hessian) * M.R == M.R - force < M.R, 'tail-only rational budget identities')
            total_routes += len(routes)
    require(modes == {'NATIVE_R34': 22, 'FRESH_R34': 6} and total_routes == 4704, 'complete later context coverage')
    for w in waves:
        require(reject(lambda: R.route(profiles, w, M.FIRST[0])), 'first-wave target admitted')
        require(reject(lambda: R.route(profiles, w, next(c for c in cases if c not in profiles[w]['case_order']))),
                'other-wave target admitted')
    return modes, total_routes


def check_admission(R, M, case):
    folder = R.E.S.PLAN / case
    binding = {'fixture': 'same-run'}
    stamp = {'case_id': case, 'pid': 173, 'bindings': binding}
    ordered = [{'kind': k, 'slot': i, 'file': f'{k}_{i:03d}.json', 'sha256': 'a'*64}
               for k, i in M.ORDER]
    plan = {'status': 'COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW',
            **stamp, 'spatial_workers_authorized': False, 'new_spatial_integrations': 0,
            'components': ordered, 'point_force_scope': M.POINT, 'whole_H_scope': M.BOX,
            'whole_box_radius_exact': str(M.R)}
    start = {'status': 'RUNNING_SELECTOR_PLAN', **stamp}
    terminal = {'status': 'FINISHED', **stamp, 'exit_code': 0, 'plan_sha256': 'b'*64}
    obs = {'status': 'PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT', 'case_id': case,
           'pid': 173, 'unified_exec_session': 29, 'observed_exit_code': 0,
           'plan_sha256': 'b'*64, 'terminal_sha256': 'c'*64, 'started_sha256': 'd'*64}

    def call(doc=obs, p=plan, s=start, t=terminal, corrupt=None):
        order = []
        expected = {'root_os_exit_observation.json': ('e'*64, doc), 'started.json': ('d'*64, s),
                    'terminal.json': ('c'*64, t), 'plan.json': ('b'*64, p)}
        fixture_leaves = {item['file'] for item in ordered}

        def load(path, h, links, parse=True):
            require(path.parent == folder, 'same case fixture folder')
            name = path.name
            order.append(name)
            if name in expected:
                pinned, result = expected[name]
                require(parse and h == pinned, 'fixture process pin')
                return result
            require(not parse and name in fixture_leaves and h == 'a'*64 and name != corrupt,
                    'fixture leaf integrity')
            return b'fixture'

        with patch.object(M, 'load', side_effect=load):
            r = R.admitted(M, case, 'e'*64, {})
        require(order[:4] == ['root_os_exit_observation.json', 'terminal.json', 'started.json', 'plan.json']
                and len(order) == 172 and r == (folder, doc, p), 'external chain precedes all 168 leaves')
        return True

    wrong_first = deepcopy(plan)
    wrong_first['components'][0]['slot'] = False
    altered_scope = {**plan, 'whole_H_scope': M.POINT}
    wrong_count = {**plan, 'components': plan['components'][:-1]}
    checks = {
        'positive_true_os_then_typed_168_leaf_hashes': call(),
        'reject_boolean_OS_and_float_pid': reject(lambda: call(doc={**obs, 'observed_exit_code': False}))
            and reject(lambda: call(doc={**obs, 'unified_exec_session': 29.0}))
            and reject(lambda: call(s={**start, 'pid': 173.0})),
        'reject_cross_case_and_mixed_start_phases': reject(lambda: call(doc={**obs, 'case_id': 'other'}))
            and reject(lambda: call(s={**start, 'bindings': {'fixture': 'other-run'}})),
        'reject_failed_terminal_and_changed_start_hash': reject(lambda: call(t={**terminal, 'exit_code': 1}))
            and reject(lambda: call(doc={**obs, 'started_sha256': '0'*64}), 'fixture process pin'),
        'reject_partial_wrong_scope_wrong_slot_and_corrupt_leaf': reject(lambda: call(p=wrong_count))
            and reject(lambda: call(p=altered_scope)) and reject(lambda: call(p=wrong_first))
            and reject(lambda: call(corrupt=ordered[-1]['file']), 'fixture leaf integrity'),
    }
    require(all(checks.values()), 'incomplete admission fixtures')
    return checks


def test_real_cold_phase_gate(R, M, catalogue, joined, profiles, pins):
    wave = R.E.WAVES[0]
    case = profiles[wave]['case_order'][0]
    args = SimpleNamespace(wave=wave, case_id=case, protocol_sha256=PROTOCOL, preflight_sha256=PREFLIGHT,
                           source_review_sha256='f'*64, observed_os_sha256='0'*64)
    approval = {'status': R.STATUS, 'checks': {'synthetic_only': True},
                'source_sha256': SOURCE, 'protocol_sha256': PROTOCOL, 'preflight_sha256': PREFLIGHT,
                'later_release_sha256': R.RELEASE, 'profiles_sha256': R.E.PROFILE_SHA,
                'new_spatial_integrations_authorized': False}
    original_load = M.load

    class ReachedPlanAdmission(Exception):
        pass

    def synthetic_load(path, want, links, parse=True):
        if path == R.APPROVAL:
            require(want == args.source_review_sha256, 'review digest forwarding')
            links[M.rel(path)] = want
            return approval
        return original_load(path, want, links, parse)

    checks = {}
    with patch.object(R, 'inputs', side_effect=lambda: (M, catalogue, joined, profiles, pins.copy())), \
            patch.object(M, 'load', side_effect=synthetic_load), \
            patch.object(R, 'admitted', side_effect=ReachedPlanAdmission):
        try:
            R.run(args)
        except ReachedPlanAdmission:
            checks['real_protocol_and_preflight_reach_admission_after_mocked_root_gate'] = True
        require(bool(checks), 'positive phase gate failed')
        approval['new_spatial_integrations_authorized'] = True
        checks['reject_root_gate_with_spatial_authorization'] = reject(lambda: R.run(args), 'approved current')
        approval['new_spatial_integrations_authorized'] = False
        args.preflight_sha256 = '0'*64
        checks['reject_changed_actual_preflight_hash'] = reject(lambda: R.run(args), 'changed checked bytes')
        args.preflight_sha256 = PREFLIGHT
        args.case_id = M.FIRST[0]
        checks['reject_first_wave_before_opening_any_plan'] = reject(lambda: R.run(args), 'explicitly assigned')
    return checks


def main():
    p, f = (PREP / name for name in ('protocol.json', 'preflight.json'))
    require([sha256(x.read_bytes()).hexdigest() for x in (PREP / 'replay.py', p, f)]
            == [SOURCE, PROTOCOL, PREFLIGHT], 'exact cold phase revision')
    import replay as R

    protocol, preflight = (json.loads(x.read_bytes()) for x in (p, f))
    M, catalogue, joined, profiles, pins = R.inputs()
    require(protocol['input_sha256'] == pins and M.yes(preflight['checks'])
            and all(M.zero(doc[key]) for doc in (protocol, preflight)
                    for key in ('actual_plan_payloads_opened', 'force_H_AD_moment_or_solve_jobs')),
            'same immutable cold authority and zero scientific jobs')
    modes, route_count = check_input_arithmetic(R, M, catalogue, joined, profiles, protocol)
    case = profiles[R.E.WAVES[0]]['case_order'][0]
    checks = {
        'actual_cold_source_and_all_frozen_input_hashes_match': True,
        'seven_profiles_partition_28_and_never_reuse_first_wave': True,
        'all_case_trial_box_contexts_and_exact_independent_tail_budgets': True,
        'reviewed_order_tail_kernels_reused_without_original_FIRST_mutation': True,
        'replay_guards_producer_evaluation_force_H_moment_and_AD_entry_points': True,
        **check_admission(R, M, case),
        **test_real_cold_phase_gate(R, M, catalogue, joined, profiles, pins),
    }
    require(all(v is True for v in checks.values()), 'all independent checks')
    receipt = {
        'status': 'PASS_INDEPENDENT_COLD_REMAINING28_SELECTOR_REPLAY_SOURCE',
        'source_sha256': SOURCE, 'protocol_sha256': PROTOCOL, 'preflight_sha256': PREFLIGHT,
        'review_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'later_release_sha256': R.RELEASE, 'profiles_sha256': R.E.PROFILE_SHA,
        'checks': checks, 'permitted_waves': list(R.E.WAVES), 'reviewed_route_count': route_count,
        'cached_input_authority_kind_counts': dict(modes),
        'actual_selector_plan_payloads_read': 0, 'actual_scientific_jobs': 0,
        'new_spatial_integrations_authorized': False,
        'scope': 'Cold source/math and synthetic process/phase fixtures. Actual 168-leaf plan replay for each '
                 'case still requires root-observed producer OS success and a separately approved root source gate.',
    }
    print(M.save(HERE / 'receipt.json', receipt))


if __name__ == '__main__':
    main()
