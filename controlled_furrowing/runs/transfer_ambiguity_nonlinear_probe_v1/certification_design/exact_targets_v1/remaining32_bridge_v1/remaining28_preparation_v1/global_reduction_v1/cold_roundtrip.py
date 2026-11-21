"""Cold wire-format fixtures for the later-wave driver. No real plans or hulls."""

from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import reduce as R


def check_native(paths):
    _, profile = paths
    wave = 'wave_02'
    case = profile[wave]['case_order'][0]
    target = profile[wave]['raw_worker_targets'][case]
    phase = {'pid': 42, 'case_id': case, 'wave_id': wave, 'source_sha256': 'frozen'}
    plan = {'plan_sha256': 'a'*64, 'independent_plan_review_sha256': 'b'*64}
    raw = {'pid': 32, 'unified_exec_session': 45, 'sha256': {str(i): 'c'*64 for i in range(715)}}
    start = {**phase, 'status': 'RUNNING'}
    terminal = {**phase, 'status': 'FINISHED', 'exit_code': 0, 'receipt_sha256': 'd'*64}
    obs = {'receipt_sha256': 'd'*64}
    receipt = {**phase, 'status': R.J.RAW_PASS, 'target': target, 'record_id': target['record_id'],
               'radius_exact': str(R.RADIUS), 'full_independent_hulls_computed': True, 'root_asserted': False,
               'overlap_or_producer_containment_used_as_proof': False,
               'admit_hull_only_global_reduction_after_separate_review': True,
               'disjoint_counts': dict.fromkeys(R.D.CHANNELS, 0), 'selector_plan_sha256': plan['plan_sha256'],
               'plan_review_sha256': plan['independent_plan_review_sha256'],
               'producer_observed_pid': raw['pid'], 'producer_observed_unified_exec_session': raw['unified_exec_session'],
               'raw_root_observed_mapping_sha256': R.S.E._digest(raw['sha256']),
               'raw_root_observed_component_files': 715,
               'case_components': [{'kind': k, 'slot': i, 'status': R.D.LEAF_PASS} for k, i in R.S.OPS]}
    R.check_native(start, terminal, receipt, obs, phase, target, plan, raw)
    reject = R.V.A.reject
    bogus = deepcopy(receipt)
    bogus['case_components'][0]['slot'] = False
    return {
        'positive_native_metadata_and_all_168_typed_slots': True,
        'reject_false_zero_mixed_phase_and_altered_target': reject(
            lambda: R.check_native(start, {**terminal, 'exit_code': False}, receipt, obs, phase, target, plan, raw))
            and reject(lambda: R.check_native({**start, 'source_sha256': 'other'}, terminal, receipt, obs,
                       phase, target, plan, raw))
            and reject(lambda: R.check_native(start, terminal, {**receipt, 'target': {}}, obs, phase, target, plan, raw)),
        'reject_boolean_slot_and_any_disjoint': reject(
            lambda: R.check_native(start, terminal, bogus, obs, phase, target, plan, raw))
            and reject(lambda: R.check_native(start, terminal, {**receipt, 'disjoint_counts': {
                           **receipt['disjoint_counts'], 'point_force': 1}}, obs, phase, target, plan, raw))}


def run_stage_gates(paths):
    p, profiles = paths
    sha = R.sha
    args = SimpleNamespace(mode='producer', wave='wave_02',
                           protocol_sha256=sha(R.HERE / 'protocol.json'), preflight_sha256=sha(R.HERE / 'preflight.json'),
                           review_sha256='a'*64, hull_admission_sha256='b'*64, producer_os_sha256='c'*64,
                           raw_wave_admission_sha256='d'*64, raw_root_observed_exit_sha256='e'*64)
    approval = {'status': R.SOURCE_PASS, 'source_sha256': sha(R.SOURCE), 'protocol_sha256': args.protocol_sha256,
                'preflight_sha256': args.preflight_sha256, 'later_waves_release_sha256': R.RELEASE,
                'profiles_sha256': R.X.PROFILE_SHA, 'permitted_waves': list(R.X.WAVES),
                'approved_modes': ['producer', 'independent'], 'checks': {'inert_fixture_only': True},
                'new_spatial_integrations_authorized': False}
    read = R.read

    class Boundary(Exception):
        pass

    def fake_read(path, pin):
        if path == R.AUTH and pin == args.review_sha256:
            return approval
        return read(path, pin)

    def replay_allowed(mode):
        args.mode = mode
        with patch.object(R, 'static', return_value=(p, profiles)), patch.object(R, 'read', side_effect=fake_read), \
                patch.object(R, 'producer_observed', return_value={'inert_fixture_only': True}) as prod, \
                patch.object(R, 'open_wave', side_effect=Boundary) as hull:
            try:
                R.run(args)
            except Boundary:
                if mode == 'independent':
                    prod.assert_called_once()
                else:
                    prod.assert_not_called()
                hull.assert_called_once()
                return True
        raise ValueError('missed metadata boundary')

    checks = {'real_protocol_and_preflight_roundtrip_producer_gate': replay_allowed('producer'),
              'independent_mode_requires_observed_producer': replay_allowed('independent')}
    with patch.object(R, 'static', return_value=(p, profiles)), patch.object(R, 'read', side_effect=fake_read), \
            patch.object(R, 'producer_observed', side_effect=AssertionError('after source guard')):
        approval['new_spatial_integrations_authorized'] = True
        checks['reject_unapproved_science_and_mode_before_actual_metadata'] = R.V.A.reject(lambda: R.run(args))
    return checks


def main():
    paths = R.static()
    checks = {**check_native(paths), **run_stage_gates(paths)}
    R.require(all(checks.values()), 'wire fixtures failed')
    path = Path(__file__).resolve()
    h = R.save(path.with_name('cold_roundtrip_receipt.json'), {
        'status': 'PASS_COLD_LATER_GLOBAL_WIRE_FIXTURES', 'source_sha256': R.sha(R.SOURCE),
        'fixture_source_sha256': R.sha(path), 'protocol_sha256': R.sha(R.HERE / 'protocol.json'),
        'preflight_sha256': R.sha(R.HERE / 'preflight.json'), 'checks': checks,
        'real_hull_leaf_reads': 0, 'real_global_reductions': 0})
    print(json.dumps({'receipt_sha256': h}))


if __name__ == '__main__':
    main()
