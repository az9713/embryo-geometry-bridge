"""Compose the independently reviewed L theorem and its current publication.

This is a component delivery check. It does not certify any of the remaining
15 roots, a continuum solution, biological calibration, or the five-gate goal.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
STAGE = 'runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/'
PINS = {
    STAGE + 'actual_independent_review.json': 'ea1ebcb2f1788eaa3fd491cca6d519e8f40568e9cdc870919a9b32f9f4adc0a9',
    STAGE + 'independent_actual/results.json': 'bdfa3b85849c47b7d3a80424c8f407711d6531ab7a43b9e5f117abd439a7a5d3',
    STAGE + 'final_mathematical_readback/receipt.json': '524a6bf025b9efff519a7d16d0df57d3901b9f35635ba57d0ff05cccf021f079',
    STAGE + 'publication/scientific_readback/receipt.json': '8b4fbd16f751da6a9e818773f4b0ff2c582cb213c3177c9b3e0a077f6ddc7001',
    STAGE + 'publication/receipt.json': 'f038f73b5101845a608bd405c66a3c30c0945fc67479f1149b5812fa4693b4a0',
    STAGE + 'publication/browser_readback.json': '56103ea64cdb209de5b218d02b814d74dbeefd232da529588b48b03e324e5be6',
    'runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json': '42d97691311c057d8152f54a6ea31a7218763e7b451df63358bdf423010e18d4',
    'PLAN.md': '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    'runs/exact_integration_cohort_v1/remaining15_integrals_preparation/manifest.json': 'e4cbc8826b763d235e24ff827d684a00526b79e18bef63b51bc74f727dcdd363',
    'runs/exact_integration_cohort_v1/remaining15_integrals_preparation/prelaunch_independent_review.json': 'e89922fe616b089f1723757b7b4f79f2d01da8b7675b32d62c81e17bbb1f0542',
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    bindings = {}

    def read(relative, expected=None):
        path = (ROOT / relative).resolve()
        assert path.is_relative_to(REPO), relative
        raw = path.read_bytes()
        actual = sha256(raw).hexdigest()
        assert expected is None or actual == expected, relative
        name = path.relative_to(REPO).as_posix()
        assert name not in bindings or bindings[name] == actual, relative
        bindings[name] = actual
        if path.suffix == '.json':
            return json.loads(raw)
        return raw.decode('utf8') if path.suffix in ('.md', '.html', '.py') else raw

    records = {p: read(p, h) for p, h in PINS.items()}
    certificate = records[STAGE + 'actual_independent_review.json']
    independent = records[STAGE + 'independent_actual/results.json']
    math = records[STAGE + 'final_mathematical_readback/receipt.json']
    science = records[STAGE + 'publication/scientific_readback/receipt.json']
    assert certificate['status'] == 'PASS_INDEPENDENT_SIGNED_ORIGINAL_CUBE_ROOT_CERTIFICATE'
    assert certificate['root_certified'] is True and certificate['strict_local_minimum_certified'] is True
    assert certificate['sufficient_test_pass'] is True and all(certificate['checks'].values())
    assert certificate['same_closed_C_sha256'] == '64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf'
    for mapping in (certificate['input_hashes'], math['sha256'], science['bound_sha256']):
        for name, expected in mapping.items():
            read(name, expected)
    eta, z, r, error, m = [Q(certificate[k]) for k in (
        'uniform_eta_I_upper_exact', 'same_box_z_exact', 'original_radius_exact',
        'root_error_upper_exact', 'uniform_curvature_lower_exact')]
    norm = Q(independent['R_infinity_norm_exact'])
    lhs = eta + z * r
    assert eta >= 0 and 0 <= z < 1 and lhs == Q(certificate['uniform_selfmap_lhs_exact']) < r
    assert error == eta / (1 - z) < r and m == (1 - z) / norm > 0
    assert len(independent['independent_full_RF_predictor_intervals']) == 611
    assert max(max(abs(Q(a)), abs(Q(b))) for a, b in independent['independent_full_RF_predictor_intervals']) == eta
    assert science['full_five_gate_goal_achieved'] is False and math['full_five_gate_goal_completed'] is False
    assert all(science['checks'].values()) and all(math['checks'].values())
    snapshot = read('evidence/readiness_snapshot_verification.json')
    assert snapshot['status'] == 'PASS' and snapshot['html_pages'] == 110 and snapshot['page_viewport_checks'] == 220
    delta = records[STAGE + 'publication/browser_readback.json']
    assert delta['status'] == 'PASS' and delta['documents_unchanged'] is True
    assert len(delta['checks']) == 12 and all(row['pass'] for row in delta['checks'])
    assert {(r['page'], r['width']) for r in delta['checks'] if r['page'].endswith('/L-CENTERED-ROOT.html')} == {
        ('controlled_furrowing/docs/L-CENTERED-ROOT.html', width) for width in (390, 1440)}
    for name, expected in snapshot['snapshot_sha256'].items():
        if name not in ('docs/L-CENTERED-ROOT.md', 'docs/L-CENTERED-ROOT.html'):
            read(name, expected)
    for name, expected in delta['document_hashes'].items():
        read('../' + name, expected)
    archive = Path(STAGE) / 'publication/prepublication_snapshot'
    before = read((archive / 'manifest.json').as_posix())
    assert len(before['files']) == 13
    for row in before['files']:
        read((archive / row['source']).as_posix(), row['sha256'])
    space = read('evidence/L_signed_root_space_update.json')
    assert space['status'] == 'PASS_SPACE_EXACT_READBACK' and all(space['checks'].values())
    assert space['sequence'] == 23 and space['current_block_count'] == 150
    assert space['independent_L_root_receipt_sha256'] == PINS[STAGE + 'actual_independent_review.json']
    assert len(space['figures']) == 2 and all(f['file_access_confirmed'] for f in space['figures'])
    for filename, figure in zip(('L_self_map_stages.png', 'L_signed_force_budget.png'), space['figures']):
        read(STAGE + 'publication/' + filename, figure['sha256'])
    assert 'ACTIVE / INCOMPLETE' in read('WORK-STATE.md')
    assert read('evidence/gates.json').get('goal_achieved') is not True
    # Immutable launch receipts record full disjoint queues, not completed roots.
    launches = ROOT / 'runs/exact_integration_cohort_v1/remaining15_integrals_preparation/launches'
    queues = [read(p.relative_to(ROOT).as_posix()) for p in sorted(launches.glob('*/started.json'))]
    jobs = [job for queue in queues for job in queue['explicit_jobs']]
    assert len(queues) == 2 and len(jobs) == len(set(jobs)) == 30
    assert all(not job.startswith(('optoGEF2_small_col_B:', 'optoGEF2_small_col_L:')) for job in jobs)
    read(Path(__file__).relative_to(ROOT).as_posix())
    receipt = {
        'status': 'PASS_L_ROOT_COMPONENT_DELIVERY_NOT_FULL_GOAL',
        'checks': {'independent_same_C_original_cube_root_and_minimum': True,
                   'exact_independent_eta_z_radius_error_and_curvature_composed': True,
                   'current_proof_and_publication_sources_rehashed': True,
                   'current_110_page_snapshot_plus_one_page_delta_and_12_targeted_views': True,
                   'thirteen_prior_files_preserved_and_B_theorem_PLAN_unchanged': True,
                   'Space_150_blocks_sequence23_and_two_native_figure_bytes': True,
                   'all30_remaining_jobs_in_disjoint_frozen_queues_not_claimed_roots': True},
        'source_and_artifact_hashes': bindings,
        'certified_case': 'optoGEF2_small_col_L', 'full_five_gate_goal_achieved': False,
        'remaining_integral_roots_certified_by_this_receipt': 0,
        'independent_bounds_exact': {'eta_I': str(eta), 'z': str(z), 'r': str(r),
                                     'selfmap_lhs': str(lhs), 'predictor_error': str(error), 'curvature': str(m)},
        'scope': 'Composes reviewed mathematical and publication premises. All15 remaining force jobs are separate, '
                 'and this receipt does not establish their roots, full physical inference, a continuum theorem, '
                 'biological validation, originality, GitHub publication, or five-gate completion.',
    }
    target = ROOT / 'evidence/L_signed_root_final_delivery.json'
    with target.open('x', encoding='utf8') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'bindings': len(bindings), 'receipt_sha256': digest(target)}))


if __name__ == '__main__':
    main()
