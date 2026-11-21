"""Compose current receipts for a completed component, never the five-gate goal."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COHORT = 'runs/exact_integration_cohort_v1/'
SCALED = COHORT + 'L_defect_diagnostic/scaled_offaxis_L_v1/'
PINS = {
    'runs/transfer_pair_pilot_v1/results.json': '6a2fcb2a0c4b7a51301d025a7725247ec81a999764b7e85a0e97b5f6e139be05',
    'runs/transfer_pair_pilot_v1/independent_actual/actual_review.json': 'bfe94d44f4634c24747755ec2183d2ef9a925f12e7098472b18fef94495b8c18',
    'runs/transfer_pair_pilot_v1/independent_actual/analysis_readback.json': '5841013f85060a27ca82a7f548583920f94e0e81d49867d66bafde919ba002a2',
    SCALED + 'actual/terminal.json': 'c768fe8cce6978ec2b6865c784dc29484121ed4e284bb84d7358df58ac4a2232',
    SCALED + 'actual/results.json': '338d7ba6806d5d9c3013d9559e0cee4333031545409fd679b799b125b6e72c95',
    SCALED + 'actual_independent_review.json': 'bcd6fd98baf2d611ce224054e51ce5cdff366e41547ec09a5f500e9b4cd74bb2',
    SCALED + 'reduction_preparation/results.json': '51d12c1ee61901f9164dbd3ab4e366b4887dc115398446d379835bdeb4df094a',
    SCALED + 'reduction_preparation/actual_independent_review.json': '54dee8659601c387712c16732ab6aed90fbf313c1e3b0f77930b3e12b19220e6',
    COHORT + 'publication_v1/scientific_readback/review.json': 'f98d102bdeee8d886ac2efb627052d0725568f57336b2dc1fc6de00269645a71',
    COHORT + 'publication_v1/scaled_completion_scientific_readback.json': '5500e4735ed2db6d67685ab16e1a28df6d41980c2e025df09131979970707c3f',
    COHORT + 'publication_v1/browser_readback.json': 'aeaf32f257282e37a2e22c680edaa6b178fdb97f93a0bd4eae5863594fc84b26',
    COHORT + 'publication_v1/scaled_completion_readback/browser_readback.json': 'bad6c30d0c0f397e64af4cecac98e7fa3c7650325f2b84c2c032e31dc9380f38',
    'runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json': '42d97691311c057d8152f54a6ea31a7218763e7b451df63358bdf423010e18d4',
    'PLAN.md': '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    bindings = {}

    def read(relative, expected=None):
        path = ROOT / relative
        assert path.resolve().is_relative_to(ROOT), relative
        actual = digest(path)
        assert expected is None or actual == expected, relative
        assert relative not in bindings or bindings[relative] == actual, relative
        bindings[relative] = actual
        raw = path.read_bytes()
        if path.suffix == '.json':
            return json.loads(raw)
        return raw.decode('utf8') if path.suffix in ('.md', '.py', '.html', '.txt') else raw

    records = {p: read(p, h) for p, h in PINS.items()}
    transfer = records['runs/transfer_pair_pilot_v1/results.json']
    assert transfer['status'] == 'COMPLETE_FOUR_CONDITIONAL_FIELDS'
    assert len(transfer['terminal_case_receipts']) == 4
    assert all(r['terminal'] is True and r['exit_code'] == 0 for r in transfer['terminal_case_receipts'])
    assert transfer['root_validated'] is False and transfer['source_compatible_witness'] is False
    assert transfer['novelty_established'] is False and transfer['source_outcomes_read_or_scored'] is False
    scaled = records[SCALED + 'reduction_preparation/results.json']
    independent = records[SCALED + 'reduction_preparation/actual_independent_review.json']
    assert independent['status'] == 'PASS_ACTUAL_SCALED_ORIGINAL_CUBE_REDUCTION_REVIEW'
    assert independent['result_sha256'] == PINS[SCALED + 'reduction_preparation/results.json']
    assert all(value is True for value in independent['checks'].values())
    assert scaled['component_count'] == 168 and scaled['force_entries'] == 611 and scaled['upper_Hessian_entries'] == 6676
    for relative, expected in scaled['source_inputs_and_packets'].items():
        read(relative, expected)
    eta0, z0, rf, b, r = [Q(scaled[k]) for k in (
        'cached_eta_upper_exact', 'cached_Z_upper_exact', 'RF_uniform_upper_exact',
        'b_upper_exact', 'normal_radius_exact')]
    assert Q(scaled['z_upper_exact']) == z0 + b < 1
    assert Q(scaled['selfmap_lhs_exact']) == eta0 + rf + (z0 + b) * r > r
    assert eta0 + (z0 + b) * r < r
    assert scaled['sufficient_conditions_pass'] is False
    assert scaled['conditional_root_error_upper_exact'] is None and scaled['conditional_curvature'] is None
    initial = read('evidence/readiness_snapshot_verification.json')
    assert initial['status'] == 'PASS' and initial['html_pages'] == 109 and initial['page_viewport_checks'] == 218
    delta = records[COHORT + 'publication_v1/scaled_completion_readback/browser_readback.json']
    assert delta['status'] == 'PASS' and delta['documents_unchanged'] is True
    assert {(row['page'], row['width']) for row in delta['checks']} == {
        ('controlled_furrowing/docs/COHORT-INTEGRATION-L.html', width) for width in (390, 1440)}
    assert all(row['pass'] is True for row in delta['checks'])
    for relative, expected in initial['snapshot_sha256'].items():
        if relative not in ('docs/COHORT-INTEGRATION-L.md', 'docs/COHORT-INTEGRATION-L.html'):
            read(relative, expected)
    for relative, expected in delta['document_hashes'].items():
        assert digest(ROOT.parent / relative) == expected
    for relative in ('docs/COHORT-INTEGRATION-L.md', 'docs/COHORT-INTEGRATION-L.html',
                     'docs/TRANSFER-PAIR-PILOT.md', 'docs/TRANSFER-PAIR-PILOT.html', 'WORK-STATE.md'):
        read(relative)
    for relative in ('evidence/cohort_L_transfer_space_update.json',
                     'evidence/cohort_L_scaled_completion_space_update.json'):
        record = read(relative)
        assert record['status'] == 'PASS_SPACE_READBACK' and all(record['checks'].values())
    assert read('evidence/cohort_L_scaled_completion_space_update.json')['sequence'] == 21
    assert 'ACTIVE / INCOMPLETE' in read('WORK-STATE.md')
    assert read('evidence/gates.json').get('goal_achieved') is not True
    bindings[Path(__file__).relative_to(ROOT).as_posix()] = digest(Path(__file__))
    target = ROOT / 'evidence/cohort_L_transfer_component_delivery.json'
    receipt = {
        'status': 'PASS_COMPLETED_COMPONENT_NOT_FULL_GOAL',
        'checks': {'four_actual_zero_exit_workers_and_no_source_or_novelty_promotion': True,
                   'same_168_component_model_and_exact_scaled_reduction_review': True,
                   'strict_selfmap_failure_and_remaining_zero_force_margin_replayed': True,
                   'current_browser_coverage_108_unchanged_pages_plus_fresh_L_delta': True,
                   'original_B_result_and_five_gate_plan_unchanged': True,
                   'Space_exact_readbacks_and_active_incomplete_goal_retained': True},
        'current_source_and_artifact_hashes': bindings,
        'scientific_outcome': 'Two conditional transfer candidates differ slightly in stored Q. '
            'Scaled L enclosures improve substantially but root existence in the original cube remains unresolved.',
        'full_goal_achieved': False,
        'scope': 'Source/receipt identity, exact scalar composition and current publication coverage. '
                 'Relies on the separately retained mathematical and saved-field reviews. '
                 'No new root/PDE/physical membership/novelty theorem or GitHub deployment.',
    }
    with target.open('x', encoding='utf8') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'bindings': len(bindings),
                      'receipt_sha256': digest(target)}))


if __name__ == '__main__':
    main()
