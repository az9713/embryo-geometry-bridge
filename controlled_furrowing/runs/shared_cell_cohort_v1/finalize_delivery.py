"""Close the shared-cell component receipt without claiming the five gates."""
from pathlib import Path
import datetime
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
truth = lambda d: isinstance(d, dict) and bool(d) and all(v is True for v in d.values())


def main():
    out = ROOT/'evidence/shared_cell_final_delivery.json'
    old = read(out)
    if old['status'] != 'PENDING_FINAL_RENDERING_AND_BINDING_CHECKS':
        raise RuntimeError('Preserve a completed delivery receipt')
    science = read(ROOT/'runs/shared_cell_cohort_v1/protocol.json')
    audit = read(ROOT/'evidence/shared_cell_cohort_independent_review.json')
    rendering = read(ROOT/'evidence/rendering_guard_independent_review.json')
    browser = read(ROOT/'evidence/readiness_snapshot_verification.json')
    independent = read(ROOT/'evidence/shared_cell_delivery_independent_review.json')
    prose = read(ROOT/'evidence/raw_prose_math_scan.json')
    docs = read(ROOT/'runs/shared_cell_cohort_v1/documentation_receipt.json')
    parent = read(ROOT/'evidence/parent_index_checkpoint_verification.json')
    space = read(ROOT/'evidence/space_shared_cell_update.json')
    expected = {p: v for p, v in browser['snapshot_sha256'].items()}
    checks = {
        'all71_scientific_inputs_current': len(science['file_hashes']) == 71 and all(sha(ROOT/p) == h for p, h in science['file_hashes'].items()),
        'independent_mathematical_combination_checks_pass': truth(audit['checks']),
        'rendering_delta_checks_pass': rendering['status'] == 'PASS' and truth(rendering['checks']),
        'final_delivery_independent_review_pass': independent['status'].startswith('PASS') and truth(independent['checks']),
        'all102_pages204_viewports_pass': browser['status'] == 'PASS' and truth(browser['checks']) and browser['html_pages'] == 102 and browser['page_viewport_checks'] == 204,
        'actual_rendered_document_snapshot_still_current': all(sha(ROOT/p) == h for p, h in expected.items()),
        'browser_report_hash_current': sha(ROOT/'evidence/report_browser_verification.json') == browser['browser_receipt_sha256'],
        'raw_prose_command_scan_pass': prose['status'] == 'PASS' and prose['html_pages'] == 102,
        'documentation_and_visual_bindings_current': docs['status'].startswith('PASS') and truth(docs['checks']) and all(sha(ROOT/p) == h for p, h in docs['bindings'].items()),
        'parent_index_actual_tested_bytes_current': parent['status'] == 'PASS' and all(c['pass'] for c in parent['checks']) and sha(ROOT.parent/'docs/index.html') == parent['page_sha256'] and sha(ROOT.parent/'scripts/build_dossier.py') == parent['builder_sha256'],
        'space_native_readback_and_figure_bindings_pass': space['status'] == 'PASS' and truth(space['checks']) and space['after_blocks'] == 91,
        'full_study_still_explicitly_incomplete': read(ROOT/'evidence/gates.json')['goal_achieved'] is False,
    }
    if not truth(checks):
        raise RuntimeError(json.dumps(checks))
    files = [
        'evidence/shared_cell_cohort_inference.json',
        'evidence/shared_cell_cohort_independent_review.json',
        'evidence/shared_cell_delivery_independent_review.json',
        'evidence/readiness_snapshot_verification.json',
        'evidence/report_browser_verification.json',
        'evidence/raw_prose_math_scan.json',
        'evidence/rendering_guard_independent_review.json',
        'evidence/documentation_rendering_corrections.json',
        'evidence/parent_index_checkpoint_verification.json',
        'evidence/space_shared_cell_update.json',
        'evidence/source_workflow_availability_refresh.json',
        'runs/shared_cell_cohort_v1/documentation_receipt.json',
        'runs/shared_cell_cohort_v1/protocol.json',
        'docs/SHARED-CELL-COHORT.md', 'docs/SHARED-CELL-COHORT.html',
        '../README.md', 'evidence/gates.json', 'scripts/build_report.py',
        'scripts/verify_report.py',
    ]
    # Validate every hash before the first write. An explicit revision-0 archive
    # was already created by the retained failed parent-README path attempt.
    hashes = {p: sha(ROOT/p) for p in files}
    archive = ROOT/'runs/shared_cell_cohort_v1/delivery_revision_0'
    archive.mkdir(exist_ok=True)
    archived = archive/'shared_cell_final_delivery.json'
    if archived.exists():
        if archived.read_bytes() != out.read_bytes() or not (archive/'failure.json').exists():
            raise RuntimeError('Explicit archived-pending recovery identity mismatch')
    else:
        archived.write_bytes(out.read_bytes())
    packet = {
        'status': 'PASS_COMPONENT_DELIVERY',
        'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'component': old['component'],
        'checks': checks,
        'file_hashes': hashes,
        'finalizer_sha256': sha(__file__),
        'rendered_snapshot_hashes': expected,
        'parent_scope': 'Parent landing page only, two viewports; not all parent scientific pages.',
        'historical_control_binding_delta': 'The 180-check mathematical review is immutable. Its old builder binding alone is superseded by the separate independently tested pipe-escaping correction. No mathematical inputs changed.',
        'full_original_five_gate_objective_complete': False,
        'confirmed_novelty': False,
        'scope': 'Complete delivery of the local conditional stationary-root/static-output component on one tiny cell. This receipt does not establish physical stability, source-time correspondence, full parameter-domain coverage, continuum accuracy, or revised larger-region prediction. Future documentation revisions must be verified as new snapshots.',
    }
    out.write_text(json.dumps(packet, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': packet['status'], 'checks': len(checks), 'sha256': sha(out)}))


if __name__ == '__main__':
    main()
