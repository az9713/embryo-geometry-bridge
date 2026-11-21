"""Bind the verified local energy-minimum component and its rendered delivery."""
from pathlib import Path
from fractions import Fraction as Q
import datetime
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'runs/cohort_energy_stability_v1'
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
truth = lambda d: isinstance(d, dict) and bool(d) and all(v is True for v in d.values())


def main():
    output = ROOT/'evidence/cohort_energy_stability_final_delivery.json'
    pending = read(output)
    if pending['status'] != 'PENDING_FINAL_PUBLICATION_BROWSER_AND_BINDING_CHECKS':
        raise RuntimeError('Preserve completed delivery receipt')
    producer = read(RUN/'protocol.json')
    review_protocol = read(ROOT/'runs/cohort_energy_stability_review_v1/protocol.json')
    review = read(ROOT/'evidence/cohort_energy_stability_independent_review.json')
    docs = read(RUN/'documentation_receipt.json')
    browser = read(ROOT/'evidence/readiness_snapshot_verification.json')
    raw = read(ROOT/'evidence/raw_prose_math_scan.json')
    parent = read(ROOT/'evidence/cohort_energy_parent_landing_verification.json')
    space = read(ROOT/'evidence/cohort_energy_space_update.json')
    inspection = read(RUN/'root_visual_inspection.json')
    publication = read(RUN/'publication_source_receipt.json')
    rows = review['record_audits']
    snapshots = browser['snapshot_sha256']
    bindings = [producer['source_hashes'], producer['input_hashes'], review_protocol['file_hashes'], docs['bindings']]
    checks = {
        'all_producer_and_review_scientific_inputs_current': all(sha(ROOT/p) == h for b in bindings for p, h in b.items()),
        'independent_review_PASS_all17_136_record_checks_and5_summary': review['status'] == 'PASS' and truth(review['checks']) and len(rows) == 17 and all(r['audit']['status'] == 'PASS' and len(r['audit']['checks']) == 8 and truth(r['audit']['checks']) and sha(ROOT/r['review_file']) == r['review_sha256'] for r in rows),
        'all17_safe_exact_q_and_m_displays': all(Q(r['audit']['independent_q_infinity_upper_exact']) < Q(1, 1250000000) and Q(r['audit']['uniform_Hessian_lower_bound_exact']) > Q(197, 500000) for r in rows),
        'publication_draft_review_checks_pass_and_exact_copy': truth(docs['checks']) and sha(ROOT/'docs/COHORT-ENERGY-STABILITY.md') == sha(RUN/'publication-draft.md'),
        'published_controls_match_reviewed_source_edits': all(sha(ROOT.parent/p) == h for p, h in publication['after_control_hashes_relative_to_repository'].items()),
        '103_pages206_viewports_fresh_PASS': browser['status'] == 'PASS' and truth(browser['checks']) and browser['html_pages'] == 103 and browser['page_viewport_checks'] == 206,
        'actual_rendered_snapshot_still_current': all(sha(ROOT/p) == h for p, h in snapshots.items()),
        'browser_receipt_current': sha(ROOT/'evidence/report_browser_verification.json') == browser['browser_receipt_sha256'],
        'raw_prose_scan_PASS103_current_HTML': raw['status'] == 'PASS' and raw['html_pages'] == 103 and all(sha(ROOT/'docs'/r['page']) == r['HTML_sha256'] for r in raw['rows']),
        'new_proof_has61_rendered_math_spans_both_widths': all(r['math'] == r['rendered'] == 61 and r['pass'] for r in read(ROOT/'evidence/report_browser_verification.json')['checks'] if r['page'] == 'COHORT-ENERGY-STABILITY.html'),
        'parent_landing_fresh2_widths_exact_bytes': parent['status'] == 'PASS' and len(parent['checks']) == 2 and all(r['pass'] for r in parent['checks']) and sha(ROOT.parent/'docs/index.html') == parent['page_sha256'] and sha(ROOT.parent/'scripts/build_dossier.py') == parent['builder_sha256'],
        'root_visual_inspection_current': truth(inspection['checks']) and sha(RUN/'cohort-energy-stability.png') == inspection['image_sha256'],
        'Space_readback98_blocks_native_figure_hash_current': space['status'] == 'PASS' and truth(space['checks']) and space['after_blocks'] == 98 and space['visual']['sha256'] == sha(RUN/'cohort-energy-stability.png'),
        'unchanged68_output_mathematical_receipt': sha(ROOT/'evidence/shared_cell_cohort_inference.json') == '1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677',
        'full_five_gate_study_and_novelty_still_incomplete': read(ROOT/'evidence/gates.json')['goal_achieved'] is False and review['original_five_gate_goal_complete'] is False and review['confirmed_novelty'] is False,
    }
    # Avoid a vacuous per-page all(): demand exactly the two expected page rows.
    target_rows = [r for r in read(ROOT/'evidence/report_browser_verification.json')['checks'] if r['page'] == 'COHORT-ENERGY-STABILITY.html']
    checks['both_new_proof_widths_explicitly_present'] = len(target_rows) == 2 and {r['width'] for r in target_rows} == {390, 1440}
    if not truth(checks):
        raise RuntimeError(json.dumps(checks))
    names = [
        'docs/COHORT-ENERGY-STABILITY.md', 'docs/COHORT-ENERGY-STABILITY.html',
        'docs/index.html', '../README.md', '../docs/index.html',
        'evidence/gates.json', 'evidence/cohort_energy_stability_independent_review.json',
        'evidence/cohort_energy_stability_theory_review.json',
        'evidence/cohort_energy_parent_landing_verification.json',
        'evidence/cohort_energy_space_update.json',
        'evidence/readiness_snapshot_verification.json',
        'evidence/report_browser_verification.json', 'evidence/raw_prose_math_scan.json',
        'runs/cohort_energy_stability_v1/results.json',
        'runs/cohort_energy_stability_v1/terminal.json',
        'runs/cohort_energy_stability_v1/documentation_receipt.json',
        'runs/cohort_energy_stability_v1/publication_source_receipt.json',
        'runs/cohort_energy_stability_v1/root_visual_inspection.json',
        'runs/cohort_energy_stability_v1/cohort-energy-stability.png',
        'runs/cohort_energy_stability_review_v1/scope_and_delivery.json',
    ]
    hashes = {p: sha(ROOT/p) for p in names}
    archived = RUN/'pending_final_delivery.json'
    if archived.exists():
        raise RuntimeError('Preserve pending receipt archive')
    archived.write_bytes(output.read_bytes())
    receipt = {
        'status': 'PASS_COMPONENT_DELIVERY', 'checks': checks,
        'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'file_hashes_relative_to_child': hashes,
        'rendered_snapshot_hashes': snapshots,
        'source_sha256': sha(__file__),
        'full_original_five_gate_goal_complete': False, 'confirmed_novelty': False,
        'scope': 'Strict local energy minima for the17 selected axisymmetric finite-quadrature root families on the same tiny parameter cell. Source-matched inference, stable physical paths, nonaxisymmetric/PDE modes, continuum error and revised larger-region validation remain incomplete. Browser tests cover103 child pages and only the parent landing. Space is updated; no public GitHub deployment claimed.',
        'prior_delivery_history': 'The previous102-page stationary-component snapshot and changed controls were archived before publication; old mathematical proofs and failures remain unchanged. Current component supersedes that documentation snapshot without upgrading the full study.',
    }
    output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': receipt['status'], 'checks': len(checks), 'sha256': sha(output)}))


if __name__ == '__main__':
    main()
