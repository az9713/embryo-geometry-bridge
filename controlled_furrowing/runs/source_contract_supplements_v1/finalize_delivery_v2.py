"""Finalize the bounded public-source audit without upgrading scientific gates."""
from pathlib import Path
import datetime
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'runs/source_contract_supplements_v1'
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
truth = lambda d: isinstance(d, dict) and bool(d) and all(v is True for v in d.values())


def main():
    out = ROOT/'evidence/source_contract_final_delivery_v2.json'
    pending = read(out)
    if pending['status'] != 'PENDING_RENDERING_AND_FINAL_SOURCE_BINDINGS':
        raise RuntimeError('Preserve completed source delivery')
    assets = read(ROOT/'evidence/public_supplement_asset_manifest.json')
    binding = read(RUN/'public_source_audit_binding.json')
    peer = read(ROOT/'evidence/source_contract_peer_review_audit.json')
    small = read(ROOT/'evidence/source_contract_small_movies_audit.json')
    others = read(ROOT/'evidence/source_contract_other_movies_audit.json')
    movie3 = read(ROOT/'evidence/movie_geometry_audit.json')
    browser = read(ROOT/'evidence/readiness_snapshot_verification.json')
    raw = read(ROOT/'evidence/raw_prose_math_scan.json')
    report = read(ROOT/'evidence/report_browser_verification.json')
    space = read(ROOT/'evidence/source_contract_space_update.json')
    visual = read(RUN/'root_visual_inspection.json')
    publication = read(RUN/'publication_source_receipt.json')
    gates = read(ROOT/'evidence/gates.json')
    decoded = sum(m['decoded_frames'] for m in small['movies']) + len(movie3['frames']) + others['total_decoded_video_frames']
    audits = ['evidence/public_supplement_asset_manifest.json', 'evidence/source_contract_workbook_package_audit.json',
              'evidence/source_contract_small_movies_audit.json', 'evidence/source_contract_other_movies_audit.json',
              'evidence/source_contract_peer_review_audit.json']
    controls = publication['after_controls_relative_to_repository']
    checks = {
        'all16_unique_assets_official_checksums_and_current_SHA': len(assets['assets']) == 16 and assets['total_bytes'] == 111690608 and truth(assets['checks']) and all(sha(ROOT/a['file']) == a['sha256'] for a in assets['assets']),
        'all_source_audit_checks_literal_true': all(truth(read(ROOT/p)['checks']) for p in audits),
        'all903_frames_from_disclosed_audits': decoded == 903,
        'full26_Movie3_decode_evidence_bound_explicitly': len(movie3['frames']) == 26 and binding['source_audit_bindings'].get('evidence/movie_geometry_audit.json') == sha(ROOT/'evidence/movie_geometry_audit.json'),
        'old_movie3_original_provenance_not_relabelled_anatomical_geometry': movie3['status'].startswith('PASS') and all(c['pass'] for c in movie3['official_checksum_checks']),
        'all_document_and_source_bindings_current': truth(binding['checks']) and all(sha(ROOT/p) == h for d in (binding['bindings'], binding['source_audit_bindings']) for p, h in d.items()),
        'new_chapter_and_request_exact_copies': sha(ROOT/'docs/PUBLIC-SUPPLEMENT-AUDIT.md') == sha(RUN/'public-source-audit-draft.md') and sha(ROOT/'docs/DRAFT-DATA-REQUEST.md') == sha(RUN/'data-request-v2.md'),
        'all4_publication_controls_current': len(controls) == 4 and all(sha(ROOT.parent/p) == h for p, h in controls.items()),
        '104_pages208_viewports_fresh_PASS': browser['status'] == 'PASS' and truth(browser['checks']) and browser['html_pages'] == 104 and browser['page_viewport_checks'] == 208 and len(report['checks']) == 208 and all(r['pass'] for r in report['checks']),
        'actual_rendered_snapshot_current': all(sha(ROOT/p) == h for p, h in browser['snapshot_sha256'].items()),
        'browser_report_hash_current': sha(ROOT/'evidence/report_browser_verification.json') == browser['browser_receipt_sha256'],
        'new_chapter_and_request_both_widths_explicit': all({r['width'] for r in report['checks'] if r['page'] == name} == {390, 1440} for name in ('PUBLIC-SUPPLEMENT-AUDIT.html', 'DRAFT-DATA-REQUEST.html')),
        'raw_prose_scan104_and_actual_HTML_hashes_current': raw['status'] == 'PASS' and raw['html_pages'] == 104 and all(sha(ROOT/'docs'/r['page']) == r['HTML_sha256'] for r in raw['rows']),
        'original_visual_root_inspection_current': truth(visual['checks']) and sha(RUN/'public-source-contract.png') == visual['image_sha256'],
        'Space106_blocks_fresh_native_readback_original_hash': space['status'] == 'PASS' and truth(space['checks']) and space['after_blocks'] == 106 and space['through_sequence'] == 15 and space['visual']['sha256'] == sha(RUN/'public-source-contract.png'),
        'Fig6_Fig2_attribution_correction_retained': read(RUN/'timing_locator_correction.json')['status'] == 'CORRECTED_BEFORE_PUBLICATION_OR_CONTACT' and peer['checks']['explicit_initiation_sentence_located_in_Fig6_not_transplanted_to_Fig2'] is True,
        'original_plan_unchanged': sha(ROOT/'PLAN.md') == '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
        'energy_minimum_and68output_math_receipts_unchanged': sha(ROOT/'evidence/shared_cell_cohort_inference.json') == '1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677' and sha(ROOT/'evidence/cohort_energy_stability_independent_review.json') == 'd4f4c4d30d16bbbf2141cabf1cee18b47fcc35ab958a2a8a203000f0a7acaf95',
        'full_goal_and_novelty_not_promoted': gates['goal_achieved'] is False and gates['latest_public_source_audit']['usable_full_cohort_geometry_recovered'] is False,
        'no_external_author_message_or_invented_sender': publication['author_request_sent'] is False and read(RUN/'public_contact_verification.json')['sender_identity_supplied'] is False,
    }
    if not truth(checks):
        raise RuntimeError(json.dumps(checks))
    paths = audits + [
        'evidence/movie_geometry_audit.json',
        'docs/PUBLIC-SUPPLEMENT-AUDIT.md', 'docs/PUBLIC-SUPPLEMENT-AUDIT.html',
        'docs/DRAFT-DATA-REQUEST.md', 'docs/DRAFT-DATA-REQUEST.html',
        'docs/index.html', '../README.md', 'evidence/gates.json',
        'evidence/readiness_snapshot_verification.json', 'evidence/report_browser_verification.json',
        'evidence/raw_prose_math_scan.json', 'evidence/source_contract_space_update.json',
        'runs/source_contract_supplements_v1/public_source_audit_binding.json',
        'runs/source_contract_supplements_v1/public-source-contract.png',
        'runs/source_contract_supplements_v1/root_visual_inspection.json',
        'runs/source_contract_supplements_v1/publication_source_receipt.json',
        'runs/source_contract_supplements_v1/public_contact_verification.json',
        'runs/source_contract_supplements_v1/timing_locator_correction.json',
        'runs/source_contract_supplements_v1/movies_4_10/parent_manifest_package_independent_review.json',
    ]
    hashes = {p: sha(ROOT/p) for p in paths}
    archive = RUN/'pending_final_delivery_v2.json'
    if archive.exists():
        raise RuntimeError('Preserve pending delivery archive')
    archive.write_bytes(out.read_bytes())
    receipt = {'status': 'PASS_BOUNDED_SOURCE_AUDIT_DELIVERY', 'checks': checks,
        'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'file_hashes_relative_to_child': hashes, 'rendered_snapshot_hashes': browser['snapshot_sha256'],
        'citation_correction': 'Initial104-page chapter citation pointed to a3-snapshot receipt; existingfull26-frame audit now explicitly linked/bound. This corrected rendering supersedes the initial104-page document delivery, without changing903total or any scientific model/source/image bytes.',
        'source_sha256': sha(__file__), 'full_original_five_gate_goal_complete': False,
        'confirmed_novelty': False, 'author_request_sent': False,
        'external_contact_authorization_question_pending': True,
        'scope': 'All16 linked supplement files audited with retained source, decode, text, package, inspection and correction receipts. This has not recovered the complete matched observation contract or completed parameter inference, physical paths, continuum error or larger-region prediction. Author contact still requires explicit authorization and sender details. No public GitHub deployment.',
        'prior_publication_snapshots': 'The earlier103-page energy-minimum delivery controls/rendering/request were archived before edits. Its numerical proofs are unchanged; its documentation snapshot is historical relative to this104-page snapshot.'}
    out.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': receipt['status'], 'checks': len(checks), 'sha256': sha(out)}))


if __name__ == '__main__':
    main()
