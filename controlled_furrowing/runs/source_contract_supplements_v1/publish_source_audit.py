"""Publish the complete bounded supplement audit and an explicitly unsent request."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'runs/source_contract_supplements_v1'
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    receipts = ['evidence/public_supplement_asset_manifest.json',
                'evidence/source_contract_workbook_package_audit.json',
                'evidence/source_contract_small_movies_audit.json',
                'evidence/source_contract_other_movies_audit.json',
                'evidence/source_contract_peer_review_audit.json']
    for name in receipts:
        d = read(ROOT/name)
        if not d.get('checks') or not all(v is True for v in d['checks'].values()):
            raise RuntimeError('Passing source audit required: '+name)
    asset = read(ROOT/receipts[0])
    if len(asset['assets']) != 16 or asset['total_bytes'] != 111690608:
        raise RuntimeError('Expected complete supplement package')
    draft = RUN/'public-source-audit-draft.md'
    request = RUN/'data-request-v2.md'
    if not draft.exists() or not request.exists() or not (RUN/'prepublication_snapshot/manifest.json').exists():
        raise RuntimeError('Prepared drafts and prior snapshot required')
    target = ROOT/'docs/PUBLIC-SUPPLEMENT-AUDIT.md'
    out = RUN/'publication_source_receipt.json'
    if target.exists() or out.exists():
        raise RuntimeError('Preserve published source checkpoint')
    files = [ROOT.parent/'README.md', ROOT/'evidence/gates.json', ROOT/'scripts/build_report.py', ROOT/'docs/DRAFT-DATA-REQUEST.md']
    before = {str(p.relative_to(ROOT.parent)): sha(p) for p in files}
    text = files[0].read_text(encoding='utf-8')
    anchor = '## Findings\n\n'
    if text.count(anchor) != 1:
        raise RuntimeError('Unique README anchor required')
    addition = ('* **All sixteen linked supplementary files are now accounted for:** '
        '[supplement audit and original diagram](controlled_furrowing/docs/PUBLIC-SUPPLEMENT-AUDIT.html). '
        'All16 linked assets (111,690,608 bytes) match official checksums; all903 movie frames were decoded, '
        'the36-page peer-review file was read, and workbook package structures were audited. '
        'No complete matched cohort of boundary tracks, reference curves and acquisition metadata was recovered. '
        'A [specific processed-data request](controlled_furrowing/docs/DRAFT-DATA-REQUEST.html) is prepared and unsent. '
        'This adds source evidence, not a new fit, uncertainty allowance or novelty claim.\n\n')
    new_readme = text.replace(anchor, anchor+addition)
    gates = read(files[1])
    g2 = next(g for g in gates['gates'] if g['gate'] == 2)
    g2['execution'] = 'ALL16 LINKED PUBLIC SUPPLEMENTS ACCOUNTED FOR; PROCESSED COHORT GEOMETRY CONTRACT INCOMPLETE'
    g2['scientific_outcome'] = ('17 calibration entries,9 scalar evaluation entries and31 larger bending labels remain available. '
        'All16 public supplement assets checksum-verified;10movies/903frames decoded and36peer-review pages read. '
        'Protocol facts are corroborated, but matched boundary tracks, embryo-specific reference curves, '
        'ROI/acquisition mappings and defensible measurement uncertainty have not been recovered from this inspected package.')
    g2['evidence'] += ['docs/PUBLIC-SUPPLEMENT-AUDIT.md', *receipts]
    gates['latest_public_source_audit'] = {
        'supplement_assets': 16, 'bytes': 111690608, 'decoded_movie_frames': 903,
        'peer_review_pages_read': 36, 'usable_full_cohort_geometry_recovered': False,
        'author_request_sent': False, 'external_contact_authorized': False,
        'source_receipt_hashes': {p: sha(ROOT/p) for p in receipts},
        'scope': 'Complete inspected article-supplement coverage, not all possible public repositories or recovered primary microscopy.'}
    gates['goal_achieved'] = False
    gates['next_executable_inference_step'] = ('Resolve source-matched reference/ROI/time/component semantics and obtain compact processed cohort geometry '
        'or authoritative clarification; a concrete request is prepared but requires explicit external-contact authorization and sender identity. '
        'The16linked supplements are now fully inventoried; nominal sampling/selection criteria are not anatomical error bounds. '
        'Full physical identified-set propagation, common loading paths, integral/continuum error and revised frozen larger-region predictions remain incomplete.')
    builder = files[2].read_text(encoding='utf-8')
    anchor = "    body+='<h2>New: all 17 selected states are local energy minima</h2>"
    if builder.count(anchor) != 1:
        raise RuntimeError('Unique child landing anchor required')
    overview = ("    body+='<h2>New: the public-source audit and the remaining data dependency</h2><p>"
        "All16 linked supplementary assets match official checksums. Ten movies contain903 decoded frames; "
        "the36-page peer-review file and workbook package structures are audited. These sources constrain "
        "the protocol but do not supply the matched cohort geometry/reference/ROI and uncertainty contract. "
        "<a href=\"PUBLIC-SUPPLEMENT-AUDIT.html\">Audit, original diagram and receipts</a>; "
        "<a href=\"DRAFT-DATA-REQUEST.html\">prepared processed-data request, not sent</a>. "
        "The existing local energy-minimum proof and residual are unchanged; the five-gate study remains incomplete.</p>'\n")
    new_builder = builder.replace(anchor, overview+anchor.replace('New:', 'Mechanics checkpoint:'))
    # All data and edit locations are checked before the first source edit.
    target.write_bytes(draft.read_bytes())
    files[0].write_text(new_readme, encoding='utf-8')
    files[1].write_text(json.dumps(gates, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    files[2].write_text(new_builder, encoding='utf-8')
    files[3].write_bytes(request.read_bytes())
    packet = {'status': 'PUBLISHED_SOURCE_AUDIT_BROWSER_PENDING',
        'before_controls_relative_to_repository': before,
        'after_controls_relative_to_repository': {str(p.relative_to(ROOT.parent)): sha(p) for p in files},
        'published_document_sha256': sha(target), 'prepared_request_sha256': sha(files[3]),
        'source_audit_hashes': {p: sha(ROOT/p) for p in receipts}, 'publisher_source_sha256': sha(__file__),
        'full_five_gate_goal_complete': False, 'author_request_sent': False, 'confirmed_novelty': False,
        'scope': 'Local documentation sources. Fresh HTML/math/assets/layout verification remains separate. No author contact or public GitHub publication.'}
    out.write_text(json.dumps(packet, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': packet['status'], 'sha256': sha(out)}))


if __name__ == '__main__':
    main()
