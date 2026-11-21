"""Publish the energy-minimum checkpoint only after independent replay passes."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'runs/cohort_energy_stability_v1'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text(encoding='utf-8'))


def main():
    review = read(ROOT/'evidence/cohort_energy_stability_independent_review.json')
    if not review['status'].startswith('PASS') or not all(v is True for v in review['checks'].values()):
        raise RuntimeError('All independent replay checks must pass before publication')
    result = read(RUN/'results.json')
    if result['status'] != 'CERTIFIED_ALL17_UNIFORM_LOCAL_ENERGY_MINIMA' or len(result['records']) != 17:
        raise RuntimeError('All17 producer results required')
    # The prior current snapshot and its exact control bytes were archived before
    # any subsequent result publication. Its mathematical packets stay intact.
    old = ROOT/'runs/checkpoint_history/shared_cell_stationary_delivery/archive_manifest.json'
    if not old.exists():
        raise RuntimeError('Prior delivery archive required')
    draft = RUN/'publication-draft.md'
    if not draft.exists():
        raise RuntimeError('Reviewed publication draft required')
    target = ROOT/'docs/COHORT-ENERGY-STABILITY.md'
    files = [ROOT.parent/'README.md', ROOT/'evidence/gates.json', ROOT/'scripts/build_report.py', ROOT.parent/'scripts/build_dossier.py']
    before = {str(p.relative_to(ROOT.parent)): sha(p) for p in files}
    if target.exists():
        raise RuntimeError('Preserve an already published checkpoint')
    text = files[0].read_text(encoding='utf-8')
    anchor = '## Findings\n\n'
    addition = ('* **The 17 selected equilibria are now certified local energy minima:** '
        '[proof, visual and independent receipts](controlled_furrowing/docs/COHORT-ENERGY-STABILITY.html). '
        'Exact-binary matrix positivity certificates, combined with the existing uniform root bounds, '
        'prove positive reduced Hessians throughout every selected root tube on the same tiny parameter cell. '
        'This establishes local energetic stability in the fixed finite-quadrature model. '
        'It supplies no common loading trajectory, continuum stability, improved fit or confirmed novelty.\n\n')
    if text.count(anchor) != 1:
        raise RuntimeError('Unique README findings anchor required')
    updated_readme = text.replace(anchor, anchor+addition)
    gates = read(files[1])
    gates['snapshot'] = ('All17 selected local stationary-root families and68 conditional static outputs '
        'have uniform certificates on ONE tiny shared parameter cell. A separate exact-preconditioner '
        'positivity calculation and independent replay now additionally prove strict local energy minima '
        'throughout all17 root tubes. Source-matched time/geometry/discrepancy, full physical identified '
        'set, common stable loading paths, continuum error and revised larger-region validation remain incomplete.')
    gates['latest_cohort_energy_stability'] = {
        'status': result['status'],
        'producer_file': 'runs/cohort_energy_stability_v1/results.json',
        'producer_sha256': sha(RUN/'results.json'),
        'independent_review_file': 'evidence/cohort_energy_stability_independent_review.json',
        'independent_review_sha256': sha(ROOT/'evidence/cohort_energy_stability_independent_review.json'),
        'selected_records': 17,
        'parameter_cell_sha256': result['parameter_cell_sha256'],
        'whole_selected_root_tubes': True,
        'common_physical_loading_path': False,
        'global_energy_minimum': False,
        'continuum_stability': False,
        'biological_fit_or_novelty_established': False,
    }
    g3 = next(g for g in gates['gates'] if g['gate'] == 3)
    g3['scientific_outcome'] += (' Subsequent independent exact-preconditioner positivity certificates '
        'now establish strict local minima of the actual17 selected cohort root families on the same C. '
        'This later component result supersedes the earlier statement that actual-cohort energy stability '
        'had not yet been certified; it does not transfer the separate old weak-load theorem or prove a common physical path.')
    g3['evidence'] += ['docs/COHORT-ENERGY-STABILITY.md', 'evidence/cohort_energy_stability_independent_review.json']
    gates['goal_achieved'] = False
    updated_gates = json.dumps(gates, indent=2, ensure_ascii=False)+'\n'
    builder = files[2].read_text(encoding='utf-8')
    anchor = "    body+='<h2>Latest: all 68 static outputs on one shared parameter cell</h2>"
    if builder.count(anchor) != 1:
        raise RuntimeError('Unique child landing checkpoint anchor required')
    addition = ("    body+='<h2>New: all 17 selected states are local energy minima</h2><p>"
        "<a href=\"COHORT-ENERGY-STABILITY.html\">Proof, visual and independent matrix certificates</a>. "
        "Validated positivity of the exact numerical preconditioners, together with the existing uniform "
        "root bounds, proves positive reduced energy Hessians throughout all seventeen selected root tubes. "
        "These are strict local minima of the fixed finite-quadrature energy on the same tiny parameter cell. "
        "A common loading history, continuum stability, source-matched biological fit and novelty remain unresolved.</p>'\n")
    updated_builder = builder.replace(anchor, addition+anchor.replace('Latest:', 'Previous checkpoint:'))
    parent_builder = files[3].read_text(encoding='utf-8')
    parent_anchor = "content+='<h2>Latest: a verified common-cell cohort, with scientific limits</h2>"
    if parent_builder.count(parent_anchor) != 1:
        raise RuntimeError('Unique parent checkpoint anchor required')
    parent_addition = ("content+='<h2>New: local energy minima for the selected cohort</h2><p>"
        "All seventeen selected controlled-furrowing root families now have independently checked "
        "positive energy-Hessian certificates throughout their tubes on the same tiny shared parameter cell. "
        "The finite-dimensional states are strict local energy minima; common physical histories, "
        "continuum stability, biological fit and novelty remain unproved. "
        "<a href=\"../controlled_furrowing/docs/COHORT-ENERGY-STABILITY.html\">Proof, figure and retained receipts</a>.</p>'\n")
    updated_parent_builder = parent_builder.replace(parent_anchor, parent_addition+parent_anchor.replace('Latest:', 'Previous checkpoint:'))
    out = RUN/'publication_source_receipt.json'
    if out.exists():
        raise RuntimeError('Preserve publication receipt')
    # Validate all sources, anchors and destinations before the first edit.
    target.write_bytes(draft.read_bytes())
    for p, content in zip(files, (updated_readme, updated_gates, updated_builder, updated_parent_builder)):
        p.write_text(content, encoding='utf-8')
    receipt = {
        'status': 'PUBLISHED_LOCAL_CHECKPOINT_SOURCE_BROWSER_PENDING',
        'before_control_hashes_relative_to_repository': before,
        'after_control_hashes_relative_to_repository': {str(p.relative_to(ROOT.parent)): sha(p) for p in files},
        'draft_sha256': sha(draft),
        'document_sha256': sha(target),
        'producer_sha256': sha(RUN/'results.json'),
        'independent_review_sha256': sha(ROOT/'evidence/cohort_energy_stability_independent_review.json'),
        'source_sha256': sha(__file__),
        'full_original_five_gate_goal_complete': False,
        'public_GitHub_publication': False,
        'scope': 'Publication source edits only. Rebuild and fresh actual browser verification required separately.',
    }
    out.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': receipt['status'], 'receipt_sha256': sha(out)}))


if __name__ == '__main__':
    main()
