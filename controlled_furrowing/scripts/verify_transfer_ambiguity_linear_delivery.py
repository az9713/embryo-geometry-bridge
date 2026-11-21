"""Verify the delivered LOCAL affine diagnostic, never the full five gates."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'runs/transfer_ambiguity_linear_v1'
PIN = {
    'actual/results.json': '68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f',
    'actual_independent/receipt.json': '8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072',
    'report_readback_v2/browser_readback.json': 'dfc491caecb436197116871c5bf5838b4fc7afbcc9531403c72c4e3d2ff7ff1f',
    'report_readback_v2/favicon_only_readback.json': '0b7f150b809704d97299d57d3c541fe1cfb07a24eb059a62165379b886d967ef',
    'report_build_receipt.json': '79ae626a72f45dda6a7b1a038b5900f27adf4e6ce96884101ebb55d4783bec8b',
    'REPORT.html': '891f5ab31e1ff18e678850e48db610b73b67ec6944b415aebae876fcfe499518',
    'local_affine_frontier.png': 'ad5f686dd1c0114ec92499c607413a3f3237bb46f914d71078a8bd06e8baa06c',
}


def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))


def main():
    hashes = {}
    def bound(path, digest):
        path = Path(path)
        assert sha(path) == digest, str(path)
        hashes[path.relative_to(ROOT.parent).as_posix()] = digest
    for name, digest in PIN.items():
        bound(OUT/name, digest)
    result = read(OUT/'actual/results.json')
    for name, digest in result['input_hashes'].items():
        bound(ROOT/name, digest)
    assert result['tangents_sha256'] == sha(OUT/'actual/tangents.npz')
    bound(OUT/'actual/tangents.npz', result['tangents_sha256'])
    review = read(OUT/'actual_independent/receipt.json')
    assert review['status'].startswith('PASS') and all(x is True for x in review['checks'].values())
    assert review['source_compatible_witness'] is False
    assert review['validated_equilibrium_or_PDE_claim'] is False and review['novelty_established'] is False
    browser = read(OUT/'report_readback_v2/browser_readback.json')
    assert browser['status'] == 'PASS_BROWSER_READBACK' and len(browser['checks']) == 2
    for name, digest in browser['document_hashes'].items():
        bound(ROOT.parent/name, digest)
    for row in browser['checks']:
        assert row['status'] == 'PASS' and row['page_errors'] == [] and row['request_failures'] == []
        assert row['math'] == row['rendered'] == 17 and row['mathErrors'] == ''
        assert all(x['status'] == 200 for x in row['links'])
        bound(OUT/'report_readback_v2'/row['screenshot'], row['screenshot_sha256'])
    assert len(result['rows']) == 10 and len(review['LP_endpoints']) == 20
    for row in result['rows']:
        endpoints = [row[k] for k in ('max_contrast_change', 'max_negative_contrast_change')]
        lo = sum(Q(e['primal_lower_exact']) for e in endpoints)
        hi = sum(Q(e['dual_upper_exact']) for e in endpoints)
        assert list(map(str, (lo, hi))) == row['contrast_range_width_um_exact'] and 0 <= lo <= hi
        assert all(e['exact_primal_feasible'] is True and e['exact_weak_duality'] is True for e in endpoints)
    near = next(x for x in result['rows'] if x['mode'] == 'small_prediction_change' and x['budget_exact'] == '1/1000')
    assert Q(near['contrast_range_width_um_exact'][0]) > Q('0.06042514')
    assert Q(near['contrast_range_width_um_exact'][1]) < Q('0.06042516')
    space = read(OUT/'space_update_receipt.json')
    assert space['status'] == 'PASS_SPACE_LOCAL_AFFINE_UPDATE' and all(space['checks'].values())
    assert space['sequence'] == 24 and space['blocks'] == 156
    for name in ('preparation/protocol.json', 'preparation/preflight.json',
                 'preparation/prelaunch_independent_review.json', 'space_update_receipt.json',
                 'report_readback/scientific_readback.json', 'REPORT.md', 'actual/terminal.json',
                 'actual/os_exit_observation.json'):
        bound(OUT/name, sha(OUT/name))
    delivery = {'status': 'PASS_LOCAL_AFFINE_COMPONENT_NOT_FIVE_GATE_COMPLETION',
                'source_sha256': sha(__file__), 'current_hashes': hashes,
                'exact_affine_optimizations': 20, 'small_constraints': 68,
                'browser_views': 2, 'equations_per_view': 17,
                'near_silent_affine_range_um_exact': near['contrast_range_width_um_exact'],
                'nonlinear_witnesses_certified': False, 'source_compatibility_established': False,
                'full_five_gate_goal_complete': False, 'novelty_established': False,
                'space_sequence_at_readback': space['sequence']}
    path = OUT/'delivery_receipt.json'
    with path.open('x', encoding='utf-8') as f:
        json.dump(delivery, f, indent=2)
        f.write('\n')
    print(json.dumps({'status': delivery['status'], 'current_hashes': len(hashes),
                      'delivery_sha256': sha(path)}))


if __name__ == '__main__':
    main()
