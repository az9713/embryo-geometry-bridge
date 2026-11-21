"""Root-side independent replay of the same closed K-L dyadic body/log cover.

Uses production exact geometry primitives, independently of segments.py's
separate P2 polynomial implementation. No trace, force, Hessian or root claim.
"""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
from math import ceil
import importlib.util
import json
import os
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT = ROOT / 'runs/exact_loading_tube_pilot_preparation_v1'

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def read(p, h):
    assert digest(p) == h, str(p)
    return json.loads(p.read_bytes())

def save(p, d):
    with p.open('x', encoding='utf8', newline='\n') as f:
        json.dump(d, f, indent=2)
        f.write('\n')
    return digest(p)

def main():
    protocol = read(PILOT/'protocol.json', 'fe64074d1a169a7672039c73dc07ca35c728cc1c299cf2a28025e0a0e7a97e43')
    plan = read(PILOT/'plan.json', 'c3b1870e5ad815bc5df5ec94d4731b0228e529ba2e0651f847102b8e0ecaa277')
    envelope = read(PILOT/'envelope.json', 'f50c04b1802176abf4f5494898c8128018d4db8f8d0bb410af7774e5a3fff925')
    native = read(HERE/'segment_actual/receipt.json', '97953457778f39011d96cc071070d2d5d342eb362f7e8e904d28d7947fd6e214')
    assert digest(HERE/'segments.py') == '9fe45835838d6da98569f32983a304a4569e86663a2794668e99c0f62f12ff69'
    observation = json.loads((HERE/'segment_actual/os_exit_observation.json').read_bytes())
    assert observation['receipt_sha256'] == digest(HERE/'segment_actual/receipt.json')
    assert observation['unified_exec_chunk_id'] == 'f19a13'
    assert type(observation['observed_exit_code']) is int and observation['observed_exit_code'] == 0
    terminal = read(HERE/'segment_actual/terminal.json', observation['terminal_sha256'])
    assert terminal['status'] == 'FINISHED' and terminal['exit_code'] == 0
    assert terminal['pid'] == observation['pid'] == native['pid']
    assert terminal['receipt_sha256'] == observation['receipt_sha256']
    for name, h in protocol['geometry_source_tree_sha256'].items():
        assert digest(ROOT/name) == h
    sys.path.insert(0, str(ROOT/'scripts'))
    import root_body_geometry as G
    assert Path(G.__file__).resolve() == ROOT/'scripts/root_body_geometry.py'
    assert digest(Path(G.__file__)) == protocol['geometry_source_tree_sha256']['scripts/root_body_geometry.py']
    anchors = [read(ROOT/name, h) for name, h in zip(plan['anchor_box_files'], plan['anchor_box_sha256'], strict=True)]
    first, last = [[Q(v) for row in a['current_centers_exact'] for v in row] for a in anchors]
    reference = [[Q(v) for v in row] for row in envelope['reference_nodes_exact']]
    free = set(envelope['free'])
    radius = Q(plan['fiber_radius_exact'])
    assert radius == Q(1,4096) == Q(native['constant_corrector_radius_exact'])
    doses = list(map(Q, plan['dose_interval_exact']))
    assert list(map(Q, native['exact_original_closed_interval'])) == doses
    out = HERE/'root_segment_actual_readback'
    out.mkdir(exist_ok=False)
    stamp = {'pid': os.getpid(), 'source_sha256': digest(Path(__file__)),
             'producer_receipt_sha256': observation['receipt_sha256']}
    save(out/'started.json', {'status':'RUNNING', **stamp})
    end = {'status':'FAILED', 'exit_code':1, **stamp}
    try:
        levels = []
        for level in native['levels']:
            count = level['slabs']
            assert type(count) is int and count in (2,4,8,16)
            rows = []
            assert len(level['checks']) == count
            for j, prior in enumerate(level['checks']):
                a, b = Q(j,count), Q(j+1,count)
                center, widths = [], []
                for i, (x,y) in enumerate(zip(first,last,strict=True)):
                    p, q = x+a*(y-x), x+b*(y-x)
                    if i not in free:
                        assert p == q == reference[i//2][i%2]
                    center.append((p+q)/2)
                    widths.append(abs(q-p)/2 + (radius if i in free else 0))
                pair = lambda values: [[values[i],values[i+1]] for i in range(0,650,2)]
                box = G.certify_box(reference, pair(center), envelope['triangles'], pair(widths))
                bh = save(out/f'body_{count}_{j}.json', box)
                badJ, badM, badRho, rhos, Js = [], [], [], [], []
                for i, t in enumerate(box['triangles']):
                    if not t['positive_J_for_entire_box']: badJ.append(i)
                    if not t['strong_monotonicity_for_entire_box']: badM.append(i)
                    if t['J_bound'] is None:
                        badRho.append(i)
                        continue
                    lo, hi = Q(t['J_bound']['lower_exact']), Q(t['J_bound']['upper_exact'])
                    Js.append((lo,hi))
                    rho = Q(ceil(max(abs(lo-1),abs(hi-1))*2**48),2**48)
                    if lo <= Q(1,10**6) or rho >= 1: badRho.append(i)
                    else: rhos.append(rho)
                floor = min(center[i]-widths[i]-reference[i//2][0]
                            -Q(float(-(1-1e-7)*float(reference[i//2][0]))) for i in free if i%2==0)
                passed = not (badJ or badM or badRho) and floor > 0 and box['axis_preserved_for_entire_box']
                assert prior['slot'] == j
                assert list(map(Q, prior['exact_dose_subinterval'])) == [doses[0]+v*(doses[1]-doses[0]) for v in (a,b)]
                assert prior['passes_body_J_monotonicity_log_rho_and_radial_bound'] is bool(passed)
                assert badJ == prior['insufficient_J_slots'] and badM == prior['insufficient_monotonicity_slots']
                assert badRho == prior['insufficient_bulk_log_rho_slots']
                assert Q(prior['radial_nodal_margin_exact']) == floor
                if len(Js) == 144:
                    assert list(map(Q, prior['J_enclosure_exact'])) == [min(v[0] for v in Js),max(v[1] for v in Js)]
                if len(rhos) == 144: assert Q(prior['maximum_bulk_log_rho_exact']) == max(rhos)
                if not badM: assert Q(prior['minimum_monotonicity_exact']) == Q(box['strong_monotonicity_lower_exact'])
                rows.append({'slot':j,'body_sha256':bh,'passes':bool(passed),
                             'J_bound':box['J_bound'],'monotonicity_lower_exact':box['strong_monotonicity_lower_exact'],
                             'max_log_rho_exact':str(max(rhos)) if len(rhos)==144 else None})
            levels.append({'slabs':count,'rows':rows,'all_pass':all(v['passes'] for v in rows)})
        chosen = next((v['slabs'] for v in levels if v['all_pass']), None)
        assert chosen == native['first_complete_tested_dyadic_level'] == 4
        h = save(out/'receipt.json', {'status':'PASS_INDEPENDENT_FULL_SPAN_KL_FOUR_SLAB_BODY_LOG_COVERAGE',
                 **stamp,'levels':levels,'first_complete_tested_dyadic_level':chosen,
                 'exact_original_closed_interval':list(map(str,doses)),'constant_radius_exact':str(radius),
                 'surface_cache_force_H_or_loading_root_proved':False})
        end.update(status='FINISHED',exit_code=0,receipt_sha256=h)
        print(json.dumps({'receipt_sha256':h,'first_complete_tested_dyadic_level':chosen}))
    finally:
        save(out/'terminal.json',end)

if __name__ == '__main__':
    main()
