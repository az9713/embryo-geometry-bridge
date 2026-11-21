"""Outcome-free, LOCAL AFFINE transfer-ambiguity diagnostic; never a physical fit.

Exact rational graph derivatives and LP weak-duality certificates enclose only
the explicitly rounded affine model. Mechanical tangents are floating and
conditional on an exact regular root. No new equilibrium or continuum claim.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import linprog
import torch

from annulus_observer import _quadratic_controls, _root_bracket
from quadratic_endpoint_observer import _state_controls, polynomial_interval, interval
from quadratic_mean_derivatives import _boundary_node_ids, BASIS
from material_tangent_inference import mechanical_tangent, check_sources, COLUMNS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/transfer_ambiguity_linear_v1'
FIELD = 'runs/transfer_pair_pilot_v1/jobs/baseline_large/checkpoint_03/field.json'
TANGENT = 'runs/material_cohort_descent_v1/results.json'
SMALL_PROTOCOL = 'runs/material_cohort_descent_v1/protocol.json'
PILOT_PROTOCOL = 'runs/transfer_pair_pilot_v1/preparation/protocol.json'
ACTIVE = (1, 2, 3, 4, 5)
LENGTH = Q('33.5')
SOURCES = ('transfer_ambiguity_linear.py', 'material_tangent_inference.py',
           'material_cohort_descent.py', 'quadratic_mean_derivatives.py',
           'quadratic_endpoint_observer.py', 'annulus_observer.py')
INPUTS = (FIELD, TANGENT, SMALL_PROTOCOL, PILOT_PROTOCOL,
          'runs/transfer_pair_pilot_v1/independent_actual/actual_review.json',
          'evidence/material_cohort_descent_independent_review.json', 'PLAN.md')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, doc):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        json.dump(doc, f, indent=2, allow_nan=False)
        f.write('\n')
    return sha(path)


def stamp():
    return datetime.now(timezone.utc).isoformat()


def frac(x):
    return Q.from_float(float(x))


def pairs(x):
    return [str(x.lo), str(x.hi)]


def add(a, b):
    return interval(a.lo + b.lo, a.hi + b.hi)


def multiply(a, b):
    v = [x*y for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
    return interval(min(v), max(v))


def divide(a, b):
    require(b.lo > 0, 'strictly positive radial derivative needed')
    return multiply(a, interval(1/b.hi, 1/b.lo))


def contrast_curve_gradient(controls, shoulder, steps=100):
    """Physical-node gradient of z(0)-z(a), with the axis radius held fixed.

    At R(t)=a, dz(a)=sum N_j dz_j - (z_t/R_t) sum N_j dr_j.
    Require one interior edge, monotone R, and unchanged axis topology.
    """
    require(type(steps) is int and steps >= 1, 'positive inverse bracket steps')
    shoulder = Q(shoulder)
    curves = _quadratic_controls(controls)
    require(curves[0][0][0][0] == 0 and shoulder > 0, 'fixed symmetry axis and positive shoulder')
    chosen = [i for i, (q, *_) in enumerate(curves) if q[0][0] < shoulder < q[2][0]]
    require(len(chosen) == 1, 'shoulder must lie strictly inside one edge, not a routing knot')
    k = chosen[0]
    q, c0, c1, c2 = curves[k]
    t = interval(*_root_bracket(c0[0], c1[0], c2[0], shoulder, steps))
    rt = polynomial_interval([c1[0], 2*c2[0]], t)
    zt = polynomial_interval([c1[1], 2*c2[1]], t)
    slope = divide(zt, rt)
    gradient = [[[interval(0), interval(0)] for _ in range(3)] for _ in curves]
    gradient[0][0][1] = interval(1)
    for j, basis in enumerate(BASIS):
        n = polynomial_interval(basis, t)
        gradient[k][j][0] = add(gradient[k][j][0], multiply(slope, n))
        gradient[k][j][1] = add(gradient[k][j][1], interval(-n.hi, -n.lo))
    z = polynomial_interval([c0[1], c1[1], c2[1]], t)
    value = interval(curves[0][0][0][1]-z.hi, curves[0][0][0][1]-z.lo)
    return gradient, {'edge_index': k, 't_exact': pairs(t), 'radial_derivative_exact': pairs(rt),
                      'slope_exact': pairs(slope), 'Q_um_exact': pairs(value),
                      'endpoint_margin_um_exact': str(min(shoulder-q[0][0], q[2][0]-shoulder))}


def state_gradient(state, shoulder):
    controls, _ = _state_controls(state, LENGTH)
    g, detail = contrast_curve_gradient(controls['top'], shoulder)
    ids = _boundary_node_ids(state)['top']
    require(len(ids) == len(g), 'same P2 trace-to-node map')
    out = [[interval(0), interval(0)] for _ in state['current_nodes']]
    for row, local in zip(ids, g):
        for node, values in zip(row, local):
            for a in (0, 1):
                out[node][a] = add(out[node][a], values[a])
    return out, detail


def max_affine(c, J, residual, delta, bounds):
    """Exact primal/weak-dual bounds for max c*d over |r+Jd|<=delta.

    The zero step must be feasible. Round solver d to dyadics then shrink it
    exactly toward zero. Any nonnegative multiplier gives a valid dual bound
    using the exact support function of the finite trust box.
    """
    c, residual, delta = list(map(Q, c)), list(map(Q, residual)), Q(delta)
    J = [list(map(Q, row)) for row in J]
    bounds = [tuple(map(Q, row)) for row in bounds]
    n = len(c)
    require(n > 0 and len(J) == len(residual) > 0 and all(len(row) == n for row in J)
            and len(bounds) == n, 'complete matching affine dimensions')
    require(delta >= max(map(abs, residual)) and all(lo <= 0 <= hi for lo, hi in bounds),
            'nonempty trust box and feasible zero anchor required')
    A = J + [[-x for x in row] for row in J]
    b = [delta-r for r in residual]+[delta+r for r in residual]
    result = linprog(-np.array(list(map(float, c))), A_ub=np.array(A, dtype=float),
                     b_ub=np.array(b, dtype=float), bounds=[tuple(map(float, x)) for x in bounds], method='highs')
    require(result.success, 'affine LP solver failed: '+result.message)
    d = [min(max(frac(x), lo), hi) for x, (lo, hi) in zip(result.x, bounds)]
    alpha = Q(1)
    for row, cap in zip(A, b):
        v = sum((a*x for a, x in zip(row, d)), Q(0))
        if v > cap:
            require(v > 0 and cap >= 0, 'radial feasibility repair needs feasible anchor')
            alpha = min(alpha, cap/v)
    d = [alpha*x for x in d]
    primal = sum((a*x for a, x in zip(c, d)), Q(0))
    lam = [max(Q(0), frac(-v)) for v in result.ineqlin.marginals]
    reduced = [cj-sum((l*row[j] for l, row in zip(lam, A)), Q(0)) for j, cj in enumerate(c)]
    dual = sum((l*v for l, v in zip(lam, b)), Q(0)) + sum(max(v*lo, v*hi) for v, (lo, hi) in zip(reduced, bounds))
    require(all(sum((a*x for a, x in zip(row, d)), Q(0)) <= cap for row, cap in zip(A, b))
            and all(lo <= x <= hi for x, (lo, hi) in zip(d, bounds)) and Q(0) <= primal <= dual,
            'exact primal feasibility / weak duality failed')
    return {'primal_lower_exact': str(primal), 'dual_upper_exact': str(dual),
            'duality_gap_exact': str(dual-primal), 'step_exact': list(map(str, d)),
            'multipliers_exact': list(map(str, lam)), 'repair_alpha_exact': str(alpha),
            'exact_primal_feasible': True, 'exact_weak_duality': True}


def synthetics():
    checks = {}
    controls = [[(Q(0), Q(0)), (Q(1,2), Q(1,4)), (Q(1), Q(1))]]
    g, detail = contrast_curve_gradient(controls, Q(1,4))
    # R=t, Z=t^2; at t=1/4, N=(3/8,3/4,-1/8), Z'/R'=1/2.
    expected = [(Q(3,16), Q(5,8)), (Q(3,8), Q(-3,4)), (Q(-1,16), Q(1,8))]
    checks['analytic_moving_radius_gradient'] = all(a.contains(x) for values, targets in zip(g[0], expected) for a,x in zip(values,targets))
    checks['analytic_contrast'] = Q(detail['Q_um_exact'][0]) <= -Q(1,16) <= Q(detail['Q_um_exact'][1])
    lp = max_affine(['1'], [['1']], ['0'], '1/2', [('-1','1')])
    checks['analytic_LP_half'] = lp['primal_lower_exact'] == lp['dual_upper_exact'] == '1/2'
    one = max_affine(['-1'], [['0']], ['0'], '0', [('0','1')])
    checks['one_sided_physical_coordinate'] = one['primal_lower_exact'] == one['dual_upper_exact'] == '0'
    for name, fn in (
        ('routing_knot_rejected', lambda: contrast_curve_gradient(controls, 1)),
        ('infeasible_anchor_rejected', lambda: max_affine(['1'], [['1']], ['2'], '1', [('-1','1')])),
        ('empty_dimension_rejected', lambda: max_affine([], [], [], '0', [])),
    ):
        try:
            fn()
            checks[name] = False
        except ValueError:
            checks[name] = True
    require(all(checks.values()), 'synthetic verification failed')
    return checks


def inputs():
    hashes = {f'scripts/{name}': sha(ROOT/'scripts'/name) for name in set(SOURCES) | set(check_sources())}
    hashes.update({name: sha(ROOT/name) for name in INPUTS})
    small = read(ROOT/SMALL_PROTOCOL)
    pilot = read(ROOT/PILOT_PROTOCOL)
    state = read(ROOT/FIELD)
    saved = read(ROOT/TANGENT)
    candidate = next(x for x in saved['candidates'] if x['variant'] == 'current_disk__spatial_max')
    return hashes, small, pilot, state, candidate


def prepare():
    checks = synthetics()
    hashes, small, pilot, state, candidate = inputs()
    # Schema guards stay explicit; archive paths are not inferred from names.
    require(len(candidate['all68_constraint_order']) == 68, 'all 68 rows needed')
    affine = candidate['affine_LP']
    require(len(affine['affine_residual_intercepts_exact']) == 68
            and len(affine['affine_derivative_coefficients_exact']) == 68
            and all(len(row) == 5 for row in affine['affine_derivative_coefficients_exact'])
            and small['baseline_variant'] == 'current_disk__spatial_max', 'same rounded small baseline')
    require(state['converged'] and state['parameters']['mu'] == 1
            and state['parameters']['bulk'] == 5000, 'same normalized mechanical model')
    baseline = pilot['calibration_candidates']['baseline']
    require(Q(baseline['gain_exact']) == Q(small['baseline_gain_exact']), 'identical exact small/large baseline gain')
    require({k:v for k,v in baseline['parameters'].items() if k != 'patch_radius'}
            == {k:v for k,v in state['parameters'].items() if k != 'patch_radius'}, 'same materials and geometry except activation radius')
    case = next(x for x in pilot['cases'] if x['case_id'] == 'baseline_large')
    require(case['parameters'] == state['parameters'], 'same frozen large-case parameters')
    grad, detail = state_gradient(state, Q('24.94698765'))
    checks.update({'all68_same_small_baseline_rows': True, 'stored_shoulder_strictly_inside_one_edge': True,
                   'no_tangent_Hessian_or_new_equilibrium_in_preparation': True})
    p = {'status': 'FROZEN_LOCAL_AFFINE_DIAGNOSTIC_REVIEW_PENDING', 'source_sha256': sha(__file__),
         'input_hashes': hashes, 'columns': [COLUMNS[i] for i in ACTIVE], 'large_field': FIELD,
         'small_affine_model': TANGENT, 'shoulder_um_exact': '498939753/20000000',
         'source_discrepancy_offsets_um_exact': ['0','1/1000','1/100','1/10','1/2','1'],
         'small_prediction_change_budgets_um_exact': ['1/10000','1/1000','1/100','1/10'],
         'trust_box_exact': affine['trust_box_exact'], 'one_sided_apical_tension': True,
         'arithmetic': 'Exact rational graph derivative; floating mechanics; rounded five-column objective; exact LP primal/weak dual bounds.',
         'scope': 'Discovery-stage LOCAL AFFINE surrogate only. No nonlinear feasible witnesses, exact integral roots, PDE error, source compatibility or novelty.',
         'heldout_response_tables_read': False, 'new_equilibrium_solves': False,
         'preparation_shoulder_details': detail, 'utc': stamp()}
    psha = save(OUT/'preparation/protocol.json', p)
    fsha = save(OUT/'preparation/preflight.json', {'status': 'PASS_PREPARATION_ONLY', 'checks': checks,
               'protocol_sha256': psha, 'source_sha256': sha(__file__)})
    print(json.dumps({'protocol_sha256': psha, 'preflight_sha256': fsha, 'source_sha256': sha(__file__)}), flush=True)


def actual(review_sha):
    ppath = OUT/'preparation/protocol.json'
    fpath = OUT/'preparation/preflight.json'
    review_path = OUT/'preparation/prelaunch_independent_review.json'
    p, pre, review = read(ppath), read(fpath), read(review_path)
    require(sha(review_path) == review_sha and review['status'] == 'PASS'
            and all(review['checks'].values()) and review['source_sha256'] == sha(__file__) == p['source_sha256']
            and review['protocol_sha256'] == sha(ppath) == pre['protocol_sha256']
            and review['preflight_sha256'] == sha(fpath), 'matching exact independent prelaunch required')
    hashes, small, pilot, state, candidate = inputs()
    require(hashes == p['input_hashes'], 'input changed after review')
    directory = OUT/'actual'
    directory.mkdir(exist_ok=False)
    start = monotonic()
    save(directory/'started.json', {'pid': os.getpid(), 'utc': stamp(), 'review_sha256': review_sha})
    try:
        torch.set_num_threads(1)
        model, u, H, B, S, mechanics = mechanical_tangent(state)
        require(all(mechanics['checks'].values()), 'floating mechanical tangent checks failed')
        gradient, detail = state_gradient(state, Q(p['shoulder_um_exact']))
        physical = [v for row in gradient for v in row]
        require(len(physical) == 2*len(state['current_nodes']), 'full physical node gradient coverage')
        projected = [multiply(physical[int(i)], interval(LENGTH)) for i in model.free]
        c_intervals = []
        for column in ACTIVE:
            value = interval(0)
            for a, v in zip(projected, S[:,column]):
                value = add(value, multiply(a, interval(frac(v))))
            c_intervals.append(value)
        c = [frac(float((v.lo+v.hi)/2)) for v in c_intervals]
        affine = candidate['affine_LP']
        J = affine['affine_derivative_coefficients_exact']
        residual = affine['affine_residual_intercepts_exact']
        bounds = affine['trust_box_exact']
        base = max(abs(Q(x)) for x in residual)
        rows = []
        for mode, budgets in (('source_residual_allowance', p['source_discrepancy_offsets_um_exact']),
                              ('small_prediction_change', p['small_prediction_change_budgets_um_exact'])):
            for budget in budgets:
                r = residual if mode == 'source_residual_allowance' else ['0']*68
                delta = base+Q(budget) if mode == 'source_residual_allowance' else Q(budget)
                upper = max_affine(c, J, r, delta, bounds)
                lower = max_affine([-x for x in c], J, r, delta, bounds)
                width = [Q(upper['primal_lower_exact'])+Q(lower['primal_lower_exact']),
                         Q(upper['dual_upper_exact'])+Q(lower['dual_upper_exact'])]
                rows.append({'mode': mode, 'budget_exact': budget, 'allowance_um_exact': str(delta),
                             'max_contrast_change': upper, 'max_negative_contrast_change': lower,
                             'contrast_range_width_um_exact': list(map(str,width)),
                             'contrast_range_width_um_display': list(map(float,width))})
        result = {'status': 'COMPLETE_LOCAL_AFFINE_DIAGNOSTIC_NOT_NONLINEAR_WITNESSES',
                  'source_sha256': sha(__file__), 'protocol_sha256': sha(ppath), 'prelaunch_review_sha256': review_sha,
                  'input_hashes': hashes, 'mechanics': mechanics, 'shoulder': detail,
                  'columns': p['columns'], 'contrast_gradient_intervals_exact': [pairs(v) for v in c_intervals],
                  'rounded_affine_contrast_gradient_exact': list(map(str,c)),
                  'rounded_affine_contrast_gradient_um_display': list(map(float,c)),
                  'source_baseline_residual_um_exact': str(base), 'rows': rows,
                  'small_J_nonsmooth_warning': 'Archived maxima columns are central graph differences; no certified differentiability or nonlinear remainder.',
                  'scope': p['scope'], 'heldout_response_tables_read': False, 'new_equilibrium_solves': False,
                  'seconds': monotonic()-start}
        np.savez_compressed(directory/'tangents.npz', S=S, B=B, u=u, free=model.free,
                            H_data=H.data, H_indices=H.indices, H_indptr=H.indptr, H_shape=H.shape)
        result['tangents_sha256'] = sha(directory/'tangents.npz')
        digest = save(directory/'results.json', result)
        save(directory/'terminal.json', {'status': 'FINISHED', 'pid': os.getpid(), 'exit_code': 0,
                                        'results_sha256': digest, 'utc': stamp()})
        print(json.dumps({'results_sha256': digest, 'rows': len(rows), 'contrast_gradient': result['rounded_affine_contrast_gradient_um_display']}), flush=True)
    except BaseException as e:
        save(directory/'terminal.json', {'status': 'FAILED', 'pid': os.getpid(), 'exit_code': 1,
                                        'error': str(e), 'utc': stamp()})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    sub.add_parser('prepare')
    run = sub.add_parser('actual')
    run.add_argument('--review-sha256', required=True)
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        actual(args.review_sha256)
