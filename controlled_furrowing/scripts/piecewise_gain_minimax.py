"""Certified minimax common gain for a declared piecewise-linear response table.

Certificates concern this surrogate only. They do not enclose nonlinear
mechanics between strength nodes or prove a physical identified set.
"""
from fractions import Fraction as Q
from itertools import combinations
import bisect
import json
import math
from pathlib import Path
import hashlib
import numpy as np
from scipy.optimize import linprog


def rat(x):
    if isinstance(x,bool):
        raise ValueError('Boolean is not a numeric response')
    try:
        return x if isinstance(x,Q) else Q(str(x))
    except (ValueError,TypeError,ZeroDivisionError) as e:
        raise ValueError('finite rational value required') from e


def envelope(lines,x):
    return max(a*x+b for a,b in lines)


def certify_cell(lines,left,right):
    """Use floating LP only to find a candidate; exact convex certificate decides."""
    lines = [(rat(a),rat(b)) for a,b in lines]+[(Q(0),Q(0))]
    left,right = rat(left),rat(right)
    if left >= right:
        raise ValueError('positive-width gain cell required')
    result = linprog([0.,1.], A_ub=np.array([[float(a),-1.] for a,b in lines]),
                     b_ub=np.array([-float(b) for a,b in lines]),
                     bounds=[(float(left),float(right)),(0,None)],method='highs')
    if not result.success:
        raise RuntimeError('candidate LP failed: '+result.message)
    approximate = float(result.x[0])
    values = [float(a)*approximate+float(b) for a,b in lines]
    top = max(values)
    active_hint = [i for i,v in enumerate(values) if top-v <= 1e-7*max(1.,abs(top))]

    def candidates(indices):
        yield left
        yield right
        for i,j in combinations(indices,2):
            a,b = lines[i]; c,d = lines[j]
            if a == c:
                continue
            x = (d-b)/(a-c)
            if left <= x <= right:
                yield x

    def certificate(x):
        loss = envelope(lines,x)
        active = [i for i,(a,b) in enumerate(lines) if a*x+b == loss]
        slopes = [lines[i][0] for i in active]
        if x == left and max(slopes) >= 0:
            condition = 'right derivative nonnegative at left boundary'
        elif x == right and min(slopes) <= 0:
            condition = 'left derivative nonpositive at right boundary'
        elif left < x < right and min(slopes) <= 0 <= max(slopes):
            condition = 'zero belongs to convex subdifferential'
        else:
            return None
        return {'gain_exact':str(x),'loss_exact':str(loss),'active_line_indices':active,
                'active_slopes_exact':list(map(str,slopes)), 'optimality_condition':condition,
                'cell_exact':[str(left),str(right)],
                'all_line_coefficients_exact':[[str(a),str(b)] for a,b in lines]}

    certified = [r for x in set(candidates(active_hint)) if (r:=certificate(x)) is not None]
    if not certified:
        certified = [r for x in set(candidates(range(len(lines)))) if (r:=certificate(x)) is not None]
    if not certified:
        raise RuntimeError('no exact cell-optimality certificate; do not promote floating LP output')
    return min(certified,key=lambda r:(Q(r['loss_exact']),Q(r['gain_exact'])))


def fit(strengths,responses,doses,observations,gain_bounds):
    """Responses shape KxP, observations NxP; all P channels count in the loss."""
    f = list(map(rat,strengths)); I = list(map(rat,doses))
    y = [[rat(x) for x in row] for row in observations]
    v = [[rat(x) for x in row] for row in responses]
    lo,hi = map(rat,gain_bounds)
    if len(f) < 2 or any(b <= a for a,b in zip(f,f[1:])) or f[0] < 0:
        raise ValueError('strictly increasing nonnegative strength nodes required')
    if not I or any(x <= 0 for x in I) or len(y) != len(I) or len(v) != len(f):
        raise ValueError('positive doses and matching observation/response rows required')
    channels = len(y[0])
    if not channels or any(len(row) != channels for row in y+v):
        raise ValueError('all channels required in every row')
    if lo < 0 or hi <= lo or min(I)*lo < f[0] or max(I)*hi > f[-1]:
        raise ValueError('gain interval must have positive width and forbid table extrapolation')
    cuts = sorted({lo,hi}|{node/dose for node in f for dose in I if lo < node/dose < hi})
    cells = []
    for left,right in zip(cuts,cuts[1:]):
        midpoint = (left+right)/2
        lines = []
        for dose,row in zip(I,y):
            j = min(bisect.bisect_right(f,midpoint*dose)-1,len(f)-2)
            for channel,observed in enumerate(row):
                slope = (v[j+1][channel]-v[j][channel])/(f[j+1]-f[j])
                a = slope*dose; b = v[j][channel]-slope*f[j]-observed
                lines.extend(((a,b),(-a,-b)))
        cells.append(certify_cell(lines,left,right))
    best = min(cells,key=lambda r:(Q(r['loss_exact']),Q(r['gain_exact'])))
    return {'status':'EXACT_SURROGATE_GAIN_OPTIMUM', 'gain_exact':best['gain_exact'],
            'loss_exact':best['loss_exact'], 'gain_display':float(Q(best['gain_exact'])),
            'loss_display':float(Q(best['loss_exact'])),
            'at_declared_gain_boundary':Q(best['gain_exact']) in (lo,hi),
            'gain_bounds_exact':[str(lo),str(hi)], 'records':len(I),'channels':channels,
            'constraints':len(I)*channels, 'certified_cells':cells,
            'scope':'Global exact minimax gain for supplied piecewise-linear table on declared gain interval only; no physical interpolation/PDE/peak error bound.'}


def verify():
    checks = {}
    # Linear response: min max(|g-1|,|2g-1|) occurs at g=2/3, loss=1/3.
    a = fit([0,4],[[0],[4]],[1,2],[[1],[1]],(0,2))
    checks['analytic_two_record_linear'] = (a['gain_exact'],a['loss_exact']) == ('2/3','1/3')
    # Two channels: max(|g-1|,|2g-5|), equal opposing residuals at g=2.
    b = fit([0,4],[[0,0],[4,8]],[1],[[1,5]],(0,4))
    checks['all_channels_affect_exact_optimum'] = (b['gain_exact'],b['loss_exact']) == ('2','1')
    c = fit([0,1,3],[[0],[2],[2]],[1],[[3]],(0,3))
    checks['plateau_smallest_gain_tiebreak'] = (c['gain_exact'],c['loss_exact']) == ('1','1')
    d = fit([0,2],[[0],[2]],[1],[[4]],(0,2))
    checks['boundary_is_flagged'] = d['at_declared_gain_boundary'] and d['gain_exact'] == '2'
    # Independent replay of every exact certificate, not the candidate LP.
    for name,result in [('a',a),('b',b),('c',c),('d',d)]:
        good = True
        for cell in result['certified_cells']:
            x,loss = Q(cell['gain_exact']),Q(cell['loss_exact'])
            lines = [(Q(m),Q(k)) for m,k in cell['all_line_coefficients_exact']]
            l,r = map(Q,cell['cell_exact'])
            active = [m for m,k in lines if m*x+k == loss]
            good = good and l <= x <= r and max(m*x+k for m,k in lines) == loss
            good = good and ((x == l and max(active) >= 0) or (x == r and min(active) <= 0)
                             or (l < x < r and min(active) <= 0 <= max(active)))
        checks['exact_certificate_replay_'+name] = good
    invalid = [lambda:fit([0,1,1],[[0],[1],[1]],[1],[[1]],(0,1)),
               lambda:fit([0,1],[[0],[1]],[0],[[1]],(0,1)),
               lambda:fit([0,1],[[0],[1]],[2],[[1]],(0,1)),
               lambda:fit([0,1],[[0],[1,2]],[1],[[1]],(0,1))]
    for i,fn in enumerate(invalid):
        try:
            fn(); checks['invalid_'+str(i)] = False
        except ValueError:
            checks['invalid_'+str(i)] = True
    root = Path(__file__).resolve().parents[1]
    receipt = {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
               'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'scope':'Exact surrogate optimization and certificate replay; not nonlinear physical inference.'}
    (root/'evidence/piecewise_gain_minimax_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    if not all(checks.values()):
        raise RuntimeError('surrogate gain verification failed')


if __name__ == '__main__':
    verify()
