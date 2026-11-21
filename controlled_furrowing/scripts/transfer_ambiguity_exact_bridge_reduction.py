"""Pure signed reduction for separately admitted singleton spatial bridge packets.

This module does not admit files, certify any physical input, or launch integrals.
Its caller must authenticate all component proofs, whole-box domains, the cached
point rows and ZQ, and the exact same binary R before using its output as evidence.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import hashlib
import json
from pathlib import Path

from flint import arb, arb_mat, ctx, fmpq

POINT = 'EXACT_CENTER_POINT_ONLY'
BOX = 'ENTIRE_CLOSED_FIXED_PARAMETER_CUBE'
RADIUS = Q(1, 2**34)
SOURCE = Path(__file__).resolve()
OUT = SOURCE.parents[1] / 'runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/bridge_reduction_v1'


def require(value, message):
    if not value:
        raise ValueError(message)


def rational(value):
    require(type(value) in (str, int, Q), 'explicit finite rational required; no float/bool conversion')
    return Q(value)


def interval(value):
    require(type(value) in (list, tuple) and len(value) == 2, 'two endpoints required')
    lo, hi = map(rational, value)
    require(lo <= hi, 'ordered interval required')
    return lo, hi


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(q, a):
    return (q*a[0], q*a[1]) if q >= 0 else (q*a[1], q*a[0])


def maxabs(a):
    return max(abs(a[0]), abs(a[1]))


def minabs(a):
    return Q(0) if a[0] <= 0 <= a[1] else min(abs(a[0]), abs(a[1]))


def assemble(packets, expected_maps, n):
    """Independent exact endpoint scatter; no producer scatter helper is used."""
    require(type(n) is int and n > 0, 'positive dimension required')
    require(len(packets) == len(expected_maps) > 0, 'complete admitted component set required')
    require(all(m['kind'] in ('bulk','surface') and type(m['slot']) is int and m['slot'] >= 0
                and type(m['global_indices']) is list and len(set(m['global_indices'])) == len(m['global_indices'])
                and all(type(i) is int and 0 <= i < n for i in m['global_indices'])
                for m in expected_maps), 'typed authenticated expected component maps required')
    order = [(m['kind'], m['slot']) for m in expected_maps]
    require(len(set(order)) == len(order), 'no duplicate expected component')
    force = [(Q(0), Q(0)) for _ in range(n)]
    hessian = {}
    for packet, projection in zip(packets, expected_maps, strict=True):
        require(type(packet['slot']) is int and packet['slot'] >= 0
                and (packet['kind'], packet['slot']) == (projection['kind'], projection['slot']),
                'ordered component identity required')
        ids = packet['global_indices']
        require(ids == projection['global_indices'] and len(set(ids)) == len(ids)
                and all(type(i) is int and 0 <= i < n for i in ids), 'authenticated free projection required')
        f, h = packet['point_force'], packet['whole_box_hessian']
        require(f['scope'] == POINT and h['scope'] == BOX and rational(h['radius_exact']) == RADIUS,
                'point force and entire larger-box Hessian must remain separate')
        require(len(f['F_I_minus_F_cache_intervals_exact']) == len(ids), 'complete local force required')
        for i, value in zip(ids, f['F_I_minus_F_cache_intervals_exact'], strict=True):
            force[i] = add(force[i], interval(value))
        upper = h['H_I_minus_H_cache_upper_exact']
        require(all(type(row) in (list,tuple) and len(row) == 4
                    and type(row[0]) is type(row[1]) is int for row in upper),
                'typed integer local upper Hessian indices required')
        require([(i,j) for i,j,*_ in upper] == list(combinations_with_replacement(range(len(ids)),2)),
                'complete ordered symmetric upper Hessian required')
        for i, j, lo, hi in upper:
            key = tuple(sorted((ids[i], ids[j])))
            hessian[key] = add(hessian.get(key, (Q(0), Q(0))), interval((lo, hi)))
    return force, hessian


def classify(rows, z, radius):
    rows = [interval(row) for row in rows]
    z, radius = rational(z), rational(radius)
    require(bool(rows) and z >= 0 and radius > 0, 'finite derivative upper bound and positive cube required')
    eta = max(map(maxabs, rows))
    lower = max(map(minabs, rows))
    lhs = eta + z*radius
    exists = z < 1 and lhs < radius
    excludes = lower > (1+z)*radius
    require(not (exists and excludes), 'inconsistent sufficient classifications')
    return {'eta_I_upper_exact':str(eta), 'RF_I_norm_lower_exact':str(lower),
            'z_upper_exact':str(z), 'radius_exact':str(radius),
            'selfmap_lhs_upper_exact':str(lhs), 'selfmap_margin_lower_exact':str(radius-lhs),
            'strict_contraction_and_selfmap':exists,
            'posterior_error_upper_exact':str(eta/(1-z)) if exists else None,
            'strict_exclusion_of_this_cube':excludes,
            'exclusion_margin_lower_exact':str(lower-(1+z)*radius),
            'classification':('SUFFICIENT_UNIQUE_LOCAL_ZERO' if exists else
                              'SUFFICIENT_NO_ZERO_IN_THIS_CUBE' if excludes else 'UNRESOLVED'),
            'minimum_or_continuum_claimed':False,
            'physical_inputs_independently_admitted_by_this_function':False}


def ball(q):
    q = rational(q)
    return arb(fmpq(q.numerator, q.denominator))


def ball_interval(value):
    lo, hi = interval(value)
    return ball(lo).union(ball(hi))


def endpoints(value):
    require(value.is_finite(), 'finite matrix result required')
    return Q(str(value.lower().fmpq())), Q(str(value.upper().fmpq()))


def contract(R, force, upper_H, cached_RF_rows, ZQ, radius=RADIUS):
    """Arb-enclose signed R*defect, add signed cached anchor, then take norms.

    R contains exact binary rationals, never shortest-decimal float spellings.
    upper_H is a sparse symmetric matrix of whole-box defect intervals.
    This function has no authority to validate a physical model or input proof.
    """
    require(ctx.prec >= 640, 'set Arb precision >=640 before constructing any input ball')
    n = len(R)
    require(n > 0 and len(force) == len(cached_RF_rows) == n
            and all(len(row) == n for row in R), 'matching complete matrix/vector dimensions required')
    zq = rational(ZQ)
    require(zq >= 0, 'nonnegative independently supplied whole-box ZQ required')
    mat = arb_mat(n,n,[ball(q) for row in R for q in row])
    fv = arb_mat(n,1,[ball_interval(row) for row in force])
    hm = arb_mat(n,n)
    for (i,j), value in upper_H.items():
        require(type(i) is type(j) is int and 0 <= i <= j < n, 'upper global Hessian index required')
        hm[i,j] = ball_interval(value)
        if i != j:
            hm[j,i] = hm[i,j]
    rf = mat*fv
    rh = mat*hm
    # Sum exact rational outward endpoints. No floating norm or midpoint decision.
    rh_rows = [sum((maxabs(endpoints(rh[i,j])) for j in range(n)), Q(0)) for i in range(n)]
    b = max(rh_rows)
    defect_rows = [endpoints(rf[i,0]) for i in range(n)]
    full_rows = [add(interval(anchor), defect) for anchor,defect in zip(cached_RF_rows,defect_rows,strict=True)]
    return {**classify(full_rows,zq+b,radius), 'b_upper_exact':str(b), 'ZQ_upper_exact':str(zq),
            'signed_RF_defect_rows_exact':[list(map(str,row)) for row in defect_rows],
            'signed_RF_I_rows_exact':[list(map(str,row)) for row in full_rows],
            'RH_absolute_row_sum_bounds_exact':list(map(str,rh_rows)), 'precision_bits':ctx.prec}


def rejected(call):
    try:
        call()
    except ValueError:
        return True
    return False


def fixtures():
    from copy import deepcopy
    ctx.prec = 640
    def packet(slot, ids, force, upper):
        return {'kind':'bulk','slot':slot,'global_indices':ids,
                'point_force':{'scope':POINT,'F_I_minus_F_cache_intervals_exact':force},
                'whole_box_hessian':{'scope':BOX,'radius_exact':str(RADIUS),
                                    'H_I_minus_H_cache_upper_exact':upper}}
    p = packet(0,[1,0],[['2','2'],['1','1']],[[0,0,'1','1'],[0,1,'2','2'],[1,1,'3','3']])
    q = packet(1,[0],[['-1','-1']],[[0,0,'-2','-2']])
    maps = [{k:v[k] for k in ('kind','slot','global_indices')} for v in (p,q)]
    F,H = assemble([p,q],maps,2)
    R = [[Q(1),Q(-1)],[Q(0),Q(1)]]
    out = contract(R,F,H,[['0','0'],['0','0']],'1/8',Q(1))
    # Direct rational reference: F=(0,2), H=((1,2),(2,1)),
    # RF=(-2,2), RH=((-1,1),(2,1)), hence ||RH||=3.
    bad = deepcopy(p)
    bad['point_force']['scope'] = BOX
    incomplete = deepcopy(p)
    incomplete['whole_box_hessian']['H_I_minus_H_cache_upper_exact'].pop()
    bool_slot = deepcopy(q)
    bool_slot['slot'] = True
    bool_H = deepcopy(p)
    bool_H['whole_box_hessian']['H_I_minus_H_cache_upper_exact'][0][0] = False
    float_H = deepcopy(p)
    float_H['whole_box_hessian']['H_I_minus_H_cache_upper_exact'][0][0] = 0.0
    bool_maps = deepcopy(maps)
    bool_maps[0]['global_indices'][0] = True
    root = classify([['-1/8','1/16']],'1/2','1')
    tie = classify([['1/2','1/2']],'1/2','1')
    exclude = classify([['-2','-13/8']],'1/2','1')
    # Nondegenerate interval contraction compared with an exact endpoint sum.
    IF = [(Q(-1,3),Q(2,7)),(Q(-3,5),Q(4,9))]
    IR = [[Q(-2,3),Q(5,7)],[Q(7,11),Q(-3,13)]]
    probe = contract(IR,IF,{},[['0','0'],['0','0']],'0')
    reference = [add(scale(row[0],IF[0]),scale(row[1],IF[1])) for row in IR]
    return {
        'signed_shared_dof_cancellation_and_permuted_symmetric_scatter':
            F == [(Q(0),Q(0)),(Q(2),Q(2))] and H == {(1,1):(Q(1),Q(1)),(0,1):(Q(2),Q(2)),(0,0):(Q(1),Q(1))},
        'signed_R_contraction_agrees_with_direct_rational_matrix':
            out['signed_RF_defect_rows_exact'] == [['-2','-2'],['2','2']] and Q(out['b_upper_exact']) == 3,
        'nondegenerate_Arb_intervals_contain_exact_endpoint_reference':all(
            interval(got)[0] <= ref[0] <= ref[1] <= interval(got)[1]
            for got,ref in zip(probe['signed_RF_defect_rows_exact'],reference,strict=True)),
        'strict_root_posterior_and_boundary':root['strict_contraction_and_selfmap']
            and Q(root['posterior_error_upper_exact']) == Q(1,4) and not tie['strict_contraction_and_selfmap'],
        'strict_cube_exclusion_is_not_global_absence':exclude['strict_exclusion_of_this_cube']
            and Q(exclude['exclusion_margin_lower_exact']) == Q(1,8),
        'reject_mixed_scope_partial_H_changed_map_duplicate_or_partial_packet':
            rejected(lambda:assemble([bad,q],maps,2)) and rejected(lambda:assemble([incomplete,q],maps,2))
            and rejected(lambda:assemble([q,p],maps,2)) and rejected(lambda:assemble([p],maps,2))
            and rejected(lambda:assemble([p,p],[maps[0],maps[0]],2)),
        'reject_float_bool_unordered_and_negative_derivative':
            rejected(lambda:rational(0.1)) and rejected(lambda:rational(True))
            and rejected(lambda:interval(['1','0'])) and rejected(lambda:classify([['0','0']],'-1','1')),
        'reject_boolean_or_float_component_and_local_index_lookalikes':
            rejected(lambda:assemble([p,bool_slot],maps,2)) and rejected(lambda:assemble([bool_H,q],maps,2))
            and rejected(lambda:assemble([float_H,q],maps,2)) and rejected(lambda:assemble([p,q],bool_maps,2)),
        'no_actual_physical_packets_or_force_jobs_read_or_run':True}


def prepare():
    folder = OUT/'preparation_v2'
    folder.mkdir(parents=True,exist_ok=False)
    checks = fixtures()
    require(all(type(v) is bool and v for v in checks.values()), 'reduction mathematical fixtures failed')
    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    (folder/'transfer_ambiguity_exact_bridge_reduction_frozen.py').write_bytes(SOURCE.read_bytes())
    protocol = {'status':'PURE_SINGLETON_REDUCTION_PREPARATION_ONLY','source_sha256':source_hash,
                'scope':'No physical input admission or actual reduction; independent source review required.',
                'formula':'eta=||R Fcache(u0)+R DF(u0)||inf; z=ZQ(r)+||R DH(B_r)||inf; eta+z*r<r, z<1.',
                'new_force_or_H_integrals':0,'physical_root_claims':0,'requires_full_168_packet_independent_admission':True}
    for name,value in (('protocol.json',protocol),('preflight.json',{'status':'PASS','checks':checks,'source_sha256':source_hash})):
        (folder/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'source_sha256':source_hash,'checks':len(checks),'scope':protocol['scope']}))


if __name__ == '__main__':
    prepare()
