"""Plan and replay separate singleton point-force and whole-box Hessian defects."""
import argparse
from contextlib import ExitStack, contextmanager
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import combinations_with_replacement
import json
from math import ceil
import os
from pathlib import Path
import subprocess
import sys
import traceback
from unittest.mock import patch

import transfer_ambiguity_cached_derivative_extension as DX
import transfer_ambiguity_exact_targets as EX
import transfer_ambiguity_observer_sweep as sweep


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
OUT = EX.OUT / 'bridge_v1'
PREP = OUT / 'preparation'
PLAN = OUT / 'plans'
ACTUAL = OUT / 'actual'
DESIGN = DX.DESIGN_RECEIPT
DERIV = EX.relative(DX.ACTUAL)
DERIV_SOURCE = 'scripts/transfer_ambiguity_cached_derivative_extension.py'
FORCE_SOURCE = 'runs/exact_integration_cohort_v1/L_defect_diagnostic/element21_force.py'
SCALED_SOURCE = 'scripts/exact_scaled_cohort_bulk.py'
PINS = {
    DERIV_SOURCE: '7ad0b7d4d42efe634123989c6092ccf768cedb92a79244cc7aaa5c10a62f370c',
    EX.relative(DX.PREP/'protocol.json'): '4af874d6e45836d11b608ae2929ff5120088fdbd549fd3615398c0068b29ebd7',
    EX.relative(DX.PREP/'preflight.json'): '899d9c9486aa0412934bf8866d24c68842c0d2085c380089e1852d906f343f4f',
    DESIGN: DX.PINS[DESIGN],
    FORCE_SOURCE: '9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329',
    SCALED_SOURCE: '6639bb2fae3c3892ec44115ff90e27242586744a2900daa16177cc5a091f1526',
}
NFREE = 611
RADIUS = Q(1,2**34)
POINT = 'EXACT_CENTER_POINT_ONLY'
BOX = 'ENTIRE_CLOSED_FIXED_PARAMETER_CUBE'
SCOPE = ('Pointwise exact spatial FE minus binary cache force at X+inj(u), and exact spatial FE minus '
         'binary cache Hessian on a separate entire r=2^-34 box, with the same frozen singleton drive. '
         'No root, local minimum, continuum or observation claim in this component phase.')
E = EX.E
T = EX.T
C = E.core
SUR = E.SUR
U = E.U
Qx = EX.exact
require, save, sha, pinned = sweep.require,sweep.save,sweep.sha,sweep.pinned


class OrderUnresolved(ValueError):
    def __init__(self, attempts):
        super().__init__('bounded surface order policy exhausted')
        self.attempts = attempts


def digest(obj):
    return E._digest(EX.doc(obj))


def boot():
    EX.bootstrap()
    require(C.ctx.prec==640 and E.core is C and Path(DX.__file__).resolve()==ROOT/DERIV_SOURCE,
            'current explicit singleton module aliases and Arb640 required')


def kernels(links):
    for path in (FORCE_SOURCE,SCALED_SOURCE):
        EX.link(links,path,PINS[path])
    # This import reuses functions only. Historical drivers have no route here.
    spec=importlib.util.spec_from_file_location('_singleton_explicit_force_parts',ROOT/FORCE_SOURCE)
    force=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(force)
    import exact_scaled_cohort_bulk as scaled
    require(Path(scaled.__file__).resolve()==ROOT/SCALED_SOURCE and scaled.E is force.E is E
            and scaled.core is force.core is C and scaled.T is force.T is T,
            'coherent imported exact kernels required')
    path=EX.relative(Path(scaled.S.__file__))
    expected=pinned(ROOT/DESIGN,PINS[DESIGN])['source_and_metadata_sha256'][path]
    EX.link(links,path,expected)
    return force,scaled


def static_inputs():
    boot()
    links,cases,manifest,data=DX.authority()
    for path,h in PINS.items():
        EX.link(links,path,h)
    design=pinned(ROOT/DESIGN,PINS[DESIGN])
    for path,h in design['source_and_metadata_sha256'].items():
        EX.link(links,path,h)
    links[EX.relative(SOURCE)]=sha(SOURCE)
    return dict(sorted(links.items())),cases,manifest,data


def authority_pins(state):
    links,cases,manifest,_=state
    return {'status':'FROZEN_EXACT_SINGLETON_BRIDGE_COLD_PREPARATION','source_sha256':sha(SOURCE),
            'scope':SCOPE,'input_sha256':links,'pilot_order':[x['pilot_id'] for x in cases],
            'target_manifest_identity_sha256':digest(manifest),'whole_H_radius_exact':str(RADIUS),'precision_bits':640,
            'point_scope':POINT,'whole_hessian_scope':BOX,
            'requires_reviewed_larger_cached_derivative_admission':True,
            'requires_independent_exact_per_component_plan_review':True,
            'point_force_route':'144 force-only quadratic/volumetric/same-Y cache corrections and 24 surface passes.',
            'whole_hessian_route':'144 scaled-off-axis or explicit-axis whole-box and 24 whole-box surface passes.',
            'surfaces_compute_unused_cross_shape_channels':True,
            'upper_H_packed_in_local_indices':True,'scatter_signed_before_any_R_application':True,
            'tail_policy':'Use measured ZQ (no quarter gate) and positive h. Min of original optional local caps and '
                          'one ten-thousandth of prospective force/H caps allocated through actual incidence and verified R norm.',
            'max_bulk_order':32,'max_surface_order':64,'surface_depth':8,
            'force_H_or_component_order_selection_run_during_preparation':False}


def exact_point_inside(point,box):
    X=box['current_centers_exact']
    radii=box['coordinate_radii_exact']
    hulls=box['coordinate_hulls_exact']
    require(len(point)==len(X)==len(radii)==len(hulls)==325
            and all(len(row)==2 for group in (point,X,radii,hulls) for row in group),
            'complete corresponding point and hull coordinates required')
    require(all(Qx(lo)==Qx(c)-Qx(r) <= Qx(value)==Qx(c) <= Qx(c)+Qx(r)==Qx(hi)
                for yr,cr,rr,hr in zip(point,X,radii,hulls,strict=True)
                for value,c,r,(lo,hi) in zip(yr,cr,rr,hr,strict=True)),
            'true exact center inclusion in independently verified whole box required')
    zero=[['0','0'] for _ in point]
    return {**box,'coordinate_radii_exact':zero,'coordinate_hulls_exact':[
                [[str(v),str(v)] for v in row] for row in point],
            'corrector_radius_exact':'0','scope':POINT}


def derivative_inputs(path,h,cases):
    review=pinned(path,h)
    require(review['status']=='PASS_REVIEWED_LARGER_CACHED_DERIVATIVE_BRIDGE_INPUTS'
            and sweep.all_true(review['checks']) and Qx(review['radius_exact'])==RADIUS,
            'positive independent whole-box cached derivative admission expected')
    proof_file=review['native_independent_review_file']
    require((ROOT/proof_file).resolve().is_relative_to(DX.OUT.resolve()),
            'native independent derivative review must reside in the explicit larger-radius archive')
    native=pinned(ROOT/proof_file,review['native_independent_review_sha256'])
    require(native['status'].startswith('PASS_INDEPENDENT') and sweep.all_true(native['checks']),
            'independent native larger-radius numeric review required')
    top=pinned(DX.ACTUAL/'results.json',review['derivative_results_sha256'])
    end=pinned(DX.ACTUAL/'terminal.json',review['derivative_terminal_sha256'])
    require(top['status']=='COMPLETE_FOUR_CACHED_DERIVATIVE_ATTEMPTS_UNREVIEWED'
            and end['status']=='FINISHED' and end['results_sha256']==review['derivative_results_sha256']
            and end['bindings']==top['bindings'] and len(end['workers'])==4
            and all(type(c['observed_exit_code']) is int and c['observed_exit_code']==0 for c in end['workers'])
            and top['bindings']['source_sha256']==PINS[DERIV_SOURCE]
            and tuple(r['pilot_id'] for r in top['case_results'])==EX.PILOTS
            and tuple(r['pilot_id'] for r in review['reviewed_cases'])==EX.PILOTS,
            'same reviewed derivative stage and four real successful subprocesses required')
    available={}
    for case,r,w,prior in zip(top['case_results'],review['reviewed_cases'],end['workers'],cases,strict=True):
        pid=case['pilot_id']
        f=DX.ACTUAL/pid
        require(w['pilot_id']==r['pilot_id']==prior['pilot_id']==pid
                and w['supervisor_pid']==end['pid'] and r['producer_result_sha256']==case['results_sha256'],
                'mixed larger derivative record or authority')
        local=pinned(f/'results.json',case['results_sha256'])
        child=pinned(f/'terminal.json',w['terminal_sha256'])
        exit_p=pinned(f/'observed_os_exit.json',w['os_exit_receipt_sha256'])
        DX.CP.admit_child(exit_p['observed_exit_code'],child,{k:w[k] for k in ('pilot_id','pid','supervisor_pid','bindings')})
        require(child['results_sha256']==case['results_sha256']
                and local['case']['box_sha256']==r['box_sha256']==prior['box_sha256']
                and local['case']['cache_packet_sha256']==r['cache_packet_sha256']==prior['cache_packet_sha256']
                and Qx(local['whole_box_cache_radius_exact'])==RADIUS
                and local['whole_box_controls_only_cached_quadrature'] is True
                and case['ZQ_upper_exact']==local['ZQ_upper_exact']
                and 0<=Qx(r['reviewed_ZQ_upper_exact']),
                'same-center whole-radius quadrature derivative enclosures required')
        z=max(Qx(r['reviewed_ZQ_upper_exact']),Qx(local['ZQ_upper_exact']))
        p=DX.spend(z,RADIUS,Qx(prior['R_inf_exact']),prior['prior_optional_tail_policy'])
        available[pid]={'ZQ_upper_exact':str(z),'planning':p,'reviewed_derivative_result_sha256':case['results_sha256'],
                        'independent_native_file':proof_file,
                        'independent_native_sha256':review['native_independent_review_sha256']}
    return available


def domain(state,pid):
    _,cases,manifest,data=state
    boot()
    ctx,rule,entry=EX.open_context(pid,manifest,data)
    case,=[x for x in cases if x['pilot_id']==pid]
    box=pinned(ROOT/entry['box_1'],entry['box_1_sha256'])
    require(box==EX.exact_node_box(ctx,RADIUS) and entry['box_1_sha256']==case['box_sha256'],
            'factory exact large box and reviewed bytes must match')
    route=pinned(ROOT/entry['routes'],entry['routes_sha256'])
    g=pinned(ROOT/case['geometry'],case['geometry_sha256'])
    b=pinned((ROOT/case['geometry']).with_name('body_geometry.json'),g['body_geometry_sha256'])
    require(g['box_sha256']==entry['box_1_sha256'] and g['routes_sha256']==entry['routes_sha256']
            and g['rule_sha256']==entry['rule_sha256'] and len(route['bulk'])==len(b['triangles'])==144
            and len(route['surface'])==len(g['trace_geometry_only'])==24
            and b['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING',
            'whole box and all proof routes must match before component generation')
    current=ctx.target.nodes()
    point_box=exact_point_inside(current,box)
    return {'ctx':ctx,'rule':rule,'entry':entry,'case':case,'box':box,'point_box':point_box,
            'point':current,'nodes':SUR.common_nodes(ctx.target,box),'routes':route,'geometry':g,'body':b}


def polynomials(fields,jet,keys):
    def norm(poly):
        return [[i,j,list(map(str,SUR.rational_hull(v)))] for (i,j),v in sorted(poly.items())]
    q,b = (fields[name] for name in keys)
    return digest({'q':[norm(x) for x in q],'B':[[norm(x) for x in row] for row in b],
                   'J_or_G':norm(jet['value']),'first':[norm(x) for x in jet['gradient']],
                   'second':[[i,j,norm(x)] for (i,j),x in sorted(jet['hessian_upper'].items())]})


def coefficientwise_inside(point, whole, keys):
    def pack(fields,jet):
        q,B=(fields[k] for k in keys)
        return [*q,*(p for row in B for p in row),jet['value'],*jet['gradient'],
                *(jet['hessian_upper'][k] for k in sorted(jet['hessian_upper']))]
    exact,box=pack(*point),pack(*whole)
    count=0
    for a,b in zip(exact,box,strict=True):
        for key in a.keys() | b.keys():
            lo,hi=SUR.rational_hull(b.get(key,0))
            qlo,qhi=SUR.rational_hull(a.get(key,0))
            require(lo<=qlo<=qhi<=hi,'point polynomial coefficient is not enclosed by its whole-box counterpart')
            count+=1
    return count


def point_maps(d,i,fields,ids):
    target,rule=d['ctx'].target,d['rule']
    view=U.RuleTarget(target,rule)
    maps=[U.C.exact_bulk_maps(view,i,k) for k in range(64)]
    Y=[T.rational(v) for j in ids for v in d['point'][j]]
    lower=[]
    for m in maps:
        require(m['global_indices']==fields['global_indices'] and m['selected_local']==fields['free_local']
                and Qx(m['auxiliary_normalized_weight_exact'])>0,'point canonical/cache local map mismatch')
        v=[]
        for star,cached in zip(m['star_full'],m['cache_full'],strict=True):
            q=sum(x*y for x,y in zip(cached,Y,strict=True))
            delta=sum((x-y)*z for x,y,z in zip(star,cached,Y,strict=True))
            v.append((q+min(0,delta),q+max(0,delta)))
        def mul(a,b):
            corners=[x*y for x in a for y in b]
            return min(corners),max(corners)
        ad,bc=mul(v[0],v[3]),mul(v[1],v[2])
        det=ad[0]-bc[1]
        require(det>0 and v[4][0]>0 and det*v[4][0]>Q(1,10**6),
                'point same-Y map path fails strict positive cached analytic domain')
        lower.append(det*v[4][0])
    return maps,Y,{'map_packets_sha256':digest(maps),'positive_line_count':len(maps),
                   'minimum_exact_convex_line_J_lower':str(min(lower))}


def build_bulk(d,i,tail,scaled):
    C.A.cache_clear()
    ctx,rule,expected=d['ctx'],d['rule'],d['routes']['bulk'][i]
    b=d['body']['triangles'][i]
    boxes={}
    polynomial_pairs=[]
    for scope,nodes in ((POINT,d['point']),(BOX,d['nodes'])):
        setup=E._bulk_setup(ctx,nodes,rule,i)
        identity,fields,_,jet,mass,grad,hess=setup
        require(all(identity[k]==expected[k] for k in ('chart','free_local','global_indices','active',
                                                         'auxiliary_element_rule_sha256')),
                'fresh bulk polynomial is on a different projection or rule')
        try:
            selector=E._select(identity,b,mass,grad,hess,*tail,32)
        except E.TailBudgetUnresolved as error:
            error.component_scope=scope
            raise
        stamp={'selector':selector,'polynomial_jets_sha256':polynomials(fields,jet,('q','B'))}
        polynomial_pairs.append((fields,jet))
        if scope==POINT:
            _,_,lines=point_maps(d,i,fields,expected['nodes'])
            stamp['same_Y_point_cache_lines']=lines
        else:
            stamp['scaled_or_axis_proof']=EX.doc(scaled.prepare_component(ctx,nodes,rule,b,selector,i))
        boxes[scope]=stamp
    return {'kind':'bulk','slot':i,'global_indices':expected['global_indices'],
            'exact_point_coefficients_enclosed_by_whole_box_count':coefficientwise_inside(*polynomial_pairs,('q','B')),
            'body_element_sha256':digest(b),'point_J_guard_scope':'Exact point is inside the independently certified same-element box.',
            **boxes}


def surface_selector(d,fields,domain,f,e,tail):
    domain=dict(domain)
    lo,hi=map(Qx,domain['G_bounds_exact'])
    q0,raw=Qx(domain['q0_exact']),Qx(domain['rho_exact'])
    rho=Q(ceil(raw*2**48),2**48)
    require(domain['status']=='CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN'
            and 0<lo<=hi and q0==(hi+lo)/2 and raw==(hi-lo)/(hi+lo) and 0<=raw<=rho<1,
            'fresh positive exact trace expansion domain required')
    domain['raw_rho_exact'],domain['rho_exact']=str(raw),str(rho)
    jet=E.JET.surface_squared_jacobian_jet(fields['q_physical_R'],fields['B_physical_R'])
    _,_,K,Tn=SUR.material_boxes(d['ctx'].row,f)
    scale=max(abs(t-k) for k in K for t in Tn)
    a,b=fields['Rlo'],fields['Rhi']
    sqrt=Q(str(C.A(q0).sqrt().lower().fmpq()))
    require(sqrt>0 and sqrt**2<=q0,'positive square-root lower endpoint required')
    mass=Q(str(C.upper(d['ctx'].target.normalization()*C.A((b*b-a*a)/2)
                       +C.A(sum(d['rule']['surface_weights'][f][e],Q(0)))).fmpq()))
    g=max(SUR.poly_norm_upper(v,a,b) for v in jet['gradient'])
    H=max(SUR.poly_norm_upper(v,a,b) for v in jet['hessian_upper'].values())
    attempts=[]
    for order in range(2,65):
        at=scale*SUR.S.binomial_tail(Q(-1,2),order,rho)/(2*sqrt)
        bt=scale*SUR.S.binomial_tail(Q(-3,2),order,rho)/(4*q0*sqrt)
        fbound,hbound=mass*at*g,mass*(at*H+bt*g*g)
        attempts.append({'order':order,'force_tail_upper_exact':str(fbound),'H_entry_tail_upper_exact':str(hbound)})
        if fbound<=tail[0] and hbound<=tail[1]:
            return order,attempts,domain,jet
    raise OrderUnresolved(attempts)


def build_surface(d,slot,tail):
    C.A.cache_clear()
    f,e=divmod(slot,12)
    ctx=d['ctx']
    expected=d['routes']['surface'][slot]
    frozen=d['geometry']['trace_geometry_only'][slot]
    require(expected['face']==f==frozen['face'] and expected['edge']==e==frozen['edge'],
            'same pinned explicit surface index required')
    choices={}
    polynomial_pairs=[]
    for mode,nodes,box in ((POINT,d['point'],d['point_box']),(BOX,d['nodes'],d['box'])):
        fields=ctx.target.surface(f,e,nodes=nodes)
        require(fields['global_indices']==expected['global_indices']==frozen['global_indices'],'surface map drift')
        E._surface_rule(ctx.target,d['rule'],f,e)
        trace=SUR.uniform_G_bounds(fields,max_depth=8)
        sites=[SUR.local_site_domain(ctx.target,box,f,e,k) for k in range(12)]
        if mode==BOX:
            require(trace==frozen['domain'] and sites==frozen['cache_site_geometry'],
                    'actual whole-box surface polynomials and cache geometries differ from independent review')
        try:
            order,attempts,trimmed,jet=surface_selector(d,fields,trace,f,e,tail)
        except OrderUnresolved as error:
            error.component_scope=mode
            raise
        polynomial_pairs.append((fields,jet))
        choices[mode]={'record_id':ctx.record_id,'face':f,'edge':e,'global_indices':fields['global_indices'],
                       'order':order,'domain':trimmed,'cache_site_domains':sites,'order_attempts':attempts,
                       'polynomial_jets_sha256':polynomials(fields,jet,('q_physical_R','B_physical_R')),
                       'node_ball_sha256':E._local_node_digest(nodes,expected['nodes'])}
    return {'kind':'surface','slot':slot,'global_indices':expected['global_indices'],
            'exact_point_coefficients_enclosed_by_whole_box_count':coefficientwise_inside(
                *polynomial_pairs,('q_physical_R','B_physical_R')),
            **choices}


@contextmanager
def planning_only(force,scaled):
    with ExitStack() as stack:
        forbidden=AssertionError('spatial moment or cached derivative calculation in selector phase')
        for obj,name in ((force,'_force_parts'),(scaled,'bulk_component'),(E,'bulk_component'),
                         (E,'surface_component'),(SUR,'run_edge'),(U,'uniform_cache'),
                         (E.P.Functional,'integral'),(E.P.Functional,'defect'),
                         (SUR.S.SurfaceFunctional,'integral'),(SUR.S.SurfaceFunctional,'defect'),
                         (C,'local_terms'),(C,'assemble_point')):
            stack.enter_context(patch.object(obj,name,side_effect=forbidden))
        yield


def packed_interval(x):
    b=DX.pair(x)
    require(all(v.denominator.bit_count()==1 for v in b),'dyadic Arb endpoints required')
    return b


def checked_F(entries,tails,n,cap):
    require(len(entries)==len(tails)==n and all(0<=DX.pair(t)[0]<=DX.pair(t)[1]<=cap for t in tails),
            'point force tail and ordered row list inconsistent')
    return [list(map(str,packed_interval(v))) for v in entries]


def checked_H(rows,n,cap):
    expected=set(combinations_with_replacement(range(n),2))
    require(len(rows)==len(expected) and {(r['i'],r['j']) for r in rows}==expected,
            'complete unique local upper Hessian required')
    require(all(0<=DX.pair(r['two_mass_tail_upper'])[0]<=DX.pair(r['two_mass_tail_upper'])[1]<=cap for r in rows),
            'local upper Hessian tails exceed frozen capsule')
    return [[r['i'],r['j'],*map(str,packed_interval(r['enclosure']))]
            for r in sorted(rows,key=lambda r:(r['i'],r['j']))]


def evaluate(d,plan,tail,force,scaled):
    boot()
    ctx,rule=d['ctx'],d['rule']
    i,n=plan['slot'],len(plan['global_indices'])
    if plan['kind']=='bulk':
        require(EX.doc(build_bulk(d,i,tail,scaled))==plan,'actual bulk point/box proof differs from frozen plan')
        p=plan[POINT]['selector']
        fields,ids=E._bulk_fields(ctx.target,d['point'],i)
        jet=E.JET.bulk_jacobian_jet(fields['q'],fields['B'])
        maps,Y,lines=point_maps(d,i,fields,ids)
        require(lines==plan[POINT]['same_Y_point_cache_lines'],'point same-Y convex cache line changed')
        functional=E._functional_type(fields)(ctx.target,fields,i,rule)
        try:
            out=force._force_parts(fields,jet,functional,list(map(C.A,Y)),C.A(ctx.drive*fields['active']),maps,p)
            point=out['normalized_force_total']
            point_tails=out['two_mass_tail_upper_exact']
            corrections=out['finite_same_Y_cache_correction']
            point_domain=out['cache_min_J_bound']
            require(DX.pair(point_domain)[0]>Q(1,10**6),'point same-Y force-domain threshold failed')
            require(len(point)==len(point_tails)==n and all(0<=v<=tail[0] for v in point_tails),
                    'exact force tail coverage required')
            packed_point=[list(map(str,packed_interval(C.pack(v)))) for v in point]
        finally:
            for name in ('integral','site_power','defect'):
                method=getattr(functional,name,None)
                if method is not None and hasattr(method,'cache_clear'):
                    method.cache_clear()
            E.M._triangle_moment.cache_clear()
            C.A.cache_clear()
        b=plan[BOX]
        whole=scaled.bulk_component(ctx,d['nodes'],rule,d['body']['triangles'][i],b['selector'],i,b['scaled_or_axis_proof'])
        cache={'point_domain':point_domain,'whole_box_cache':whole['cache_consistency']}
    else:
        require(EX.doc(build_surface(d,i,tail))==plan,'actual surface point/box proof differs from frozen plan')
        p=E.surface_component(ctx,d['point'],rule,plan[POINT])
        whole=E.surface_component(ctx,d['nodes'],rule,plan[BOX])
        packed_point=checked_F(p['force_integral_minus_cache'],p['force_two_mass_tail_upper'],n,tail[0])
        point_tails=[DX.pair(x)[1] for x in p['force_two_mass_tail_upper']]
        corrections=[C.A(DX.pair(v)[0]).union(C.A(DX.pair(v)[1])) for v in p['force_auxiliary_minus_cache']]
        cache={'point_site_count':len(p['actual_cache_sites']),'whole_box_site_count':len(whole['actual_cache_sites'])}
    require(whole['global_indices']==plan['global_indices'],'computed derivative indexes drifted')
    return {'kind':plan['kind'],'slot':i,'record_id':ctx.record_id,'precision_bits':640,
            'global_indices':plan['global_indices'],'proof_sha256':digest(plan),
            'point_force':{'scope':POINT,'F_I_minus_F_cache_intervals_exact':packed_point,
                           'force_tail_upper_exact':list(map(str,point_tails)),
                           'Qstar_minus_cache_correction_exact':[list(map(str,SUR.rational_hull(v))) for v in corrections]},
            'whole_box_hessian':{'scope':BOX,'radius_exact':str(RADIUS),
                                 'H_I_minus_H_cache_upper_exact':checked_H(whole['H_integral_minus_cache_upper'],n,tail[1])},
            'cache_domain_evidence':cache,'root_asserted':False}


def scatter(packets,expected_order,n):
    """Signed FE scatter, before any multiplication by R. Caller admits packet authority."""
    require([(p['kind'],p['slot']) for p in packets]==expected_order and len(set(expected_order))==len(expected_order),
            'complete ordered component packet coverage required')
    F=[(Q(0),Q(0)) for _ in range(n)]
    H={}
    add=lambda x,y:(x[0]+y[0],x[1]+y[1])
    for p in packets:
        f,h=p['point_force'],p['whole_box_hessian']
        ids=p['global_indices']
        require(f['scope']==POINT and h['scope']==BOX and Qx(h['radius_exact'])==RADIUS
                and len(set(ids))==len(ids) and all(type(j) is int and 0<=j<n for j in ids),
                'separate point/whole-box scopes and unique free projection required')
        for i,pair in zip(ids,f['F_I_minus_F_cache_intervals_exact'],strict=True):
            v=tuple(map(Qx,pair))
            require(len(v)==2 and v[0]<=v[1],'ordered signed force interval required')
            F[i]=add(F[i],v)
        upper=h['H_I_minus_H_cache_upper_exact']
        require([(i,j) for i,j,*_ in upper] == list(combinations_with_replacement(range(len(ids)),2)),
                'missing, duplicated or unordered upper Hessian monomial')
        for i,j,lo,hi in upper:
            key=tuple(sorted((ids[i],ids[j])))
            value=(Qx(lo),Qx(hi))
            require(value[0]<=value[1],'ordered signed Hessian interval required')
            H[key]=add(H.get(key,(Q(0),Q(0))),value)
    return F,H


def plan_gate(review_path,review_sha):
    review=pinned(review_path,review_sha)
    bindings={'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),
              'preflight_sha256':sha(PREP/'preflight.json'),'geometry_review_sha256':DX.PINS[DX.GEOMETRY_REVIEW],
              'cached_root_review_sha256':EX.PINS[EX.CACHE_REVIEW]}
    require(review['status']=='PASS' and sweep.all_true(review['checks'])
            and all(review[k]==v for k,v in bindings.items())
            and sha(PREP/'transfer_ambiguity_exact_bridge_frozen.py')==bindings['source_sha256'],
            'exact singleton bridge source prelaunch review required')
    p=pinned(PREP/'preflight.json',bindings['preflight_sha256'])
    require(p['status']=='PASS' and sweep.all_true(p['checks']) and p['selector_plans_run']==0
            and p['protocol_sha256']==bindings['protocol_sha256'] and p['source_sha256']==bindings['source_sha256'],
            'frozen cold preflight changed')
    state=static_inputs()
    require(authority_pins(state)==pinned(PREP/'protocol.json',bindings['protocol_sha256']),
            'exact singleton bridge protocol changed')
    return state,{**bindings,'source_prelaunch_review_sha256':review_sha}


def build_plan(args):
    state,b=plan_gate(args.review,args.review_sha256)
    derivative=derivative_inputs(args.derivative_admission,args.derivative_admission_sha256,state[1])[args.pilot_id]
    d=domain(state,args.pilot_id)
    folder=PLAN/args.pilot_id
    folder.mkdir(parents=True,exist_ok=False)
    sizing=derivative['planning']
    bindings={**b,'derivative_admission_sha256':args.derivative_admission_sha256,
              'independent_native_derivative_review_sha256':derivative['independent_native_sha256']}
    output_pid=os.getpid()
    save(folder/'started.json',{'pilot_id':args.pilot_id,'bindings':bindings,'pid':output_pid})
    if sizing['strict_positive_headroom'] is not True:
        h=save(folder/'unresolved.json',{'status':'UNRESOLVED_NONPOSITIVE_CACHED_HEADROOM','derivative':derivative,
                                        'bindings':bindings,'selector_or_component_launched':False})
        save(folder/'terminal.json',{'status':'UNRESOLVED_NONPOSITIVE_CACHED_HEADROOM','pid':output_pid,
                                     'exit_code':0,'pilot_id':args.pilot_id,'unresolved_sha256':h})
        return
    tail=tuple(Qx(sizing[k]) for k in ('per_component_force_entry_tail_cap_exact','per_component_H_entry_tail_cap_exact'))
    force,scaled=kernels({})
    manifest=[]
    try:
        with planning_only(force,scaled):
            for kind,count in (('bulk',144),('surface',24)):
                for slot in range(count):
                    rel=f'{kind}_{slot:03d}.json'
                    try:
                        value=(build_bulk(d,slot,tail,scaled) if kind=='bulk' else build_surface(d,slot,tail))
                        h=save(folder/rel,EX.doc(value))
                        manifest.append({'kind':kind,'slot':slot,'file':rel,'sha256':h})
                    except (E.TailBudgetUnresolved,OrderUnresolved) as error:
                        save(folder/f'{kind}_{slot:03d}_unresolved.json',{'kind':kind,'slot':slot,'attempts':error.attempts,
                                    'component_scope':error.component_scope,
                                    'bindings':bindings,'status':'UNRESOLVED_SELECTOR_ORDER'})
                        raise
        plan={'status':'COMPLETE_POINT_FORCE_WHOLE_H_SELECTOR_PLAN_UNREVIEWED','pilot_id':args.pilot_id,
              'record_id':d['ctx'].record_id,'scope':SCOPE,'point_membership':POINT,'H_membership':BOX,
              'whole_box_radius_exact':str(RADIUS),'context_box_sha256':d['entry']['box_1_sha256'],
              'point_box_sha256':digest(d['point_box']),'tail_caps_exact':list(map(str,tail)),
              'derivative':derivative,'components':manifest,'bindings':bindings,
              'spatial_component_launch_authorized':False}
        sha_plan=save(folder/'plan.json',plan)
        save(folder/'terminal.json',{'status':'FINISHED','pid':output_pid,'exit_code':0,
                                     'plan_sha256':sha_plan,'pilot_id':args.pilot_id})
        print(json.dumps({'plan_sha256':sha_plan,'components':len(manifest)}),flush=True)
    except BaseException as error:
        save(folder/'failure.json',{'error':str(error),'traceback':traceback.format_exc(),'completed':manifest})
        save(folder/'terminal.json',{'status':'FAILED','pid':output_pid,'exit_code':1,'pilot_id':args.pilot_id})
        raise


def launch_gate(args):
    state,b=plan_gate(args.review,args.review_sha256)
    dr=derivative_inputs(args.derivative_admission,args.derivative_admission_sha256,state[1])[args.pilot_id]
    review=pinned(args.plan_review,args.plan_review_sha256)
    links={**b,'derivative_admission_sha256':args.derivative_admission_sha256,
           'independent_native_derivative_review_sha256':dr['independent_native_sha256']}
    folder=PLAN/args.pilot_id
    plan_path=folder/'plan.json'
    actual_plan_sha=sha(plan_path)
    terminal_sha=sha(folder/'terminal.json')
    require(review['status']=='PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_PLAN' and sweep.all_true(review['checks'])
            and review['pilot_id']==args.pilot_id and review['plan_sha256']==actual_plan_sha
            and review['source_sha256']==b['source_sha256'] and review['plan_terminal_sha256']==terminal_sha
            and review['derivative_admission_sha256']==args.derivative_admission_sha256
            and type(review['plan_observed_os_exit_code']) is int and review['plan_observed_os_exit_code']==0,
            'current own independently reviewed and OS-complete actual plan required')
    end=pinned(folder/'terminal.json',terminal_sha)
    plan=pinned(plan_path,actual_plan_sha)
    order=[(k,i) for k,n in (('bulk',144),('surface',24)) for i in range(n)]
    require(end['status']=='FINISHED' and end['exit_code']==0 and end['pid']==review['plan_observed_pid']
            and end['plan_sha256']==actual_plan_sha and end['pilot_id']==args.pilot_id
            and plan['status']=='COMPLETE_POINT_FORCE_WHOLE_H_SELECTOR_PLAN_UNREVIEWED'
            and plan['bindings']==links and plan['derivative']==dr and plan['pilot_id']==args.pilot_id
            and [(p['kind'],p['slot']) for p in plan['components']]==order
            and all(p['file']==f"{p['kind']}_{p['slot']:03d}.json" for p in plan['components']),
            'changed or partial 168-component proof manifest')
    require(all(sha(folder/r['file'])==r['sha256'] for r in plan['components']),
            'component plan bytes must be authenticated before payload reads')
    return state,plan,{**links,'plan_sha256':actual_plan_sha,'plan_review_sha256':args.plan_review_sha256}


def worker(args):
    state,plan,links=launch_gate(args)
    run=ACTUAL/args.pilot_id
    parent=sweep.read(run/'started.json')
    require(parent['pid']==args.supervisor_pid==os.getppid() and parent['bindings']==links,'authorized live supervisor required')
    folder=run/args.part
    own={'pid':os.getpid(),'pilot_id':args.pilot_id,'part':args.part,'bindings':links,'supervisor_pid':args.supervisor_pid}
    save(folder/'started.json',own)
    terminal={**own,'status':'FAILED'}
    files=[]
    try:
        d=domain(state,args.pilot_id)
        force,scaled=kernels({})
        tail=tuple(map(Qx,plan['tail_caps_exact']))
        for p in plan['components']:
            if p['kind']!=args.part:
                continue
            value=pinned(PLAN/args.pilot_id/p['file'],p['sha256'])
            require(value['kind']==args.part and value['slot']==p['slot'],'single explicit component index required')
            result=evaluate(d,value,tail,force,scaled)
            h=save(folder/p['file'],{**result,'bindings':links,'plan_file_sha256':p['sha256']})
            files.append({'file':EX.relative(folder/p['file']),'sha256':h,'kind':p['kind'],'slot':p['slot']})
            print(json.dumps({'part':args.part,'completed':len(files),'slot':p['slot']}),flush=True)
        require(len(files)=={'bulk':144,'surface':24}[args.part],'incomplete component set')
        terminal['manifest_sha256']=save(folder/'manifest.json',{'files':files,**own})
        terminal['status']='FINISHED'
    except BaseException as e:
        save(folder/'failure.json',{'error':str(e),'traceback':traceback.format_exc(),'completed':files})
        raise
    finally:
        save(folder/'terminal.json',terminal)


def launch(args):
    _,_,links=launch_gate(args)
    run=ACTUAL/args.pilot_id
    run.mkdir(parents=True,exist_ok=False)
    pid=os.getpid()
    save(run/'started.json',{'pid':pid,'pilot_id':args.pilot_id,'bindings':links})
    terminal={'pid':pid,'pilot_id':args.pilot_id,'status':'FAILED','bindings':links,'workers':[]}
    files=[]
    try:
        for part in ('bulk','surface'):
            folder=run/part
            folder.mkdir(exist_ok=False)
            comm=[sys.executable,'-B',str(SOURCE),'worker','--pilot-id',args.pilot_id,'--part',part,'--supervisor-pid',str(pid)]
            for name in ('review','review_sha256','derivative_admission','derivative_admission_sha256',
                         'plan_review','plan_review_sha256'):
                comm+=['--'+name.replace('_','-'),str(getattr(args,name))]
            with (folder/'stdout.log').open('x',encoding='utf-8') as log:
                child=subprocess.Popen(comm,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                save(folder/'supervisor_child_started.json',{'pid':child.pid,'supervisor_pid':pid,'command':comm})
                code=child.wait()
            path=folder/'terminal.json'
            observed={'pid':child.pid,'supervisor_pid':pid,'observed_exit_code':code,'part':part,
                      'terminal_sha256':sha(path) if path.is_file() else None}
            h=save(folder/'observed_os_exit.json',observed)
            terminal['workers'].append({**observed,'os_exit_receipt_sha256':h})
            require(observed['terminal_sha256'] is not None,'child terminal missing')
            ended=pinned(path,observed['terminal_sha256'])
            identity={'pid':child.pid,'supervisor_pid':pid,'pilot_id':args.pilot_id,'part':part,'bindings':links}
            DX.CP.admit_child(code,ended,identity)
            require(sweep.read(folder/'started.json')==identity,'actual worker start identity drifted')
            manifest=pinned(folder/'manifest.json',ended['manifest_sha256'])
            require(all(manifest[k]==v for k,v in identity.items()),'mixed component part manifest')
            files.extend(manifest['files'])
        require([(x['kind'],x['slot']) for x in files]
                ==[(k,i) for k,n in (('bulk',144),('surface',24)) for i in range(n)]
                and all(sha(ROOT/r['file'])==r['sha256'] for r in files),'complete immutable 168-component manifest required')
        terminal['manifest_sha256']=save(run/'manifest.json',{'status':'COMPLETE_168_SINGLETON_BRIDGE_PACKETS_UNREVIEWED',
                            'bindings':links,'pilot_id':args.pilot_id,'files':files,
                            'independent_actual_component_review_required':True,'root_claimed':False})
        terminal['status']='FINISHED'
        print(json.dumps({'manifest_sha256':terminal['manifest_sha256'],'packets':168}),flush=True)
    except BaseException as error:
        terminal['error']=str(error)
        raise
    finally:
        save(run/'terminal.json',terminal)


def fixtures(state):
    fail=AssertionError('selector construction is forbidden in cold preparation')
    force,scaled=kernels({})
    checks=[]
    with (planning_only(force,scaled), patch.object(E,'_select',side_effect=fail),
          patch.object(SUR,'uniform_G_bounds',side_effect=fail)):
        for pid in EX.PILOTS:
            d=domain(state,pid)
            require(d['ctx'].drive==Qx(d['case']['actual_exact_drive']) and len(d['point'])==325,
                    'bad exact point/box route in cold preflight')
            checks.append(d)
    # Two matched rows cancel on purpose. Point and whole-box roles are enforced by the scatter.
    def packet(i,x):
        return {'kind':'bulk','slot':i,'global_indices':[0],
                'point_force':{'scope':POINT,'F_I_minus_F_cache_intervals_exact':[[str(x),str(x)]]},
                'whole_box_hessian':{'scope':BOX,'radius_exact':str(RADIUS),'H_I_minus_H_cache_upper_exact':[[0,0,'1','1']]}}
    a,b=packet(0,2),packet(1,-1)
    ordered=[('bulk',0),('bulk',1)]
    F,H=scatter([a,b],ordered,1)
    wrong={**a,'point_force':{**a['point_force'],'scope':BOX}}
    bad=list(checks[0]['point'])
    bad[0]=[Q(7),Q(7)]
    point_polys={'q':[{(0,0):Q(1)}],'B':[[{(0,0):Q(1)}]]}
    jet={'value':{(0,0):Q(1)},'gradient':[{}],'hessian_upper':{(0,0):{}}}
    whole=dict(point_polys,q=[{(0,0):C.A(0).union(C.A(2))}])
    disjoint=dict(point_polys,q=[{(0,0):Q(2)}])
    return {'four_separate_actual_point_and_box_bindings_verified':len(checks)==4,
            'center_must_be_inside_whole_box_and_equal_exact_frozen_center':
                EX.reject(lambda:exact_point_inside(bad,checks[0]['box'])),
            'signed_scatter_before_norm_or_R_with_symmetric_upper_H':F==[(Q(1),Q(1))] and H=={(0,0):(Q(2),Q(2))},
            'actual_polynomial_coefficient_containment_requires_true_interval_inclusion':
                (coefficientwise_inside((point_polys,jet),(whole,jet),('q','B'))==3
                 and EX.reject(lambda:coefficientwise_inside((point_polys,jet),(disjoint,jet),('q','B')))),
            'reject_partial_duplicate_or_mixed_force_scope':(EX.reject(lambda:scatter([a],ordered,1))
                 and EX.reject(lambda:scatter([a,a],ordered,1)) and EX.reject(lambda:scatter([wrong,b],ordered,1))),
            'no_actual_selectors_or_spatial_jobs':not PLAN.exists() and not ACTUAL.exists()}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/'transfer_ambiguity_exact_bridge_frozen.py').open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        state=static_inputs()
        checks=fixtures(state)
        require(sweep.all_true(checks),'cold singleton bridge preparation failed')
        ph=save(PREP/'protocol.json',authority_pins(state))
        pf=save(PREP/'preflight.json',{'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,
                   'checks':checks,'selector_plans_run':0,'new_spatial_or_cached_AD_jobs':0})
        print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':pf}),flush=True)
    except BaseException as error:
        save(PREP/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='stage',required=True)
    sub.add_parser('prepare')
    for mode in ('plan','launch','worker'):
        p=sub.add_parser(mode)
        p.add_argument('--pilot-id',choices=EX.PILOTS,required=True)
        names=['review','derivative-admission']+(['plan-review'] if mode!='plan' else [])
        for name in names:
            p.add_argument('--'+name,required=True,type=lambda s: Path(s).resolve())
            p.add_argument('--'+name+'-sha256',required=True)
        if mode=='worker':
            p.add_argument('--part',choices=('bulk','surface'),required=True)
            p.add_argument('--supervisor-pid',type=int,required=True)
    args=cli.parse_args()
    {'prepare':lambda _:prepare(),'plan':build_plan,'launch':launch,'worker':worker}[args.stage](args)

