"""Extend four frozen cached-quadrature derivative enclosures to r=2^-34."""
import argparse
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic
import traceback
from unittest.mock import patch

import numpy as np

import algebraic_root_feasibility as core
import transfer_ambiguity_cached_pilots as CP
import transfer_ambiguity_exact_targets as EX
import transfer_ambiguity_observer_sweep as sweep


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
BASE = EX.relative(EX.OUT)
GEOM = f'{BASE}/geometry'
REVIEW_DIR = f'{BASE}/independent_actual/replay'
DESIGN = f'{BASE}/force_bridge_design'
GEOMETRY_REVIEW = f'{REVIEW_DIR}/receipt.json'
GEOMETRY_END = f'{REVIEW_DIR}/terminal.json'
GEOMETRY_OS = f'{REVIEW_DIR}/os_exit_observation.json'
DESIGN_RECEIPT = f'{DESIGN}/receipt.json'
OUT = EX.OUT / 'cached_derivative_extension_v1'
PREP = OUT / 'preparation'
ACTUAL = OUT / 'actual'
RADIUS = Q(1, 2**34)
OLD_RADIUS = Q(1, 2**36)
PRECISION = 256
NFREE = 611
PINS = {
    'scripts/algebraic_root_feasibility.py': '7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
    'scripts/transfer_ambiguity_cached_pilots.py': '1b22ebc307c0cee51bf4193ff3999242b1cfe278d07c10cc63619f2795489c58',
    'scripts/transfer_ambiguity_exact_targets.py': '53ff92a87d164871e04323759b27e015fc8e19d039dc0dd15c2d6080acda7cc8',
    f'{BASE}/preparation/protocol.json': '28efd2467bbeae4da32797e0405b9a9595cdbb03ded9ce94221f18592305687c',
    f'{BASE}/preparation/target_manifest.json': 'cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3',
    f'{BASE}/preparation/preflight.json': 'a4dacbab0dafb730db12321892d741979fea4654c436996b964cc9cfde1b0808',
    f'{BASE}/preparation/prelaunch_independent_review.json': '17c0e8feb5f91898787994ba84bd1afd65e4be5e8ff52f01b0f5a48853f3ce6b',
    f'{GEOM}/results.json': 'fabf695e7f83189b2a088570723f1c103cafa7c5e973ac704f9c7adc25b9bdc9',
    f'{GEOM}/terminal.json': '060ffa111367f95bbb7eccf17084e888ccf6e0036edae6f8f4246ebd47bccad1',
    GEOMETRY_REVIEW: '884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    GEOMETRY_END: '5011dc016786a4c3fd2a184cd4c620ad10e9d374d25d8d7ae02aac2d06d08c14',
    GEOMETRY_OS: 'c2be4be8b403e170e20ec54c34f17a69b4aec4299e6f9f1ab7ecb3d4afa337d9',
    DESIGN_RECEIPT: 'c1a11d4d54fb33f7a0675448eb9e1c3f4ff4a3b772a60d83e656918aa912b1e3',
    f'{DESIGN}/README.md': '8092734ed018004fbbeef23c9f1d3c5110029ee5282440d5eb0801664600bf54',
    EX.CACHE_REVIEW: EX.PINS[EX.CACHE_REVIEW], EX.CACHE_ADMISSION: EX.PINS[EX.CACHE_ADMISSION],
}
SCOPE = ('Fresh whole-radius derivative defect of the exact-real binary quadrature potential at four '
         'different fixed physical parameter tuples. Same per-point center, preconditioner, quadrature arrays and '
         'requested rational drive as the old cached proof. This phase supplies neither an exact spatial force/Hessian '
         'nor its true defect or bias. Planning caps concern a future exact-spatial bridge.')
require, save, read, sha, pinned = sweep.require, sweep.save, sweep.read, sweep.sha, sweep.pinned


def bound(path, links, digest, parse=True):
    EX.link(links,path,digest)
    return pinned(ROOT/path,digest) if parse else None


def pair(packed):
    lo, hi = (EX.exact(packed[k]) for k in ('lower_exact','upper_exact'))
    require(lo <= hi,'ordered finite exact bound required')
    return lo, hi


def max_abs(packed):
    return max(map(abs,pair(packed)))


def validate_center(proof, case):
    eta = EX.exact(proof['independent_eta_upper_exact'])
    z0 = EX.exact(proof['independent_Z0_upper_exact'])
    rr = EX.exact(proof['preconditioner_infinity_norm_exact'])
    require(proof['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
            and proof['pilot_id']==case['pilot_id'] and proof['field_sha256']==case['field_sha256']
            and proof['cache_packet_sha256']==case['cache_packet_sha256']
            and proof['cache_array_contents_sha256']==case['embedded_cache_sha256']
            and EX.exact(proof['exact_strength'])==EX.exact(case['actual_exact_drive'])
            and proof['rounded_strength_used_in_proof'] is False
            and EX.exact(proof['radius_exact'])==OLD_RADIUS
            and len(proof['signed_center_force_intervals'])==NFREE
            and 0 <= max(map(max_abs,proof['signed_center_force_intervals'])) <= eta
            and 0 <= z0 < 1 and rr > 0,
            'certified same-center scalar and signed cached rows required')
    return eta,z0,rr


def authority():
    links={}
    for path,h in PINS.items():
        bound(path,links,h,parse=False)
    require(Path(core.__file__).resolve() == ROOT/'scripts/algebraic_root_feasibility.py'
            and EX.T.core is core, 'same explicit cache algebra required')
    initial, manifest, contexts = EX.prepared()
    # This replays the 23 terminal and four independent cached proof admissions,
    # plus the native nonlinear/budget and singleton authorities, before field loads.
    for path,h in initial['input_hashes'].items():
        EX.link(links,path,h)
    review = pinned(ROOT/GEOMETRY_REVIEW,PINS[GEOMETRY_REVIEW])
    terminal = pinned(ROOT/GEOMETRY_END,PINS[GEOMETRY_END])
    os_receipt = pinned(ROOT/GEOMETRY_OS,PINS[GEOMETRY_OS])
    produced = pinned(ROOT/f'{GEOM}/results.json',PINS[f'{GEOM}/results.json'])
    producer_end = pinned(ROOT/f'{GEOM}/terminal.json',PINS[f'{GEOM}/terminal.json'])
    require(review['status']=='PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES'
            and sweep.all_true(review['checks']) and len(review['geometry_records'])==8
            and review['exact_spatial_root_count']==review['spatial_derivative_or_tail_bound_count']==0
            and terminal['status']=='FINISHED' and type(terminal['exit_code']) is int and terminal['exit_code']==0
            and os_receipt['status']=='PASS_OBSERVED_INDEPENDENT_GEOMETRY_REPLAY_TRUE_OS_EXIT'
            and type(os_receipt['observed_exit_code']) is int and os_receipt['observed_exit_code']==0
            and terminal['receipt_sha256']==os_receipt['receipt_sha256']==PINS[GEOMETRY_REVIEW]
            and terminal['pid']==os_receipt['pid']==review['source_and_phase']['pid']
            and os_receipt['terminal_sha256']==PINS[GEOMETRY_END]
            and producer_end['status']=='FINISHED' and type(producer_end['exit_code']) is int
            and producer_end['exit_code']==0 and review['observed_producer_os']['exit_code']==0
            and producer_end['pid']==review['observed_producer_os']['pid']
            and producer_end['results_sha256']==PINS[f'{GEOM}/results.json']
            and produced['force_or_H_jobs_authorized'] is False,
            'reviewed true terminal fresh geometry and precise geometry-only scope required')
    for path,h in review['checked_file_hashes'].items():
        EX.link(links,path,h)
    memo = pinned(ROOT/DESIGN_RECEIPT,PINS[DESIGN_RECEIPT])
    require(memo['status']=='READ_ONLY_FORCE_BRIDGE_DESIGN_COMPLETE'
            and memo['memo_sha256']==PINS[f'{DESIGN}/README.md']
            and EX.exact(memo['proposed_r_exact'])==RADIUS
            and EX.exact(memo['prior_cached_derivative_radius_exact'])==OLD_RADIUS
            and memo['larger_radius_requires_fresh_cached_derivative_enclosure'] is True,
            'separate planning-only force bridge design required')
    caps = {x['pilot_id']:x for x in memo['optional_tail_planning_only']}
    require(tuple(caps)==EX.PILOTS and tuple(x['pilot_id'] for x in manifest['targets'])==EX.PILOTS,
            'same ordered four current singleton identities required')
    accepted=[]
    for entry in manifest['targets']:
        pid = entry['pilot_id']
        p, result, independent, proof_path = contexts[pid]
        case_path = f'{GEOM}/{pid}/box_1/geometry.json'
        g, = [x for x in produced['geometry_cases'] if x['file']==case_path]
        independent_geo, = [x for x in review['geometry_records'] if x['pilot_id']==pid
                             and EX.exact(x['radius_exact'])==RADIUS]
        geometry = bound(case_path,links,g['sha256'])
        geo_proof = bound(independent_geo['receipt'],links,independent_geo['sha256'])
        proposed_box = bound(entry['box_1'],links,entry['box_1_sha256'])
        predict = bound(entry['model_contract'],links,entry['model_contract_sha256'])
        routes = bound(entry['routes'],links,entry['routes_sha256'])
        require(geometry['status']=='PASS_FRESH_GEOMETRY_ONLY_UNREVIEWED'
                and geo_proof['status']=='PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY'
                and geometry['pilot_id']==geo_proof['pilot_id']==pid
                and g['pilot_id']==pid
                and geometry['box_sha256']==independent_geo['box_sha256']==geo_proof['box_sha256']==entry['box_1_sha256']
                and geometry['body_geometry_sha256']==geo_proof['body_geometry_sha256']
                and geometry['model_contract_sha256']==entry['model_contract_sha256']
                and geo_proof['proof_center_field_sha256']==p['field_sha256']==predict['source_field_sha256']
                and EX.exact(geometry['radius_exact'])==EX.exact(proposed_box['corrector_radius_exact'])==RADIUS
                and geometry['new_cached_AD_assembly'] is False
                and geo_proof['larger_cached_derivative_proof_reused'] is False
                and predict['larger_radius_cached_derivative_asserted'] is False
                and proposed_box['root_membership_assumed'] is False,
                'exact current-center larger box required without extrapolated derivative premise')
        f=[0]*NFREE
        hh=[0]*NFREE
        for component in routes['bulk']+routes['surface']:
            projection = component['global_indices']
            for i in projection:
                f[i]+=1
                hh[i]+=len(projection)
        sizing = caps[pid]
        require((max(f),max(hh))==(sizing['max_force_incidence'],sizing['max_H_entry_row_incidence'])==(6,72),
                'same actual 168-route force and Hessian incidence counts required')
        row = {'pilot_id':pid,'record_id':predict['record_id'],'field':p['field'],'field_sha256':p['field_sha256'],
               'cache_packet':result['cache_packet'],'cache_packet_sha256':result['cache_packet_sha256'],
               'embedded_cache_sha256':predict['frozen_binary_cache_sha256'],
               'model_contract_sha256':entry['model_contract_sha256'],
               'preconditioner_binary_sha256':predict['binary_preconditioner_sha256'],
               'geometry':case_path,'geometry_sha256':g['sha256'],
               'whole_box_geometry_review':independent_geo['receipt'],
               'whole_box_geometry_review_sha256':independent_geo['sha256'],
               'box':entry['box_1'],'box_sha256':entry['box_1_sha256'],
               'verified_point':proof_path,'verified_point_sha256':independent['receipt_sha256'],
               'actual_exact_drive':p['exact_drive']['activation_strength_exact'],
               'rounded_solver_strength':p['exact_drive']['activation_strength_binary_exact'],
               'cached_proof_radius_exact':str(OLD_RADIUS),'new_derivative_radius_exact':str(RADIUS),
               'observer_goal_exact':p['e_goal_dimensionless_exact'],'prior_optional_tail_policy':sizing}
        point=bound(proof_path,links,independent['receipt_sha256'])
        eta,z0,rr=validate_center(point,row)
        require(OLD_RADIUS<RADIUS<EX.exact(row['observer_goal_exact'])
                and EX.exact(row['actual_exact_drive'])!=EX.exact(row['rounded_solver_strength'])
                and EX.exact(sizing['frozen_R_infinity_norm_upper_exact'])>=rr,
                'strict radius science budget, exact drive or R norm coverage changed')
        row.update({'point_eta_upper_exact':str(eta),'center_Z0_upper_exact':str(z0),'R_inf_exact':str(rr)})
        accepted.append(row)
    return dict(sorted(links.items())),accepted,manifest,contexts


def protocol(record):
    bindings, cases, _, _ = record
    return {'status':'FROZEN_FOUR_CACHED_DERIVATIVE_EXTENSIONS','scope':SCOPE,
            'source_sha256':sha(SOURCE),'input_sha256':bindings,'ordered_cases':cases,
            'old_cached_certificate_radius_exact':str(OLD_RADIUS),'larger_radius_exact':str(RADIUS),
            'precision_bits':PRECISION,'threads':1,'whole_box_site_counts':CP.TERMS,
            'optional_engineering_ZQ_threshold_exact':'1/4','threshold_is_gate':False,
            'center_and_R_preserved':True,'exact_requested_drive_required':True,
            'planning_rule':'h=1-ZQ; if h>0 use b*=h/2 and signed total force allowance (h-b*)r/2. '
                            'Local tail sizing: one ten-thousandth of these prospective caps through the verified incidence and R norm, '
                            'at most the legacy optional local caps. No exact FE force defect or bias is inferred.',
            'execution':'Four serial fresh subprocesses, with write-once true OS terminal and output receipts.',
            'new_spatial_integrals':0,'new_equilibrium_solves':0,'new_preconditioners':0}


def spend(zq,radius,norm,previous):
    require(type(zq) is type(radius) is type(norm) is Q and zq>=0 and radius>0 and norm>0
            and (previous['max_force_incidence'],previous['max_H_entry_row_incidence'])==(6,72),
            'exact valid bound, radius, R norm and verified incidence required')
    h=1-zq
    row={'headroom_exact':str(h),'below_optional_quarter_reference':zq<=Q(1,4),
         'strict_positive_headroom':h>0,'quarter_reference_is_gate':False}
    if h<=0:
        return {**row,'planning_status':'NO_POSITIVE_HEADROOM','b_star_exact':None,
                'preconditioned_signed_total_force_allowance_exact':None}
    b=h/2
    f=(h-b)*radius/2
    share=Q(1,10000)
    gf=share*f/(norm*previous['max_force_incidence'])
    gh=share*b/(norm*previous['max_H_entry_row_incidence'])
    return {**row,'planning_status':'POSITIVE_CONDITIONAL_HEADROOM',
            'b_star_exact':str(b),'preconditioned_signed_total_force_allowance_exact':str(f),
            'unspent_strict_selfmap_margin_at_both_planning_caps_exact':str((h-b)*radius-f),
            'max_force_incidence':previous['max_force_incidence'],
            'max_H_entry_row_incidence':previous['max_H_entry_row_incidence'],
            'per_component_force_entry_tail_cap_exact':str(min(gf, EX.exact(previous['tail_force_budget_exact']))),
            'per_component_H_entry_tail_cap_exact':str(min(gh, EX.exact(previous['tail_H_entry_budget_exact']))),
            'optional_sizing_only':True}


def abs_matrix(inverse):
    require(inverse.ndim==2 and inverse.shape[0]==inverse.shape[1] and np.isfinite(inverse).all(),
            'square finite selected binary R required')
    n=inverse.shape[0]
    return core.arb_mat(n,n,[core.A(abs(float(x))) for x in inverse.ravel()])


def propagate(abs_R,variation,n):
    require(abs_R.nrows()==abs_R.ncols()==variation.nrows()==n and variation.ncols()==1
            and all(pair(core.pack(variation[i,0]))[0]>=0 for i in range(n)),
            'preconditioned variation requires a nonnegative whole-row majorant')
    worst=abs_R*variation
    ends=[max_abs(core.pack(worst[i,0])) for i in range(n)]
    return max(ends), ends


def calculate(case,manifest,contexts,folder):
    import torch
    torch.set_num_threads(1)
    ctx, _, entry = EX.open_context(case['pilot_id'],manifest,contexts)
    box=pinned(ROOT/case['box'],case['box_sha256'])
    require(box==EX.exact_node_box(ctx,RADIUS) and entry['model_contract_sha256']==case['model_contract_sha256']
            and EX.array_sha(ctx.R)==case['preconditioner_binary_sha256']
            and ctx.drive==EX.exact(case['actual_exact_drive']),
            'explicit target center, requested strength or fixed R drift')
    packet=pinned(ROOT/case['cache_packet'],case['cache_packet_sha256'])
    core.ctx.prec=PRECISION
    core.A.cache_clear()
    require(core.ctx.prec==PRECISION and torch.get_num_threads()==1, 'one CPU thread and Arb256 required')
    restored,u,inverse=core.restore_cache(packet)
    require(np.array_equal(u,ctx.u) and np.array_equal(inverse,ctx.R)
            and np.array_equal(restored.free,ctx.free) and packet['frozen_binary_cache']['sha256']==case['embedded_cache_sha256'],
            'restored fresh cache does not match the declared singleton')
    point=pinned(ROOT/case['verified_point'],case['verified_point_sha256'])
    eta,z0,rr=validate_center(point,case)
    strength=EX.exact(case['actual_exact_drive'])
    (_, _, cache, min_point_J), counts=CP.point_assembly(core,restored,u,strength,NFREE,CP.TERMS)
    require(pair(core.pack(min_point_J))[0]>Q(1,10**6), 'positive cached point J required')
    evidence={'site_counts':counts,'point_lower_quadrature_J':core.pack(min_point_J),
               'proof_drive_exact':str(strength),'reused_center_eta_exact':str(eta),
               'reused_center_Z0_upper_exact':str(z0),'exact_R_inf':str(rr),
               'independent_signed_anchor_file':case['verified_point'],
               'independent_signed_anchor_sha256':case['verified_point_sha256'],
               'center_cache_rebuilt_at_precision_bits':PRECISION}
    save(folder/'center_and_analytic_cache.json',evidence)
    try:
        variation,min_j,domain=core.variation_rows(cache,NFREE,core.A(RADIUS))
    except ValueError as error:
        return {'status':'UNRESOLVED_FRESH_CACHED_ANALYTIC_DOMAIN','case':case,'scope':SCOPE,
                'center':evidence,'reason':str(error),'ZQ_upper_exact':None,'planning':None}
    require(domain['term_counts']==CP.TERMS and set(domain['strict_lower_margins'])=={
        'bulk_det_min','bulk_hoop_min','surface_radius_min','surface_tangent_squared_min','normalized_weight_min'},
        'all actual cached analytic domains required')
    if pair(core.pack(min_j))[0]<=Q(1,10**6) or any(pair(v)[0]<=0 for v in domain['strict_lower_margins'].values()):
        return {'status':'UNRESOLVED_FRESH_CACHED_ANALYTIC_DOMAIN','case':case,'scope':SCOPE,
                'center':evidence,'domain':domain,'ZQ_upper_exact':None,'planning':None}
    var, row_bounds=propagate(abs_matrix(inverse),variation,NFREE)
    zq=z0+var
    verdict=CP.contraction(eta,zq,RADIUS,EX.exact(case['observer_goal_exact']))
    return {'status':'FRESH_LARGER_RADIUS_CACHED_DERIVATIVE_UNREVIEWED','case':case,'scope':SCOPE,
            'center':evidence,'quadrature_positive_domain':domain,'minimum_box_quadrature_J':core.pack(min_j),
            'whole_box_cache_radius_exact':str(RADIUS),'preconditioned_variation_upper_exact':str(var),
            'whole_row_variation_upper_exact':[str(max_abs(core.pack(variation[i,0]))) for i in range(NFREE)],
            'preconditioned_row_variation_upper_exact':list(map(str,row_bounds)),
            'ZQ_upper_exact':str(zq),'whole_box_controls_only_cached_quadrature':True,
            'auxiliary_cached_selfmap_check':verdict,'planning':spend(zq,RADIUS,rr,case['prior_optional_tail_policy']),
            'exact_spatial_derivatives_computed':False,'exact_spatial_root_certified':False,
            'independent_derivative_review_required':True}


def fixtures(cases):
    check=cases[0]['prior_optional_tail_policy']
    z=Q(1,3)
    budget=spend(z,Q(1),Q(1),check)
    left=spend(Q(1),Q(1),Q(1),check)
    fail=spend(Q(5,4),Q(1),Q(1),check)
    core.ctx.prec=PRECISION
    core.A.cache_clear()
    inv = np.array([[-2.,1.],[3.,0.]])
    var=core.arb_mat(2,1,[core.A(Q(1,8)),core.A(Q(1,2))])
    mag,rows=propagate(abs_matrix(inv),var,2)
    tie=CP.contraction(Q(1,4),Q(3,4),Q(1),Q(2))
    center = pinned(ROOT/cases[0]['verified_point'],cases[0]['verified_point_sha256'])
    other = cases[-1]
    return {'all_four_own_center_R_and_large_geometry_admitted':len(cases)==4 and
               all(EX.exact(c['new_derivative_radius_exact'])==RADIUS for c in cases),
            'unsigned_R_propagates_nonnegative_whole_row_bounds':mag==Q(3,4) and rows==[Q(3,4),Q(3,8)],
            'headroom_plan_uses_actual_bound_even_above_quarter':(
               budget['b_star_exact']=='1/3' and budget['preconditioned_signed_total_force_allowance_exact']=='1/6'
               and budget['below_optional_quarter_reference'] is False),
            'zero_and_negative_headroom_retain_value_without_planning':
               left['headroom_exact']=='0' and fail['headroom_exact']=='-1/4' and
               left['b_star_exact'] is fail['b_star_exact'] is None,
            'strict_cached_selfmap_rejects_equality':not tie['strict_selfmap_contraction_and_observer_radius'],
            'reject_mixed_anchor_and_nonnegative_bound_failure':(
               EX.reject(lambda: validate_center(center,other))
               and EX.reject(lambda: propagate(abs_matrix(inv),-var,2))
               and EX.reject(lambda: spend(z,1.,Q(1),check))),
            'no_actual_spatial_or_cached_AD_jobs_in_preflight':not ACTUAL.exists()}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/'transfer_ambiguity_cached_derivative_extension_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    try:
        forbid=AssertionError('actual cached AD or derivative assembly is not permitted in preparation')
        with (patch.object(core,'restore_cache',side_effect=forbid),
              patch.object(core,'local_terms',side_effect=forbid),
              patch.object(core,'assemble_point',side_effect=forbid),
              patch.object(core,'variation_rows',side_effect=forbid)):
            a=authority()
            checks=fixtures(a[1])
        require(sweep.all_true(checks),'cached derivative extension preflight failed')
        h=save(PREP/'protocol.json',protocol(a))
        pf=save(PREP/'preflight.json',{'status':'PASS','checks':checks,'source_sha256':sha(SOURCE),
                                     'protocol_sha256':h,'actual_site_AD_assemblies':0,'actual_whole_box_variations':0})
        print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':h,'preflight_sha256':pf}),flush=True)
    except BaseException as e:
        save(PREP/'failure.json',{'error':str(e),'traceback':traceback.format_exc()})
        raise


def admissible(review_path,review_sha):
    links={'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),
           'preflight_sha256':sha(PREP/'preflight.json'),
           'actual_geometry_review_sha256':PINS[GEOMETRY_REVIEW],
           'cached_root_independent_review_sha256':PINS[EX.CACHE_REVIEW]}
    review=pinned(review_path,review_sha)
    require(review['status']=='PASS' and sweep.all_true(review['checks'])
            and all(review[k]==v for k,v in links.items()), 'independent larger cached derivative prelaunch gate required')
    pf=pinned(PREP/'preflight.json',links['preflight_sha256'])
    require(pf['status']=='PASS' and sweep.all_true(pf['checks'])
            and pf['source_sha256']==links['source_sha256']==sha(PREP/'transfer_ambiguity_cached_derivative_extension_frozen.py')
            and pf['protocol_sha256']==links['protocol_sha256']
            and pf['actual_site_AD_assemblies']==pf['actual_whole_box_variations']==0,
            'changed frozen preparation')
    state=authority()
    require(protocol(state)==pinned(PREP/'protocol.json',links['protocol_sha256']),
            'changed derivative scope or input identities')
    return state,{**links,'prelaunch_review_sha256':review_sha}


def worker(review_path,review_sha,pid,parent):
    (_,cases,manifest,contexts),bindings=admissible(review_path,review_sha)
    launch=read(ACTUAL/'started.json')
    require(parent==launch['pid']==os.getppid() and launch['bindings']==bindings
            and pid in launch['pilot_order'],'fresh authorized serial process required')
    folder=ACTUAL/pid
    case,=[r for r in cases if r['pilot_id']==pid]
    identity={'pilot_id':pid,'pid':os.getpid(),'supervisor_pid':parent,'bindings':bindings}
    save(folder/'started.json',identity)
    terminal={**identity,'status':'FAILED'}
    start=monotonic()
    try:
        result=calculate(case,manifest,contexts,folder)
        terminal['results_sha256']=save(folder/'results.json',{**result,'bindings':bindings})
        terminal['status']='FINISHED'
    except BaseException as e:
        terminal['error']=str(e)
        save(folder/'failure.json',{'error':str(e),'traceback':traceback.format_exc()})
        raise
    finally:
        save(folder/'terminal.json',{**terminal,'seconds':monotonic()-start})


def launch(review_path,review_sha):
    _,bindings=admissible(review_path,review_sha)
    ACTUAL.mkdir(exist_ok=False)
    me=os.getpid()
    save(ACTUAL/'started.json',{'pid':me,'bindings':bindings,'pilot_order':list(EX.PILOTS)})
    terminal={'status':'FAILED','pid':me,'bindings':bindings,'workers':[]}
    outputs=[]
    try:
        for pid in EX.PILOTS:
            folder=ACTUAL/pid
            folder.mkdir(exist_ok=False)
            command=[sys.executable,'-B',str(SOURCE),'worker','--review',str(review_path),
                     '--review-sha256',review_sha,'--pilot-id',pid,'--supervisor-pid',str(me)]
            with (folder/'stdout.log').open('x',encoding='utf-8') as log:
                child=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                identity={'pilot_id':pid,'pid':child.pid,'supervisor_pid':me,'bindings':bindings}
                save(folder/'supervisor_child_started.json',{**identity,'command':command})
                code=child.wait()
            child_terminal=folder/'terminal.json'
            observed={**identity,'observed_exit_code':code,
                      'terminal_sha256':sha(child_terminal) if child_terminal.is_file() else None}
            observation_sha=save(folder/'observed_os_exit.json',observed)
            terminal['workers'].append({**observed,'os_exit_receipt_sha256':observation_sha})
            require(observed['terminal_sha256'] is not None,'child terminal missing')
            finished=pinned(child_terminal,observed['terminal_sha256'])
            CP.admit_child(code,finished,identity)
            require(read(folder/'started.json')==identity,'true observed child identity changed')
            result=pinned(folder/'results.json',finished['results_sha256'])
            require(result['case']['pilot_id']==pid and result['bindings']==bindings,'mixed result')
            outputs.append({'pilot_id':pid,'status':result['status'],'ZQ_upper_exact':result['ZQ_upper_exact'],
                            'results_sha256':finished['results_sha256']})
        terminal['results_sha256']=save(ACTUAL/'results.json',{'status':'COMPLETE_FOUR_CACHED_DERIVATIVE_ATTEMPTS_UNREVIEWED',
                                          'case_results':outputs,'bindings':bindings,'scope':SCOPE,
                                          'exact_spatial_bridge_admitted':False})
        terminal['status']='FINISHED'
    except BaseException as error:
        terminal['error']=str(error)
        raise
    finally:
        save(ACTUAL/'terminal.json',terminal)
        print(json.dumps(terminal),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='mode',required=True)
    sub.add_parser('prepare')
    for mode in ('launch','worker'):
        stage=sub.add_parser(mode)
        stage.add_argument('--review',required=True,type=Path)
        stage.add_argument('--review-sha256',required=True)
        if mode=='worker':
            stage.add_argument('--pilot-id',required=True)
            stage.add_argument('--supervisor-pid',required=True,type=int)
    args=cli.parse_args()
    if args.mode=='prepare':
        prepare()
    elif args.mode=='launch':
        launch(args.review.resolve(),args.review_sha256)
    else:
        worker(args.review.resolve(),args.review_sha256,args.pilot_id,args.supervisor_pid)

