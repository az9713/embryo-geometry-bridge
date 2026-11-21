"""Apply same-target cached roots, derivative domains and positive symmetric R parts."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import re
from time import perf_counter
import traceback


HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve()
EXACT=HERE.parent
DESIGN=EXACT.parent
ROOT=EXACT.parents[3]
OBS=DESIGN/'full_outcome_coverage_v1/observer_application_v1'
OLD=EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay'
NEW=EXACT/'remaining32_preconditioner_spd_v1/independent_actual_preparation/actual_replay'
OLD_TARGET=EXACT/'preparation/target_manifest.json'
NEW_TARGET=EXACT/'remaining32_targets_v1/preparation/manifest.json'
OLD_GEO=EXACT/'independent_actual/replay/receipt.json'
NEW_GEO=EXACT/'remaining32_geometry_v1/independent_actual_preparation/actual_replay/receipt.json'
NEW_TRIAL=EXACT/'remaining32_geometry_v1/preparation/protocol.json'
VAR=EXACT/'remaining32_derivative_extension_v1/independent_actual_preparation/actual_replay/receipt.json'
JOIN=EXACT/'remaining32_bridge_v1/admission/bridge_input_admission.json'
OLD_LARGE=EXACT/'cached_derivative_extension_v1/independent_actual_preparation'
METHOD=HERE/'METHOD.md'
ACTUAL=HERE/'actual'
GATE_STATUS='PASS_ROOT_ALL36_CACHED_MINIMUM_APPLICATION_SOURCE'
PILOT_MIN=EXACT/'cached_minimum_application_v1/independent_readback.json'
EXCEPTION='max_contrast_change_small_02'
R34=Q(1,2**34)
R32=Q(1,2**32)
R36=Q(1,2**36)
SITE_COUNTS={'bulk':9216,'surface':288}
MARGINS={'bulk_det_min','bulk_hoop_min','bulk_J_min','normalized_weight_min','surface_radius_min',
         'surface_tangent_squared_min','surface_J_min'}
SCOPE=('Local strict and cube-unique minima of 36 normalized fixed binary cached-quadrature energies. '
       'Same original parameter point, exact optical drive, free coordinate map, cache and proof center. '
       'No exact-spatial minimum, dynamic loading, continuum or experimental compatibility conclusion.')
PINS={
    OBS/'results.json':'1df2bc372edd12687fad17edf63013509994e6bbe41c02150cefb25952eab9f6',
    OBS/'observed_os_exit.json':'16af0dbf81345dba9bec8b18e9c02c22d6321d102a59f035d47981b8aaf1894f',
    OBS/'independent_readback/receipt.json':'a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d',
    OLD/'receipt.json':'769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
    OLD/'os_exit_observation.json':'b5a4a1ab78dc56c7a1147bb75b657dc141d3c2ed4a08ea31074bb5c5a2a0a3a0',
    NEW/'receipt.json':'de31a30f3270829c7fd89ebc75db58e45d2d9a911cd02dfc6e564493748d8615',
    NEW/'os_exit_observation.json':'fd8ac1cb3b920b57096649544bba105442404fe19576aae99c53f849fe15a8d5',
    OLD_TARGET:'cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3',
    NEW_TARGET:'dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242',
    OLD_GEO:'884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    NEW_GEO:'74179dde74571f792b6f3949a06d6e26076b8c746e6700c355264d0de3debf04',
    NEW_GEO.with_name('os_exit_observation.json'):'d39b2d5e97ad6d15bfb68f3f8360e0214c45d11b37dbbd3eaae48313a2e8fc0c',
    NEW_TRIAL:'754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850',
    VAR:'2a1be56b3050215e05da8e11b4785a65cdc8a229210fa82ebee37f1aa8236d05',
    VAR.with_name('os_exit_observation.json'):'69ededa4a104e07ac85e251e6a019b93fc11016cdd83954aaee9f7a987f83325',
    JOIN:'2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb',
    OLD_LARGE/'bridge_admission.json':'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
    OLD_LARGE/'actual_replay_v2/receipt.json':'f766b8b41142ce1f82f3e672f423c32485ae436487d9418c42d907b6f27e8c74',
    OLD_LARGE/'actual_replay_v2/os_exit_observation.json':'be63043a9bafc6369f89ac2dcfe83674541e753cf221e923af18813ca2db3bc0',
    PILOT_MIN:'a7bfee53506f54837c6a0c04f2c67b4ae0b15deca76e133c6a3f16101a25fe10',
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def rel(path):
    return path.relative_to(ROOT).as_posix()


def canon(name):
    path=(ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and rel(path)==name,'canonical study-relative authority required')
    return path


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path,h,links):
    raw=path.read_bytes()
    require(type(h) is str and re.fullmatch('[0-9a-f]{64}',h) is not None and sha256(raw).hexdigest()==h
            and links.setdefault(rel(path),h)==h,'changed or conflicting immutable input')
    return json.loads(raw)


def dig(obj):
    return sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(path,obj):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(obj,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def exact(value):
    require(type(value) in (int,str,Q),'literal exact rational required')
    return Q(value)


def literal(x,expected=None):
    return type(x) is int and (x>0 if expected is None else x==expected)


def yes(x):
    return type(x) is dict and bool(x) and all(v is True for v in x.values())


def selected(rows,k,v):
    row,=(row for row in rows if row[k]==v)
    return row


def posterior(eta,z,r,lhs,error):
    eta,z,r,lhs,error=map(exact,(eta,z,r,lhs,error))
    require(r>0 and eta>=0 and 0<=z<1 and eta+z*r==lhs<r and 0<=error==eta/(1-z)<r,
            'native strict closed-cube contraction and posterior required')
    return eta,z,r,error


def determinant_domains(domain):
    boxes=domain['strict_lower_margins']
    require(domain['counts']==SITE_COUNTS,'complete cached quadrature site catalogue required')
    if 'bulk_J_raw' in domain:
        named={**boxes,'bulk_J_min':domain['bulk_J_raw'],'surface_J_min':domain['surface_J_lower']}
        require(set(boxes)==MARGINS-{'bulk_J_min','surface_J_min'},'recentered analytic-domain schema changed')
    else:
        named=boxes
    require(set(named)==MARGINS and all(0<exact(b['lower_exact'])<=exact(b['upper_exact']) for b in named.values())
            and exact(named['bulk_J_min']['lower_exact'])>Q(1,10**6),
            'same-cube positive cached sites, smooth log/radius/tangent and inactive low-J penalty required')
    return {k:v['lower_exact'] for k,v in named.items()}


def native_root(case,leaf,contract):
    pid=case['case_id']; root=case['root']; state=leaf
    status=leaf['status']; s=Q(0); xi=exact(root['archived_xi']); r=exact(root['radius'])
    if 'own_radius_attempts' in leaf:
        require(status=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT' and leaf['chosen_radius_exact']==root['radius']
                and leaf['rounded_strength_used_in_proof'] is False,'valid wave native root authority required')
        state=selected(leaf['own_radius_attempts'],'radius_exact',leaf['chosen_radius_exact'])
        require(state['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
                and state['whole_cube_AD_recomputed_for_this_radius'] is True,
                'same-radius wave derivative was not independently evaluated')
        fields=('eta_upper_exact','ZQ_upper_exact','radius_exact','strict_selfmap_lhs_upper_exact','posterior_error_upper_exact')
        field,drive,cached=leaf['field_sha256'],leaf['exact_requested_strength'],leaf['cache_packet_sha256']
        goal=leaf['observer_radius_budget_exact']; domain=state['analytic_domain']
    elif status=='PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET':
        fields=('eta_upper_exact','ZQ_upper_exact','radius_exact','strict_selfmap_lhs_upper_exact','posterior_error_upper_exact')
        field,drive,cached=leaf['archived_field_sha256'],leaf['exact_requested_strength'],leaf['recentered_cache_packet_sha256']
        goal=leaf['archived_observer_goal_exact']; domain=leaf['positive_whole_cached_domain']
        s=exact(leaf['actual_binary_center_shift_sigma_exact'])
        require(s>0 and xi==exact(leaf['archived_addition_gap_xi0_exact'])
                and leaf['fresh_center_and_entire_cube_recomputed'] is True
                and leaf['strict_old_observer_goal_pass'] is True and leaf['original_field_or_scores_replaced'] is False
                and leaf['replayed_u1_binary_sha256']==contract['center_binary_sha256']
                and leaf['original_binary_R_sha256']==contract['binary_R_sha256'],
                'recentered root must retain exact archived sigma, xi, fields and unchanged R')
    elif status=='PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT':
        fields=('eta_upper_exact','independent_ZQ_upper_exact','radius_exact','strict_selfmap_lhs_upper_exact','posterior_error_upper_exact')
        field,drive,cached=leaf['field_sha256'],leaf['exact_strength'],leaf['cache_packet_sha256']
        goal=leaf['observer_goal_exact']; domain=leaf['independent_whole_cube_domains']
        require(leaf['original_unresolved_test_preserved'] is True,'historical unsuccessful M proof must remain')
    else:
        require(status in ('PASS_INDEPENDENT_CACHED_LOCAL_ROOT','PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT')
                and leaf['rounded_strength_used_in_proof'] is False,'unrecognized native root proof class')
        is_pilot=status=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
        fields=('independent_eta_upper_exact','independent_z_upper_exact' if is_pilot else 'independent_ZQ_upper_exact',
                'radius_exact','independent_selfmap_lhs_upper_exact' if is_pilot else 'strict_selfmap_lhs_upper_exact',
                'independent_posterior_error_upper_exact')
        field,drive,cached=leaf['field_sha256'],leaf['exact_strength'],leaf['cache_packet_sha256']
        goal=leaf['science_radius_budget_exact' if is_pilot else 'observer_goal_exact']; domain=leaf['analytic_domain']
    eta,z,native_r,e=posterior(*(state[k] for k in fields))
    require(pid==leaf.get('case_id',leaf.get('pilot_id')) and field==case['field_sha256']
            and leaf['record_id']==case['record_id']
            and native_r==r and r in (R32,R34,R36) and e==exact(root['posterior'])
            and s==exact(root['sigma']) and s+e<exact(goal) and xi>=0
            and exact(root['epsilon_um'])==Q(67,2)*(s+e+xi),
            'not the original independently scored local root and observer budget')
    return {'eta':eta,'z':z,'radius':r,'e':e,'sigma':s,'xi':xi,'goal':exact(goal),'field_sha256':field,
            'cache_packet_sha256':cached,'drive':exact(drive),'analytic_lower':determinant_domains(domain),
            'native_status':status}


def norm_and_q(spd,pilot):
    if pilot:
        require(spd['status']=='PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD'
                and spd['exact_target_symmetric_part_spd'] is True,'earlier actual same-case positive symmetric R part required')
        q=exact(spd['max_reviewed_q_upper_exact'])
        require(q==max(exact(spd[k]) for k in ('fresh_congruence_q_upper_exact','producer_congruence_q_upper_exact')),
                'earlier independent congruence operands differ')
    else:
        require(spd['status']=='PASS_INDEPENDENT_REMAINING32_EXACT_SYMMETRIC_PART_SPD'
                and spd['strict_congruence_bound'] is True and spd['root_or_minimum_claimed'] is False,
                'same-case new independent positive symmetric R part required')
        q=exact(spd['independent_q_infinity_upper_exact'])
    n=exact(spd['exact_R_infinity_norm'])
    require(0<=q<1 and n>0,'strict independently verified positive symmetric part with exact norm required')
    return n,q


def os_chain(folder,receipt,os_exit,pid,expected,links):
    receipt_h=PINS[receipt]
    terminal=load(folder/'terminal.json',os_exit['terminal_sha256'],links)
    require(os_exit['status']==expected and literal(os_exit['observed_exit_code'],0) and literal(terminal['exit_code'],0)
            and terminal['status']=='FINISHED' and literal(pid) and literal(os_exit['pid'],pid)
            and literal(terminal['pid'],pid) and os_exit['receipt_sha256']==terminal['receipt_sha256']==receipt_h,
            'externally observed same-process terminal receipt required')


def verify_box(box,identity,r):
    idx=box['free']; active=set(idx)
    require(box['record_id']==identity['record_id'] and exact(box['corrector_radius_exact'])==r
            and box['parameter_cell_sha256']==identity['parameter_cell_sha256']
            and exact(box['dose_exact'])==exact(identity['dose_exact'])
            and box['model_observation_contract_sha256']==identity['model_contract_sha256']
            and box['point_cache_sha256']==identity['cache_packet_sha256']
            and len(idx)==len(active)==611 and all(type(i) is int and 0<=i<650 for i in idx)
            and len(box['reference_nodes_exact'])==325 and len(box['triangles'])==144
            and box['root_membership_assumed'] is False and box['any_spatial_derivative_bound_assumed'] is False
            and len(box['predictor_hulls_exact'])==611 and all(tuple(map(exact,h))==(0,0) for h in box['predictor_hulls_exact']),
            'fixed-input whole-body cube does not match the preconditioner and scalar energy model')
    n=0
    for i,(x,width,h) in enumerate(zip(box['current_centers_exact'],box['coordinate_radii_exact'],box['coordinate_hulls_exact'],strict=True)):
        require(len(x)==len(width)==len(h)==2,'fixed two-coordinate node hull required')
        for d,(c,w,ends) in enumerate(zip(x,width,h,strict=True)):
            c,w=exact(c),exact(w)
            require(w==(r if 2*i+d in active else 0) and tuple(map(exact,ends))==(c-w,c+w),
                    'declared r34 cube lacks exact coordinate hull coverage')
            n+=1
    require(n==650,'complete 650-entry matched geometry cube required')


def same_centers(a,b):
    shared=('record_id','dose_exact','parameter_cell_sha256','model_observation_contract_sha256','point_cache_sha256',
            'predictor_sha256','reference_nodes_exact','current_centers_exact','triangles','free','predictor_hulls_exact')
    require(all(a[k]==b[k] for k in shared),'cannot use proof and geometry cubes around different exact centers')


def match_model(c,e,identity,packet,r,geo,pilot):
    field_key='source_field_sha256' if pilot else 'archived_field_sha256'
    Rkey='binary_preconditioner_sha256' if pilot else 'binary_R_sha256'
    drive=c['exact_requested_drive']; params=c['parameter_cell']
    dose,gain,strength=(exact(drive[k]) for k in ('optical_dose_exact','common_gain_exact','activation_strength_exact'))
    require(c['pilot_id']==identity['pilot_id'] and c['record_id']==identity['record_id']==e['record_id']
            and r['field_sha256']==c[field_key]==identity['field_sha256']==packet['input_state_sha256']
            and r['cache_packet_sha256']==c['point_cache_sha256']==identity['cache_packet_sha256']
            and packet['frozen_binary_cache']['sha256']==c['frozen_binary_cache_sha256']==identity['embedded_cache_sha256']
            and c[Rkey]==identity['R_binary_sha256'] and e['model_contract_sha256']==identity['model_contract_sha256']
            and dig(params)==c['parameter_cell_sha256']==identity['parameter_cell_sha256']
            and tuple(map(exact,params['common_gain']))==(gain,gain)
            and all(tuple(map(exact,v))==(Q(packet['parameters'][k]),)*2 for k,v in params.items() if k!='common_gain')
            and (dose,gain,strength)==tuple(map(exact,(identity['dose_exact'],identity['gain_exact'],identity['drive_exact'])))
            and 0<r['drive']==strength==dose*gain==exact(packet['activation_strength_exact'])
            and dose==exact(packet['optical_dose_exact']) and gain==exact(packet['common_gain_exact'])
            and strength!=exact(packet['activation_strength_binary_exact'])
            and 0<exact(c['energy_scale_binary_exact'])==exact(identity['energy_scale_binary_exact'])
            ==exact(packet['normalization_binary_exact']) and literal(packet['free_dof'],611)
            and packet['material_parameter_storage_exact']==c['material_parameters_binary_exact'],
            'operator mismatch in R, exact input, archived field, material, cache, normalization or cell')
    if pilot:
        require(geo['status']=='PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY'
                and geo['proof_center_field_sha256']==r['field_sha256'] and r['sigma']==0 and r['radius']==R36,
                'same original prior-pilot geometry center required')
    else:
        require(geo['status']=='PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY'
                and geo['center_binary_sha256']==c['center_binary_sha256']==identity['center_binary_sha256'],
                'same current 32-case geometry proof center required')
        require(geo['point_cache_sha256']==r['cache_packet_sha256'] and geo['R_binary_sha256']==identity['R_binary_sha256']
                and geo['archived_field_sha256']==r['field_sha256'] and exact(geo['exact_drive'])==strength
                and exact(geo['archived_shift_exact'])==exact(c['archived_shift_exact'])==r['sigma']
                and exact(geo['archived_addition_gap_exact'])==exact(c['archived_addition_gap_exact'])==r['xi']
                and geo['closed_full_nodal_coordinates_rebuilt']==650 and geo['root_asserted'] is False,
                'signed original observer archive or exact same geometry changed')
    require(geo['pilot_id']==identity['pilot_id'] and geo['record_id']==identity['record_id']
            and exact(geo['radius_exact'])==R34 and exact(geo['strong_monotonicity_lower_exact'])>0
            and 0<exact(geo['body_J_bounds_exact'][0])<=exact(geo['body_J_bounds_exact'][1])
            and geo['all_144_body_and_24_trace_domains'] is True and literal(geo['raw_surface_cache_geometry_lines'],288),
            'accepted same-center whole-box rather than cached-site-only geometry required')


def pilot_derivative(d,mark,sp,r,box_h,n):
    identity=sp['target_identity']
    z=max(exact(d['reviewed_ZQ_upper_exact']),exact(d['producer_ZQ_upper_exact']))
    require(d['status']=='PASS_INDEPENDENT_LARGER_RADIUS_CACHED_DERIVATIVE'
            and d['pilot_id']==mark['pilot_id']==identity['pilot_id']
            and d['box_sha256']==mark['box_sha256']==box_h
            and d['cache_packet_sha256']==mark['cache_packet_sha256']==identity['cache_packet_sha256']
            and exact(d['radius_exact'])==exact(mark['radius_exact'])==R34
            and exact(d['prior_point_eta_upper_exact'])==r['eta'] and exact(d['proof_drive_exact'])==r['drive']
            and exact(d['preconditioner_infinity_norm_exact'])==n
            and d['reviewed_ZQ_upper_exact']==mark['reviewed_ZQ_upper_exact']
            and literal(d['certified_site_counts']['bulk'],9216) and literal(d['certified_site_counts']['surface'],288)
            and 0<=min(map(exact,(d['reviewed_ZQ_upper_exact'],d['producer_ZQ_upper_exact'])))<=z<1,
            'prior four need their own independently verified r34 cached-derivative geometry match')
    return z,determinant_domains(d['whole_box_domains'])


def fresh_derivative(d,mark,j,identity,r,root_h,n):
    z=max(exact(d['independent_ZQ_upper_exact']),exact(d['producer_ZQ_upper_exact']))
    require(d['status']==mark['status']=='PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE'
            and d['case_id']==mark['case_id']==j['case_id']==identity['pilot_id']
            and d['trial_box_sha256']==j['trial_box_sha256'] and exact(d['radius_exact'])==R34
            and d['model_contract_sha256']==identity['model_contract_sha256']
            and d['center_binary_sha256']==identity['center_binary_sha256']
            and d['R_binary_sha256']==identity['R_binary_sha256']
            and d['cache_packet_sha256']==identity['cache_packet_sha256']
            and d['native_cached_anchor_sha256']==root_h
            and exact(d['proof_drive_exact'])==r['drive'] and exact(d['native_eta_upper_exact'])==r['eta']
            and exact(d['archived_center_shift_exact'])==r['sigma'] and exact(d['archived_xi_exact'])==r['xi']
            and exact(d['R_inf_exact'])==n==exact(j['R_inf_exact'])
            and z==exact(d['reviewed_ZQ_upper_exact'])==exact(j['reviewed_ZQ_upper_exact'])
            and 0<=min(map(exact,(d['independent_ZQ_upper_exact'],d['producer_ZQ_upper_exact'])))<=z<1
            and d['root_asserted'] is False and d['exact_spatial_derivative_claimed'] is False,
            'fresh r34 derivative bound does not describe the same signed anchor and normalized cached energy')
    return z,determinant_domains(d['analytic_domain'])


def localize(pid,r,z,aux):
    if r['e']<R34:
        return r['e'],'native_posterior_inside_r34',None
    require(pid==EXCEPTION and r['radius']==R32 and aux is not None
            and aux['diagnostic_cached_r34_selfmap_pass'] is True,'native root lacks an admitted interior r34 witness')
    eta,aux_z,small,e=posterior(aux['native_eta_upper_exact'],aux['reviewed_ZQ_upper_exact'],aux['radius_exact'],
                               aux['diagnostic_cached_r34_selfmap_lhs_upper_exact'],aux['diagnostic_cached_r34_posterior_upper_exact'])
    require(eta==r['eta'] and aux_z==z and small==R34<r['radius'] and r['sigma']+e<r['goal'],
            'small cached self-map must locate the same unique native root on the same centered cube')
    return e,'same_root_via_strict_enclosed_small_cube_selfmap',eta+z*R34


def input_rows():
    links={rel(SOURCE):sha(SOURCE),rel(METHOD):sha(METHOD)}
    p={path:load(path,h,links) for path,h in PINS.items()}
    obs=p[OBS/'results.json']; verified=p[OBS/'independent_readback/receipt.json']; os_obs=p[OBS/'observed_os_exit.json']
    old=p[OLD/'receipt.json']; new=p[NEW/'receipt.json']; join=p[JOIN]
    old_geo=p[OLD_GEO]; new_geo=p[NEW_GEO]; old_der=p[OLD_LARGE/'actual_replay_v2/receipt.json']
    old_wrap=p[OLD_LARGE/'bridge_admission.json']; variation=p[VAR]
    require(obs['status']=='CACHED_QUADRATURE_PRIMARY_OBSERVER_APPLICATION_AWAITING_INDEPENDENT_READBACK'
            and verified['status']=='PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS'
            and verified['cached_fixed_parameter_roots']==obs['cached_root_fields']==36
            and verified['primary_signed_residuals']==obs['primary_calibration_comparisons']==136
            and verified['scope']==obs['primary_observer']=='current_disk__spatial_max'
            and obs['exact_spatial_root_observer_results']==verified['exact_spatial_observer_applications_from_this_readback']==0
            and obs['full_five_gate_goal_achieved'] is verified['full_five_gate_goal_achieved'] is False
            and os_obs['status']=='PASS_ROOT_OBSERVED_CACHED_36_OBSERVER_APPLICATION_EXIT' and literal(os_obs['observed_exit_code'],0)
            and os_obs['results_sha256']==PINS[OBS/'results.json']
            and verified['current_input_sha256'][rel(OBS/'results.json')]==PINS[OBS/'results.json'],
            'same original all36 native roots and independently audited observer archive required')
    require(old['status']=='PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD'
            and old['fresh_actual_matrix_congruences']==4 and old['root_or_minimum_certified'] is False
            and new['status']=='PASS_INDEPENDENT_ACTUAL_REMAINING32_EXACT_SYMMETRIC_PARTS'
            and new['independently_verified_symmetric_parts_SPD']==32 and new['remaining_unresolved_sufficient_bounds']==0
            and new['root_or_minimum_claimed'] is False and new['preconditioner_or_archived_field_changed'] is False
            and p[PILOT_MIN]['status']=='PASS_INDEPENDENT_FOUR_CACHED_MINIMUM_APPLICATION',
            '4 prior and 32 new whole recorded symmetric-part SPD premises required')
    require(old_geo['status']=='PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES'
            and new_geo['status']=='PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY'
            and join['status']=='PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS' and yes(join['checks'])
            and join['target_manifest_sha256']==PINS[NEW_TARGET] and join['geometry_receipt_sha256']==PINS[NEW_GEO]
            and join['variation_receipt_sha256']==PINS[VAR] and join['spatial_workers_authorized'] is False
            and join['geometry_observed_exit_sha256']==PINS[NEW_GEO.with_name('os_exit_observation.json')]
            and join['variation_observed_exit_sha256']==PINS[VAR.with_name('os_exit_observation.json')]
            and variation['status']=='PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES'
            and old_wrap['status']=='PASS_REVIEWED_LARGER_CACHED_DERIVATIVE_BRIDGE_INPUTS'
            and old_der['status']=='PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES'
            and old_wrap['native_independent_review_sha256']==PINS[OLD_LARGE/'actual_replay_v2/receipt.json']
            and old_wrap['native_independent_os_exit_sha256']==PINS[OLD_LARGE/'actual_replay_v2/os_exit_observation.json']
            and old_wrap['reviewed_cases']==old_der['reviewed_cases'],'correct geometry and r34 cached derivative proofs required')
    os_chain(OLD,OLD/'receipt.json',p[OLD/'os_exit_observation.json'],old['pid'],'PASS_OBSERVED_INDEPENDENT_ACTUAL_SPD_TRUE_OS_EXIT',links)
    os_chain(NEW,NEW/'receipt.json',p[NEW/'os_exit_observation.json'],new['source_and_phase']['pid'],
             'PASS_ROOT_OBSERVED_REMAINING32_SPD_REPLAY_EXIT',links)
    os_chain(NEW_GEO.parent,NEW_GEO,p[NEW_GEO.with_name('os_exit_observation.json')],new_geo['source_and_phase']['pid'],
             'PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_REPLAY_EXIT',links)
    os_chain(VAR.parent,VAR,p[VAR.with_name('os_exit_observation.json')],variation['pid'],
             'PASS_OBSERVED_INDEPENDENT_REMAINING32_VARIATION_TRUE_OS_EXIT',links)
    os_chain(OLD_LARGE/'actual_replay_v2',OLD_LARGE/'actual_replay_v2/receipt.json',p[OLD_LARGE/'actual_replay_v2/os_exit_observation.json'],
             old_der['source_and_phase']['pid'],'PASS_OBSERVED_INDEPENDENT_DERIVATIVE_REPLAY_TRUE_OS_EXIT',links)
    a4={v['pilot_id']:v for v in old['ordered_cases']}; a32={v['case_id']:v for v in new['case_results']}
    ids=[r['case_id'] for r in obs['cases']]
    check_observer={r['case_id']:r for r in verified['observer_transfer_terms']}
    require(len(ids)==len(set(ids))==36 and len(a4)==4 and len(a32)==32 and not set(a4)&set(a32)
            and set(ids)==set(a4)|set(a32)==set(check_observer)
            and tuple(a32)==tuple(r['case_id'] for r in join['cases'])
            and len(variation['case_results'])==10,'complete unique 4+32 case and independent rooted observer coverage required')
    by_hash={h:name for name,h in obs['source_and_actual_sha256'].items()}
    bundles=[]
    for case in obs['cases']:
        pid=case['case_id']; rt=case['root']; pilot=pid in a4
        def take(file,h):
            return load(canon(file),h,links)
        e=selected(p[OLD_TARGET if pilot else NEW_TARGET]['targets'],'pilot_id' if pilot else 'case_id',pid)
        contract=take(e['model_contract'],e['model_contract_sha256'])
        packet=take(contract['point_cache'],contract['point_cache_sha256'])
        native=take(by_hash[rt['cached_proof_leaf_sha256']],rt['cached_proof_leaf_sha256'])
        own=native_root(case,native,contract)
        score=check_observer[pid]
        require(tuple(map(exact,(score['posterior_exact'],score['sigma_exact'],score['archived_xi_exact'],score['epsilon_um_exact'])))
                ==(own['e'],own['sigma'],own['xi'],Q(67,2)*(own['e']+own['sigma']+own['xi'])),
                'application cannot revise old primary scores or reset recenter shifts')
        mark=a4[pid] if pilot else a32[pid]
        spd_h=mark['case_receipt_sha256' if pilot else 'sha256']
        sp=take(mark['case_receipt_file' if pilot else 'receipt'],spd_h)
        ident=sp['target_identity']
        n,q=norm_and_q(sp,pilot)
        require(ident==mark['target_identity'] and ident['pilot_id']==pid and exact(mark['exact_R_infinity_norm'])==n,
                'native signed congruence and listed matrix digest/norm do not match')
        tensors=native.get('per_tensor_sha256',native.get('cache_tensor_content_sha256'))
        if tensors is not None:
            require(tensors['preconditioner']==ident['R_binary_sha256']
                    and (pilot or tensors['u']==ident['center_binary_sha256']), 'root proof binary R or center differs')
        if 'preconditioner_infinity_norm_exact' in native:
            require(exact(native['preconditioner_infinity_norm_exact'])==n,'native rooted preconditioner norm differs')
        if 'preconditioner_binary_sha256' in native:
            require(native['preconditioner_binary_sha256']==ident['R_binary_sha256'],'historical narrow M R identity differs')
        if pilot:
            gm,=(g for g in old_geo['geometry_records'] if g['pilot_id']==pid and exact(g['radius_exact'])==R34)
            geo=take(gm['receipt'],gm['sha256']); box=take(e['box_1'],e['box_1_sha256'])
            short=take(e['box_0'],e['box_0_sha256']); verify_box(short,ident,R36); same_centers(short,box)
            dm=selected(old_der['reviewed_cases'],'pilot_id',pid)
            dh=dm['native_case_sha256']; d=take(dm['native_case_file'],dh)
            box_h=e['box_1_sha256']
            require(gm['box_sha256']==geo['box_sha256']==box_h,'pilot r34 geometry box mismatch')
            z,domain=pilot_derivative(d,dm,sp,own,box_h,n)
            mode='prior_four_independently_extended_r34_cached_derivatives'
            aux=None
        else:
            gm=selected(new_geo['geometry_records'],'pilot_id',pid)
            geo=take(gm['receipt'],gm['sha256']); gplan=selected(p[NEW_TRIAL]['trials'],'case_id',pid)
            box_h=gplan['trial_box_sha256']; box=take(gplan['trial_box'],box_h)
            j=selected(join['cases'],'case_id',pid)
            require(j['trial_box_sha256']==geo['box_sha256']==gm['box_sha256']==sp['trial_box_sha256']==box_h
                    and j['geometry_replay_case_sha256']==gm['sha256']==sp['native_geometry_case_sha256']
                    and j['cached_derivative_case_sha256']==sp['native_cached_derivative_case_sha256']
                    and sp['R_file']==contract['point_cache'] and sp['R_file_sha256']==contract['point_cache_sha256']
                    and exact(j['R_inf_exact'])==n and contract['native_case_sha256']==rt['cached_proof_leaf_sha256'],
                    'new native geometry, R/cached derivative bound or root anchor mismatch')
            dh=j['cached_derivative_case_sha256']
            if j['cached_input_kind']=='FRESH_R34':
                dm=selected(variation['case_results'],'case_id',pid)
                require(dm['native_sha256']==dh,'fresh cached derivative proof changed')
                d=take(dm['native_file'],dh)
                z,domain=fresh_derivative(d,dm,j,ident,own,rt['cached_proof_leaf_sha256'],n)
                mode='fresh_independent_r34_cached_derivative'; aux=d
            else:
                require(j['cached_input_kind']=='NATIVE_R34' and own['radius']==R34
                        and dh==rt['cached_proof_leaf_sha256'] and own['z']==exact(j['reviewed_ZQ_upper_exact']),
                        'same-r34 native certified root derivative does not match joint scope')
                z=own['z']; domain=own['analytic_lower']; mode='native_independent_r34_cached_derivative'; aux=None
        match_model(contract,{**e,'record_id':box['record_id']},ident,packet,own,geo,pilot)
        verify_box(box,ident,R34)
        located,method,lhs=localize(pid,own,z,aux)
        bundles.append({'pilot_id':pid,'target_identity':ident,'native':own,'sympart_q':q,'R_norm':n,
                        'z_r34':z,'analytic_lower_r34':domain,'geometry_sha256':gm['sha256'],
                        'geometry_box_sha256':box_h,'root_proof_sha256':rt['cached_proof_leaf_sha256'],
                        'spd_proof_sha256':spd_h,'cached_derivative_proof_sha256':dh,'cached_derivative_r34_kind':mode,
                        'r34_error':located,'localization_method':method,'r34_auxiliary_selfmap_lhs':lhs})
    require(len(bundles)==36 and Counter(x['cached_derivative_r34_kind'] for x in bundles)=={
                'prior_four_independently_extended_r34_cached_derivatives':4,
                'fresh_independent_r34_cached_derivative':10,'native_independent_r34_cached_derivative':22}
            and [b['pilot_id'] for b in bundles if b['r34_auxiliary_selfmap_lhs'] is not None]==[EXCEPTION]
            and {b['pilot_id'] for b in bundles if b['native']['sigma']>0}=={
                'max_contrast_change_small_13','max_negative_contrast_change_small_16'},
            '36-case r34 curvature, sole smaller localization, or recenter scope mismatch')
    return bundles,dict(sorted(links.items()))


def curvature(z,n,q,e):
    z,n,q,e=map(exact,(z,n,q,e))
    require(0<=z<1 and n>0 and 0<=q<1 and 0<=e<R34,
            'strict cached Hessian, SPD or interior premises failed')
    m=(1-z)/n
    return m,m/2,R34-e


def stringify(obj):
    if isinstance(obj,Q):
        return str(obj)
    if type(obj) is dict:
        return {k:stringify(v) for k,v in obj.items()}
    if type(obj) is list:
        return [stringify(v) for v in obj]
    return obj


def result_row(bundle):
    r=bundle['native']; z,n,q,e=(bundle[k] for k in ('z_r34','R_norm','sympart_q','r34_error'))
    m,half,margin=curvature(z,n,q,e)
    require(0<margin<=R34 and r['sigma']+e<r['goal'] and r['xi']>=0,
            'no reset of archived score or observer state budget in minimum application')
    return stringify({'status':'CANDIDATE_STRICT_CACHED_LOCAL_MINIMUM_AWAITING_INDEPENDENT_READBACK',
        'pilot_id':bundle['pilot_id'],'target_identity':bundle['target_identity'],
        'proof_sha256':{k:v for k,v in bundle.items() if k in ('geometry_sha256','geometry_box_sha256',
            'root_proof_sha256','spd_proof_sha256','cached_derivative_proof_sha256')},
        'native_root_radius_exact':r['radius'],'native_posterior_exact':r['e'],'native_ZQ_upper_exact':r['z'],
        'fixed_r34_cached_ZQ_upper_exact':z,'fixed_r34_cached_derivative_source':bundle['cached_derivative_r34_kind'],
        'cached_site_analytic_domain_lower_exact':bundle['analytic_lower_r34'],
        'same_center_interior_error_upper_exact':e,'same_center_interior_margin_lower_exact':margin,
        'same_root_localization_method':bundle['localization_method'],
        'auxiliary_smaller_cube_selfmap_lhs_upper_exact':bundle['r34_auxiliary_selfmap_lhs'],
        'radius_of_geometry_and_curvature_cube_exact':R34,'exact_R_infinity_norm':n,'symmetric_R_congruence_q_upper_exact':q,
        'normalized_cached_H_min_eigenvalue_lower_exact':m,'normalized_cached_energy_gap_coefficient_exact':half,
        'archived_sigma_exact':r['sigma'],'archived_xi_exact':r['xi'],
        'archived_scoring_native_posterior_exact':r['e'],
        'archived_observer_score_or_saved_field_replaced':False,
        'conditional_unique_minimum_on_entire_closed_r34_cube':True,
        'exact_spatial_energy_or_loading_claimed':False})


def protocol(bundles,links):
    return {'status':'FROZEN_36_CACHED_ENERGY_MINIMUM_APPLICATION','source_sha256':sha(SOURCE),
            'method_sha256':sha(METHOD),'source_and_input_sha256':links,
            'observer_authority_sha256':PINS[OBS/'independent_readback/receipt.json'],
            'four_spd_authority_sha256':PINS[OLD/'receipt.json'],'remaining32_spd_authority_sha256':PINS[NEW/'receipt.json'],
            'joined_r34_authority_sha256':PINS[JOIN],'four_r34_derivative_authority_sha256':PINS[OLD_LARGE/'bridge_admission.json'],
            'ordered_case_ids':[r['pilot_id'] for r in bundles],'fixed_geometry_and_cached_derivative_radius_exact':str(R34),
            'frozen_input_cases':stringify(bundles),'auxiliary_small_cached_selfmap_only_case_id':EXCEPTION,
            'native_root_posteriors_preserved':True,'archived_observer_inputs_or_scores_changed':False,
            'source_review_required_before_application':True,'independent_actual_scalar_readback_required':True,
            'no_new_mechanical_jobs':True,'exact_spatial_minimum_asserted':False,'scope':SCOPE}


def preflight(ph,checks):
    return {'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,'checks':checks,
            'actual_minimum_applications':0,'new_matrix_force_Hessian_AD_integral_or_solve_jobs':0}


def binding(prot,ph,fh):
    names=('observer_authority_sha256','four_spd_authority_sha256','remaining32_spd_authority_sha256',
           'joined_r34_authority_sha256','four_r34_derivative_authority_sha256')
    return {'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':fh,**{k:prot[k] for k in names}}


def phase_gate(review,pf,p,b):
    require(review['status']==GATE_STATUS and yes(review['checks']) and all(review[k]==v for k,v in b.items())
            and review['ordered_case_ids']==p['ordered_case_ids'] and review['exact_spatial_minimum_asserted'] is False
            and pf['status']=='PASS' and yes(pf['checks'])
            and literal(pf['actual_minimum_applications'],0)
            and literal(pf['new_matrix_force_Hessian_AD_integral_or_solve_jobs'],0)
            and all(pf[k]==b[k] for k in ('source_sha256','protocol_sha256')) and p['source_sha256']==b['source_sha256'],
            'external readback of the current cached-minimum scalar phase required')


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def fixtures(prot,ph):
    # For R=[[2,1/2],[-1/2,2]] and H=I/2 the exact defect is 1/4 and R norm is 5/2.
    m,half,rem=curvature(Q(1,4),Q(5,2),Q(1,10),R34/2)
    require(m==Q(3,10)<Q(1,2) and half==Q(3,20) and rem==R34/2,
            'strict independently adequate curvature and energy-gap fixture differs')
    require(all(reject(lambda:curvature(*v)) for v in ((1,1,0,0),(0,1,1,0),(0,1,0,R34),
                    (0,0,0,0),(0,1,-1,0),(0,1,False,0))) and
            reject(lambda:norm_and_q({'status':'PASS_INDEPENDENT_REMAINING32_EXACT_SYMMETRIC_PART_SPD',
                                     'strict_congruence_bound':False,'root_or_minimum_claimed':False},False)),
            'contraction alone or equality must never certify coercivity')
    e=R34/6; z=Q(1,4)
    r={'e':R34,'radius':R32,'eta':R34/8,'sigma':Q(0),'goal':Q(1)}
    small={'native_eta_upper_exact':str(R34/8),'reviewed_ZQ_upper_exact':str(z),'radius_exact':str(R34),
           'diagnostic_cached_r34_selfmap_lhs_upper_exact':str(3*R34/8),
           'diagnostic_cached_r34_posterior_upper_exact':str(e),'diagnostic_cached_r34_selfmap_pass':True}
    require(localize(EXCEPTION,r,z,small)==(e,'same_root_via_strict_enclosed_small_cube_selfmap',3*R34/8)
            and localize('any',{**r,'e':R34/2},z,None)==(R34/2,'native_posterior_inside_r34',None)
            and reject(lambda:localize('other',r,z,small))
            and reject(lambda:localize(EXCEPTION,{**r,'radius':R34},z,small))
            and reject(lambda:localize(EXCEPTION,r,z,{**small,'diagnostic_cached_r34_selfmap_pass':False})),
            'same-root nested localization and unresolved narrow case guard failed')
    boxes={k:{'lower_exact':'1','upper_exact':'2'} for k in MARGINS}
    require(set(determinant_domains({'strict_lower_margins':boxes,'counts':SITE_COUNTS}))==MARGINS
            and reject(lambda:determinant_domains({'strict_lower_margins':{**boxes,'bulk_J_min':{
                    'lower_exact':'1/1000000','upper_exact':'1'}},'counts':SITE_COUNTS})),
            'body geometry cannot replace the cached analytic site premises')
    h='e'*64; b=binding(prot,ph,h); pf=preflight(ph,{'fixture':True})
    review={'status':GATE_STATUS,**b,'checks':{'fixture':True},'ordered_case_ids':prot['ordered_case_ids'],
            'exact_spatial_minimum_asserted':False}
    phase_gate(review,pf,prot,b)
    require(all(reject(lambda:phase_gate(v,pf,prot,b)) for v in (
        {**review,'preflight_sha256':'f'*64},{**review,'exact_spatial_minimum_asserted':True},
        {**review,'ordered_case_ids':prot['ordered_case_ids'][1:]},{**review,'checks':{'truthy':1}}))
        and reject(lambda:phase_gate(review,{**pf,'actual_minimum_applications':False},prot,b)),
        'independent live phase and malformed hash/scope fixtures failed')
    return {'exact_nonsymmetric_R_curvature_and_energy_gap_fixture':True,
            'reject_contraction_alone_and_strict_equality_without_SymR_proof':True,
            'same_root_nested_localization_and_wrong_scope_rejections':True,
            'cached_analytic_site_premises_not_replaced_by_geometry':True,
            'external_source_gate_real_schema_without_preflight_self_hash':True,
            '36_unchanged_same_r34_geometry_and_derivative_authority_checks':True,
            '35_native_localizations_and_one_separate_plus02_small_selfmap_metadata':True}


def prepare():
    require(not ACTUAL.exists(),'full cached minimum application outputs are write once')
    all_cases,links=input_rows()
    p=protocol(all_cases,links)
    ph=sha256((json.dumps(p,indent=2,allow_nan=False)+'\n').encode('utf-8')).hexdigest()
    checks=fixtures(p,ph)
    with (HERE/'apply_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    require(save(HERE/'protocol.json',p)==ph,'frozen protocol serialization drift')
    fh=save(HERE/'preflight.json',preflight(ph,checks))
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':fh,
                      'bound_existing_proof_cases':len(all_cases),'actual_minimum_applications':0}),flush=True)


def actual(args):
    p=load(HERE/'protocol.json',args.protocol_sha256,{})
    pf=load(HERE/'preflight.json',args.preflight_sha256,{})
    require(args.source_review.resolve()==HERE/'root_source_review.json','canonical external source readback required')
    review=load(args.source_review,args.source_review_sha256,{})
    b=binding(p,args.protocol_sha256,args.preflight_sha256)
    phase_gate(review,pf,p,b)
    require(sha(HERE/'apply_frozen.py')==b['source_sha256'],'frozen application source changed')
    bundles,links=input_rows()
    require(p==protocol(bundles,links),'current exact proof context changes the authorized scalar application')
    phase={'pid':os.getpid(),**b,'root_source_review_sha256':args.source_review_sha256}
    ACTUAL.mkdir(exist_ok=False)
    save(ACTUAL/'started.json',{'status':'RUNNING',**phase})
    terminal={'status':'FAILED','exit_code':1,**phase}
    start=perf_counter()
    try:
        rows=[result_row(bundle) for bundle in bundles]
        require([r['pilot_id'] for r in rows]==p['ordered_case_ids'] and len(rows)==36,'complete strict geometry state coverage')
        for file,h in links.items():
            require(sha(canon(file))==h,'immutable preconditioner, observer or target source changed while applying')
        output={'status':'CANDIDATE_36_CACHED_ENERGY_LOCAL_MINIMA_AWAITING_INDEPENDENT_ACTUAL_READBACK',
                'source_and_phase':phase,'scope':SCOPE,'cases':rows,'source_and_input_sha256':links,
                'cached_local_minimum_candidate_count':36,'same_geometry_and_cached_derivative_radius_exact':str(R34),
                'same_native_cached_roots':36,'independent_actual_replay_required_before_promotion':True,
                'one_auxiliary_cached_localization_case_id':EXCEPTION,'archived_score_or_saved_field_changed':False,
                'exact_spatial_minima_asserted':0,'full_five_gate_goal_achieved':False,
                'new_matrix_force_Hessian_AD_integral_or_solve_jobs':0,'seconds':perf_counter()-start}
        h=save(ACTUAL/'results.json',output)
        terminal.update(status='FINISHED',exit_code=0,results_sha256=h)
        print(json.dumps({'status':output['status'],'results_sha256':h,'cached_minimum_candidates':36}),flush=True)
    except BaseException as e:
        save(ACTUAL/'failure.json',{'error':str(e),'traceback':traceback.format_exc(),**phase})
        raise
    finally:
        save(ACTUAL/'terminal.json',terminal)


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    modes=cli.add_subparsers(dest='mode',required=True)
    modes.add_parser('prepare')
    launch=modes.add_parser('apply')
    for k in ('protocol-sha256','preflight-sha256','source-review-sha256'):
        launch.add_argument('--'+k,required=True)
    launch.add_argument('--source-review',type=lambda x:Path(x).resolve(),required=True)
    args=cli.parse_args()
    prepare() if args.mode=='prepare' else actual(args)


if __name__=='__main__':
    main()
