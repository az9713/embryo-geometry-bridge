"""Independent exact metadata, cube and local-map checks for the remaining 32 stages."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path


HERE=Path(__file__).resolve().parent
GOUT=HERE.parent
EXACT=GOUT.parent
DOUT=EXACT/'remaining32_derivative_extension_v1'
ROOT=EXACT.parents[3]
TROOT=EXACT/'remaining32_targets_v1'
TREC=TROOT/'independent_prelaunch/review.json'
TMAN=TROOT/'preparation/manifest.json'
TSHA='dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242'
TRSHA='67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9'
HALF=Q(1,2**36)
TRIAL=Q(1,2**34)
LARGE=Q(1,2**32)
POLICIES={TRIAL:'EXISTING_MATCHED_WHOLE_CUBE_CACHED_ZQ',
          HALF:'EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD',
          LARGE:'FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH'}
STAGES={'geometry':GOUT,'variation':DOUT}
ZERO_WORK=('new_geometry_certificates','new_AD_derivatives','new_spatial_integrations','new_solves')


def check(ok,s):
    if not ok:
        raise ValueError(s)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path,h):
    raw=path.read_bytes()
    check(sha256(raw).hexdigest()==h,f'changed authority {path.name}')
    return json.loads(raw)


def canon(s):
    p=(ROOT/s).resolve()
    check(p.is_relative_to(ROOT) and p.relative_to(ROOT).as_posix()==s,'canonical local path required')
    return p


def exact(s):
    check(type(s) in (str,int,Q),'exact rational text required')
    return Q(s)


def typed_ints(seq,n):
    return type(seq) is list and all(type(i) is int and 0<=i<n for i in seq)


def box_hulls(box,r):
    free=box['free']
    check(typed_ints(free,650) and len(set(free))==len(free)==611
          and len(box['reference_nodes_exact'])==325 and len(box['triangles'])==144
          and all(typed_ints(tri,325) and len(tri)==6 for tri in box['triangles'])
          and exact(box['corrector_radius_exact'])==r
          and box['root_membership_assumed'] is False and box['any_spatial_derivative_bound_assumed'] is False,
          'invalid frozen root-independent hull provenance or free map')
    lookup=set(free)
    count=0
    for n,group in enumerate(zip(box['current_centers_exact'],box['coordinate_radii_exact'],box['coordinate_hulls_exact'],strict=True)):
        check(all(len(row)==2 for row in group),'two coordinate components required')
        for component,(v,spread,ends) in enumerate(zip(*group,strict=True)):
            j=2*n+component; value=exact(v); radius=r if j in lookup else Q(0)
            check(exact(spread)==radius and tuple(map(exact,ends))==(value-radius,value+radius),
                  'incorrect closed coordinate hull endpoints')
            count+=1
    check(count==650,'complete 325-node coordinate cube required')


def routes(box,route):
    check(len(route['bulk'])==144 and len(route['surface'])==24 and route['positive_deformed_geometry_proved'] is False,
          'complete unpromoted route catalogues required')
    inverse={j:i for i,j in enumerate(box['free'])}
    for row,k,n in [(z,'bulk',6) for z in route['bulk']]+[(z,'surface',3) for z in route['surface']]:
        nodes=row['nodes']
        check(typed_ints(nodes,325) and len(nodes)==len(set(nodes))==n,'local element/trace node count mismatch')
        local=[2*k+d for k in nodes for d in (0,1)]
        entries=[(i,inverse[j]) for i,j in enumerate(local) if j in inverse]
        check(typed_ints(row['free_local'],2*n) and typed_ints(row['global_indices'],611)
              and entries==list(zip(row['free_local'],row['global_indices'],strict=True)),
              'wrong local-to-global free-coordinate order')
        if k=='bulk':
            e=row['element']
            check(type(e) is int and 0<=e<144 and set(nodes)==set(box['triangles'][e]),
                  'wrong P2 element topology')
        else:
            check(type(row['face']) is int and type(row['edge']) is int and exact(row['Rlo_exact'])<exact(row['Rhi_exact']),
                  'wrong trace route or reference radial interval')
    check([r['element'] for r in route['bulk']]==list(range(144))
          and [(r['face'],r['edge']) for r in route['surface']]==[(f,e) for f in (0,1) for e in range(12)],
          'complete ordered chart coverage required')


def geometry_variation_inputs(expected):
    manifest=pinned(TMAN,TSHA)
    native_review=pinned(TREC,TRSHA)
    check(native_review['status']=='PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS'
          and manifest['geometry_or_variation_actual_launch_authorized'] is False
          and manifest['spatial_bridge_actual_launch_authorized'] is False
          and native_review['manifest_sha256']==TSHA and native_review['actual_mechanics_launch_authorized'] is False
          and all(v is True for v in native_review['checks'].values()),'same existing independent target authority required')
    docs={}
    for stage,folder in STAGES.items():
        proto=pinned(folder/'preparation/protocol.json',expected[f'{stage}_protocol_sha256'])
        pre=pinned(folder/'preparation/preflight.json',expected[f'{stage}_preflight_sha256'])
        check(pre['status']=='PASS' and all(v is True for v in pre['checks'].values())
              and all(type(x[k]) is int and x[k]==0 for k in ZERO_WORK for x in (proto,pre))
              and pre['protocol_sha256']==expected[f'{stage}_protocol_sha256']
              and proto['stage']==pre['stage']==stage and proto['status']=='FROZEN_EXPLICIT_REMAINING32_STAGE'
              and proto['max_concurrent_workers_per_supervisor']==1
              and proto['exact_spatial_bridge_authorized'] is False
              and proto['stage_actual_independent_review_required'] is True
              and proto['geometry_or_variation_actual_requires_independent_review'] is True
              and all(x['source_sha256']==expected['source_sha256'] and x['target_manifest_sha256']==TSHA for x in (proto,pre))
              and exact(proto['closed_radius_exact'])==TRIAL
              and sha(folder/'preparation/launcher_frozen.py')==expected['source_sha256'],
              'current write-once stage preparation and same-source preflight required')
        for name,digest in proto['kernel_source_sha256'].items():
            check(sha(canon(name))==digest,'changed cold transitive kernel')
        docs[stage]=proto
    geo,der=docs['geometry'],docs['variation']
    check(len(geo['trials'])==len(manifest['targets'])==len(native_review['per_case'])==32
          and geo['case_order']==native_review['ordered_case_ids']==[r['case_id'] for r in geo['trials']]
          and len(set(geo['case_order']))==32,'canonical stage catalogs required')
    policies=Counter(); new_boxes=[]; shifted=[]; checked=[]
    for trial,old,proof in zip(geo['trials'],manifest['targets'],native_review['per_case'],strict=True):
        pid=trial['case_id']; c=pinned(canon(trial['model_contract']),trial['model_contract_sha256'])
        check(all(trial[k]==v for k,v in old.items()) and pid==old['case_id']==proof['case_id']==c['pilot_id']
              and all(trial[k]==proof[k] for k in ('record_id','model_contract_sha256','box_sha256','routes_sha256','rule_sha256')),
              'mutated reviewed model, route or old box proposal')
        native_radius=exact(c['source_cached_radius_exact'])
        z=exact(c['native_cached_ZQ_upper_exact'])
        policy=POLICIES[native_radius]; policies[policy]+=1
        matched=policy==POLICIES[TRIAL]
        original_matched=native_radius!=LARGE
        center=c['archived_shift_exact']; xi=c['archived_addition_gap_exact']; goal=exact(c['observer_goal_exact'])
        e=exact(proof['native_cached_posterior_upper_exact']); s=exact(center)
        check(trial['derivative_policy']==proof['derivative_policy']==c['derivative_policy']==policy
              and exact(proof['native_cached_radius_exact'])==native_radius and exact(proof['native_ZQ_upper_exact'])==z
              and 0<=z<1 and 0<=e<native_radius and 0<=s and 0<=exact(xi) and s+e<goal
              and exact(trial['trial_radius_exact'])==TRIAL
              and exact(trial['planned_geometry_radius_exact'])==(HALF if native_radius==HALF else TRIAL)
              and trial['native_cached_bound_usable_at_trial_radius'] is matched
              and proof['planned_radius_uses_native_matched_ZQ'] is original_matched
              and c['matched_cached_ZQ_available'] is original_matched
              and trial['native_geometry_claimed'] is False
              and trial['archived_center_shift_exact']==center==proof['observer_shift_sigma_exact']
              and xi==proof['observer_archived_xi_exact']
              and trial['center_binary_sha256']==c['center_binary_sha256']==proof['center_sha256']
              and trial['R_binary_sha256']==c['binary_R_sha256']==proof['binary_R_sha256']
              and trial['point_cache_sha256']==c['point_cache_sha256']
              and trial['energy_scale_binary_exact']==c['energy_scale_binary_exact']
              and trial['exact_drive']==c['exact_requested_drive']['activation_strength_exact']==proof['exact_drive']
              and exact(trial['exact_drive'])==exact(c['exact_requested_drive']['optical_dose_exact'])
                               *exact(c['exact_requested_drive']['common_gain_exact']),
              'mismatched radius authority, point, observer or exact drive')
        orig=pinned(canon(trial['box']),trial['box_sha256'])
        box=pinned(canon(trial['trial_box']),trial['trial_box_sha256'])
        box_hulls(box,TRIAL)
        if native_radius==HALF:
            box_hulls(orig,HALF)
            new_boxes.append(pid)
            changes={'coordinate_radii_exact','coordinate_hulls_exact','corrector_radius_exact'}
            check(set(box)==set(orig) and {k for k in box if box[k]!=orig[k]}==changes
                  and trial['trial_box_sha256']!=trial['box_sha256'] and not matched,
                  'three enlarged trials must keep exactly the same operator and center')
        else:
            check(box==orig and trial['trial_box_sha256']==trial['box_sha256'] and trial['trial_box']==trial['box'],
                  'unedited existing exact r34 plan required')
        route=pinned(canon(trial['routes']),trial['routes_sha256'])
        routes(box,route)
        check(sha(canon(trial['rule']))==trial['rule_sha256'],'derived auxiliary quadrature rule changed')
        if s:
            shifted.append(pid)
            native=pinned(canon(c['native_case']),c['native_case_sha256'])
            exact_step=canon(c['point_cache']).with_name('exact_binary_step.json')
            step=pinned(exact_step,c['recenter_exact_step_sha256'])
            check(c['verification_center_kind']=='VERIFIED_RECENTERED_U1'
                  and native['status']=='PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET'
                  and native['archived_field_sha256']==c['archived_field_sha256']
                  and native['archived_observer_goal_exact']==c['observer_goal_exact']
                  and native['archived_addition_gap_xi0_exact']==xi
                  and native['replayed_u1_binary_sha256']==step['u1_sha256']==trial['center_binary_sha256']
                  and c['point_cache_sha256']==native['recentered_cache_packet_sha256']
                  and native['actual_binary_center_shift_sigma_exact']==step['sigma_exact']==center
                  and exact(native['shift_plus_posterior_exact'])==s+e and e==exact(native['posterior_error_upper_exact']),
                  'released recenter archive must retain exact sigma, xi, field and goals')
        else:
            check(c['recenter_exact_step_sha256'] is None,'recenter step cannot masquerade as original')
        checked.append({'case_id':pid,'policy':policy,'shift_exact':center,'xi_exact':xi,
                        'trial_radius_exact':trial['trial_radius_exact'],'trial_box_sha256':trial['trial_box_sha256']})
    chosen=[row for row in geo['trials'] if not row['native_cached_bound_usable_at_trial_radius']]
    check(der['trials']==chosen and der['case_order']==[r['case_id'] for r in chosen]
          and policies==Counter({POLICIES[TRIAL]:22,POLICIES[LARGE]:7,POLICIES[HALF]:3})
          and set(shifted)=={'max_contrast_change_small_13','max_negative_contrast_change_small_16'}
          and len(new_boxes)==3 and len(chosen)==10,'22 matched, seven narrower, three enlarged coverage required')
    return {'orders':{s:p['case_order'] for s,p in docs.items()},'new_trial_boxes':new_boxes,
            'preserved_recentered_shifts':shifted,'cases':checked,'policy_counts':dict(policies),
            'trial_boxes_checked':len(checked),'separate_preserved_smaller_boxes_checked':len(new_boxes),
            'total_distinct_boxes_checked':len(checked)+len(new_boxes),
            'coordinates_per_trial':650,'bulk_routes_per_trial':144,'surface_routes_per_trial':24}
