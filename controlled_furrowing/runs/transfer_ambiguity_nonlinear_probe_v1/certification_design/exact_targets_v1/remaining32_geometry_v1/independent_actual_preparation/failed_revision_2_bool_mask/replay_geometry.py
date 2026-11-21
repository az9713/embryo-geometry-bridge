"""Rebuild, replay and independently check 32 fixed-input geometry certificates."""
import argparse
from collections import Counter
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import ceil
import os
from pathlib import Path
import re
import sys
from time import perf_counter
import traceback
from unittest.mock import patch


HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve()
GOUT=HERE.parent
EXACT=GOUT.parent
ROOT=EXACT.parents[3]
ACTUAL=GOUT/'actual'
REPLAY=HERE/'actual_replay'
TARGET=EXACT/'remaining32_targets_v1'
GEOMETRY=GOUT/'preparation'
DVAR=EXACT/'remaining32_derivative_extension_v1/preparation'
GATE=GOUT/'independent_prelaunch/review.json'
T_GATE=TARGET/'independent_prelaunch/review.json'
G_SOURCE=ROOT/'scripts/transfer_ambiguity_remaining32_geometry_variation.py'
T_SOURCE=ROOT/'scripts/transfer_ambiguity_remaining32_exact_targets.py'
OLD_AUDIT=EXACT/'independent_actual/audit.py'
REJECTED=HERE/'cold_revision_1_activation_mask_hold'
ADMISSION=ACTUAL/'root_terminal_admission.json'
ADMISSION_SHA='50ce32f50b5e330b9c78fc1863527c49eb6c5e4394bdb1c4e4313cea8628efe0'
STAGE='geometry'
SESSION=17811
RADIUS=Q(1,2**34)
ADMIT_STATUS='PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_TERMINALS'
ROOT_STATUS='PASS_ROOT_REMAINING32_GEOMETRY_REPLAY_SOURCE'
RESULT_STATUS='PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY'
CASE_PASS='PASS_CANDIDATE_WHOLE_BOX_GEOMETRY_AWAITS_INDEPENDENT_REVIEW'
CASE_OTHER=('UNRESOLVED_WHOLE_BODY_GEOMETRY','UNRESOLVED_WHOLE_BOX_ANALYTIC_SUBDOMAIN')
SCOPE=('Recompute exact whole-box geometry at 32 independently admitted fixed binary centers. Reuse the '
       'reviewed Bernstein body and trace kernels, verify exact box endpoints, geometry arithmetic, axis charts '
       'and the finite auxiliary site domains independently. No moment, integration, force, Hessian, AD, '
       'solve, observer or spatial root/stability claim.')
PINS={
    G_SOURCE:'68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a',
    T_SOURCE:'d1a9c9e4839ff85a9926b6fd5ecc11e4cd8e0831c7e6263ecfbe068ef22ea3f4',
    GEOMETRY/'protocol.json':'754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850',
    GEOMETRY/'preflight.json':'d70a7217c4dbbf68f789bb59193321252ea4192576c3f41790cc7daea38fb15a',
    GEOMETRY/'launcher_frozen.py':'68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a',
    DVAR/'protocol.json':'6d04fc82779430e66c2d27b5ab7aa508f522f7ffac03b87f8f2e8e2587e1d128',
    DVAR/'preflight.json':'abbbeb7c0c303db2d2675398f8ea09155a19fd969b121e7936d7f4a1a468e264',
    T_GATE:'67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9',
    TARGET/'preparation/protocol.json':'a12669fd1e50cbd7435d016ba60223791c8a70119f622048c28e24ef5a9f456b',
    TARGET/'preparation/preflight.json':'86ddb06c75c2deedea29d62bb7e03da6142d55b41a0291e1567c9e4dba6aab4f',
    TARGET/'preparation/manifest.json':'dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242',
    GATE:'0c30f413a4796703f6f4a8248924587fcb2152a2dfe7449464fb69c99d763511',
    GOUT/'independent_prelaunch/review_stage_inputs.py':'b8c2963b39da989fa44245fe3e2d1a7e3d2720cfe75ecfd80c9d7b574957806c',
    OLD_AUDIT:'2d49bd52fb489b2b9b5766329e88792eefea21d7205f13856f5c7e5c9e3b46f6',
    EXACT/'independent_actual/protocol.json':'868fca02dc73aa666a348629286850acbf95e7b05dfdddf242129c026ee09b76',
    EXACT/'independent_actual/preflight.json':'f526e50a9b89b71bea3179c6d67617758ceaadb423ed39424db379bff667f87a',
    EXACT/'independent_actual/replay/receipt.json':'884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    REJECTED/'replay_geometry_frozen.py':'0b479982baf86d4248b4d4c1eeec9548a9a9920151f706c618b9369651966775',
    REJECTED/'protocol.json':'4b8a4e79a3c917f13d73cc53c446ea2b8b26ebaa2f275cc1067af3ffd4ac22cd',
    REJECTED/'preflight.json':'4ebf2ea52b7c6edae00f359dd010d858eccde056ea069f6220b81d022ebb19de',
}
PARENT_FILES=('started.json','terminal.json','manifest.json')
CHILD_META=('launched.json','started.json','observed_os_exit.json','terminal.json')
CHILD_PAYLOAD=('results.json','body_geometry.json')
CHARTS=('off_axis','axis_vertex','axis_edge_physical')


def check(ok,msg):
    if not ok:
        raise ValueError(msg)


def rel(path):
    return path.relative_to(ROOT).as_posix()


def canon(name):
    path=(ROOT/name).resolve()
    check(path.is_relative_to(ROOT) and rel(path)==name,'canonical study-relative path required')
    return path


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path,expected,*,parse=True):
    check(type(expected) is str and re.fullmatch('[0-9a-f]{64}',expected) is not None,'explicit SHA-256 required')
    raw=path.read_bytes()
    check(sha256(raw).hexdigest()==expected,f'changed authority {path}')
    return json.loads(raw) if parse else raw


def read_pre(path):
    return pinned(path,PINS[path])


def save(path,contents):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(contents,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def digest(obj):
    return sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def integer(value,wanted=None):
    return type(value) is int and (value>0 if wanted is None else value==wanted)


def exact(value):
    check(type(value) in (str,int,Q),'finite exact rational required')
    return Q(value)


def yes(value):
    return type(value) is dict and bool(value) and all(x is True for x in value.values())


def stage_binding():
    return {'source_sha256':PINS[G_SOURCE],'target_manifest_sha256':PINS[TARGET/'preparation/manifest.json'],
            'target_context_review_sha256':PINS[T_GATE],'launcher_review_sha256':PINS[GATE],
            'geometry_protocol_sha256':PINS[GEOMETRY/'protocol.json'],'geometry_preflight_sha256':PINS[GEOMETRY/'preflight.json'],
            'variation_protocol_sha256':PINS[DVAR/'protocol.json'],'variation_preflight_sha256':PINS[DVAR/'preflight.json']}


def plan():
    for path,h in PINS.items():
        pinned(path,h,parse=False)
    g=read_pre(GEOMETRY/'protocol.json'); p=read_pre(GEOMETRY/'preflight.json')
    review=read_pre(GATE); target=read_pre(T_GATE); previous=read_pre(EXACT/'independent_actual/replay/receipt.json')
    expected=stage_binding(); order=g['case_order']
    check(review['status']=='PASS_INDEPENDENT_REMAINING32_GEOMETRY_VARIATION_PRELAUNCH' and yes(review['checks'])
          and all(review[k]==v for k,v in expected.items() if k!='launcher_review_sha256')
          and review['exact_spatial_bridge_authorized'] is False
          and tuple(order)==tuple(review['permitted_case_ids'][STAGE])==tuple(target['ordered_case_ids'])
          ==tuple(r['case_id'] for r in g['trials']) and len(set(order))==len(order)==32
          and review['exact_metadata_reconstruction']['total_distinct_boxes_checked']==35
          and target['status']=='PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS' and yes(target['checks'])
          and target['actual_mechanics_launch_authorized'] is False
          and p['status']=='PASS' and yes(p['checks']) and p['protocol_sha256']==PINS[GEOMETRY/'protocol.json']
          and g['stage']==STAGE and g['exact_spatial_bridge_authorized'] is False and exact(g['closed_radius_exact'])==RADIUS
          and previous['status']=='PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES'
          and previous['exact_spatial_root_count']==0,'matched independent source authorities required')
    for name,h in g['kernel_source_sha256'].items():
        pinned(canon(name),h,parse=False)
    return {'status':'FROZEN_INDEPENDENT_ACTUAL_REMAINING32_GEOMETRY_REPLAY',
            'source_sha256':sha(SOURCE),'frozen_input_sha256':{rel(p):h for p,h in PINS.items()},
            'producer_binding':expected,'kernel_closure_sha256':g['kernel_source_sha256'],
            'required_root_observation_status':ADMIT_STATUS,'required_root_source_review_status':ROOT_STATUS,
            'ordered_ids':order,'launch_group':digest(order),'unified_exec_session':SESSION,
            'radius_exact':str(RADIUS),'producer_root_admission_sha256':ADMISSION_SHA,
            'proof_engine_reuse':'Eight-case independent endpoint algebra plus its reviewed Bernstein body/trace kernels. Own all-32 routing and exact X+u reconstruction, including separately audited archived xi and sigma.',
            'scope':SCOPE,'before_actual_geometry_payload_opens':0,
            'rejected_cold_revision':'The source-only first draft tested membership in a 0/1 mask. '
                'This version independently rebuilds the geometric indicator by element position.'}


def validate_os_docs(start,end,manifest,pid,group,order,bindings,manifest_hash):
    identity={'pid':pid,'stage':STAGE,'case_order':order,'launch_group':group,'bindings':bindings}
    check(all(integer(d['pid'],pid) and all(d[k]==v for k,v in identity.items()) for d in (start,end,manifest))
          and start['status']=='RUNNING' and end['status']=='FINISHED' and integer(end['exit_code'],0)
          and end['manifest_sha256']==manifest_hash
          and manifest['status']=='COMPLETE_ACTUAL_WORK_AWAITS_INDEPENDENT_REPLAY' and manifest['root_asserted'] is False
          and end['workers']==manifest['workers'] and [w['case_id'] for w in manifest['workers']]==order,
          'actual producer full terminal/manifest or typed process chain invalid')


def validate_worker(case,item,docs,pins,pid,group,bindings,supervisor_start_hash):
    launch,start,os_record,end=docs
    child=launch['pid']
    ident={'pid':child,'supervisor_pid':pid,'stage':STAGE,'case_id':case,'launch_group':group,'bindings':bindings}
    all_rows=(launch,start,os_record,end,item)
    check(integer(child) and child!=pid and all(integer(d['pid'],child) and integer(d['supervisor_pid'],pid)
           and all(d[k]==v for k,v in ident.items()) for d in all_rows)
          and start['status']=='RUNNING' and end['status']=='FINISHED'
          and all(integer(d[k],0) for d,k in ((os_record,'observed_exit_code'),(end,'exit_code'),(item,'observed_exit_code')))
          and pins['terminal.json']==os_record['terminal_sha256']
          and end['results_sha256']==pins['results.json']==item['results_sha256']
          and item['candidate_status'] in (CASE_PASS,*CASE_OTHER)
          and item=={**os_record,'os_exit_sha256':pins['observed_os_exit.json'],
                    'results_sha256':pins['results.json'],'candidate_status':item['candidate_status']},
          'incomplete, mixed or unsuccessful typed worker terminal')
    cmd=launch['command']
    tail=['worker','--stage',STAGE,'--case-id',case,'--supervisor-pid',str(pid),'--group',group,
          '--supervisor-started-sha256',supervisor_start_hash,'--review',str(GATE),'--review-sha256',PINS[GATE],
          '--target-review',str(T_GATE),'--target-review-sha256',PINS[T_GATE]]
    check(type(cmd) is list and len(cmd)==len(tail)+3 and cmd[1]=='-B' and Path(cmd[2]).resolve()==G_SOURCE
          and cmd[3:]==tail,'different executed reviewed worker command')
    return ident


def observed_inputs(admission_path,admission_sha,p):
    check(admission_path.resolve()==ADMISSION,'canonical root-observation receipt required')
    check(admission_sha==p['producer_root_admission_sha256']==ADMISSION_SHA,'matching frozen root observation required')
    admitted=pinned(admission_path,admission_sha)
    pid=admitted['supervisor_pid']; group=p['launch_group']; order=p['ordered_ids']; binding=p['producer_binding']
    check(admitted['status']==ADMIT_STATUS and admitted['stage']==STAGE and integer(pid)
          and integer(admitted['unified_exec_session'],SESSION) and integer(admitted['observed_exit_code'],0)
          and admitted['launch_group']==group,'same root-observed external stage exit required')
    supervisor=ACTUAL/'launches'/group
    folders={case:ACTUAL/'cases'/case for case in order}
    paths={rel(supervisor/n) for n in PARENT_FILES}|{rel(f/k) for f in folders.values() for k in (*CHILD_META,*CHILD_PAYLOAD)}
    input_hash=admitted['sha256']
    check(len(paths)==195 and set(input_hash)==paths,'required complete 32 geometry terminal and payload path set')
    read=lambda f,k:pinned(f/k,input_hash[rel(f/k)])
    start,end,m=(read(supervisor,k) for k in PARENT_FILES)
    validate_os_docs(start,end,m,pid,group,order,binding,input_hash[rel(supervisor/'manifest.json')])
    cases=[]
    for row in m['workers']:
        case=row['case_id']; folder=folders[case]
        pins={name:input_hash[rel(folder/name)] for name in (*CHILD_META,*CHILD_PAYLOAD)}
        current=[read(folder,n) for n in CHILD_META]
        ident=validate_worker(case,row,current,pins,pid,group,binding,input_hash[rel(supervisor/'started.json')])
        cases.append((case,folder,pins,ident,row['candidate_status']))
    return admitted,cases


def load_exact(path,name):
    module_spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    check(Path(module.__file__).resolve()==path and sha(path)==PINS[path],'changed reviewed audit helper')
    return module


def source_gate(r,p,pf,binding):
    check(r['status']==ROOT_STATUS and yes(r['checks']) and all(r[k]==v for k,v in binding.items())
          and r['exact_spatial_bridge_authorized'] is False
          and pf['status']=='PASS' and yes(pf['checks'])
          and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
          and p['source_sha256']==binding['source_sha256']
          and integer(pf['actual_geometry_payloads_opened'],0) and integer(pf['actual_numeric_or_geometry_jobs'],0),
          'own frozen replay phase and external root authorization required')


def history_box(ctx,trial,box,contract,archived):
    a=ctx.target.a
    reference=[list(map(Q,x)) for x in a['reference'].tolist()]
    old_ref=[list(map(Q,x)) for x in archived['reference_nodes']]
    free=a['free'].tolist()
    curr=list(map(Q,a['u'].tolist()))
    old=list(map(Q,archived['free_displacement']))
    triangles=a['triangles'].tolist()
    rim=Q(archived['parameters']['outer_radius'])
    free_expected=[2*i+d for i,(r,z) in enumerate(reference) for d in (0,1) if r!=rim and (r!=0 or d==1)]
    check(reference==old_ref and len(reference)==325
          and free==free_expected==box['free'] and len(set(free))==len(free)==len(old)==len(curr)==611
          and all(type(k) is int for k in free) and triangles==archived['triangles']==box['triangles'] and len(triangles)==144
          and sha256(a['u'].tobytes(order='C')).hexdigest()==contract['center_binary_sha256']==trial['center_binary_sha256']
          and sha256(a['preconditioner'].tobytes(order='C')).hexdigest()==contract['binary_R_sha256']==trial['R_binary_sha256'],
          'identical exact reference, free injection, binary proof center and R required')
    radii=[[Q(0),Q(0)] for r in reference]
    center=[list(row) for row in reference]
    base=[list(row) for row in reference]
    for i,q,u in zip(free,old,curr,strict=True):
        node,d=divmod(i,2)
        base[node][d]+=q
        center[node][d]+=u
        radii[node][d]=RADIUS
    stored=[list(map(Q,row)) for row in archived['current_nodes']]
    sigma=max(abs(a-b) for a,b in zip(curr,old,strict=True))
    xi=max(abs(v-w) for row,c in zip(base,stored,strict=True) for v,w in zip(row,c,strict=True))
    check(all([Q(float(v)) for v in row]==s for row,s in zip(base,stored,strict=True))
          and sigma==exact(trial['archived_center_shift_exact'])==exact(contract['archived_shift_exact'])
          and xi==exact(contract['archived_addition_gap_exact'])
          and reference==[list(map(exact,r)) for r in box['reference_nodes_exact']]
          and center==[list(map(exact,r)) for r in box['current_centers_exact']]
          and radii==[list(map(exact,r)) for r in box['coordinate_radii_exact']]
          and [[[v-d,v+d] for v,d in zip(row,rad,strict=True)] for row,rad in zip(center,radii,strict=True)]
               ==[[list(map(exact,ends)) for ends in row] for row in box['coordinate_hulls_exact']]
          and exact(box['corrector_radius_exact'])==RADIUS
          and box['root_membership_assumed'] is False and box['any_spatial_derivative_bound_assumed'] is False
          and len(box['predictor_hulls_exact'])==611 and all(tuple(map(exact,ends))==(0,0) for ends in box['predictor_hulls_exact']),
          'one of 650 actual whole-box hulls or archived sigma/xi identities differs')
    check((sigma>0)==(contract['verification_center_kind']=='VERIFIED_RECENTERED_U1'),
          'recentered target may not replace the scored field')
    return reference,center,radii,triangles,{d:i for i,d in enumerate(free)}


def activity_mask(archived,cached,reference,triangles):
    params=archived['parameters']
    patch,zlo,zhi=map(Q,(params['patch_radius'],params['activation_z_min'],params['activation_z_max']))
    expected=[]
    for nodes in triangles:
        coords=[reference[i] for i in nodes[:3]]
        r,z=(sum(p[k] for p in coords)/3 for k in (0,1))
        expected.append(int(r<patch and zlo<z<zhi))
    saved=archived['activation_elements']
    check(len(saved)==len(expected)==len(cached)==144 and all(type(x) is int and x in (0,1) for x in saved)
          and expected==saved==[Q(row[0]) for row in cached] and sum(expected)==24,
          'same indexed geometric 24-element activation mask required')
    return expected


def replay_case(a,bmath,trial,producer,body,status,crosslinks):
    case=trial['case_id']
    ctx,rule,entry=a.open_context(case,PINS[TARGET/'preparation/manifest.json'])
    check(all(trial[k]==v for k,v in entry.items()) and producer['trial']==trial,
          'actual geometry must use exactly the frozen trial')
    def authority(path,h):
        check(crosslinks.setdefault(path,h)==h,'conflicting per-case immutable link')
        return pinned(canon(path),h)
    box=authority(trial['trial_box'],trial['trial_box_sha256'])
    route=authority(trial['routes'],trial['routes_sha256'])
    contract=authority(trial['model_contract'],trial['model_contract_sha256'])
    for key in ('rule','predictor','row','box'):
        authority(trial[key],trial[key+'_sha256'])
    archived=authority(contract['archived_field'],contract['archived_field_sha256'])
    authority(contract['archived_verification'],contract['archived_verification_sha256'])
    check(ctx.pilot_id==case==contract['pilot_id'] and trial['record_id']==ctx.record_id==contract['record_id']
          and exact(trial['exact_drive'])==ctx.drive==ctx.dose*ctx.base['common_gain']
          and trial['point_cache_sha256']==ctx.row['point_cache_sha256']==contract['point_cache_sha256']
          and exact(trial['energy_scale_binary_exact'])==ctx.target.energy_scale>0
          and trial['model_contract_sha256']==ctx.row['model_observation_contract_sha256'],
          'same source-driven exact singleton context')
    ref,u,rad,tri,free=history_box(ctx,trial,box,contract,archived)
    active=activity_mask(archived,ctx.target.a['active'].tolist(),ref,tri)
    body_new=a.X.G.certify_box(ref,u,tri,rad)
    check(body_new==body,'fresh whole-box Bernstein body replay disagrees')
    summary={'pilot_id':case,'record_id':ctx.record_id,'radius_exact':str(RADIUS),'box_sha256':trial['trial_box_sha256'],
             'body_geometry_sha256':producer['whole_body_sha256'],'model_contract_sha256':trial['model_contract_sha256'],
             'point_cache_sha256':trial['point_cache_sha256'],'center_binary_sha256':trial['center_binary_sha256'],
             'R_binary_sha256':trial['R_binary_sha256'],'exact_drive':trial['exact_drive'],
             'energy_scale_binary_exact':trial['energy_scale_binary_exact'],'archived_field_sha256':contract['archived_field_sha256'],
             'archived_shift_exact':contract['archived_shift_exact'],'archived_addition_gap_exact':contract['archived_addition_gap_exact'],
             'archived_observer_goal_exact':contract['observer_goal_exact'],'derivative_policy':trial['derivative_policy'],
             'root_asserted':False,'scope':SCOPE}
    if body_new['status']!='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING':
        check(status==CASE_OTHER[0] and not producer['bulk_geometry_only'] and not producer['trace_geometry_only'],
              'unresolved embedding cannot return positive domain claims')
        return {**summary,'status':'UNRESOLVED_INDEPENDENT_WHOLE_BODY_GEOMETRY','reason':body_new['status']}
    check(len(body_new['triangles'])==len(route['bulk'])==len(producer['bulk_geometry_only'])==144,
          '144 uniquely matched bulk routes required')
    s=a.E.SUR
    nodes=s.common_nodes(ctx.target,box)
    charts=Counter(); bulk_rhos=[]; domains=[]; failures=[]
    for i,(f,out,ids,z) in enumerate(zip(route['bulk'],producer['bulk_geometry_only'],tri,body_new['triangles'],strict=True)):
        calc,geom_ids=a.E._bulk_fields(ctx.target,nodes,i)
        chart=CHARTS[sum(ref[j][0]==0 for j in ids[:3])]
        charts[chart]+=1
        check(type(f['element']) is int and f['element']==i and f['nodes']==geom_ids==ids
              and f['chart']==calc['chart']['kind']==chart
              and f['auxiliary_element_rule_sha256']==a.E._bulk_rule(ctx.target,rule,i)
              and (f['free_local'],f['global_indices'])==bmath.proj(ids,free)==(calc['free_local'],calc['global_indices'])
              and exact(f['active'])==calc['active']==active[i]
              and all(out[k]==v for k,v in f.items()),'bulk rule, chart or free projection drift')
        j=bmath.pair((z['J_bound']['lower_exact'],z['J_bound']['upper_exact']))
        domains.append(j)
        raw=max(abs(v-1) for v in j)
        rounding=2**48
        rho=Q(ceil(raw*rounding),rounding)
        if not (j[0]>Q(1,10**6) and rho<1):
            check(out['status']=='UNRESOLVED_BULK_DOMAIN','nonanalytic bulk enclosure reported complete')
            failures.append(f'bulk:{i}')
        else:
            check(bmath.pair(out['J_bounds_exact'])==j and exact(out['raw_rho_exact'])==raw
                  and raw<=exact(out['rho_exact'])==rho and rho<raw+Q(1,rounding),
                  'body J and outward bulk rho differ')
            bulk_rhos.append(rho)
    check(charts==dict(zip(CHARTS,(132,6,6),strict=True)),'all three reviewed P2 chart classes required')
    check(len(route['surface'])==len(producer['trace_geometry_only'])==24,'24 matched traces required')
    surface_rhos=[]; gmins=[]; hoops=[]; tangents=[]; nleaves=0
    for slot,(row,out) in enumerate(zip(route['surface'],producer['trace_geometry_only'],strict=True)):
        face,edge=divmod(slot,12)
        ids=ctx.target.a[f'surface_{face}_edges'][edge].tolist()
        calc=ctx.target.surface(face,edge,nodes=nodes)
        a.E._surface_rule(ctx.target,rule,face,edge)
        check(row['nodes']==ids and integer(row['face'],face) and integer(row['edge'],edge)
              and exact(row['Rlo_exact'])==calc['Rlo']==ref[ids[0]][0]
              and exact(row['Rhi_exact'])==calc['Rhi']==ref[ids[-1]][0]
              and row['trace_rule_sha256']==digest(a.doc([rule['surface_s'][face],rule['surface_weights'][face][edge]]))
              and (row['free_local'],row['global_indices'])==bmath.proj(ids,free)==(calc['free_local'],calc['global_indices'])
              and all(out[k]==v for k,v in row.items()),'surface trace chart, auxiliary rule or free map differs')
        try:
            d=s.uniform_G_bounds(calc,max_depth=8)
            if d['status']!='CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN':
                raise ValueError(d['status'])
            # The normal branch reconstructs independently from frozen exact site maps.
            if 'domain' not in out:
                for iq in range(12):
                    s.local_site_domain(ctx.target,box,face,edge,iq)
        except ValueError:
            check(out['status']=='UNRESOLVED_SURFACE_DOMAIN','nonanalytic trace marked complete')
            failures.append(f'surface:{face}:{edge}')
            continue
        check(d==out['domain'] and len(out['cache_site_geometry'])==12,'fresh trace leaves or 12 site domains differ')
        g,_,rho,hoop,n=bmath.surface_bounds(d)
        gmins.append(g); surface_rhos.append(rho); hoops.append(hoop); nleaves+=n
        for iq,site in enumerate(out['cache_site_geometry']):
            tangents.append(bmath.independent_site(s,ctx.target,box,face,edge,iq,ids,row['global_indices'],site))
    if failures:
        check(status==CASE_OTHER[1],'terminal status loses an unresolved analytic component')
        return {**summary,'status':'UNRESOLVED_INDEPENDENT_WHOLE_BOX_ANALYTIC_SUBDOMAIN','unresolved_components':failures}
    independent_j,j,margin=bmath.check_body(body_new)
    check(status==CASE_PASS and independent_j==domains and len(bulk_rhos)==144
          and len(surface_rhos)==24 and len(tangents)==288,'full geometry candidate coverage required')
    return {**summary,'status':'PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY',
            'all_144_body_and_24_trace_domains':True,'raw_surface_cache_geometry_lines':288,
            'body_J_bounds_exact':list(map(str,j)),'strong_monotonicity_lower_exact':str(margin),
            'max_bulk_rho_outward_exact':str(max(bulk_rhos)),'max_surface_rho_raw_exact':str(max(surface_rhos)),
            'min_surface_G_exact':str(min(gmins)),'min_trace_hoop_exact':str(min(hoops)),
            'min_site_tangent_squared_exact':str(min(tangents)),'surface_Bernstein_leaves':nleaves,
            'chart_counts':dict(charts),'closed_full_nodal_coordinates_rebuilt':650}


def preflight_doc(ph,checks):
    return {'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,'checks':checks,
            'actual_geometry_payloads_opened':0,'actual_numeric_or_geometry_jobs':0}


def rejects(action):
    try:
        action()
    except (ValueError,AssertionError):
        return True
    return False


def center_fixture():
    from types import SimpleNamespace

    import numpy as np
    original=np.asarray([(i/24,j/12) for j in range(13) for i in range(25)],dtype='<f8')
    width=Q(1)
    ref=[list(map(Q,p)) for p in original.tolist()]
    free=[2*i+d for i,(r,z) in enumerate(ref) for d in (0,1) if r!=width and (r!=0 or d==1)]
    check(len(free)==611,'synthetic axis/rim map has wrong size')
    original_displacements=np.zeros(611,dtype='<f8')
    proposed=original_displacements.copy(); proposed[0]=2.0**-50
    small_R=np.ones((1,1),dtype='<f8')
    tri=[[0,1,2,3,4,5]]*144
    a={'reference':original,'free':np.asarray(free),'u':proposed,'triangles':np.asarray(tri),'preconditioner':small_R}
    ctx=SimpleNamespace(target=SimpleNamespace(a=a))
    c={'archived_shift_exact':str(Q(1,2**50)),'archived_addition_gap_exact':'0',
       'verification_center_kind':'VERIFIED_RECENTERED_U1',
       'center_binary_sha256':sha256(proposed.tobytes()).hexdigest(),
       'binary_R_sha256':sha256(small_R.tobytes()).hexdigest()}
    t={'archived_center_shift_exact':c['archived_shift_exact'],'center_binary_sha256':c['center_binary_sha256'],
       'R_binary_sha256':c['binary_R_sha256']}
    archived={'reference_nodes':original.tolist(),'current_nodes':original.tolist(),
              'triangles':tri,'free_displacement':original_displacements.tolist(),'parameters':{'outer_radius':1.0}}
    center=[list(row) for row in ref]
    n,k=divmod(free[0],2); center[n][k]+=Q(1,2**50)
    radii=[[RADIUS if 2*i+j in free else Q(0) for j in (0,1)] for i in range(325)]
    box={'free':free,'triangles':tri,'reference_nodes_exact':ref,'current_centers_exact':center,
         'coordinate_radii_exact':radii,'coordinate_hulls_exact':[[[q-d,q+d] for q,d in zip(row,r,strict=True)]
                              for row,r in zip(center,radii,strict=True)],'corrector_radius_exact':str(RADIUS),
         'root_membership_assumed':False,'any_spatial_derivative_bound_assumed':False,'predictor_hulls_exact':[[0,0]]*611}
    history_box(ctx,t,box,c,archived)
    check(rejects(lambda:history_box(ctx,{**t,'archived_center_shift_exact':'0'},box,c,archived))
          and rejects(lambda:history_box(ctx,t,{**box,'coordinate_hulls_exact':box['coordinate_hulls_exact'][:-1]},c,archived)),
          'recenter reset or incomplete box kept false proof authority')


def cold_fixtures(p,producer_cases,ph):
    math=load_exact(OLD_AUDIT,'prior_geometry_endpoint_audit_cold')
    check(math.product((-2,4),(3,5))==(-10,20) and math.square((-2,4))==(0,16)
          and math.combination((Q(1,2),-2),((-2,4),(3,5)))==(-11,-4),
          'retained independent exact endpoint arithmetic changed')
    center_fixture()
    active_nodes=[[Q(1,4),Q(1,2)]]*6
    inactive_nodes=[[Q(3,4),Q(1,2)]]*6
    mask=[int(i%6==2) for i in range(144)]
    mask_archived={'parameters':{'patch_radius':Q(1,2),'activation_z_min':Q(0),'activation_z_max':Q(1)},
                   'activation_elements':mask}
    triangles=[list(range(6)) if k else list(range(6,12)) for k in mask]
    binary_mask=[[float(x)] for x in mask]
    check(activity_mask(mask_archived,binary_mask,active_nodes+inactive_nodes,triangles)==mask
          and rejects(lambda:activity_mask({**mask_archived,'activation_elements':list(range(24))},
                       binary_mask,active_nodes+inactive_nodes,triangles))
          and rejects(lambda:activity_mask({**mask_archived,'activation_elements':mask[1:]+mask[:1]},
                       binary_mask,active_nodes+inactive_nodes,triangles)),
          'indexed indicator cannot be confused with active IDs or a shifted indicator')
    dummy_pf='f'*64
    binding={'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':dummy_pf,
             'producer_root_admission_sha256':ADMISSION_SHA}
    pf=preflight_doc(ph,{'exact_fixture':True})
    r={'status':ROOT_STATUS,**binding,'checks':{'mock_root_gate':True},'exact_spatial_bridge_authorized':False}
    source_gate(r,p,pf,binding)
    for bad in ({**r,'preflight_sha256':'e'*64},{**r,'checks':{'truthy':1}},{**r,'exact_spatial_bridge_authorized':True}):
        check(rejects(lambda:source_gate(bad,p,pf,binding)),'wrong actual root phase allowed')
    first=producer_cases[0]
    case,folder,pins,ident,status=first
    worker=[pinned(folder/f,pins[f]) for f in CHILD_META]
    os_row=worker[2]
    terminal_index={**os_row,'os_exit_sha256':pins['observed_os_exit.json'],
                    'results_sha256':pins['results.json'],'candidate_status':status}
    group=p['launch_group']; binding=p['producer_binding']; parent=ident['supervisor_pid']
    start_h=worker[0]['command'][13]
    validate_worker(case,terminal_index,worker,pins,parent,group,binding,start_h)
    for replaced,index in (({'pid':float(ident['pid'])},0),({'supervisor_pid':float(parent)},1),
                           ({'observed_exit_code':False},2),({'exit_code':0.0},3),({'case_id':'wrong'},3)):
        changed=list(worker); changed[index]={**changed[index],**replaced}
        check(rejects(lambda:validate_worker(case,terminal_index,changed,pins,parent,group,binding,start_h)),
              'mixed PID, terminal or loose numeric alias accepted')
    return {'exact_interval_and_chart_rounding_reviewed':True,
            'synthetic_325_node_recenter_and_frozen_archive_gap_rebuilt':True,
            'indexed_source_144_activation_mask_and_distinct_24_active_ids_reconstructed':True,
            'root_source_gate_positive_negative_actual_schema_roundtrip':True,
            'all_32_real_OS_and_exact_worker_command_observations_match_before_payload_access':True,
            'changed_child_phase_float_PID_boolean_exit_and_wrong_case_rejected':True}


def prepare():
    check(not REPLAY.exists(),'cold phase cannot touch prior actual replay')
    p=plan()
    admission,cases=observed_inputs(ADMISSION,ADMISSION_SHA,p)
    # Figure the exact future protocol hash without writing it or any worker payload.
    raw=(json.dumps(p,indent=2,allow_nan=False)+'\n').encode('utf-8')
    ph=sha256(raw).hexdigest()
    checks=cold_fixtures(p,cases,ph)
    with (HERE/'replay_geometry_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    check(save(HERE/'protocol.json',p)==ph,'cold serialization hash differs')
    pf=save(HERE/'preflight.json',{**preflight_doc(ph,checks),'observed_metadata_cases':len(cases),
                                 'observed_producer_supervisor_pid':admission['supervisor_pid']})
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':pf,
                      'producer_cases_metadata_bound':len(cases),'payloads_opened':0}),flush=True)


def run(args):
    check(sys.flags.optimize==0,'prior exact independent assertions must remain enabled')
    p=pinned(HERE/'protocol.json',args.protocol_sha256)
    pf=pinned(HERE/'preflight.json',args.preflight_sha256)
    check(args.source_review.resolve()==HERE/'root_source_review.json','canonical root source review path required')
    root_source=pinned(args.source_review,args.source_review_sha256)
    binding={'source_sha256':sha(SOURCE),'protocol_sha256':args.protocol_sha256,'preflight_sha256':args.preflight_sha256,
             'producer_root_admission_sha256':args.admission_sha256}
    source_gate(root_source,p,pf,binding)
    check(p==plan() and sha(HERE/'replay_geometry_frozen.py')==binding['source_sha256'],'changed cold replay')
    admission,cases=observed_inputs(args.admission,args.admission_sha256,p)
    bmath=load_exact(OLD_AUDIT,'prior_geometry_endpoint_audit_actual')
    sys.path.insert(0,str(ROOT/'scripts'))
    import algebraic_root_feasibility as core
    import review_algebraic_root as maps
    import transfer_ambiguity_remaining32_exact_targets as a
    check(Path(a.__file__).resolve()==T_SOURCE and sha(T_SOURCE)==PINS[T_SOURCE]
          and all(p['kernel_closure_sha256'][rel(Path(obj.__file__).resolve())]==sha(Path(obj.__file__).resolve())
                  for obj in (core,maps,a.X.G,a.E,a.E.SUR,a.E.AX)),
          'same reviewed source, exact target and geometry or AD guard aliases required')
    checked={**p['frozen_input_sha256'],**p['kernel_closure_sha256'],**admission['sha256'],
             rel(ADMISSION):args.admission_sha256,rel(args.source_review):args.source_review_sha256}
    phase={'pid':os.getpid(),**binding,'source_review_sha256':args.source_review_sha256,
           'producer_supervisor_pid':admission['supervisor_pid'],'producer_launch_group':admission['launch_group']}
    REPLAY.mkdir(exist_ok=False)
    save(REPLAY/'started.json',{'status':'RUNNING',**phase})
    terminal={'status':'FAILED','exit_code':1,**phase}
    t0=perf_counter()
    try:
        records=[]
        with ExitStack() as guard:
            guard.enter_context(a.X.no_integrals())
            for module,key in ((core,'variation_rows'),(maps,'terms_from_arrays'),(core,'bulk_local'),(core,'surface_local')):
                guard.enter_context(patch.object(module,key,side_effect=AssertionError('fresh cache, AD or spatial force prohibited')))
            trials=read_pre(GEOMETRY/'protocol.json')['trials']
            for (case,folder,pins,ident,status),trial in zip(cases,trials,strict=True):
                producer=pinned(folder/'results.json',pins['results.json'])
                check(all(producer[k]==v for k,v in ident.items()) and integer(producer['pid'],ident['pid'])
                      and integer(producer['supervisor_pid'],ident['supervisor_pid']) and case==trial['case_id']
                      and producer['status']==status and producer['whole_body_sha256']==pins['body_geometry.json']
                      and integer(producer['new_AD_derivatives'],0) and integer(producer['new_spatial_integrations'],0)
                      and producer['root_asserted'] is False,'current exact geometry-only payload authority required')
                body=pinned(folder/'body_geometry.json',pins['body_geometry.json'])
                evidence=replay_case(a,bmath,trial,producer,body,status,checked)
                native={**evidence,'producer_case_results_sha256':pins['results.json'],
                        'producer_case_body_file':rel(folder/'body_geometry.json'),
                        'producer_case_results_file':rel(folder/'results.json'),
                        'producer_worker_terminal_sha256':pins['terminal.json'],
                        'producer_worker_observed_os_exit_sha256':pins['observed_os_exit.json'],
                        'replay_phase':phase}
                path=REPLAY/f'{case}.json'
                h=save(path,native)
                records.append({k:native[k] for k in ('pilot_id','record_id','radius_exact','box_sha256','body_geometry_sha256','status')}
                               |{'receipt':rel(path),'sha256':h})
                print(json.dumps({'case_id':case,'status':native['status'],'geometry_sha256':h}),flush=True)
        completed=sum(r['status']=='PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY' for r in records)
        check(len(records)==32 and [r['pilot_id'] for r in records]==p['ordered_ids'],'complete replay required')
        for name,h in checked.items():
            check(sha(canon(name))==h,'producer, frozen context or mathematical source changed during replay')
        result={'status':RESULT_STATUS if completed==32 else 'UNRESOLVED_INDEPENDENT_REMAINING32_GEOMETRY_REPLAY',
                'source_and_phase':phase,'scope':SCOPE,'geometry_records':records,
                'successfully_replayed_geometry_cases':completed,'total_geometry_cases':len(records),
                'independent_all_32_geometry':completed==32,'checked_input_sha256':checked,
                'observed_producer_os':{'pid':admission['supervisor_pid'],'session':SESSION,'exit_code':0,'observer':'root'},
                'exact_spatial_bridge_authorized':False,'exact_spatial_roots_proved':0,'root_asserted':False,
                'new_spatial_integral_force_Hessian_AD_or_solve_jobs':0,'seconds':perf_counter()-t0}
        rh=save(REPLAY/'receipt.json',result)
        terminal.update(status='FINISHED',exit_code=0,receipt_sha256=rh)
        print(json.dumps({'status':result['status'],'receipt_sha256':rh,'completed_cases':completed}),flush=True)
    except BaseException as error:
        save(REPLAY/'failure.json',{'error':str(error),'traceback':traceback.format_exc(),**phase})
        raise
    finally:
        save(REPLAY/'terminal.json',terminal)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_subparsers(dest='stage',required=True)
    actions.add_parser('prepare')
    launcher=actions.add_parser('run')
    for key in ('protocol-sha256','preflight-sha256','source-review-sha256','admission-sha256'):
        launcher.add_argument('--'+key,required=True)
    launcher.add_argument('--source-review',required=True,type=lambda q:Path(q).resolve())
    launcher.add_argument('--admission',required=True,type=lambda q:Path(q).resolve())
    args=parser.parse_args()
    if args.stage=='prepare':
        prepare()
    else:
        run(args)


if __name__=='__main__':
    main()
