"""Independently replay exact symmetric-part congruence estimates for 32 frozen matrices."""
import argparse
from base64 import b64decode
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import re
import sys
from time import perf_counter
import traceback

import numpy as np


HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve()
OUT=HERE.parent
EXACT=OUT.parent
ROOT=EXACT.parents[3]
PRODUCED=OUT/'actual'
REPLAY=HERE/'actual_replay'
TARGET=EXACT/'remaining32_targets_v1'
JOIN=EXACT/'remaining32_bridge_v1/admission'
CALLER=ROOT/'scripts/transfer_ambiguity_remaining32_preconditioner_spd.py'
MATH=ROOT/'scripts/review_cohort_energy_stability.py'
ADMISSION=PRODUCED/'root_terminal_admission.json'
EXEC_SESSION=42727
ADMIT_STATUS='PASS_ROOT_OBSERVED_REMAINING32_SPD_TERMINAL'
REVIEW_STATUS='PASS_ROOT_REMAINING32_SPD_REPLAY_SOURCE'
PROD_PASS='CANDIDATE_EXACT_SYMMETRIC_PART_SPD_AWAITING_INDEPENDENT_REVIEW'
PROD_UNRESOLVED='UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_TEST'
CASE_PASS='PASS_INDEPENDENT_REMAINING32_EXACT_SYMMETRIC_PART_SPD'
RESULT_PASS='PASS_INDEPENDENT_ACTUAL_REMAINING32_EXACT_SYMMETRIC_PARTS'
RANK=611
PRECISION=768
PRODUCER_PRECISION=640
SCOPE=('A separate exact-binary symmetric-part SPD premise for each fixed preconditioner. '
       'Reuses the previously reviewed alternate endpoint association with independently decoded R and K. '
       'No physical force, Hessian, AD, equilibrium, exact-spatial root or minimum calculation.')
PINS={
    CALLER:'b28a1daae0bec20b00423b218a2eb2efd759d91afdd52eb013c04405da55d3d3',
    OUT/'preparation/producer_frozen.py':'b28a1daae0bec20b00423b218a2eb2efd759d91afdd52eb013c04405da55d3d3',
    OUT/'preparation/protocol.json':'310fca3da31ba8e5c8bdf4b21cc13c51559bb316ee7cb60f7b6573216b92355c',
    OUT/'preparation/preflight.json':'07adb3539dfa9e3062680e4efb713e9657eb64b6e38445ac683c395016e13270',
    OUT/'preparation/independent_prelaunch/review.json':'88975bf6ce01b3c20c5cd30f8334a9e5831745c8b42b4d8212a0d1df5b68d0fe',
    MATH:'4c2acbc96ce3f19954b55d14fde46a8532536e2d39c2250e0123a17b78663c73',
    TARGET/'preparation/manifest.json':'dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242',
    TARGET/'independent_prelaunch/review.json':'67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9',
    JOIN/'bridge_input_admission.json':'2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb',
    EXACT/'preconditioner_spd_v1/independent_actual_preparation/preparation_v2/replay_frozen.py':
        '4cefc71d02ae7a735d1e5abef28b2386ad4922652bed108e4e5d244658b87699',
    EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json':
        '769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
}
PARENT=('started.json','terminal.json','manifest.json')


def require(condition,message):
    if not condition:
        raise ValueError(message)


def rel(path):
    return path.relative_to(ROOT).as_posix()


def canon(name):
    path=(ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and rel(path)==name,'canonical in-study authority link required')
    return path


def sha(path):
    h=sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def is_hash(s):
    return type(s) is str and re.fullmatch('[0-9a-f]{64}',s) is not None


def read(path,expected,*,parse=True):
    require(is_hash(expected),f'explicit SHA required: {path}')
    if not parse:
        require(sha(path)==expected,f'changed exact input: {path}')
        return None
    raw=path.read_bytes()
    require(sha256(raw).hexdigest()==expected,f'changed exact input: {path}')
    return json.loads(raw)


def literal(x,n=None):
    return type(x) is int and (x>0 if n is None else x==n)


def yes(x):
    return type(x) is dict and bool(x) and all(v is True for v in x.values())


def exact(x):
    require(type(x) in (str,int,Q),'finite exact rational required')
    return Q(x)


def digest(x):
    return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def save(path,data):
    with path.open('x',encoding='utf8',newline='\n') as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def binding_for_producer():
    return {'source_sha256':PINS[CALLER],'protocol_sha256':PINS[OUT/'preparation/protocol.json'],
            'preflight_sha256':PINS[OUT/'preparation/preflight.json'],
            'target_manifest_sha256':PINS[TARGET/'preparation/manifest.json'],
            'target_context_review_sha256':PINS[TARGET/'independent_prelaunch/review.json'],
            'joined_bridge_input_sha256':PINS[JOIN/'bridge_input_admission.json'],
            'independent_prelaunch_review_sha256':PINS[OUT/'preparation/independent_prelaunch/review.json']}


def static_plan():
    for path,h in PINS.items():
        read(path,h,parse=False)
    docs=lambda path:read(path,PINS[path])
    protocol=docs(OUT/'preparation/protocol.json'); preflight=docs(OUT/'preparation/preflight.json')
    approval=docs(OUT/'preparation/independent_prelaunch/review.json'); target=docs(TARGET/'preparation/manifest.json')
    tr=docs(TARGET/'independent_prelaunch/review.json'); join=docs(JOIN/'bridge_input_admission.json')
    old=docs(EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json')
    bindings=binding_for_producer()
    ids=protocol['ordered_case_ids']
    require(protocol['status']=='FROZEN_REMAINING32_EXACT_SYMMETRIC_PRECONDITIONER_PART_PROPOSALS'
            and approval['status']=='PASS_INDEPENDENT_REMAINING32_PRECONDITIONER_SPD_PRELAUNCH' and yes(approval['checks'])
            and all(approval[k]==h for k,h in bindings.items() if k!='independent_prelaunch_review_sha256')
            and all(protocol[k]==bindings[k] for k in ('source_sha256','target_manifest_sha256','target_context_review_sha256','joined_bridge_input_sha256'))
            and approval['root_or_minimum_claimed'] is False and protocol['root_or_minimum_claimed'] is False
            and protocol['single_serial_process'] is True and protocol['separate_independent_actual_replay_required'] is True
            and literal(protocol['rank'],RANK) and literal(protocol['precision_bits'],PRODUCER_PRECISION)
            and preflight['status']=='PASS' and yes(preflight['checks']) and literal(preflight['new_production_matrix_jobs'],0)
            and all(preflight[k]==bindings[k] for k in ('source_sha256','protocol_sha256'))
            and tr['status']=='PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS' and yes(tr['checks'])
            and join['status']=='PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS'
            and join['spatial_workers_authorized'] is False and yes(join['checks'])
            and len(set(ids))==len(ids)==32 and ids==approval['permitted_case_ids']==tr['ordered_case_ids']
            ==[e['case_id'] for e in target['targets']]==[e['case_id'] for e in join['cases']]
            ==[r['pilot_id'] for r in protocol['ordered_targets']]
            and old['status']=='PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD',
            'same complete frozen SPD producer and accepted proof-class authority required')
    return {'status':'FROZEN_INDEPENDENT_ACTUAL_REMAINING32_SPD_REPLAY','scope':SCOPE,'source_sha256':sha(SOURCE),
            'input_sha256':{rel(path):h for path,h in PINS.items()},'producer_bindings':bindings,
            'rank':RANK,'precision_bits':PRECISION,'ordered_case_ids':ids,'external_producer_observed_session':EXEC_SESSION,
            'required_producer_admission_status':ADMIT_STATUS,'required_root_source_review_status':REVIEW_STATUS,
            'actual_producer_result_or_matrix_payloads_opened_in_preparation':0,
            'actual_matrix_jobs_in_preparation':0}


def verify_parent(admission,starter,terminal,manifest,manifest_h,plan):
    pid=admission['pid']
    require(admission['status']==ADMIT_STATUS and literal(pid) and literal(admission['observed_exit_code'],0)
            and literal(admission['unified_exec_session'],plan['external_producer_observed_session'])
            and all(literal(row['pid'],pid) and row['bindings']==plan['producer_bindings']
                    and row['ordered_case_ids']==plan['ordered_case_ids'] for row in (starter,terminal,manifest))
            and starter['status']=='RUNNING' and terminal['status']=='FINISHED' and literal(terminal['exit_code'],0)
            and terminal['manifest_sha256']==manifest_h==admission['manifest_sha256']
            and manifest['status']=='COMPLETE_REMAINING32_SYMMETRIC_PART_ATTEMPTS_AWAITING_INDEPENDENT_REPLAY'
            and manifest['externally_observed_exit_required'] is True and manifest['external_independent_actual_review_required'] is True
            and manifest['root_or_minimum_claimed'] is False
            and [r['case_id'] for r in manifest['cases']]==plan['ordered_case_ids'],
            'same complete serial native process and genuine externally observed success required')


def admissible_paths(rows):
    all_paths={rel(PRODUCED/n) for n in PARENT}
    for row in rows:
        planned=PRODUCED/row['case_id']
        f=rel(planned/'result.json')
        require(row['file']==f and is_hash(row['sha256']) and row['status'] in (PROD_PASS,PROD_UNRESOLVED),
                'complete exact bounded result catalogue required')
        all_paths.add(f)
        k=row['K_file_sha256']
        if k is None:
            require(row['status']==PROD_UNRESOLVED,'missing proposal may not be a certified candidate')
        else:
            require(is_hash(k),'invalid proposal hash'); all_paths.add(rel(planned/'congruence_K.npy'))
    return all_paths


def admitted(path,h,plan):
    require(path.resolve()==ADMISSION,'canonical root terminal observation required')
    a=read(path,h)
    links=a['sha256']
    require(type(links) is dict and all(is_hash(v) for v in links.values())
            and links[rel(PRODUCED/'terminal.json')]==a['terminal_sha256']
            and links[rel(PRODUCED/'manifest.json')]==a['manifest_sha256'], 'paired native terminal and manifest hashes required')
    starter,terminal,manifest=(read(PRODUCED/k,links[rel(PRODUCED/k)]) for k in PARENT)
    verify_parent(a,starter,terminal,manifest,links[rel(PRODUCED/'manifest.json')],plan)
    rows=manifest['cases']
    require(len(rows)==32 and set(links)==admissible_paths(rows) and len({r['case_id'] for r in rows})==32
            and all(links[r['file']]==r['sha256'] and (r['K_file_sha256'] is None or
                    links[rel(PRODUCED/r['case_id']/'congruence_K.npy')]==r['K_file_sha256']) for r in rows)
            and literal(manifest['candidate_spd_count'],sum(r['status']==PROD_PASS for r in rows))
            and literal(manifest['unresolved_count'],sum(r['status']==PROD_UNRESOLVED for r in rows)),
            'exact 32-case archive set and statuses required before opening any result, R or K')
    return a,rows


def matrix_bytes(a,n):
    require(type(a) is np.ndarray and a.shape==(n,n) and a.dtype.str=='<f8' and np.isfinite(a).all(),
            'finite square exact binary64 array required')
    return a.tobytes(order='C')


def direct_R_norm(R):
    return max(sum(abs(Q.from_float(float(x))) for x in row) for row in R)


def payload(packet):
    p=packet['frozen_binary_cache']
    raw=b64decode(p['data'],validate=True)
    require(literal(p['bytes'],len(raw)) and sha256(raw).hexdigest()==p['sha256'],'authenticated cache payload mismatch')
    with np.load(BytesIO(raw),allow_pickle=False) as files:
        require(len(set(files.files))==len(files.files) and all(k in files.files for k in ('preconditioner','u','free','reference','triangles')),
                'duplicate or missing cache tensors')
        tensors={k:files[k] for k in files.files}
    require({k:list(v.shape) for k,v in tensors.items()}==p['array_shapes']
            and all(a.dtype.kind in 'fibu' and np.isfinite(a).all() for a in tensors.values()),
            'typed archive tensor shapes and finite entries required')
    return tensors


def byte_equal(a,b):
    return a.dtype==b.dtype and a.shape==b.shape and a.tobytes(order='C')==b.tobytes(order='C')


def matched_tensor_except_u(old,new):
    return set(old)==set(new) and all(old[k].dtype==new[k].dtype and old[k].shape==new[k].shape
                                      and (k=='u' or byte_equal(old[k],new[k])) for k in old)


def identity_and_R(entry,joined,row,local_links):
    def linked(name,h):
        require(local_links.setdefault(name,h)==h,'immutable input pin conflict')
        return read(canon(name),h)
    pid=row['pilot_id']; target=row['identity']
    c=linked(entry['model_contract'],entry['model_contract_sha256'])
    old=linked(c['original_point_cache'],c['original_point_cache_sha256'])
    field=linked(c['archived_field'],c['archived_field_sha256'])
    linked(c['archived_verification'],c['archived_verification_sha256'])
    native=linked(c['native_case'],c['native_case_sha256'])
    archived=payload(old)
    displacement_hash=sha256(archived['u'].tobytes(order='C')).hexdigest()
    s=exact(c['archived_shift_exact']); xi=exact(c['archived_addition_gap_exact'])
    require(c['pilot_id']==pid==entry['case_id']==joined['case_id'] and s>=0 and xi>=0
            and target['record_id']==entry['record_id']==c['record_id']==joined['record_id']
            and row['archived_shift_exact']==c['archived_shift_exact'] and row['archived_addition_gap_exact']==c['archived_addition_gap_exact']
            and old['status']=='FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE'
            and c['archived_field']==old['input_state']
            and c['archived_field_sha256']==old['input_state_sha256']
            and c['archived_verification_sha256']==old['field_verification_sha256']
            and byte_equal(archived['u'],np.asarray(field['free_displacement'],dtype='<f8'))
            and byte_equal(archived['reference'],np.asarray(field['reference_nodes'],dtype='<f8'))
            and archived['triangles'].tolist()==field['triangles']
            and tuple(archived['reference'].shape)==(325,2) and len(archived['triangles'])==144
            and old['parameters']==field['parameters'] and literal(old['free_dof'],RANK),
            'same source-driven original archive field, mesh and parameters required')
    rim=Q(field['parameters']['outer_radius'])
    expected_free=[2*i+k for i,(r,z) in enumerate(archived['reference'].tolist()) for k in (0,1)
                   if Q(r)!=rim and (Q(r)!=0 or k==1)]
    require(archived['free'].shape==(RANK,) and archived['free'].dtype.kind in 'iu'
            and archived['free'].tolist()==expected_free and len(expected_free)==RANK,
            'same exact rank and independent ordered free injection required')
    if s:
        packet=linked(c['point_cache'],c['point_cache_sha256'])
        state=payload(packet)
        step=linked(rel(canon(c['point_cache']).with_name('exact_binary_step.json')),c['recenter_exact_step_sha256'])
        actual_sigma=max(abs(Q.from_float(a)-Q.from_float(b)) for a,b in zip(
                         state['u'].tolist(),archived['u'].tolist(),strict=True))
        require(packet['status']=='FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER'
                and set(packet)==set(old)|{'recenter'} and all(packet[k]==v for k,v in old.items() if k not in ('status','frozen_binary_cache'))
                and matched_tensor_except_u(archived,state)
                and c['verification_center_kind']=='VERIFIED_RECENTERED_U1'
                and step['case_id']==pid and step['u0_sha256']==displacement_hash
                and step['u1_sha256']==c['center_binary_sha256'] and exact(step['sigma_exact'])==s==actual_sigma
                and packet['recenter']['old_field_is_not_new_proof_center'] is True
                and native['status']=='PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET'
                and native['archived_field_sha256']==c['archived_field_sha256']
                and native['actual_binary_center_shift_sigma_exact']==c['archived_shift_exact']
                and native['archived_addition_gap_xi0_exact']==c['archived_addition_gap_exact'],
                'recentered proof must retain the same R, fields, original shift and every operator tensor other than u')
    else:
        packet,state=old,archived
        require(c['point_cache']==c['original_point_cache'] and c['point_cache_sha256']==c['original_point_cache_sha256']
                and c['recenter_exact_step_sha256'] is None,'unmoved original proof retains its exact cache')
    R=state['preconditioner']; rb=sha256(matrix_bytes(R,RANK)).hexdigest()
    center=sha256(state['u'].tobytes(order='C')).hexdigest()
    dose,gain,drive,norm=(exact(packet[k]) for k in ('optical_dose_exact','common_gain_exact','activation_strength_exact','normalization_binary_exact'))
    cell=c['parameter_cell']
    require(0<drive==dose*gain==exact(c['exact_requested_drive']['activation_strength_exact'])
            and dose==exact(c['exact_requested_drive']['optical_dose_exact'])
            and gain==exact(c['exact_requested_drive']['common_gain_exact'])
            and drive!=exact(packet['activation_strength_binary_exact'])==Q(field['activation_strength'])
            and all(tuple(map(exact,v))==(gain,gain) if k=='common_gain' else
                    tuple(map(exact,v))==(Q(packet['parameters'][k]),)*2 for k,v in cell.items())
            and 0<norm==exact(c['energy_scale_binary_exact'])==Q(float(state['energy_scale'][0]))
            and digest(cell)==c['parameter_cell_sha256']
            and row['R_source_file']==c['point_cache'] and row['R_source_file_sha256']==c['point_cache_sha256']
            and c['frozen_binary_cache_sha256']==packet['frozen_binary_cache']['sha256']
            and rb==c['binary_R_sha256'] and center==c['center_binary_sha256']
            and joined['cache_packet_sha256']==c['point_cache_sha256']
            and joined['trial_box_sha256']==row['trial_box_sha256']
            and joined['model_contract_sha256']==entry['model_contract_sha256']
            and joined['cached_derivative_case_sha256']==row['cached_derivative_case_sha256']
            and joined['geometry_replay_case_sha256']==row['geometry_replay_case_sha256'],
            'same exact binary matrix, archived field, prescribed optical input and physical operator required')
    independently={'pilot_id':pid,'record_id':entry['record_id'],'field_sha256':c['archived_field_sha256'],
                   'cache_packet_sha256':c['point_cache_sha256'],'embedded_cache_sha256':packet['frozen_binary_cache']['sha256'],
                   'R_binary_sha256':rb,'center_binary_sha256':center,'model_contract_sha256':entry['model_contract_sha256'],
                   'parameter_cell_sha256':c['parameter_cell_sha256'],'dose_exact':str(dose),'gain_exact':str(gain),
                   'drive_exact':str(drive),'energy_scale_binary_exact':str(norm)}
    require(independently==target,'canonical exact target fingerprint drift')
    return R,independently


def check_producer_rows(q,rows,n):
    require(type(rows) is list and len(rows)==n and all(v>=0 for v in rows)
            and q==max(rows),'complete ordered nonnegative congruence row bounds required')


def numeric_review(m,reported,rrow,target,R,k_file,links):
    pid=target['pilot_id']; identity=target['identity']
    direct=direct_R_norm(R); integer=m.norm_R_integer(R)
    require(0<direct==integer==exact(reported['exact_R_infinity_norm'])==exact(target['reviewed_R_norm_exact'])
            and reported['producer_congruence'].get('R_binary_sha256',identity['R_binary_sha256'])==identity['R_binary_sha256'],
            'independent exact dyadic norm or frozen R identity differs')
    fields={'status':'UNRESOLVED_INDEPENDENT_REMAINING32_SPD_SUFFICIENT_BOUND','pilot_id':pid,
            'target_identity':identity,'trial_box_sha256':target['trial_box_sha256'],
            'native_cached_derivative_case_sha256':target['cached_derivative_case_sha256'],
            'native_geometry_case_sha256':target['geometry_replay_case_sha256'],
            'producer_result_sha256':rrow['sha256'],'exact_R_infinity_norm':str(direct),
            'R_file':target['R_source_file'],'R_file_sha256':target['R_source_file_sha256'],
            'archived_shift_exact':target['archived_shift_exact'],'archived_addition_gap_exact':target['archived_addition_gap_exact'],
            'precision_bits':PRECISION,'preconditioner_changed':False,'root_or_minimum_claimed':False,'scope':SCOPE}
    if rrow['K_file_sha256'] is None:
        require(reported['status']==PROD_UNRESOLVED and k_file is None and 'K_file' not in reported
                and reported['second_association_readback'] is None
                and reported['producer_congruence']['status']=='UNRESOLVED_FLOATING_PROPOSAL',
                'missing proposal must retain explicit unresolved outcome')
        return {**fields,'failure':'No authenticated proposal K.'}
    ksha=rrow['K_file_sha256']
    require(reported['K_file']==rel(k_file) and reported['K_file_sha256']==ksha==links[rel(k_file)],
            'unexpected exact candidate archive link')
    raw=k_file.read_bytes()
    require(sha256(raw).hexdigest()==ksha,'immutable proposal K file bytes differ')
    K=np.load(BytesIO(raw),allow_pickle=False)
    require(sha256(matrix_bytes(K,RANK)).hexdigest()==reported['K_binary_sha256'], 'immutable K entry bytes differ')
    p=reported['producer_congruence']; alternate=reported['second_association_readback']
    require(alternate['source_sha256']==PINS[MATH] and exact(alternate['exact_R_infinity_norm'])==direct,
            'producer second-association premise drift')
    q0,q1=(exact(x['q_infinity_upper_exact']) for x in (p,alternate))
    for obj,q in ((p,q0),(alternate,q1)):
        check_producer_rows(q,list(map(exact,obj['row_absolute_sum_upper_exact'])),RANK)
    old_pass=p['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD' and q0<1 and q1<1
    require(old_pass==(reported['status']==PROD_PASS),'producer status and strict sufficient bound differ')
    q,rows=m.congruence_bounds(R,K,PRECISION)
    check_producer_rows(q,rows,RANK)
    require(literal(m.ctx.prec,PRECISION) and sha256(matrix_bytes(R,RANK)).hexdigest()==identity['R_binary_sha256']
            and sha256(matrix_bytes(K,RANK)).hexdigest()==reported['K_binary_sha256'] and sha(k_file)==ksha,
            'R/K cannot drift while re-evaluating directed exact endpoints')
    ok=old_pass and q<1
    return {**fields,'status':CASE_PASS if ok else fields['status'],
            'K_file':rel(k_file),'K_file_sha256':ksha,'K_binary_sha256':reported['K_binary_sha256'],
            'independent_q_infinity_upper_exact':str(q),'independent_row_absolute_sum_upper_exact':list(map(str,rows)),
            'producer_q_upper_exact':str(q0),'producer_second_q_upper_exact':str(q1),
            'strict_congruence_bound':ok,'certificate_form':'K*((R+R^T)/2*K^T) evaluated in outward exact-endpoint Arb'}


def phase_binding(protocol_sha,preflight_sha,admission_sha):
    return {'source_sha256':sha(SOURCE),'protocol_sha256':protocol_sha,'preflight_sha256':preflight_sha,
            'producer_root_admission_sha256':admission_sha}


def preflight_record(ph,checks):
    return {'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,'checks':checks,
            'actual_producer_result_or_matrix_payloads_opened':0,'actual_matrix_jobs':0}


def reviewed_phase(p,pf,r,binding):
    require(r['status']==REVIEW_STATUS and yes(r['checks']) and all(r[k]==v for k,v in binding.items())
            and r['root_or_minimum_claimed'] is False and pf['status']=='PASS' and yes(pf['checks'])
            and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and literal(pf['actual_producer_result_or_matrix_payloads_opened'],0) and literal(pf['actual_matrix_jobs'],0)
            and p['source_sha256']==binding['source_sha256'], 'future native serial replay requires external source gate')


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def endpoint_module():
    read(MATH,PINS[MATH],parse=False)
    sys.path.insert(0,str(ROOT/'scripts'))
    import review_cohort_energy_stability as bounds
    require(Path(bounds.__file__).resolve()==MATH and sha(MATH)==PINS[MATH],'current reviewed endpoint arithmetic required')
    return bounds


def fixtures(plan,ph):
    binding=phase_binding(ph,'c'*64,'d'*64)
    pf=preflight_record(ph,{'check':True})
    gate={'status':REVIEW_STATUS,**binding,'checks':{'bound_io':True},'root_or_minimum_claimed':False}
    reviewed_phase(plan,pf,gate,binding)
    for bad in ({**gate,'preflight_sha256':'e'*64},{**gate,'checks':{'literal':1}},{**gate,'root_or_minimum_claimed':True}):
        require(reject(lambda:reviewed_phase(plan,pf,bad,binding)),'mixed source/phase/claim stage admitted')
    require(reject(lambda:reviewed_phase(plan,{**pf,'actual_matrix_jobs':False},gate,binding)),
            'Boolean masquerading as cold literal zero')
    group={'pid':9234,'bindings':plan['producer_bindings'],'ordered_case_ids':plan['ordered_case_ids'][:2]}
    prod={**plan,'ordered_case_ids':group['ordered_case_ids']}
    adm={'pid':9234,'status':ADMIT_STATUS,'observed_exit_code':0,'unified_exec_session':EXEC_SESSION,
         'manifest_sha256':'f'*64}
    start={**group,'status':'RUNNING'}
    term={**group,'status':'FINISHED','exit_code':0,'manifest_sha256':'f'*64}
    chosen=[{'case_id':k,'file':rel(PRODUCED/k/'result.json'),'sha256':'2'*64,'status':status,'K_file_sha256':kh}
            for k,status,kh in zip(prod['ordered_case_ids'],(PROD_PASS,PROD_UNRESOLVED),('3'*64,None),strict=True)]
    man={**group,'status':'COMPLETE_REMAINING32_SYMMETRIC_PART_ATTEMPTS_AWAITING_INDEPENDENT_REPLAY',
         'cases':chosen,'root_or_minimum_claimed':False,'externally_observed_exit_required':True,
         'external_independent_actual_review_required':True}
    verify_parent(adm,start,term,man,'f'*64,prod)
    require(len(admissible_paths(chosen))==6
            and reject(lambda:admissible_paths([{**chosen[0],'K_file_sha256':None},chosen[1]]))
            and reject(lambda:verify_parent({**adm,'observed_exit_code':False},start,term,man,'f'*64,prod))
            and reject(lambda:verify_parent(adm,{**start,'pid':9234.0},term,man,'f'*64,prod))
            and reject(lambda:verify_parent(adm,start,{**term,'exit_code':0.0},man,'f'*64,prod))
            and reject(lambda:verify_parent(adm,start,term,{**man,'cases':chosen[::-1]},'f'*64,prod)),
            'foreign native authority, incomplete K catalogue or ambiguous OS values accepted')
    from base64 import b64encode

    mat=np.array([[2.,5.],[-5.,2.]],dtype='<f8')
    center=np.zeros(2,dtype='<f8')
    data={'preconditioner':mat,'u':center,'free':np.arange(2),'reference':np.zeros((1,2),dtype='<f8'),'triangles':np.zeros((1,6),dtype=int)}
    buffer=BytesIO(); np.savez(buffer,**data); blob=buffer.getvalue()
    packed={'frozen_binary_cache':{'bytes':len(blob),'data':b64encode(blob).decode(),'sha256':sha256(blob).hexdigest(),
                                  'array_shapes':{k:list(v.shape) for k,v in data.items()}}}
    decoded=payload(packed)
    require(all(byte_equal(data[k],v) for k,v in decoded.items())
            and matched_tensor_except_u(decoded,{**decoded,'u':np.array([1.,0.])})
            and not matched_tensor_except_u(decoded,{**decoded,'preconditioner':mat.T})
            and reject(lambda:matrix_bytes(np.eye(2,dtype=bool),2))
            and reject(lambda:matrix_bytes(np.array([[np.nan,0],[0,1.]]),2)),
            'binary tensor serialization or proof center separation regression')
    m=endpoint_module()
    I=np.eye(2,dtype='<f8')
    g,gr=m.congruence_bounds(mat,I/2,PRECISION)
    s,sr=m.congruence_bounds(I,np.diag([1.,0.]),PRECISION)
    n,nr=m.congruence_bounds(np.diag([1.,-1.]),I,PRECISION)
    for q,rows in ((g,gr),(s,sr),(n,nr)):
        check_producer_rows(q,rows,2)
    require(g==Q(1,2) and s==1 and n==2 and direct_R_norm(decoded['preconditioner'])==m.norm_R_integer(mat)==7
            and reject(lambda:check_producer_rows(g,gr[:1],2)),
            'direct norm or strict congruence boundary fixtures differ')
    return {'actual_phase_roundtrip_positive_negative_without_self_referential_hash':True,
            'external_authority_order_literal_exit_and_missing_K_fail_closed':True,
            'synthetic_exact_binary_cache_and_center_only_tensor_transition':True,
            'non_symmetric_SPD_indefinite_and_exact_q_one_endpoint_fixtures':True,
            'exact_row_norm_and_nonfinite_or_Boolean_matrix_rejected':True}


def prepare():
    require(not REPLAY.exists(),'replay artifacts write once')
    p=static_plan()
    encoded=(json.dumps(p,indent=2,allow_nan=False)+'\n').encode('utf-8')
    ph=sha256(encoded).hexdigest()
    checks=fixtures(p,ph)
    with (HERE/'replay_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    require(save(HERE/'protocol.json',p)==ph,'preparation serialization drift')
    pf=save(HERE/'preflight.json',preflight_record(ph,checks))
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':pf,
                      'cold_target_count':len(p['ordered_case_ids']),'actual_results_or_arrays_opened':0}),flush=True)


def run(args):
    require(args.source_review.resolve()==HERE/'root_source_review.json','canonical external root source review required')
    p=read(HERE/'protocol.json',args.protocol_sha256)
    pf=read(HERE/'preflight.json',args.preflight_sha256)
    review=read(args.source_review,args.source_review_sha256)
    binding=phase_binding(args.protocol_sha256,args.preflight_sha256,args.admission_sha256)
    reviewed_phase(p,pf,review,binding)
    require(p==static_plan() and sha(HERE/'replay_frozen.py')==binding['source_sha256'], 'changed frozen replay phase')
    adm,observed=admitted(args.admission,args.admission_sha256,p)
    m=endpoint_module()
    source=read(OUT/'preparation/protocol.json',PINS[OUT/'preparation/protocol.json'])
    merged={**p['input_sha256'],**adm['sha256'],rel(ADMISSION):args.admission_sha256,
            rel(args.source_review):args.source_review_sha256}
    for name,h in source['source_and_input_sha256'].items():
        require(merged.setdefault(name,h)==h,'conflicting independent source closure')
        if name.endswith('.py'):
            read(canon(name),h,parse=False)
    prepared=read(TARGET/'preparation/manifest.json',PINS[TARGET/'preparation/manifest.json'])['targets']
    joined=read(JOIN/'bridge_input_admission.json',PINS[JOIN/'bridge_input_admission.json'])['cases']
    REPLAY.mkdir(exist_ok=False)
    started=perf_counter()
    phase={'pid':os.getpid(),**binding,'source_review_sha256':args.source_review_sha256,
           'producer_pid':adm['pid'],'producer_exec_session':adm['unified_exec_session'],
           'producer_manifest_sha256':adm['manifest_sha256'],'producer_terminal_sha256':adm['terminal_sha256']}
    save(REPLAY/'started.json',{'status':'RUNNING',**phase})
    terminal={'status':'FAILED','exit_code':1,**phase}
    rows=[]
    try:
        for e,j,target,rrow in zip(prepared,joined,source['ordered_targets'],observed,strict=True):
            case=target['pilot_id']; path=canon(rrow['file'])
            reported=read(path,rrow['sha256'])
            require(reported['case_id']==rrow['case_id']==case and literal(reported['pid'],adm['pid'])
                    and reported['bindings']==p['producer_bindings'] and reported['status']==rrow['status']
                    and literal(reported['precision_bits'],PRODUCER_PRECISION) and reported['target']==target
                    and reported['external_independent_actual_review_required'] is True
                    and reported['preconditioner_changed'] is False and reported['root_or_minimum_claimed'] is False
                    and literal(reported['new_force_Hessian_AD_integration_or_solve_jobs'],0)
                    and rrow['R_binary_sha256']==target['identity']['R_binary_sha256'],
                    'same-named frozen producer candidate and no-field-change scope required')
            R,ident=identity_and_R(e,j,target,merged)
            require(ident['R_binary_sha256']==rrow['R_binary_sha256'],'actual raw exact binary R differs from manifest')
            K=path.with_name('congruence_K.npy') if rrow['K_file_sha256'] is not None else None
            rr=numeric_review(m,reported,rrow,target,R,K,merged)
            receipt=REPLAY/f'{case}.json'
            dh=save(receipt,{**rr,'case_id':case,'replay_phase':phase})
            rows.append({'case_id':case,'status':rr['status'],'receipt':rel(receipt),'sha256':dh,
                         'target_identity':ident,'trial_box_sha256':target['trial_box_sha256'],
                         'exact_R_infinity_norm':rr['exact_R_infinity_norm'],
                         'independent_q_infinity_upper_exact':rr.get('independent_q_infinity_upper_exact')})
            print(json.dumps({'case_id':case,'status':rr['status'],'native_case_sha256':dh}),flush=True)
        require(len(rows)==32 and [r['case_id'] for r in rows]==p['ordered_case_ids'],'all32 actual native cases required')
        for name,h in merged.items():
            read(canon(name),h,parse=False)
        count=sum(r['status']==CASE_PASS for r in rows)
        result={'status':RESULT_PASS if count==32 else 'REMAINING32_SPD_REPLAY_RETAINS_UNRESOLVED_CASES',
                'source_and_phase':phase,'scope':SCOPE,'case_results':rows,'producer_and_input_sha256':merged,
                'independently_verified_symmetric_parts_SPD':count,'remaining_unresolved_sufficient_bounds':32-count,
                'root_or_minimum_claimed':False,'preconditioner_or_archived_field_changed':False,
                'new_force_Hessian_AD_integration_or_solve_jobs':0,'seconds':perf_counter()-started}
        rh=save(REPLAY/'receipt.json',result)
        terminal.update(status='FINISHED',exit_code=0,receipt_sha256=rh)
        print(json.dumps({'status':result['status'],'receipt_sha256':rh,'certified':count,'unresolved':32-count}),flush=True)
    except BaseException as error:
        save(REPLAY/'failure.json',{'error':str(error),'traceback':traceback.format_exc(),'completed':rows,**phase})
        raise
    finally:
        save(REPLAY/'terminal.json',terminal)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest='mode',required=True)
    modes.add_parser('prepare')
    actual=modes.add_parser('run')
    for k in ('protocol-sha256','preflight-sha256','admission-sha256','source-review-sha256'):
        actual.add_argument('--'+k,required=True)
    for k in ('admission','source-review'):
        actual.add_argument('--'+k,type=lambda v:Path(v).resolve(),required=True)
    args=parser.parse_args()
    if args.mode=='prepare':
        prepare()
    else:
        run(args)


if __name__=='__main__':
    main()
