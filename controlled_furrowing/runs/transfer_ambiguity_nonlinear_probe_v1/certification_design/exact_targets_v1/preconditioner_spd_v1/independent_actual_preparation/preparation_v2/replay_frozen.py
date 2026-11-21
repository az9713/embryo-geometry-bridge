"""Independently replay the four completed exact symmetric-part congruences."""

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
import traceback


HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
SPD = HERE.parent
EXACT = SPD.parent
ROOT = EXACT.parents[3]
PREP = SPD/'preparation'
PRODUCED = SPD/'actual'
OUT = HERE/'actual_replay'
REPLAY_PREPARATION = HERE/'preparation_v2'
CALLER = ROOT/'scripts/transfer_ambiguity_preconditioner_spd.py'
CALLER_PROTOCOL = PREP/'protocol.json'
CALLER_PREFLIGHT = PREP/'preflight.json'
CALLER_REVIEW = PREP/'independent_prelaunch/receipt.json'
MATH = ROOT/'scripts/review_cohort_energy_stability.py'
DERIV = EXACT/'cached_derivative_extension_v1/independent_actual_preparation/bridge_admission.json'
RESULT_SHA = 'f70baf99b77e7748c2996140c67acae2492b3ba719b6d8b7eda6d2f3ad0dc319'
TERMINAL_SHA = 'bc716e69dfe9e2431c4096e3c87517c4da490990e005b60e86d2be75e5d3696f'
PIDS = tuple(f'{endpoint}_{suffix}' for endpoint in ('max_contrast_change','max_negative_contrast_change')
             for suffix in ('small_17','large_03'))
SUPERVISOR_PID = 6012
EXEC_SESSION = 93286
N = 611
PRECISION = 640
ADMISSION_STATUS = 'PASS_ROOT_OBSERVED_SINGLETON_PRECONDITIONER_SPD_TERMINALS'
SOURCE_REVIEW_STATUS = 'PASS_ROOT_SPD_INDEPENDENT_ACTUAL_REPLAY_SOURCE'
SCOPE = ('Four independently replayed exact symmetric parts of the same frozen binary preconditioners. '
         'A separate pure congruence arithmetic association is evaluated afresh with decoded R and K. '
         'No spatial force, Hessian, local minimum, PDE, source or new equilibrium result.')
PINS = {
    CALLER:'67ef2be40fc423d721390fd360f6244ecfc529923aa5c41a04981edec4cd20d7',
    PREP/'transfer_ambiguity_preconditioner_spd_frozen.py':
        '67ef2be40fc423d721390fd360f6244ecfc529923aa5c41a04981edec4cd20d7',
    CALLER_PROTOCOL:'1b6369802375d358abf14d8db75ae1b0e3d418be8705c2040395ce2cc434c905',
    CALLER_PREFLIGHT:'489cce78693ac12b6a0851c48b39895a0a65022611de356367b0b5275a2846cc',
    CALLER_REVIEW:'e295d53e69cf1cc924dae538a1988f1487040c94b2199e893ec5c6c115887279',
    DERIV:'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
}


def rel(path):
    return path.relative_to(ROOT).as_posix()


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def digest_text(s):
    return type(s) is str and re.fullmatch('[0-9a-f]{64}',s) is not None


def pinned(path, h, *, parse=True):
    require(digest_text(h),'explicit SHA-256 required')
    path=path.resolve()
    require(path.is_relative_to(ROOT),f'outside study: {path}')
    raw=path.read_bytes()
    require(sha256(raw).hexdigest()==h,f'changed immutable bytes: {path}')
    return json.loads(raw) if parse else raw


def link(name):
    path=(ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and rel(path)==name,'canonical study-relative authority path required')
    return path


def save(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def positive_int(x):
    return type(x) is int and x>0


def strict_int_zero(x):
    return type(x) is int and x==0


def yes(x):
    return type(x) is dict and bool(x) and all(v is True for v in x.values())


def exact(value):
    require(type(value) in (str,int,Q),'finite exact rational without float or Boolean required')
    return Q(value)


def static_protocol():
    for path,h in PINS.items():
        pinned(path,h,parse=False)
    proc=pinned(CALLER_PROTOCOL,PINS[CALLER_PROTOCOL])
    frozen=pinned(CALLER_PREFLIGHT,PINS[CALLER_PREFLIGHT])
    cold=pinned(CALLER_REVIEW,PINS[CALLER_REVIEW])
    expected={'source_sha256':PINS[CALLER],'protocol_sha256':PINS[CALLER_PROTOCOL],
              'preflight_sha256':PINS[CALLER_PREFLIGHT],
              **{k:proc[k] for k in ('minimum_design_sha256','cached_root_review_sha256','target_manifest_sha256')}}
    require(cold['status']=='PASS_INDEPENDENT_SINGLETON_PRECONDITIONER_SPD_PRELAUNCH'
            and yes(cold['checks']) and all(cold[k]==v for k,v in expected.items())
            and proc['status']=='FROZEN_FOUR_SINGLETON_SYMMETRIC_PRECONDITIONER_PART_PROOFS'
            and proc['precision_bits']==PRECISION and proc['rank']==N
            and tuple(r['pilot_id'] for r in proc['ordered_targets'])==PIDS
            and proc['larger_derivative_admission_sha256']==PINS[DERIV]
            and frozen['status']=='PASS' and yes(frozen['checks'])
            and frozen['source_sha256']==PINS[CALLER] and frozen['protocol_sha256']==PINS[CALLER_PROTOCOL]
            and strict_int_zero(frozen['actual_production_matrix_jobs']),
            'frozen producer source, cold approval and all four cases required')
    inherited=proc['source_and_input_sha256']
    require(len(inherited)==227 and rel(MATH) in inherited,'missing reviewed module and source closure')
    for name,h in inherited.items():
        path=link(name)
        if name.endswith('.py'):
            pinned(path,h,parse=False)
    return {'status':'FROZEN_INDEPENDENT_ACTUAL_SPD_REPLAY','scope':SCOPE,'source_sha256':sha(SOURCE),
            'producer_bindings':{**expected,'prelaunch_review_sha256':PINS[CALLER_REVIEW]},
            'source_and_static_sha256':{rel(k):v for k,v in PINS.items()},
            'source_and_case_input_sha256':inherited,
            'observed_supervisor_pid':SUPERVISOR_PID,'observed_exec_session':EXEC_SESSION,
            'expected_result_sha256':RESULT_SHA,'expected_terminal_sha256':TERMINAL_SHA,
            'required_metadata_status':ADMISSION_STATUS,'required_source_review_status':SOURCE_REVIEW_STATUS,
            'pilot_order':list(PIDS),'rank':N,'precision_bits':PRECISION,
            'method':'Decode physical binary cache and supplied K separately. Reopen same approved singleton. '
                     'Evaluate the reviewed K*(P*K^T) association again on these exact inputs at 640 bits. '
                     'Independently recompute the raw binary64 R row-sum norm.',
            'production_metadata_body_or_numeric_arrays_read_during_preparation':0,
            'actual_production_matrix_jobs_during_preparation':0}


def admitted_parent(a,start,final,results,bindings):
    require(a['status']==ADMISSION_STATUS and strict_int_zero(a['observed_exit_code'])
            and positive_int(a['supervisor_pid']) and a['supervisor_pid']==SUPERVISOR_PID
            and positive_int(a['unified_exec_session']) and a['unified_exec_session']==EXEC_SESSION
            and all(positive_int(d['pid']) and d['pid']==SUPERVISOR_PID for d in (start,final))
            and start['status']=='RUNNING' and final['status']=='FINISHED'
            and strict_int_zero(final['exit_code']) and final['results_sha256']==RESULT_SHA
            and start['bindings']==final['bindings']==results['bindings']==bindings
            and tuple(start['ordered_pilots'])==PIDS
            and tuple(r['pilot_id'] for r in final['workers'])==tuple(r['pilot_id'] for r in results['cases'])==PIDS
            and results['status']=='COMPLETE_FOUR_PRECONDITIONER_SPD_ATTEMPTS_AWAITING_REVIEW'
            and results['no_F_H_AD_solve_or_minimum'] is True,
            'recorded actual successful supervisor and exact case order required')


def admit_child(pilot,meta,root_row,worker,hashes,bindings):
    start,launched,observed,terminal,result=meta
    pid=launched['pid']
    require(positive_int(pid) and positive_int(launched['supervisor_pid']) and pid!=SUPERVISOR_PID
            and all(positive_int(d['pid']) and d['pid']==pid
                    and positive_int(d['supervisor_pid']) and d['supervisor_pid']==SUPERVISOR_PID
                    and d['pilot_id']==pilot for d in meta)
            and start['bindings']==terminal['bindings']==result['bindings']==bindings
            and start['status']=='RUNNING' and terminal['status']=='FINISHED'
            and strict_int_zero(observed['observed_exit_code']) and strict_int_zero(terminal['exit_code'])
            and terminal['result_sha256']==root_row['result_sha256']==hashes['result.json']
            and observed['terminal_sha256']==root_row['terminal_sha256']==hashes['terminal.json']
            and root_row['observed_os_exit_sha256']==hashes['observed_os_exit.json']
            and positive_int(worker['pid']) and positive_int(worker['supervisor_pid'])
            and strict_int_zero(worker['observed_exit_code'])
            and worker=={**observed,'os_exit_sha256':hashes['observed_os_exit.json']},
            'matched literal-int worker/parent PIDs, terminal, and true OS exits required')
    claimed='CANDIDATE_EXACT_SYMMETRIC_PART_SPD_AWAITING_REVIEW'
    require(result['status']==root_row['status']==claimed and result['external_independent_actual_review_required'] is True
            and result['preconditioner_changed'] is False and result['root_or_minimum_claimed'] is False
            and strict_int_zero(result['actual_spatial_H_or_force_jobs'])
            and strict_int_zero(result['actual_AD_or_equilibrium_solves'])
            and type(result['precision_bits']) is int and result['precision_bits']==PRECISION,
            'completed producer candidate with no promoted root, minimum, AD or spatial claim required')
    return pid


def admitted(path,digest,protocol):
    a=pinned(path,digest)
    require(a['status']==ADMISSION_STATUS and strict_int_zero(a['observed_exit_code'])
            and positive_int(a['supervisor_pid']) and a['supervisor_pid']==SUPERVISOR_PID
            and positive_int(a['unified_exec_session']) and a['unified_exec_session']==EXEC_SESSION,
            'external observed true process exit required before reading actual results')
    parent_names=('started.json','terminal.json','results.json')
    child_json=('started.json','supervisor_child_started.json','observed_os_exit.json','terminal.json','result.json')
    kp='congruence_K.npy'
    expected={rel(PRODUCED/k) for k in parent_names}
    expected.update(rel(PRODUCED/p/k) for p in PIDS for k in (*child_json,kp))
    hashes=a['sha256']
    require(set(hashes)==expected and all(digest_text(h) for h in hashes.values())
            and hashes[rel(PRODUCED/'results.json')]==RESULT_SHA and hashes[rel(PRODUCED/'terminal.json')]==TERMINAL_SHA,
            'exact complete terminal and K manifest required')
    start,final,results=[pinned(PRODUCED/s,hashes[rel(PRODUCED/s)]) for s in parent_names]
    bindings=protocol['producer_bindings']
    admitted_parent(a,start,final,results,bindings)
    candidates=[]
    seen=[]
    for pilot,row,wd in zip(PIDS,results['cases'],final['workers'],strict=True):
        folder=PRODUCED/pilot
        fh={k:hashes[rel(folder/k)] for k in (*child_json,kp)}
        meta=[pinned(folder/k,fh[k]) for k in child_json]
        pid=admit_child(pilot,meta,row,wd,fh,bindings)
        cmd=meta[1]['command']
        require(type(cmd) is list and len(cmd)==12 and Path(cmd[2]).resolve()==CALLER
                and cmd[1:2]==['-B'] and cmd[3:9]==['worker','--pilot-id',pilot,'--supervisor-pid',str(SUPERVISOR_PID),
                                                    '--review'] and Path(cmd[9]).resolve()==CALLER_REVIEW
                and cmd[10:]==['--review-sha256',PINS[CALLER_REVIEW]],
                'supervisor command did not bind exact independently approved worker')
        result=meta[-1]
        require(result['K_file']==rel(folder/kp) and result['K_file_sha256']==fh[kp],
                'wrong candidate K artifact and byte authority')
        seen.append(pid)
        candidates.append((pilot,result,folder/kp,fh[kp]))
    require(len(set(seen))==len(seen)==4,'distinct observed child PIDs required')
    # All process metadata was checked before opening any binary K bytes.
    require(all(sha(k)==h for _,_,k,h in candidates),'an actual K candidate no longer matches its admitted bytes')
    return a,candidates,seen


def matrix_bytes(R,n):
    import numpy as np
    require(type(R) is np.ndarray and R.shape==(n,n) and R.dtype.str=='<f8'
            and np.isfinite(R).all(), 'finite binary64 matrix, exact rank and byte order required')
    return R.tobytes(order='C')


def direct_R_norm(R):
    return max(sum(abs(Q.from_float(float(x))) for x in row) for row in R)


def decode_R(packet):
    import numpy as np
    data=packet['frozen_binary_cache']
    raw=b64decode(data['data'],validate=True)
    require(sha256(raw).hexdigest()==data['sha256'],'actual cache payload differs')
    with np.load(BytesIO(raw),allow_pickle=False) as arrays:
        require(all(k in arrays.files for k in ('preconditioner','u','free')), 'cache arrays incomplete')
        return arrays['preconditioner'].copy(),arrays['u'].copy(),arrays['free'].copy()


def identity_and_matrix(pilot,reported,target,cases,manifest,data,ex,dx):
    case,=[c for c in cases if c['pilot_id']==pilot]
    target_manifest_entry,=[r for r in manifest['targets'] if r['pilot_id']==pilot]
    ex.bootstrap()
    ctx,_,entry=ex.open_context(pilot,manifest,data)
    require(entry==target_manifest_entry and ctx.pilot_id==target['pilot_id']==case['pilot_id']==pilot
            and ex.T.core.ctx.prec==PRECISION, 'fresh exact same-target context required')
    contract=pinned(link(entry['model_contract']),entry['model_contract_sha256'])
    packet=pinned(link(case['cache_packet']),case['cache_packet_sha256'])
    proof=pinned(link(case['verified_point']),case['verified_point_sha256'])
    _,_,reference_norm=dx.validate_center(proof,case)
    R,u,free=decode_R(packet)
    rb=matrix_bytes(R,N)
    rh=sha256(rb).hexdigest()
    require(ctx.R.flags.writeable is False and rb==matrix_bytes(ctx.R,N)
            and ctx.u.shape==u.shape==ctx.free.shape==free.shape==(N,)
            and ctx.u.dtype==u.dtype and ctx.free.dtype==free.dtype
            and ctx.u.tobytes()==u.tobytes() and ctx.free.tobytes()==free.tobytes(),
            'separate decoded binary R and saved coordinate order differ from reviewed factory')
    gain=exact(packet['common_gain_exact'])
    dose=exact(packet['optical_dose_exact'])
    strength=exact(packet['activation_strength_exact'])
    scale=exact(packet['normalization_binary_exact'])
    cell=ctx.row['parameter_cell']
    cell_hash=sha256(json.dumps(cell,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    require(0<strength==gain*dose==ctx.drive==exact(case['actual_exact_drive'])==exact(proof['exact_strength'])
            and strength!=exact(packet['activation_strength_binary_exact'])
            and (gain,dose)==(ctx.base['common_gain'],ctx.dose)
            and all(lo==hi==ctx.base[k] for k,(lo,hi) in ctx.cell.items())
            and 0<scale==ctx.target.energy_scale==exact(contract['energy_scale_binary_exact'])
            and cell==contract['parameter_cell']
            and cell_hash==ctx.row['parameter_cell_sha256']==contract['parameter_cell_sha256']
            and case['cache_packet_sha256']==contract['point_cache_sha256']==ctx.row['point_cache_sha256']
            and case['embedded_cache_sha256']==contract['frozen_binary_cache_sha256']==packet['frozen_binary_cache']['sha256']
            and case['field_sha256']==contract['source_field_sha256']==packet['input_state_sha256']
            and rh==case['preconditioner_binary_sha256']==contract['binary_preconditioner_sha256'],
            'different exact drive, cell, cache, normalization or matrix identity')
    check={'pilot_id':pilot,'record_id':ctx.record_id,'field_sha256':case['field_sha256'],
           'cache_packet_sha256':case['cache_packet_sha256'],'embedded_cache_sha256':case['embedded_cache_sha256'],
           'R_binary_sha256':rh,'model_contract_sha256':entry['model_contract_sha256'],'parameter_cell_sha256':cell_hash,
           'dose_exact':str(dose),'gain_exact':str(gain),'drive_exact':str(strength),'energy_scale_binary_exact':str(scale)}
    require(reported['target']==target and target['identity']==check
            and target['point_proof_sha256']==case['verified_point_sha256']
            and exact(target['cached_reference_R_norm_upper_exact'])==reference_norm
            and target['larger_radius_box_sha256']==entry['box_1_sha256']==case['box_sha256'],
            'producer target does not match independent runtime cell and archive')
    return R,check,reference_norm


def accepted_rows(q,rows,n):
    require(type(rows) is list and len(rows)==n and all(x>=0 for x in rows)
            and 0<=q==max(rows)<1,'strict q < 1 and full nonnegative outward rows required')


def numeric_review(pilot,result,Kfile,kh,target,phase):
    import numpy as np
    cases,manifest,data,ex,dx,math=phase
    R,identity,reference_norm=identity_and_matrix(pilot,result,target,cases,manifest,data,ex,dx)
    raw=pinned(Kfile,kh,parse=False)
    K=np.load(BytesIO(raw),allow_pickle=False)
    require(sha256(matrix_bytes(K,N)).hexdigest()==result['K_binary_sha256'],'binary K and archive SHA differ')
    ex.bootstrap()
    q,rows=math.congruence_bounds(R,K,PRECISION)
    require(math.ctx.prec==PRECISION and ex.T.core.ctx is math.ctx,'640-bit fresh congruence arithmetic required')
    accepted_rows(q,rows,N)
    norm=direct_R_norm(R)
    require(norm==math.norm_R_integer(R)==reference_norm and norm>0,
            'direct exact binary row norm differs from reviewed prior or independent integer computation')
    p=result['producer_congruence']
    peer=result['second_association_readback']
    qp=exact(p['q_infinity_upper_exact'])
    q2=exact(peer['q_infinity_upper_exact'])
    accepted_rows(qp,list(map(exact,p['row_absolute_sum_upper_exact'])),N)
    accepted_rows(q2,list(map(exact,peer['row_absolute_sum_upper_exact'])),N)
    require(p['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD'
            and type(p['dimension']) is int and p['dimension']==N
            and type(p['precision_bits']) is int and p['precision_bits']==PRECISION
            and p['R_binary_sha256']==identity['R_binary_sha256']
            and p['floating_Cholesky_and_K_are_untrusted_proposals_only'] is True
            and p['R_symmetry_assumed'] is False
            and q==q2 and rows==list(map(exact,peer['row_absolute_sum_upper_exact']))
            and norm==exact(peer['exact_R_infinity_norm']) and peer['source_sha256']==sha(MATH)
            and sha(Kfile)==kh and sha256(matrix_bytes(R,N)).hexdigest()==identity['R_binary_sha256'],
            'independent fresh congruence and exact norm must match same-input outward evidence')
    return {'status':'PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD','scope':SCOPE,
            'pilot_id':pilot,'target_identity':identity,'record_id':identity['record_id'],'precision_bits':PRECISION,
            'cache_packet_sha256':identity['cache_packet_sha256'],
            'producer_result_sha256':sha(PRODUCED/pilot/'result.json'),'K_file':rel(Kfile),'K_file_sha256':kh,
            'K_binary_sha256':result['K_binary_sha256'],'cached_point_proof_sha256':target['point_proof_sha256'],
            'larger_radius_box_sha256':target['larger_radius_box_sha256'],
            'derivative_case_sha256':target['native_cached_derivative_case_sha256'],
            'alternate_association':'K*(P*K^T)','fresh_congruence_q_upper_exact':str(q),
            'producer_congruence_q_upper_exact':str(qp),'max_reviewed_q_upper_exact':str(max(q,qp)),
            'fresh_absolute_row_sum_upper_exact':list(map(str,rows)),
            'exact_R_infinity_norm':str(norm),'strict_congruence_margin_lower_exact':str(1-max(q,qp)),
            'exact_target_symmetric_part_spd':True,'R_symmetry_assumed':False,'root_or_minimum_claimed':False,
            'spatial_force_Hessian_AD_or_solve_jobs':0}


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def preflight_record(source_sha,protocol_sha,checks):
    # A receipt cannot contain its own digest. The launch verifies that externally.
    return {'status':'PASS','source_sha256':source_sha,'protocol_sha256':protocol_sha,'checks':checks,
            'actual_producer_metadata_or_arrays_opened':0,'fresh_actual_production_congruences':0,
            'new_F_H_AD_or_solve_jobs':0}


def reviewed_phase(protocol,preflight,review,binding):
    require(review['status']==SOURCE_REVIEW_STATUS and yes(review['checks'])
            and all(review[k]==v for k,v in binding.items())
            and protocol['source_sha256']==binding['source_sha256']
            and preflight['status']=='PASS' and yes(preflight['checks'])
            and all(preflight[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and strict_int_zero(preflight['actual_producer_metadata_or_arrays_opened'])
            and strict_int_zero(preflight['fresh_actual_production_congruences']),
            'reviewed unchanged actual replay source and phase required')


def fixtures(protocol):
    import numpy as np
    x=np.array([[2.,1.],[2.**-53,2.]],dtype='<f8')
    first=PIDS[0]
    groups={'source':'fixture'}
    pid=7
    token={'pilot_id':first,'pid':pid,'supervisor_pid':SUPERVISOR_PID}
    records=[
        {**token,'status':'RUNNING','bindings':groups},
        {**token,'command':['fixture']},
        {**token,'observed_exit_code':0,'terminal_sha256':'a'*64},
        {**token,'exit_code':0,'status':'FINISHED','result_sha256':'b'*64,'bindings':groups},
        {**token,'status':'CANDIDATE_EXACT_SYMMETRIC_PART_SPD_AWAITING_REVIEW','bindings':groups,
         'external_independent_actual_review_required':True,'preconditioner_changed':False,
         'root_or_minimum_claimed':False,'actual_spatial_H_or_force_jobs':0,'actual_AD_or_equilibrium_solves':0,
         'precision_bits':PRECISION},
    ]
    hashes={'terminal.json':'a'*64,'result.json':'b'*64,'observed_os_exit.json':'c'*64}
    row={'result_sha256':'b'*64,'terminal_sha256':'a'*64,
         'observed_os_exit_sha256':'c'*64,'status':records[-1]['status']}
    worker={**records[2],'os_exit_sha256':'c'*64}
    require(admit_child(first,records,row,worker,hashes,groups)==pid,'synthetic worker identity')
    def edited(i,key,value):
        copied=[dict(x) for x in records]
        copied[i][key]=value
        return rejects(lambda:admit_child(first,copied,row,worker,hashes,groups))
    accepted_rows(Q(1,2),[Q(1,2),Q(1,3)],2)
    binding={'source_sha256':protocol['source_sha256'],'protocol_sha256':'a'*64,
             'preflight_sha256':'b'*64,'producer_root_admission_sha256':'c'*64}
    pf=preflight_record(binding['source_sha256'],binding['protocol_sha256'],{'synthetic':True})
    review={'status':SOURCE_REVIEW_STATUS,'checks':{'synthetic':True},**binding}
    reviewed_phase(protocol,pf,review,binding)
    return {'exact_binary_identity_and_direct_R_norm_oracle':len(matrix_bytes(x,2))==32 and direct_R_norm(x)==3
                 and rejects(lambda:matrix_bytes(np.zeros((2,2),dtype=int),2)),
            'all_process_PIDs_require_literal_int_including_terminal_and_started':all(
                edited(i,key,float(records[i][key])) for i in range(len(records)) for key in ('pid','supervisor_pid'))
                and rejects(lambda:admit_child(first,records,row,{**worker,'pid':float(pid)},hashes,groups)),
            'reject_boolean_worker_exits_or_stale_result_hash':(
                edited(2,'observed_exit_code',False) and edited(3,'exit_code',False)
                and edited(3,'result_sha256','d'*64)),
            'strict_q_equality_and_negative_row_rejected':(
                rejects(lambda:accepted_rows(Q(1),[Q(1),Q(0)],2))
                and rejects(lambda:accepted_rows(Q(1,2),[Q(1,2),Q(-1,2)],2))),
            'gate_uses_real_generated_preflight_schema_without_impossible_self_digest':(
                'preflight_sha256' not in pf and rejects(lambda:reviewed_phase(protocol,pf,review,
                    {**binding,'preflight_sha256':'d'*64}))
                and rejects(lambda:reviewed_phase(protocol,{**pf,'protocol_sha256':'d'*64},review,binding))),
            'no_independent_actual_outputs':not OUT.exists()}


def prepare():
    REPLAY_PREPARATION.mkdir(exist_ok=False)
    with (REPLAY_PREPARATION/'replay_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    protocol=static_protocol()
    checks=fixtures(protocol)
    require(yes(checks),'cold independent SPD replay fixtures failed')
    ph=save(REPLAY_PREPARATION/'protocol.json',protocol)
    fh=save(REPLAY_PREPARATION/'preflight.json',preflight_record(sha(SOURCE),ph,checks))
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':fh}),flush=True)


def run(args):
    protocol=pinned(REPLAY_PREPARATION/'protocol.json',args.protocol_sha256)
    pf=pinned(REPLAY_PREPARATION/'preflight.json',args.preflight_sha256)
    r=pinned(args.source_review,args.source_review_sha256)
    binding={'source_sha256':sha(SOURCE),'protocol_sha256':args.protocol_sha256,'preflight_sha256':args.preflight_sha256,
             'producer_root_admission_sha256':args.admission_sha256}
    reviewed_phase(protocol,pf,r,binding)
    require(sha(REPLAY_PREPARATION/'replay_frozen.py')==binding['source_sha256']
            and static_protocol()==protocol,'frozen replay source, authority or phase changed')
    authority,candidates,children=admitted(args.admission,args.admission_sha256,protocol)
    inputs=protocol['source_and_case_input_sha256']
    for name,digest in inputs.items():
        require(sha(link(name))==digest,'same actual case input closure required')
    sys.path.insert(0,str(ROOT/'scripts'))
    import transfer_ambiguity_exact_targets as ex
    import transfer_ambiguity_cached_derivative_extension as dx
    import review_cohort_energy_stability as math
    for mod in (ex,dx,math):
        path=Path(mod.__file__).resolve()
        require(sha(path)==inputs[rel(path)],'runtime source alias drift')
    require(Path(math.__file__).resolve()==MATH and dx.EX is ex,'same factory and independent arithmetic module required')
    ex.bootstrap()
    _,cases,manifest,data=dx.authority()
    wrapper=pinned(DERIV,PINS[DERIV])
    derivative={r['pilot_id']:r for r in wrapper['reviewed_cases']}
    producer=pinned(CALLER_PROTOCOL,PINS[CALLER_PROTOCOL])
    require(tuple(derivative)==tuple(row['pilot_id'] for row in producer['ordered_targets'])==PIDS
            and all(row['native_cached_derivative_case_sha256']==derivative[row['pilot_id']]['native_case_sha256']
                    and row['larger_radius_box_sha256']==derivative[row['pilot_id']]['box_sha256']
                    for row in producer['ordered_targets']),'same larger-domain and case association required')
    OUT.mkdir(exist_ok=False)
    mark={'pid':os.getpid(),**binding,'source_review_sha256':args.source_review_sha256}
    save(OUT/'started.json',{'status':'RUNNING',**mark})
    terminal={'status':'FAILED','exit_code':1,**mark}
    try:
        rows=[]
        phase=(cases,manifest,data,ex,dx,math)
        for (pilot,reported,k,h),target in zip(candidates,producer['ordered_targets'],strict=True):
            case=numeric_review(pilot,reported,k,h,target,phase)
            path=OUT/(pilot+'.json')
            digest=save(path,{**case,'replay_phase':mark})
            rows.append({q:case[q] for q in ('pilot_id','status','target_identity','cache_packet_sha256',
                                             'producer_result_sha256','K_file_sha256','K_binary_sha256','fresh_congruence_q_upper_exact',
                                             'max_reviewed_q_upper_exact','exact_R_infinity_norm')}
                        |{'case_receipt_file':rel(path),'case_receipt_sha256':digest})
            print(json.dumps({'pilot':pilot,'q_display':float(exact(case['fresh_congruence_q_upper_exact'])),
                              'case_receipt_sha256':digest}),flush=True)
        require(tuple(row['pilot_id'] for row in rows)==PIDS,'all four matrix proofs required')
        receipt={'status':'PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD',**mark,
                 'scope':SCOPE,'checks':{'literal_supervisor_and_distinct_worker_OS_exits_verified':True,
                     'exact_source_and_cache_R_K_energy_drive_cell_associations':True,
                     'all_four_fresh_Arb640_congruence_products_strict':True,
                     'independently_recomputed_exact_binary_R_norms':True},
                 'producer_observed_exec_session':authority['unified_exec_session'],
                 'producer_observed_pid':authority['supervisor_pid'],'producer_worker_pids':children,
                 'producer_metadata_and_K_file_sha256':authority['sha256'],
                 'producer_results_sha256':RESULT_SHA,'producer_terminal_sha256':TERMINAL_SHA,
                 'ordered_cases':rows,'fresh_actual_matrix_congruences':4,'spatial_force_Hessian_AD_or_solve_jobs':0,
                 'root_or_minimum_certified':False}
        h=save(OUT/'receipt.json',receipt)
        terminal.update(status='FINISHED',exit_code=0,receipt_sha256=h)
        print(json.dumps({'status':receipt['status'],'receipt_sha256':h}),flush=True)
    except BaseException as e:
        save(OUT/'failure.json',{'error':str(e),'traceback':traceback.format_exc(),**mark})
        raise
    finally:
        save(OUT/'terminal.json',terminal)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='stage',required=True)
    sub.add_parser('prepare')
    actual=sub.add_parser('run')
    for key in ('protocol','preflight','admission','source-review'):
        actual.add_argument('--'+key+'-sha256',required=True)
        if key in ('admission','source-review'):
            actual.add_argument('--'+key,type=lambda s:Path(s).resolve(),required=True)
    args=parser.parse_args()
    if args.stage=='prepare':
        prepare()
    else:
        run(args)
