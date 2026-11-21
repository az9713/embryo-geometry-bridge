"""Cold approval audit for the two exact binary cached recenter trials."""

import base64
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np


HERE=Path(__file__).resolve().parent
PREP=HERE.parent
RECENTER=PREP.parent
ROOT=RECENTER.parents[5]
SOURCE=ROOT/'scripts/transfer_ambiguity_recenter_cached.py'
MAIN='9b149c5889c8dc6a77b9407674200d9e6bcbbee8e9875f12a13b97f7962de7d4'
PROTOCOL='bb5f0aef19bf5a7d5a6f954d7b41f03219bd19056be4e4eb4672454a588bec63'
PREFLIGHT='9d38373a9f9895880ff1d03d6ba5c1165a790cbff8123787cee9a0ddf6fe0a35'
DESIGN='d9621eef3cc7fec129b2f597dd4ab584393a6ed0b610051498734fdd2cad1d04'
ORDER=('max_negative_contrast_change_small_16','max_contrast_change_small_13')


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def pin(p,h):
    raw=p.read_bytes()
    require(sha256(raw).hexdigest()==h,f'altered bytes: {p}')
    return json.loads(raw)


def rejected(function):
    try:
        function()
    except ValueError:
        return True
    return False


def gate_audit(c,p):
    proposal=HERE/'synthetic_review_unwritten.json'
    current={'source_sha256':MAIN,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,
             'recenter_design_sha256':DESIGN}
    review={'status':c.REVIEW,'checks':{'synthetic_v1':True},**current}
    digest='1'*64
    raw=c.read
    file_read=Path.read_bytes
    blocked_names={'cache_packet.json','field.json','verification.json'}
    def metadata_only(path):
        require(path.name not in blocked_names,'actual cache or original field opened during cold gate')
        return file_read(path)
    def read(path,h):
        if path==proposal:
            require(h==digest,'changed synthetic approval SHA')
            return review
        return raw(path,h)
    stop=AssertionError('actual seed or proof during gate')
    with patch.object(Path,'read_bytes',metadata_only),patch.object(c,'read',side_effect=read),\
         patch.object(c,'extract_actual',side_effect=stop),patch.object(c,'new_seed',side_effect=stop),\
         patch.object(c,'decode',side_effect=stop),patch.object(c.ad,'assemble',side_effect=stop),\
         patch.object(c.ad,'variation',side_effect=stop):
        protocol,bindings=c.gate(proposal,digest)
        require(protocol==p and bindings=={**current,'prelaunch_review_sha256':digest},
                'real roundtrip with frozen preflight and hashes failed')
        for field in ('protocol_sha256','preflight_sha256','recenter_design_sha256','source_sha256','status'):
            old=review[field]
            review[field]='0'*64
            require(rejected(lambda:c.gate(proposal,digest)),f'mixed launch phase accepted: {field}')
            review[field]=old
    return current


def old_scope(c,protocol):
    out=[]
    for expected,item in zip(ORDER,protocol['cases'],strict=True):
        r=item['case']
        wave,record=c.CASES[expected]
        require(r['case_id']==expected and r['record_id']==record and len(item['original_failed_radii'])>0,
                'two intended unresolved old cases and preserved failed radii required')
        basis=c.path_for(item['native_proof'])
        n=pin(basis,item['native_proof_sha256'])
        require(n['status']=='UNRESOLVED_INDEPENDENT_CACHED_ROOT'
                and n['chosen_radius_exact'] is None and n['posterior_error_upper_exact'] is None
                and n['case_id']==expected and n['field_sha256']==r['field_sha256']
                and n['cache_packet_sha256']==item['cache_packet_sha256']
                and n['per_tensor_sha256']==item['per_tensor_sha256']
                and n['embedded_cache_sha256']==item['cache_payload_sha256']
                and Q(n['observer_radius_budget_exact'])==Q(r['e_goal_dimensionless_exact']),
                'native original case or identical measurement scope drift')
        attempt=c.select(n['own_radius_attempts'],radius_exact=str(c.RADIUS))
        bands=[(Q(x['lower_exact']),Q(x['upper_exact'])) for x in attempt['signed_center_force_intervals']]
        require(len(bands)==611 and all(a<=b for a,b in bands)
                and 0<max(abs(a+b)/2 for a,b in bands)<=Q(attempt['eta_upper_exact']),
                'old signed RFQ, rather than a new seed, is the proposal authority')
        phase=c.NATIVE/'actual_replay'/wave
        os=pin(phase/'os_exit_observation.json',item['native_os_sha256'])
        top=pin(phase/'receipt.json',os['receipt_sha256'])
        end=pin(phase/'terminal.json',os['terminal_sha256'])
        prod=c.OLD/'actual'/wave
        a=pin(prod/'root_terminal_admission.json',top['root_terminal_admission_sha256'])
        t=prod/'terminal.json'
        tr=pin(t,a['sha256'][c.filename(t)])
        child=prod/expected/'terminal.json'
        exit=prod/expected/'observed_os_exit.json'
        ct=pin(child,a['sha256'][c.filename(child)])
        co=pin(exit,a['sha256'][c.filename(exit)])
        require(all(type(x) is int and x==0 for x in (a['observed_exit_code'],os['observed_exit_code'],end['exit_code'],
                                                       co['observed_exit_code']))
                and all(type(x) is int and x==os['pid']>0 for x in (os['pid'],top['pid'],end['pid']))
                and all(type(x) is int and x==a['supervisor_pid']>0 for x in (a['supervisor_pid'],tr['pid'],co['supervisor_pid'],ct['supervisor_pid']))
                and all(type(x) is int and x==co['pid']>0 for x in (co['pid'],ct['pid']))
                and co['terminal_sha256']==a['sha256'][c.filename(child)],
                'actual archived parent/worker process identities must be literal integers')
        out.append({'case_id':expected,'wave_id':wave,'record_id':record,'native_phase_sha256':os['receipt_sha256'],
                    'native_case_sha256':item['native_proof_sha256'],'old_cache_packet_sha256':item['cache_packet_sha256'],
                    'old_frozen_R_sha256':item['per_tensor_sha256']['preconditioner'],
                    'observer_goal_exact':r['e_goal_dimensionless_exact'],
                    'signed_midpoint_norm_preview_exact':str(max(abs(a+b)/2 for a,b in bands)),
                    'no_real_seed_generated':True})
    require(tuple(row['case_id'] for row in out)==ORDER and not c.ACTUAL.exists(),'cold review requires two exact planned cases')
    return out


def rounding_and_budget(c):
    require(Q.from_float(c.round_binary(1+Q(1,2**53)))==1
            and Q.from_float(c.round_binary(1+Q(3,2**53)))==1+Q(1,2**51)
            and Q.from_float(c.round_binary(-1-Q(3,2**53)))==-1-Q(1,2**51),
            'ties must round to the actual even adjacent binary value on both sides')
    # The posterior passes in isolation but exactly consumes the old budget after adding shift.
    eta,z0,dz,sigma,goal,xi,r=map(Q,('1/8','1/4','1/4','1','5/4','1/16','1'))
    fail=c.classify(eta,z0,dz,sigma,goal,xi,r)
    good=c.classify(eta,z0,dz,sigma,goal+Q(1,10),xi,r)
    require(fail['strict_contraction_selfmap'] is True and fail['strict_old_observer_goal_pass'] is False
            and Q(fail['posterior_error_upper_exact'])==Q(1,4)<goal
            and Q(fail['shift_plus_posterior_exact'])==goal
            and good['strict_old_observer_goal_pass'] is True
            and Q(good['physical_nodal_coordinate_error_um_exact'])==Q(67,2)*(sigma+Q(1,4)+xi),
            'strict old budget or once-accounted addition gap failed')
    require(c.classify(Q(1,2),Q(1,4),Q(1,4),Q(0),Q(2),Q(0),Q(1))['strict_contraction_selfmap'] is False
            and c.classify(Q(0),Q(1),Q(0),Q(0),Q(2),Q(0),Q(1))['strict_contraction_selfmap'] is False,
            'strict self-map at equality or derivative equality should fail')


def full_synthetic_trial(c):
    # This routes the actual case_trial function through a small nonphysical cache and
    # an inert fresh-bound callback. No 611-vector or mechanics routine can run.
    a={'u':np.array([1.,-.5],dtype='<f8'),'reference':np.array([.125,.25],dtype='<f8'),
       'free':np.array([0,1],dtype='<i8'),'volume_weights':np.array([.25,.5],dtype='<f8'),
       'preconditioner':np.eye(2,dtype='<f8'),'energy_scale':np.array([1.5],dtype='<f8')}
    buf=BytesIO()
    np.savez_compressed(buf,**a)
    raw=buf.getvalue()
    old={'status':'FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE','parameters':{'fixture':'unaltered'},'header_marker':'fixed',
         'frozen_binary_cache':{'data':base64.b64encode(raw).decode(),'sha256':c.sha(raw),'bytes':len(raw),
                                'array_shapes':{k:list(v.shape) for k,v in a.items()}}}
    old_rows=[{'lower_exact':'1/8','upper_exact':'1/8'},{'lower_exact':'-1/16','upper_exact':'0'}]
    wanted=np.array([7/8,-15/32],dtype='<f8')
    item={'case':{'case_id':'fixture','record_id':'fixture','field_sha256':'a'*64,
                  'observer_variant':'current_disk__spatial_max','frozen_conditional_allowance_um_exact':'2',
                  'e_goal_dimensionless_exact':'1','xi_dimensionless_exact':'1/100'},
          'per_tensor_sha256':{k:c.array_sha(v) for k,v in a.items()},'native_proof_sha256':'b'*64,
          'cache_packet_sha256':'c'*64,'cache_payload_sha256':c.sha(raw),'old_result_sha256':'d'*64,
          'original_failed_radii':[{'status':'UNRESOLVED'}]}
    archived={}
    evaluated=[]
    def save(path,doc):
        require(path.name not in archived,'duplicate synthetic output')
        archived[path.name]=doc
        return c.sha(json.dumps(doc,sort_keys=True).encode())
    def prove(array,params,drive):
        require(params==old['parameters'] and drive==Q(3,10)
                and np.array_equal(array['u'],wanted) and not np.array_equal(array['u'],a['u'])
                and all(c.array_sha(v)==c.array_sha(a[k]) for k,v in array.items() if k!='u'),
                'fresh point and whole-cube callback did not receive the right changed center and same operator')
        evaluated.append(c.array_sha(array['u']))
        return c.RADIUS/8,Q(1,8),Q(1,8),{'fixture':'fresh'},[]
    with patch.object(c,'extract_actual',return_value=(old,a,Q(3,10),old_rows)),\
         patch.object(c,'save',side_effect=save),patch.object(c,'prove_after_recenter',side_effect=prove):
        rec=c.case_trial(item,HERE/'synthetic_no_outputs',{'fixture':'phase'})
    require(list(archived)==['exact_binary_step.json','recentered_cache_packet.json'] and len(evaluated)==1
            and rec['status']=='CANDIDATE_RECENTERED_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW'
            and rec['root_certified'] is False and rec['fresh_point_and_whole_cube_AD'] is True
            and rec['fresh_derivative_center_u1_sha256']==evaluated[0]
            and rec['old_failed_radius_attempts']==item['original_failed_radii']
            and rec['original_observer_field_replaced'] is False and rec['original_scores_or_allowance_changed'] is False,
            'synthetic producer callback and frozen observer proof did not complete end-to-end')
    recache=archived['recentered_cache_packet.json']
    bins=c.decode(recache)
    require(recache['recenter']['updated_tensor_keys']==['u'] and recache['header_marker']==old['header_marker']
            and recache['parameters']==old['parameters'],'unexpected change to archived header or model')
    for key in ('reference','free','volume_weights','preconditioner','energy_scale'):
        changed={**bins,key:bins[key].copy()}
        changed[key].flat[0]+=1
        require(rejected(lambda:c.verify_tensors(a,changed,item['per_tensor_sha256'],wanted)),
                f'changed operator tensor should be rejected: {key}')


def main():
    require(sha(SOURCE)==sha(PREP/'transfer_ambiguity_recenter_cached_frozen.py')==MAIN,'source freeze differs')
    protocol=pin(PREP/'protocol.json',PROTOCOL)
    pf=pin(PREP/'preflight.json',PREFLIGHT)
    require(pf['status']=='PASS' and pf['source_sha256']==MAIN and pf['protocol_sha256']==PROTOCOL
            and all(v is True for v in pf['checks'].values()) and pf['actual_production_jobs']==0
            and tuple(protocol['case_order'])==ORDER and Q(protocol['radius_exact'])==Q(1,2**36)
            and Q(protocol['length_um_exact'])==Q(67,2),'frozen scope or bounds drift')
    sys.path.insert(0,str(ROOT/'scripts'))
    import transfer_ambiguity_recenter_cached as c
    require(Path(c.__file__).resolve()==SOURCE and Path(c.ad.__file__).resolve()==ROOT/'scripts/review_algebraic_root.py',
            'wrong source or AD-kernel import')
    gate=gate_audit(c,protocol)
    metadata=old_scope(c,protocol)
    rounding_and_budget(c)
    full_synthetic_trial(c)
    f={'case_id':'synthetic','pid':17,'supervisor_pid':16}
    t={**f,'status':'FINISHED','exit_code':0,'results_sha256':'e'*64}
    c.child_admit(0,t,f,'e'*64)
    require(rejected(lambda:c.child_admit(False,t,f,'e'*64))
            and rejected(lambda:c.child_admit(0,{**t,'case_id':'wrong'},f,'e'*64))
            and rejected(lambda:c.child_admit(0,t,f,'f'*64)),
            'nonliteral OS exit, foreign record or wrong bytes accepted')
    alias=not rejected(lambda:c.child_admit(0,{**t,'pid':17.0},f,'e'*64))
    rec={'status':c.REVIEW,**gate,
         'checks':{'cold_exact_gate_roundtrips_frozen_preflight_and_accepts_only_current_review':True,
                   'two_exact_original_unresolved_case_scopes_and_archived_literal_OS_chains_verified':True,
                   'one_exact_round_to_nearest_ties_even_and_sign_verified_synthetically':True,
                   'unchanged_R_reference_maps_weights_normalization_and_material_verified_synthetically':True,
                   'fresh_new_center_point_and_whole_cube_AD_routed_without_reusing_old_Z0':True,
                   'strict_shift_plus_posterior_compared_with_old_goal_and_old_gap_counted_once':True,
                   'root_and_observer_boundary_equalities_rejected':True,
                   'actual_worker_uses_wait_and_serial_exclusive_output_gates':True,
                   'producer_is_candidate_only_and_old_artifacts_retained':True},
         'review_source_sha256':sha(Path(__file__)),'same_radius_exact':str(c.RADIUS),'archived_metadata_readback':metadata,
         'schema_note':{'terminal_pid_numeric_alias_can_equal_typed_expected_pid':alias,
                        'actual_review_required':'Require literal integer supervisor/worker PID in every new JSON side of the observed OS chain. Frozen source emits os.getpid() and observed return codes.'},
         'actual_seed_cache_model_AD_or_solver_jobs':0,
         'scope':'Independent cold approval for two frozen recenter trials. Synthetic cache/seed fixture does not read actual source fields or caches and substitutes the AD bound callback. No root, minimum, continuum or observer fit claim.'}
    path=HERE/'receipt.json'
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(rec,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':rec['status'],'receipt_sha256':sha(path),'actual_jobs':0}),flush=True)


if __name__=='__main__':
    main()
