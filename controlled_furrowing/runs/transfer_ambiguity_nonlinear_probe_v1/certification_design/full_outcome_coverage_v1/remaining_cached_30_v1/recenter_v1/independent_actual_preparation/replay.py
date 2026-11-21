"""Independently replay the two frozen cached-energy recenter candidates."""

import argparse
from base64 import b64decode
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from math import isfinite, nextafter, inf
import os
from pathlib import Path
import re
import struct
import sys
import traceback
from types import SimpleNamespace
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve()
RECENTER = HERE.parent
ROOT = RECENTER.parents[5]
PREP = RECENTER/'preparation'
PRODUCED = RECENTER/'actual'
OUT = HERE/'actual_replay'
CALLER = ROOT/'scripts/transfer_ambiguity_recenter_cached.py'
CALLER_PROTOCOL = PREP/'protocol.json'
CALLER_PREFLIGHT = PREP/'preflight.json'
CALLER_REVIEW = PREP/'independent_prelaunch/receipt.json'
MAPS = ROOT/'scripts/review_algebraic_root.py'
MATH = ROOT/'scripts/algebraic_root_feasibility.py'
DESIGN = RECENTER.parent/'negative_K_recenter_v1/receipt.json'
RESULT_SHA = 'c63b1259067a22e4418f42ad0e1e6307036c57435239c7d5266b017d29beeef7'
TERMINAL_SHA = '1a9aeae9de12b0be1e949253e699215965ec50aa345060978280ced0280406ea'
CASES = ('max_negative_contrast_change_small_16', 'max_contrast_change_small_13')
SUPERVISOR = 21924
EXEC_SESSION = 38756
RADIUS = Q(1,2**36)
LENGTH = Q(67,2)
N = 611
PRECISION = 256
COUNTS = {'bulk':9216, 'surface':288}
MARGINS = {'bulk_det_min','bulk_hoop_min','surface_radius_min','surface_tangent_squared_min','normalized_weight_min'}
ADMISSION = 'PASS_ROOT_OBSERVED_CACHED_RECENTER_TERMINALS'
REVIEW = 'PASS_ROOT_RECENTER_INDEPENDENT_ACTUAL_REPLAY_SOURCE'
SCOPE = ('Two fresh closed-form cached force and Hessian replays at the recorded one-step exact binary centers. '
         'Same original frozen scalar quadrature operator, original observers and budget. No exact spatial '
         'force/integration, whole-body geometry, stable path, minimum, PDE or refit.')
PINS = {
    CALLER:'9b149c5889c8dc6a77b9407674200d9e6bcbbee8e9875f12a13b97f7962de7d4',
    PREP/'transfer_ambiguity_recenter_cached_frozen.py':'9b149c5889c8dc6a77b9407674200d9e6bcbbee8e9875f12a13b97f7962de7d4',
    CALLER_PROTOCOL:'bb5f0aef19bf5a7d5a6f954d7b41f03219bd19056be4e4eb4672454a588bec63',
    CALLER_PREFLIGHT:'9d38373a9f9895880ff1d03d6ba5c1165a790cbff8123787cee9a0ddf6fe0a35',
    CALLER_REVIEW:'6caa412e20bb6c90f9be6a87c802d1b2d6f3e1b5a697b2e098ea1795e887a463',
    DESIGN:'d9621eef3cc7fec129b2f597dd4ab584393a6ed0b610051498734fdd2cad1d04',
    MAPS:'6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a',
    MATH:'7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
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
    p=pinned(CALLER_PROTOCOL,PINS[CALLER_PROTOCOL])
    f=pinned(CALLER_PREFLIGHT,PINS[CALLER_PREFLIGHT])
    c=pinned(CALLER_REVIEW,PINS[CALLER_REVIEW])
    binding={'source_sha256':PINS[CALLER],'protocol_sha256':PINS[CALLER_PROTOCOL],
             'preflight_sha256':PINS[CALLER_PREFLIGHT],'recenter_design_sha256':PINS[DESIGN]}
    require(c['status']=='PASS_INDEPENDENT_COLD_CACHED_RECENTER_PRELAUNCH' and yes(c['checks'])
            and all(c[k]==v for k,v in binding.items())
            and f['status']=='PASS' and yes(f['checks']) and strict_int_zero(f['actual_production_jobs'])
            and all(f[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and p['status']=='FROZEN_TWO_OLD_CACHED_OPERATOR_RECENTER_TRIALS'
            and tuple(p['case_order'])==tuple(x['case']['case_id'] for x in p['cases'])==CASES
            and p['source_sha256']==PINS[CALLER] and exact(p['radius_exact'])==RADIUS
            and exact(p['length_um_exact'])==LENGTH and p['precision_bits']==PRECISION
            and strict_int_zero(p['actual_seed_cache_model_AD_proposal_or_solve_jobs']),
            'frozen original recenter plan and reviewed scope required')
    pins=p['source_and_input_sha256']
    require(pins[rel(DESIGN)]==PINS[DESIGN] and pins[rel(MAPS)]==PINS[MAPS],'original inputs changed')
    for s,h in pins.items():
        if s.endswith('.py'):
            pinned(link(s),h,parse=False)
    return {'status':'FROZEN_INDEPENDENT_ACTUAL_RECENTER_REPLAY','scope':SCOPE,'source_sha256':sha(SOURCE),
            'producer_binding':{**binding,'prelaunch_review_sha256':PINS[CALLER_REVIEW]},
            'pinned_sources':{rel(k):h for k,h in PINS.items()},'original_inputs':pins,'ordered_cases':list(CASES),
            'observed_supervisor_pid':SUPERVISOR,'observed_exec_session':EXEC_SESSION,
            'expected_producer_results_sha256':RESULT_SHA,'expected_producer_terminal_sha256':TERMINAL_SHA,
            'required_root_admission_status':ADMISSION,'source_review_status':REVIEW,
            'radius_exact':str(RADIUS),'precision_bits':PRECISION,'rank':N,
            'method':('Independently reproduce nearest-even rounding against both adjacent binary numbers. '
                      'Decode both caches, compare every tensor except u and all original physical header fields. '
                      'Share only the exact frozen-array local map generator. Use the separate analytic bulk_local '
                      'and surface_local derivatives and closed-form whole-ball variation enclosure at u1. '
                      'Fuse the verified posterior and exact shift with the original observer gap.'),
            'actual_seed_cache_or_AD_before_gate':0}


def verify_child(case,child,worker,row,digests,binding,mark):
    start,launched,observed,end,result=child
    pid=launched['pid']
    require(positive_int(pid) and pid!=SUPERVISOR
            and all(positive_int(m['pid']) and m['pid']==pid
                    and positive_int(m['supervisor_pid']) and m['supervisor_pid']==SUPERVISOR
                    and m['case_id']==case for m in (*child,worker,row))
            and strict_int_zero(worker['observed_exit_code']) and strict_int_zero(row['observed_exit_code'])
            and strict_int_zero(observed['observed_exit_code']) and strict_int_zero(end['exit_code'])
            and end['status']=='FINISHED' and start['status']=='RUNNING'
            and start['bindings']==end['bindings']==result['bindings']==binding
            and worker=={**observed,'os_exit_sha256':digests['observed_os_exit.json']}
            and row=={**observed,'results_sha256':digests['results.json'],'status':result['status']}
            and observed['terminal_sha256']==digests['terminal.json']
            and end['results_sha256']==digests['results.json'],
            'literal PIDs and genuine bounded child success chain required')
    require(result['status']=='CANDIDATE_RECENTERED_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW'
            and result['root_certified'] is False and result['fresh_point_and_whole_cube_AD'] is True
            and result['independent_actual_recentered_AD_replay_required'] is True
            and result['original_observer_field_replaced'] is False and result['original_scores_or_allowance_changed'] is False
            and result['exact_spatial_or_minimum_claimed'] is False
            and strict_int_zero(result['new_equilibrium_solves_refits_or_preconditioners'])
            and exact(result['trial_radius_exact'])==RADIUS and type(result['AD_precision_bits']) is int
            and result['AD_precision_bits']==PRECISION
            and result['recentered_cache_packet_sha256']==digests['recentered_cache_packet.json']
            and result['exact_step_sha256']==digests['exact_binary_step.json']
            and result['archived_field_sha256']==mark['case']['field_sha256']
            and result['old_failed_result_sha256']==mark['old_result_sha256']
            and result['old_failed_radius_attempts']==mark['original_failed_radii']
            and result['old_native_anchor_sha256']==mark['native_proof_sha256'],
            'same exact producer candidate, old center and failed attempt authority required')
    return pid


def admitted(path,h,p):
    a=pinned(path,h)
    require(a['status']==ADMISSION and positive_int(a['supervisor_pid']) and a['supervisor_pid']==SUPERVISOR
            and positive_int(a['unified_exec_session']) and a['unified_exec_session']==EXEC_SESSION
            and strict_int_zero(a['observed_exit_code']),'externally observed root process exit required')
    parent=('started.json','terminal.json','results.json')
    child=('started.json','supervisor_child_started.json','observed_os_exit.json','terminal.json','results.json')
    deferred=('exact_binary_step.json','recentered_cache_packet.json')
    wanted={rel(PRODUCED/k) for k in parent}
    wanted.update(rel(PRODUCED/c/k) for c in CASES for k in (*child,*deferred))
    observed_pins=a['sha256']
    require(set(observed_pins)==wanted and all(digest_text(h) for h in observed_pins.values())
            and observed_pins[rel(PRODUCED/'results.json')]==RESULT_SHA
            and observed_pins[rel(PRODUCED/'terminal.json')]==TERMINAL_SHA,'complete two-case actual artifact manifest required')
    start,end,results=(pinned(PRODUCED/k,observed_pins[rel(PRODUCED/k)]) for k in parent)
    b=p['producer_binding']
    require(all(positive_int(x['pid']) and x['pid']==SUPERVISOR for x in (start,end))
            and start['status']=='RUNNING' and end['status']=='FINISHED' and strict_int_zero(end['exit_code'])
            and end['results_sha256']==RESULT_SHA and start['bindings']==end['bindings']==results['bindings']==b
            and tuple(start['case_order'])==tuple(x['case_id'] for x in end['workers'])
            ==tuple(x['case_id'] for x in results['cases'])==CASES
            and results['status']=='COMPLETE_TWO_RECENTER_ATTEMPTS_AWAITING_INDEPENDENT_REVIEW'
            and results['root_certified'] is False,'completed same producer phase and ordered cases required')
    plan=pinned(CALLER_PROTOCOL,PINS[CALLER_PROTOCOL])
    records=[]
    pids=[]
    for case,w,row,item in zip(CASES,end['workers'],results['cases'],plan['cases'],strict=True):
        folder=PRODUCED/case
        digests={k:observed_pins[rel(folder/k)] for k in (*child,*deferred)}
        data=[pinned(folder/k,digests[k]) for k in child]
        pid=verify_child(case,data,w,row,digests,b,item)
        cmd=data[1]['command']
        require(type(cmd) is list and len(cmd)==12 and cmd[1]=='-B' and Path(cmd[2]).resolve()==CALLER
                and cmd[3:9]==['worker','--case-id',case,'--supervisor-pid',str(SUPERVISOR),'--review']
                and Path(cmd[9]).resolve()==CALLER_REVIEW and cmd[10:]==['--review-sha256',PINS[CALLER_REVIEW]],
                'matching observed independently approved producer command required')
        records.append((item,data[-1],folder,digests))
        pids.append(pid)
    require(len(set(pids))==2,'two distinct real child process identifiers required')
    return a,records,pids


def independent_round(target):
    value=float(target)
    alternatives=(nextafter(value,-inf),value,nextafter(value,inf))
    require(all(isfinite(x) for x in alternatives),'finite binary64 candidate and neighbors required')
    def distance(x):
        exact_error=abs(Q.from_float(x)-target)
        odd=struct.unpack('<Q',struct.pack('<d',x))[0]&1
        return exact_error,odd
    return min(alternatives,key=distance)


def signed_pair(x):
    lo,hi=exact(x['lower_exact']),exact(x['upper_exact'])
    require(lo<=hi,'ordered finite signed bounds required')
    return lo,hi


def exact_seed(u,forces):
    import numpy as np
    require(u.dtype.str=='<f8' and u.ndim==1 and len(u)==len(forces),'complete typed frozen seed input required')
    values=[]
    witness=[]
    for binary,interval in zip(u,forces,strict=True):
        old=Q.from_float(float(binary))
        lo,hi=signed_pair(interval)
        midpoint=(lo+hi)/2
        target=old-midpoint
        rounded=independent_round(target)
        delta=Q.from_float(rounded)-old
        witness.append(list(map(str,(old,midpoint,target,Q.from_float(rounded),delta))))
        values.append(rounded)
    return np.asarray(values,dtype='<f8'),witness,max(abs(exact(row[-1])) for row in witness)


def raw_array_sha(x):
    return sha256(x.tobytes(order='C')).hexdigest()


def arrays(packet):
    import numpy as np
    packed=packet['frozen_binary_cache']
    blob=b64decode(packed['data'],validate=True)
    require(len(blob)==packed['bytes'] and sha256(blob).hexdigest()==packed['sha256'],'changed encoded packet')
    with np.load(BytesIO(blob),allow_pickle=False) as data:
        require(len(set(data.files))==len(data.files),'duplicate NPZ keys')
        out={k:data[k].copy() for k in data.files}
    require({k:list(v.shape) for k,v in out.items()}==packed['array_shapes']
            and all(v.dtype.kind in 'iufb' and np.isfinite(v).all() for v in out.values()),'finite exact cache tensors required')
    return out


def only_u_changes(old,new,prior):
    require(old.keys()==new.keys()==prior.keys() and 'u' in prior,'complete same frozen operator tensor keys')
    for k in prior:
        require(old[k].shape==new[k].shape and old[k].dtype==new[k].dtype and raw_array_sha(old[k])==prior[k],
                'frozen original map/type/shape/hash changed')
        if k!='u':
            require(raw_array_sha(new[k])==prior[k],'same weights, maps, R and global normalization required')
    require(raw_array_sha(new['u'])!=prior['u'],'fresh proposed evaluation center must differ')


def recorded_center(item,result,folder,digests,protocol):
    import numpy as np
    c=item['case']
    links=protocol['original_inputs']
    for s,h in ((item['native_proof'],item['native_proof_sha256']),
                (item['old_result'],item['old_result_sha256']),(item['cache_packet'],item['cache_packet_sha256']),
                (c['field'],c['field_sha256']),(c['verification'],c['verification_sha256'])):
        require(links[s]==h,'original independent input link changed')
    old=pinned(link(item['cache_packet']),item['cache_packet_sha256'])
    raw_proof=pinned(link(item['native_proof']),item['native_proof_sha256'])
    old_r=pinned(link(item['old_result']),item['old_result_sha256'])
    original=pinned(link(c['field']),c['field_sha256'])
    pinned(link(c['verification']),c['verification_sha256'],parse=False)
    new=pinned(folder/'recentered_cache_packet.json',digests['recentered_cache_packet.json'])
    before,after=arrays(old),arrays(new)
    delta_keys={'status','frozen_binary_cache','recenter'}
    require(set(new)==set(old)|{'recenter'} and all(new[k]==v for k,v in old.items() if k not in delta_keys)
            and new['status']=='FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER'
            and old['status']=='FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE'
            and old['case_id']==raw_proof['case_id']==result['case_id']==c['case_id']
            and old['record_id']==raw_proof['record_id']==result['record_id']==c['record_id']
            and old['input_state']==c['field'] and old['input_state_sha256']==raw_proof['field_sha256']==c['field_sha256']
            and old['field_verification_sha256']==c['verification_sha256']
            and old['parameters']==original['parameters']==c['material']['parameters']
            and old['closed_parameter_cell_exact']==c['closed_parameter_cell_exact']
            and old['parameter_order']==c['parameter_order']
            and old['material_parameter_storage_exact']==c['material']['parameter_storage_exact']
            and all(old[k]==v for k,v in c['drive'].items())
            and all(new['frozen_binary_cache'][k]==v for k,v in old['frozen_binary_cache'].items()
                    if k not in ('data','sha256','bytes')),
            'unexpected original target, cached tensor encoding or header mutation')
    dose=exact(old['optical_dose_exact'])
    gain=exact(old['common_gain_exact'])
    drive=exact(old['activation_strength_exact'])
    rcell=old['closed_parameter_cell_exact']
    require(0<drive==dose*gain==exact(raw_proof['exact_requested_strength'])==exact(result['requested_strength_exact'])
            and drive!=exact(old['activation_strength_binary_exact'])==Q.from_float(float(original['activation_strength']))
            and (dose,gain)==(exact(c['step']['dose_exact']),exact(c['material']['gain_exact']))
            and all(exact(lo)==exact(hi) for lo,hi in rcell.values())
            and tuple(map(exact,rcell['common_gain']))==(gain,gain)
            and all(Q.from_float(float(old['parameters'][k]))==exact(rcell[k][0]) for k in rcell if k!='common_gain'),
            'unrounded physical drive and frozen singleton required')
    stable=item['per_tensor_sha256']
    require(stable==raw_proof['per_tensor_sha256'] and (before['u'].shape,before['preconditioner'].shape)==((N,),(N,N))
            and before['preconditioner'].dtype.str==before['u'].dtype.str=='<f8'
            and before['free'].shape==(N,) and before['free'].dtype.kind in 'iu'
            and list(before['free'])==c['free_map'] and len(set(c['free_map']))==N
            and np.array_equal(before['u'],original['free_displacement'])
            and np.array_equal(before['reference'],original['reference_nodes'])
            and np.array_equal(before['triangles'],original['triangles'])
            and Q.from_float(float(before['energy_scale'][0]))==exact(old['normalization_binary_exact'])>0
            and old['frozen_binary_cache']['sha256']==item['cache_payload_sha256']==raw_proof['embedded_cache_sha256']
            and raw_proof['cache_packet_sha256']==item['cache_packet_sha256']==result['original_cache_sha256'],
            'exact original binary cache and saved position required')
    only_u_changes(before,after,stable)
    require(old_r['status']=='UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET'
            and old_r['accepted_radius'] is None and old_r['radius_attempts']==item['original_failed_radii']
            and raw_proof['status']=='UNRESOLVED_INDEPENDENT_CACHED_ROOT'
            and raw_proof['chosen_radius_exact'] is None and raw_proof['posterior_error_upper_exact'] is None
            and raw_proof['rounded_strength_used_in_proof'] is False,
            'original unresolved numerical authority must be precisely retained')
    old_attempt,=[r for r in raw_proof['own_radius_attempts'] if exact(r['radius_exact'])==RADIUS]
    require(old_attempt['status']=='UNRESOLVED_INDEPENDENT_CACHED_CONTRACTION'
            and old_attempt['whole_cube_AD_recomputed_for_this_radius'] is True
            and old_attempt['posterior_error_upper_exact'] is None,'old signed RF anchor must be original independent result')
    seed,rows,sigma=exact_seed(before['u'],old_attempt['signed_center_force_intervals'])
    uh=raw_array_sha(seed)
    evidence=pinned(folder/'exact_binary_step.json',digests['exact_binary_step.json'])
    require(uh==raw_array_sha(after['u'])==result['u1_binary_sha256']==result['fresh_derivative_center_u1_sha256']
            and seed.tobytes()==after['u'].tobytes()
            and evidence['case_id']==c['case_id'] and evidence['ordered_rows']['values']==rows
            and evidence['ordered_rows']['columns']==['old_u','native_signed_RF_midpoint','exact_unrounded_u1','rounded_binary_u1',
                                                      'exact_binary_shift']
            and 0<sigma==exact(evidence['sigma_exact'])==exact(result['actual_shift_sigma_exact'])
            and evidence['native_proof_sha256']==item['native_proof_sha256']
            and evidence['round_once_ties_to_even'] is True and evidence['u0_sha256']==stable['u']
            and evidence['u1_sha256']==uh and result['new_per_tensor_sha256']=={k:raw_array_sha(v) for k,v in after.items()}
            and result['original_R_binary_sha256']==stable['preconditioner'],'actual one-round binary seed differs')
    marker=new['recenter']
    require(marker['original_point_cache_sha256']==item['cache_packet_sha256']
            and marker['original_payload_sha256']==item['cache_payload_sha256']
            and marker['original_field_retains_scoring_authority'] is True
            and marker['old_field_is_not_new_proof_center'] is True and marker['updated_tensor_keys']==['u']
            and marker['u1_sha256']==uh
            and marker['exact_step_sha256']==digests['exact_binary_step.json'] and marker['bindings']==result['bindings'],
            'same old prediction and new binary proof center authority required')
    goal,xi=map(exact,(c['e_goal_dimensionless_exact'],c['xi_dimensionless_exact']))
    require(0<RADIUS<goal==exact(raw_proof['observer_radius_budget_exact'])==exact(result['original_e_goal_exact'])
            and 0<=xi==exact(result['original_binary_addition_gap_exact'])
            and c['observer_variant']==result['same_observer_variant']=='current_disk__spatial_max'
            and exact(c['frozen_conditional_allowance_um_exact'])==exact(result['same_frozen_allowance_um_exact']),
            'same original independently reviewed conditional observer required')
    return after,old['parameters'],drive,sigma,goal,xi,uh


def classify(eta,z,sigma,goal,xi,r):
    require(eta>=0 and z>=0 and sigma>=0 and xi>=0 and 0<r<goal,'exact remaining budget and positive radius required')
    lhs=eta+z*r
    local=z<1 and lhs<r
    distance=eta/(1-z) if local else None
    if local:
        require(0<=distance<r,'posterior and strict self-map disagree')
    transfer=local and sigma+distance<goal
    physical=LENGTH*(xi+sigma+distance) if transfer else None
    return {'strict_contraction_selfmap':local,'strict_old_observer_goal_pass':transfer,
            'strict_selfmap_lhs_upper_exact':str(lhs),'posterior_error_upper_exact':None if distance is None else str(distance),
            'shift_plus_posterior_exact':str(sigma+distance) if local else None,
            'physical_nodal_coordinate_error_um_exact':None if physical is None else str(physical)}


def evaluate_closed_form(core,maps,a,params,drive):
    # The exact linear site map generator is shared. All site derivatives, point
    # assembly and ball variation below use the separate closed-form routines.
    core.ctx.prec=PRECISION
    core.A.cache_clear()
    maps.A.cache_clear()
    count=dict.fromkeys(COUNTS,0)
    def terms():
        for site in maps.terms_from_arrays(a,params,drive):
            count[site['kind']]+=1
            yield site
    g,H,sites,point_J=core.assemble_point(terms(),N)
    require(count==COUNTS and len(sites)==sum(COUNTS.values())
            and point_J is not None and exact(core.pack(point_J)['lower_exact'])>Q(1,10**6),
            'fresh complete positive center analytic domain required')
    rij=[Q.from_float(x) for row in a['preconditioner'].tolist() for x in row]
    R=core.arb_mat(N,N,list(map(core.A,rij)))
    abs_R=core.arb_mat(N,N,[core.A(abs(v)) for v in rij])
    identity=core.arb_mat(N,N)
    for i in range(N):
        identity[i,i]=1
    rf=R*g
    point=[core.pack(rf[i,0]) for i in range(N)]
    eta=max(max(map(abs,signed_pair(v))) for v in point)
    def row_bounds(matrix):
        return [exact(core.pack(sum((core.upper(matrix[i,j]) for j in range(matrix.ncols())),
                                    core.arb(0)))['upper_exact']) for i in range(matrix.nrows())]
    Z0=row_bounds(identity-R*H)
    v,bulk_J,domain=core.variation_rows(sites,N,core.A(RADIUS))
    lower=domain['strict_lower_margins']
    require(domain['term_counts']==COUNTS and set(lower)==MARGINS
            and all(signed_pair(b)[0]>0 for b in lower.values()) and bulk_J is not None
            and signed_pair(core.pack(bulk_J))[0]>Q(1,10**6),'fresh complete positive whole-cube domain required')
    # Surface J keeps a positive lower bound as scale * R * |t|. The exact
    # map generator fixes positive reference R0 and tangent length in scale.
    scale=min(signed_pair(core.pack(t['parameters'][2]))[0] for t in sites if t['kind']=='surface')
    min_R=signed_pair(lower['surface_radius_min'])[0]
    min_t2=signed_pair(lower['surface_tangent_squared_min'])[0]
    surface_J=core.pack(core.A(scale)*core.A(min_R)*core.A(min_t2).sqrt())
    require(signed_pair(surface_J)[0]>0,'fresh surface Jacobian lower bound required')
    DZ=row_bounds(abs_R*v)
    return eta,max(Z0),max(DZ),point,{'counts':COUNTS,'bulk_J_raw':core.pack(bulk_J),
                'surface_J_lower':surface_J,'point_bulk_J_raw':core.pack(point_J),'strict_lower_margins':lower,
                'mapped_whole_variation_rows_exact':list(map(str,DZ))},list(map(str,Z0))


def numeric_case(item,p,folder,digests,protocol,core,maps):
    a,params,drive,sigma,goal,xi,uh=recorded_center(item,p,folder,digests,protocol)
    result={'case_id':p['case_id'],'record_id':p['record_id'],'archived_field_sha256':p['archived_field_sha256'],
            'old_failed_radius_attempts_retained':item['original_failed_radii'],
            'original_cache_packet_sha256':item['cache_packet_sha256'],
            'original_binary_R_sha256':item['per_tensor_sha256']['preconditioner'],
            'recentered_cache_packet_sha256':digests['recentered_cache_packet.json'],
            'exact_step_sha256':digests['exact_binary_step.json'],'producer_case_sha256':digests['results.json'],
            'old_native_signed_point_sha256':item['native_proof_sha256'],'replayed_u1_binary_sha256':uh,
            'exact_requested_strength':str(drive),'archived_observer_goal_exact':str(goal),
            'actual_binary_center_shift_sigma_exact':str(sigma),'archived_addition_gap_xi0_exact':str(xi),
            'radius_exact':str(RADIUS),'precision_bits':PRECISION,'original_field_or_scores_replaced':False,
            'no_new_preconditioner_or_equilibrium_solves':True,'exact_spatial_or_minimum_claimed':False}
    try:
        eta,z0,dz,rf,domain,center_rows=evaluate_closed_form(core,maps,a,params,drive)
    except ValueError as error:
        return {**result,'status':'UNRESOLVED_INDEPENDENT_RECENTERED_CACHED_ANALYTIC_DOMAIN','failure':str(error)}
    require(len(p['signed_center_force_intervals'])==len(rf)==N and all(
        lo<=other_hi and other_lo<=hi
        for (lo,hi),(other_lo,other_hi) in zip(map(signed_pair,rf),map(signed_pair,p['signed_center_force_intervals']),strict=True)),
        'independent closed-form and producer Jet gradients enclosures disagree')
    dec=classify(eta,z0+dz,sigma,goal,xi,RADIUS)
    status=('PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET'
            if dec['strict_old_observer_goal_pass'] else 'UNRESOLVED_INDEPENDENT_RECENTERED_CACHED_CONTRACTION_OR_SHIFT')
    return {**result,**dec,'status':status,'signed_center_force_intervals':rf,'eta_upper_exact':str(eta),
            'center_Z0_upper_exact':str(z0),'variation_upper_exact':str(dz),'ZQ_upper_exact':str(z0+dz),
            'point_defect_rows_exact':center_rows,'positive_whole_cached_domain':domain,
            'fresh_site_counts':COUNTS,'fresh_center_and_entire_cube_recomputed':True}


def preflight_document(source_hash,protocol_hash,checks):
    return {'status':'PASS','source_sha256':source_hash,'protocol_sha256':protocol_hash,'checks':checks,
            'actual_producer_packets_opened':0,'actual_new_center_arrays_read_or_written':0,'actual_mechanics_jobs':0}


def reviewed_phase(p,pf,r,binding):
    # The external reader checks the preflight SHA. Its content binds two other documents.
    require(r['status']==REVIEW and yes(r['checks']) and all(r[k]==v for k,v in binding.items())
            and p['source_sha256']==binding['source_sha256'] and pf['status']=='PASS' and yes(pf['checks'])
            and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and all(strict_int_zero(pf[k]) for k in ('actual_producer_packets_opened',
                                                     'actual_new_center_arrays_read_or_written','actual_mechanics_jobs')),
            'externally reviewed exact frozen recenter replay preflight required')


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def tiny_closed_form_fixture():
    import numpy as np
    sys.path.insert(0,str(ROOT/'scripts'))
    import algebraic_root_feasibility as core
    require(Path(core.__file__).resolve()==MATH and sha(MATH)==PINS[MATH],'synthetic math module drift')
    def tiny_sites(_,__,___):
        A=core.A
        yield {'kind':'bulk','q':list(map(A,(1,0,0,1,1))),'B':[[A(v)] for v in (0,0,0,0,1)],
               'indices':[0],'weight':A(1),'parameters':tuple(map(A,(1,2,3)))}
        yield {'kind':'surface','q':list(map(A,(2,1,0))),'B':[[A(v)] for v in (1,0,0)],
               'indices':[0],'weight':A(1),'parameters':tuple(map(A,(2,0,1)))}
    mock=SimpleNamespace(A=core.A,terms_from_arrays=tiny_sites)
    with patch.object(sys.modules[__name__],'N',1),patch.object(sys.modules[__name__],'COUNTS',{'bulk':1,'surface':1}):
        eta,z0,dz,rf,domains,rows=evaluate_closed_form(core,mock,{'preconditioner':np.array([[.125]])},{},Q(1,3))
    lo,hi=signed_pair(rf[0])
    return lo<=Q(1,4)<=hi and eta==Q(1,4)==z0 and dz>0 and list(map(Q,rows))==[Q(1,4)] \
        and domains['counts']=={'bulk':1,'surface':1} and signed_pair(domains['surface_J_lower'])[0]>0


def fixtures(p):
    import numpy as np
    binding={'source_sha256':p['source_sha256'],'protocol_sha256':'1'*64,'preflight_sha256':'2'*64,
             'producer_root_admission_sha256':'3'*64}
    pf=preflight_document(binding['source_sha256'],binding['protocol_sha256'],{'synthetic':True})
    r={'status':REVIEW,'checks':{'synthetic':True},**binding}
    reviewed_phase(p,pf,r,binding)
    u=np.array([1.,-.5],dtype='<f8')
    forces=[{'lower_exact':'1/8','upper_exact':'1/8'},{'lower_exact':'-1/16','upper_exact':'0'}]
    shifted,rows,s=exact_seed(u,forces)
    before={'u':u,'R':np.eye(2,dtype='<f8'),'map':np.array([1,2]),'weight':np.array([.25,.5])}
    after={**before,'u':shifted}
    hashes={k:raw_array_sha(v) for k,v in before.items()}
    only_u_changes(before,after,hashes)
    bounds=classify(Q(1,8),Q(1,2),Q(1),Q(5,4),Q(1,16),Q(1))
    good=classify(Q(1,8),Q(1,2),Q(1),Q(6,4),Q(1,16),Q(1))
    test=classify(Q(1,2),Q(1,2),Q(0),Q(2),Q(0),Q(1))
    return {'gate_roundtrip_with_actual_preflight_schema_and_wrong_external_pins':(
             'preflight_sha256' not in pf and rejects(lambda:reviewed_phase(p,pf,r,{**binding,'preflight_sha256':'4'*64}))
             and rejects(lambda:reviewed_phase(p,{**pf,'protocol_sha256':'4'*64},r,binding))),
            'exact_nearest_even_at_competing_positive_and_negative_ties':(
             Q.from_float(independent_round(1+Q(1,2**53)))==1
             and Q.from_float(independent_round(1+Q(3,2**53)))==1+Q(1,2**51)
             and Q.from_float(independent_round(-1-Q(3,2**53)))==-1-Q(1,2**51)),
            'independent_exact_signed_seed_oracle_and_array_mutation_rejected':(
             np.array_equal(shifted,[7/8,-15/32]) and rows[1][-1]=='1/32' and s==Q(1,8)
             and rejects(lambda:only_u_changes(before,{**after,'R':2*after['R']},hashes))
             and rejects(lambda:only_u_changes(before,{**after,'map':np.array([2,1])},hashes))
             and rejects(lambda:only_u_changes(before,{**after,'weight':np.array([.3,.5])},hashes))),
            'strict_contraction_and_shift_inclusive_budget_equalities_guarded':(
             bounds['strict_contraction_selfmap'] is True and bounds['strict_old_observer_goal_pass'] is False
             and good['strict_old_observer_goal_pass'] is True
             and exact(good['physical_nodal_coordinate_error_um_exact'])==LENGTH*Q(21,16)
             and test['strict_contraction_selfmap'] is False),
            'synthetic_one_dof_closed_form_bulk_and_surface_assembly_and_variation':tiny_closed_form_fixture(),
            'reject_float_boolean_process_ids_or_exit_and_actual_outputs_absent':(
             positive_int(17) and not positive_int(17.) and not positive_int(True)
             and strict_int_zero(0) and not strict_int_zero(False) and not strict_int_zero(0.)
             and not OUT.exists())}


def prepare():
    p=static_protocol()
    checks=fixtures(p)
    require(yes(checks),'cold recenter replay fixture failure')
    with (HERE/'replay_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    ph=save(HERE/'protocol.json',p)
    fh=save(HERE/'preflight.json',preflight_document(sha(SOURCE),ph,checks))
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':fh}),flush=True)


def run(args):
    p=pinned(HERE/'protocol.json',args.protocol_sha256)
    pf=pinned(HERE/'preflight.json',args.preflight_sha256)
    review=pinned(args.source_review,args.source_review_sha256)
    binding={'source_sha256':sha(SOURCE),'protocol_sha256':args.protocol_sha256,'preflight_sha256':args.preflight_sha256,
             'producer_root_admission_sha256':args.admission_sha256}
    reviewed_phase(p,pf,review,binding)
    require(sha(HERE/'replay_frozen.py')==binding['source_sha256'] and static_protocol()==p,'changed frozen replay phase')
    a,current,pids=admitted(args.admission,args.admission_sha256,p)
    for s,h in p['original_inputs'].items():
        pinned(link(s),h,parse=False)
    sys.path.insert(0,str(ROOT/'scripts'))
    import algebraic_root_feasibility as core
    import review_algebraic_root as maps
    require(Path(core.__file__).resolve()==MATH and Path(maps.__file__).resolve()==MAPS
            and sha(MATH)==PINS[MATH] and sha(MAPS)==PINS[MAPS] and core.ctx is maps.ctx,
            'exact reviewed core and frozen-array map identities required')
    OUT.mkdir(exist_ok=False)
    mark={'pid':os.getpid(),**binding,'source_review_sha256':args.source_review_sha256}
    save(OUT/'started.json',{'status':'RUNNING',**mark})
    terminal={**mark,'status':'FAILED','exit_code':1}
    try:
        rows=[]
        valid=0
        for item,reported,folder,digests in current:
            record=numeric_case(item,reported,folder,digests,p,core,maps)
            path=OUT/(record['case_id']+'.json')
            dh=save(path,{**record,'replay_phase':mark})
            pass_case=record['status']=='PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET'
            valid+=int(pass_case)
            rows.append({k:record.get(k) for k in ('case_id','record_id','status','archived_field_sha256','replayed_u1_binary_sha256',
                    'recentered_cache_packet_sha256','original_binary_R_sha256','archived_observer_goal_exact',
                    'actual_binary_center_shift_sigma_exact','posterior_error_upper_exact','shift_plus_posterior_exact',
                    'physical_nodal_coordinate_error_um_exact','ZQ_upper_exact','strict_selfmap_lhs_upper_exact')}
                    |{'native_file':rel(path),'native_sha256':dh})
            print(json.dumps({'case_id':record['case_id'],'status':record['status'],'native_sha256':dh}),flush=True)
        require(tuple(r['case_id'] for r in rows)==CASES,'two same-target replays required')
        status=('PASS_INDEPENDENT_TWO_RECENTERED_CACHED_ROOTS_AND_OLD_OBSERVER_BUDGETS'
                if valid==2 else 'INDEPENDENT_RECENTER_REPLAY_RETAINED_UNRESOLVED_CASE')
        receipt={'status':status,**mark,'scope':SCOPE,'producer_observed_exec_session':a['unified_exec_session'],
                 'producer_observed_pid':a['supervisor_pid'],'producer_observed_child_pids':pids,
                 'producer_results_sha256':RESULT_SHA,'producer_terminal_sha256':TERMINAL_SHA,
                 'producer_metadata_and_packet_sha256':a['sha256'],'case_results':rows,'cached_local_roots_certified':valid,
                 'exact_spatial_root_count':0,'new_refits_equilibrium_solves_or_R':0,
                 'archived_observers_predictions_or_scores_replaced':False}
        rh=save(OUT/'receipt.json',receipt)
        terminal.update(status='FINISHED',exit_code=0,receipt_sha256=rh)
        print(json.dumps({'status':status,'receipt_sha256':rh}),flush=True)
    except BaseException as error:
        save(OUT/'failure.json',{'error':str(error),'traceback':traceback.format_exc(),**mark})
        raise
    finally:
        save(OUT/'terminal.json',terminal)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest='mode',required=True)
    sub.add_parser('prepare')
    actual=sub.add_parser('run')
    for k in ('protocol','preflight','admission','source-review'):
        actual.add_argument('--'+k+'-sha256',required=True)
        if k in ('admission','source-review'):
            actual.add_argument('--'+k,required=True,type=lambda s:Path(s).resolve())
    args=cli.parse_args()
    if args.mode=='prepare':
        prepare()
    else:
        run(args)


