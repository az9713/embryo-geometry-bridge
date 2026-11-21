"""Prepare and separately run exact-binary Sym(R) proposals for the remaining 32 targets."""
import argparse
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory
from time import perf_counter
import traceback
from unittest.mock import patch

from flint import ctx as arb_context
import numpy as np

import algebraic_root_feasibility as core
import cohort_energy_stability as certifier
import review_algebraic_root as maps
import review_cohort_energy_stability as endpoint
import transfer_ambiguity_remaining32_exact_targets as adapter


SOURCE=Path(__file__).resolve()
ROOT=SOURCE.parents[1]
EXACT=adapter.X.OUT
TARGET=EXACT/'remaining32_targets_v1'
GEOMETRY=EXACT/'remaining32_geometry_v1'
JOIN=EXACT/'remaining32_bridge_v1/admission'
OUT=EXACT/'remaining32_preconditioner_spd_v1'
PREP=OUT/'preparation'
ACTUAL=OUT/'actual'
METHOD=OUT/'METHOD.md'
RANK=611
PRECISION=640
RADIUS=Q(1,2**34)
GATE_STATUS='PASS_INDEPENDENT_REMAINING32_PRECONDITIONER_SPD_PRELAUNCH'
SCOPE=('Exact SPD premise for the symmetric part of each of 32 fixed binary preconditioners. '
       'No change to the archived fields or exact physical inputs. No spatial force, Hessian, root, '
       'minimum, continuum, observer, or loading-path claim.')
PINS={
    ROOT/'scripts/cohort_energy_stability.py':'a3fe42c9d79907103c0e02af35db3bdc4c2f89aee59755c006e9ab09e322ee7d',
    ROOT/'scripts/review_cohort_energy_stability.py':'4c2acbc96ce3f19954b55d14fde46a8532536e2d39c2250e0123a17b78663c73',
    ROOT/'scripts/algebraic_root_feasibility.py':'7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf',
    ROOT/'scripts/review_algebraic_root.py':'6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a',
    ROOT/'scripts/transfer_ambiguity_remaining32_exact_targets.py':'d1a9c9e4839ff85a9926b6fd5ecc11e4cd8e0831c7e6263ecfbe068ef22ea3f4',
    TARGET/'preparation/protocol.json':'a12669fd1e50cbd7435d016ba60223791c8a70119f622048c28e24ef5a9f456b',
    TARGET/'preparation/preflight.json':'86ddb06c75c2deedea29d62bb7e03da6142d55b41a0291e1567c9e4dba6aab4f',
    TARGET/'preparation/manifest.json':'dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242',
    TARGET/'independent_prelaunch/review.json':'67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9',
    GEOMETRY/'preparation/protocol.json':'754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850',
    GEOMETRY/'independent_actual_preparation/actual_replay/receipt.json':'74179dde74571f792b6f3949a06d6e26076b8c746e6700c355264d0de3debf04',
    GEOMETRY/'independent_actual_preparation/actual_replay/os_exit_observation.json':'d39b2d5e97ad6d15bfb68f3f8360e0214c45d11b37dbbd3eaae48313a2e8fc0c',
    JOIN/'bridge_input_admission.json':'2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb',
    JOIN/'os_exit_observation.json':'636d312911742e22722377ee540e13ed725c116de6d1853f137349d03b9dee96',
    EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json':'769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
}
CANDIDATE='CANDIDATE_EXACT_SYMMETRIC_PART_SPD_AWAITING_INDEPENDENT_REVIEW'
UNRESOLVED='UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_TEST'


def require(value,message):
    if not value:
        raise ValueError(message)


def rel(path):
    return path.relative_to(ROOT).as_posix()


def canon(name):
    check=(ROOT/name).resolve()
    require(check.is_relative_to(ROOT) and rel(check)==name,'canonical study-relative authority required')
    return check


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def checked(path,expected):
    require(type(expected) is str and re.fullmatch('[0-9a-f]{64}',expected) is not None and sha(path)==expected,
            f'changed authoritative source or exact data: {path.name}')


def read(path,expected):
    checked(path,expected)
    return json.loads(path.read_bytes())


def read_pre(path):
    return read(path,PINS[path])


def save(path,value):
    with path.open('x',encoding='utf8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')
    return sha(path)


def yes(checks):
    return type(checks) is dict and bool(checks) and all(v is True for v in checks.values())


def iszero(x):
    return type(x) is int and x==0


def exact(value):
    require(type(value) in (str,int,Q),'exact finite rational without Boolean/float alias required')
    return Q(value)


def digest(value):
    return sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def boot():
    for module in (adapter,certifier,endpoint,core,maps):
        path=Path(module.__file__).resolve()
        checked(path,PINS[path])
    adapter.X.bootstrap()
    core.A.cache_clear()
    maps.A.cache_clear()
    require(certifier.ctx is endpoint.ctx is core.ctx is maps.ctx is arb_context
            and arb_context.prec==PRECISION,'explicit matched Arb640 and exact source identities required')


def array_hash(a,rank=RANK):
    require(type(a) is np.ndarray and type(rank) is int and rank>0 and a.shape==(rank,rank)
            and a.dtype.str=='<f8' and np.isfinite(a).all(),
            'finite declared square binary64 matrix required')
    return sha256(a.tobytes(order='C')).hexdigest()


def pinned_link(links,name,expected,parse=False):
    require(links.setdefault(name,expected)==expected,'conflicting immutable input pin')
    path=canon(name)
    return read(path,expected) if parse else checked(path,expected)


def identity(entry,joined,geo,contract):
    pid=entry['case_id']; d=contract['exact_requested_drive']; norm=exact(joined['R_inf_exact'])
    dose,gain,drive=(exact(d[k]) for k in ('optical_dose_exact','common_gain_exact','activation_strength_exact'))
    require(pid==joined['case_id']==geo['pilot_id']==contract['pilot_id']
            and entry['record_id']==joined['record_id']==geo['record_id']==contract['record_id']
            and entry['model_contract_sha256']==joined['model_contract_sha256']==geo['model_contract_sha256']
            and joined['trial_box_sha256']==geo['box_sha256'] and exact(geo['radius_exact'])==RADIUS
            and joined['cache_packet_sha256']==geo['point_cache_sha256']==contract['point_cache_sha256']
            and geo['R_binary_sha256']==contract['binary_R_sha256']
            and geo['center_binary_sha256']==contract['center_binary_sha256']
            and geo['archived_field_sha256']==contract['archived_field_sha256']
            and geo['archived_shift_exact']==contract['archived_shift_exact']
            and geo['archived_addition_gap_exact']==contract['archived_addition_gap_exact']
            and exact(geo['energy_scale_binary_exact'])==exact(contract['energy_scale_binary_exact'])>0
            and 0<drive==dose*gain==exact(geo['exact_drive']) and drive!=exact(d['activation_strength_binary_exact'])
            and all(exact(a)==exact(b) for a,b in contract['parameter_cell'].values())
            and tuple(map(exact,contract['parameter_cell']['common_gain']))==(gain,gain)
            and digest(contract['parameter_cell'])==contract['parameter_cell_sha256']
            and norm>0 and 0<=exact(joined['reviewed_ZQ_upper_exact'])<1 and geo['root_asserted'] is False,
            'scope or same preconditioner/cache/observer/gain/normalization identity drift')
    return {'pilot_id':pid,'record_id':entry['record_id'],'field_sha256':contract['archived_field_sha256'],
            'cache_packet_sha256':contract['point_cache_sha256'],'embedded_cache_sha256':contract['frozen_binary_cache_sha256'],
            'R_binary_sha256':contract['binary_R_sha256'],'center_binary_sha256':contract['center_binary_sha256'],
            'model_contract_sha256':entry['model_contract_sha256'],'parameter_cell_sha256':contract['parameter_cell_sha256'],
            'dose_exact':str(dose),'gain_exact':str(gain),'drive_exact':str(drive),
            'energy_scale_binary_exact':str(exact(contract['energy_scale_binary_exact']))}


def authorities():
    links={rel(path):h for path,h in PINS.items()}
    for path,h in PINS.items():
        checked(path,h)
    links[rel(SOURCE)]=sha(SOURCE)
    links[rel(METHOD)]=sha(METHOD)
    my_target=TARGET/'preparation/manifest.json'; m=read_pre(my_target)
    target_proto=read_pre(TARGET/'preparation/protocol.json')
    tg=read_pre(TARGET/'independent_prelaunch/review.json')
    verified=read_pre(GEOMETRY/'independent_actual_preparation/actual_replay/receipt.json')
    geom_os=read_pre(GEOMETRY/'independent_actual_preparation/actual_replay/os_exit_observation.json')
    joined=read_pre(JOIN/'bridge_input_admission.json')
    joined_os=read_pre(JOIN/'os_exit_observation.json')
    geom_plan=read_pre(GEOMETRY/'preparation/protocol.json')
    require(tg['status']=='PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS' and yes(tg['checks'])
            and tg['manifest_sha256']==PINS[my_target]==joined['target_manifest_sha256']
            and m['geometry_or_variation_actual_launch_authorized'] is False
            and verified['status']=='PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY'
            and verified['successfully_replayed_geometry_cases']==32 and verified['exact_spatial_bridge_authorized'] is False
            and geom_os['status']=='PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_REPLAY_EXIT' and iszero(geom_os['observed_exit_code'])
            and geom_os['receipt_sha256']==PINS[GEOMETRY/'independent_actual_preparation/actual_replay/receipt.json']==joined['geometry_receipt_sha256']
            and joined['geometry_observed_exit_sha256']==PINS[GEOMETRY/'independent_actual_preparation/actual_replay/os_exit_observation.json']
            and joined['status']=='PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS' and yes(joined['checks'])
            and joined['spatial_workers_authorized'] is False
            and joined_os['status']=='PASS_OBSERVED_ALL32_JOINED_INPUT_READBACK_TERMINAL' and iszero(joined_os['observed_exit_code'])
            and joined_os['joined_admission_sha256']==PINS[JOIN/'bridge_input_admission.json'],
            'same independent 32-cache and whole-box evidence required')
    for group in (target_proto['local_source_closure_sha256'],joined['source_and_input_sha256']):
        for name,h in group.items():
            pinned_link(links,name,h)
    order=[x['case_id'] for x in m['targets']]
    require(len(set(order))==len(order)==32 and tuple(order)==adapter.IDS==tuple(tg['ordered_case_ids'])
            ==tuple(x['case_id'] for x in geom_plan['trials'])==tuple(x['case_id'] for x in joined['cases'])
            ==tuple(x['pilot_id'] for x in verified['geometry_records']), 'canonical full 32-case order required')
    rows=[]
    for entry,j,g,trial in zip(m['targets'],joined['cases'],verified['geometry_records'],geom_plan['trials'],strict=True):
        c=pinned_link(links,entry['model_contract'],entry['model_contract_sha256'],parse=True)
        geo=pinned_link(links,g['receipt'],g['sha256'],parse=True)
        require(geo['status']==g['status']=='PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY'
                and g['sha256']==j['geometry_replay_case_sha256']
                and j['trial_box_sha256']==trial['trial_box_sha256'], 'unmatched whole-box geometry inputs')
        phase=identity(entry,j,geo,c)
        rows.append({'pilot_id':entry['case_id'],'identity':phase,'R_source_file':c['point_cache'],
                     'R_source_file_sha256':c['point_cache_sha256'],'trial_box_sha256':trial['trial_box_sha256'],
                     'geometry_replay_case_sha256':g['sha256'],'cached_derivative_case_sha256':j['cached_derivative_case_sha256'],
                     'reviewed_R_norm_exact':str(exact(j['R_inf_exact'])),
                     'archived_shift_exact':c['archived_shift_exact'],'archived_addition_gap_exact':c['archived_addition_gap_exact']})
    return rows,dict(sorted(links.items()))


def protocol(rows,links):
    return {'status':'FROZEN_REMAINING32_EXACT_SYMMETRIC_PRECONDITIONER_PART_PROPOSALS',
            'source_sha256':sha(SOURCE),'target_manifest_sha256':PINS[TARGET/'preparation/manifest.json'],
            'target_context_review_sha256':PINS[TARGET/'independent_prelaunch/review.json'],
            'joined_bridge_input_sha256':PINS[JOIN/'bridge_input_admission.json'],
            'source_and_input_sha256':links,'ordered_case_ids':[r['pilot_id'] for r in rows],'ordered_targets':rows,
            'rank':RANK,'precision_bits':PRECISION,'single_serial_process':True,
            'backend':{key:importlib.metadata.version(key) for key in ('numpy','scipy','python-flint')},
            'python_version':sys.version.split()[0],
            'target':'P=(R+R^T)/2 from exact frozen binary64 entries, no rounded P as target',
            'criterion':'Both outward q upper bounds strictly less than 1. Floating inverse Cholesky supplies K only.',
            'root_or_minimum_claimed':False,'separate_independent_actual_replay_required':True,
            'matrix_stage_launch_requires_independent_cold_review':True,
            'new_production_matrix_jobs_in_preparation':0,'scope':SCOPE}


def check_phase(review,pf,p,binding):
    require(review['status']==GATE_STATUS and yes(review['checks']) and all(review[k]==h for k,h in binding.items())
            and review['permitted_case_ids']==p['ordered_case_ids'] and review['root_or_minimum_claimed'] is False
            and pf['status']=='PASS' and yes(pf['checks']) and iszero(pf['new_production_matrix_jobs'])
            and all(pf[k]==binding[k] for k in ('source_sha256','protocol_sha256'))
            and all(p[k]==binding[k] for k in ('source_sha256','target_manifest_sha256','target_context_review_sha256',
                                             'joined_bridge_input_sha256')), 'same external cold SPD phase review required')


def args_binding(p,protocol_sha,preflight_sha):
    return {'source_sha256':sha(SOURCE),'protocol_sha256':protocol_sha,'preflight_sha256':preflight_sha,
            **{k:p[k] for k in ('target_manifest_sha256','target_context_review_sha256','joined_bridge_input_sha256')}}


def require_identity(expected,observed):
    require(expected==observed,'same fixed target canonical field, center, drive, cache and R required')


def context_matrix(row):
    boot()
    pid=row['pilot_id']
    ctx,_,e=adapter.open_context(pid,PINS[TARGET/'preparation/manifest.json'])
    R=ctx.R
    contract=read(canon(e['model_contract']),e['model_contract_sha256'])
    require(arb_context.prec==PRECISION and R.flags.writeable is False and ctx.u.shape==ctx.free.shape==(RANK,)
            and ctx.target.packet['input_state_sha256']==contract['archived_field_sha256']
            and ctx.row['point_cache_file']==row['R_source_file']==contract['point_cache']
            and ctx.row['point_cache_sha256']==row['R_source_file_sha256']
            and contract['archived_shift_exact']==row['archived_shift_exact']
            and contract['archived_addition_gap_exact']==row['archived_addition_gap_exact'],
            'immutable exact R in the originally archived scalar model required')
    observed={'pilot_id':ctx.pilot_id,'record_id':ctx.record_id,'field_sha256':ctx.target.packet['input_state_sha256'],
              'cache_packet_sha256':ctx.row['point_cache_sha256'],'embedded_cache_sha256':ctx.row['frozen_binary_cache_sha256'],
              'R_binary_sha256':array_hash(R),'center_binary_sha256':sha256(ctx.u.tobytes(order='C')).hexdigest(),
              'model_contract_sha256':e['model_contract_sha256'],'parameter_cell_sha256':ctx.row['parameter_cell_sha256'],
              'dose_exact':str(ctx.dose),'gain_exact':str(ctx.base['common_gain']),'drive_exact':str(ctx.drive),
              'energy_scale_binary_exact':str(ctx.target.energy_scale)}
    require_identity(row['identity'],observed)
    return R


def matrix_candidate(R,row,folder,rank=RANK):
    rsha=array_hash(R,rank)
    require(row['identity']['R_binary_sha256']==rsha,'producer R bytes disagree with target receipt')
    result,K=certifier.certify_sympart(R,precision_bits=PRECISION)
    n=endpoint.norm_R_integer(R)
    require(array_hash(R,rank)==rsha and result.get('R_binary_sha256',rsha)==rsha and arb_context.prec==PRECISION
            and n==exact(row['reviewed_R_norm_exact'])>0,'matrix, precision or exact R row norm changed')
    second=None; proposal={}; accepted=False
    if K is not None:
        ksha=array_hash(K,rank)
        kp=folder/'congruence_K.npy'
        with kp.open('xb') as f:
            np.save(f,K,allow_pickle=False)
        kh=sha(kp); frozen=np.load(kp,allow_pickle=False)
        require(array_hash(frozen,rank)==ksha,'K proposal exact entry byte serialization changed')
        q,rows=endpoint.congruence_bounds(R,frozen,PRECISION)
        original=list(map(exact,result['row_absolute_sum_upper_exact']))
        upper=exact(result['q_infinity_upper_exact'])
        require(len(original)==len(rows)==rank and min(original)>=0 and min(rows)>=0
                and upper==max(original) and q==max(rows) and array_hash(R,rank)==rsha and sha(kp)==kh,
                'complete separate outward congruence bounds required')
        accepted=result['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD' and 0<=upper<1 and 0<=q<1
        second={'association':'K*(P*K^T)','q_infinity_upper_exact':str(q),'row_absolute_sum_upper_exact':list(map(str,rows)),
                'exact_R_infinity_norm':str(n),'source_sha256':sha(Path(endpoint.__file__))}
        proposal={'K_file':rel(kp),'K_file_sha256':kh,'K_binary_sha256':ksha}
    return {'status':CANDIDATE if accepted else UNRESOLVED,'target':row,'precision_bits':PRECISION,
            'exact_R_infinity_norm':str(n),'producer_congruence':result,'second_association_readback':second,**proposal,
            'external_independent_actual_review_required':True,'preconditioner_changed':False,
            'new_force_Hessian_AD_integration_or_solve_jobs':0,'root_or_minimum_claimed':False,'scope':SCOPE}


def preflight(ph,checks):
    return {'status':'PASS','source_sha256':sha(SOURCE),'protocol_sha256':ph,'checks':checks,
            'new_production_matrix_jobs':0,'binary_case_R_arrays_decoded':0,
            'new_force_Hessian_AD_integration_or_solve_jobs':0}


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def tiny_fixtures(p,ph):
    boot()
    I=np.eye(2,dtype='<f8')
    asymmetric=np.array([[2.,5.],[-5.,2.]],dtype='<f8')
    indefinite=np.array([[1.,4.],[0.,1.]],dtype='<f8')
    good,_=certifier.certify_sympart(asymmetric,K=I/2,precision_bits=PRECISION)
    bad,_=certifier.certify_sympart(indefinite,K=I,precision_bits=PRECISION)
    boundary=np.diag([1.,0.])
    at_one,_=certifier.certify_sympart(I,K=boundary,precision_bits=PRECISION)
    q,_=endpoint.congruence_bounds(asymmetric,I/2,PRECISION)
    bq,_=endpoint.congruence_bounds(indefinite,I,PRECISION)
    sq,_=endpoint.congruence_bounds(I,boundary,PRECISION)
    tricky=np.array([[2.,1.],[2.**-53,2.]],dtype='<f8')
    P=endpoint.exact_P(tricky)
    exact_entry=Q(1,2)+Q(1,2**54)
    rounded=Q.from_float(float(((tricky+tricky.T)/2)[0,1]))
    require(P[0,1].contains(endpoint.A(exact_entry)) and not P[0,1].contains(endpoint.A(rounded)),
            'exact binary Sym(R) differentiates from the floating proposal')
    require(good['status']=='CERTIFIED_EXACT_SYMMETRIC_PART_SPD' and exact(good['q_infinity_upper_exact'])==q==Q(1,2)
            and exact(bad['q_infinity_upper_exact'])==bq==2 and bad['status']=='UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_BOUND'
            and exact(at_one['q_infinity_upper_exact'])==sq==1 and at_one['status']=='UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_BOUND'
            and endpoint.norm_R_integer(asymmetric)==7,'strict SPD criterion or independent row norm changed')
    tiny={'identity':{'R_binary_sha256':array_hash(asymmetric,2)},'reviewed_R_norm_exact':'7'}
    with TemporaryDirectory(dir=OUT) as d:
        d=Path(d).resolve()
        require(d.is_relative_to(OUT),'scoped synthetic K archive required')
        t=matrix_candidate(asymmetric,tiny,d,2)
        require(t['status']==CANDIDATE and sha(d/'congruence_K.npy')==t['K_file_sha256']
                and array_hash(np.load(d/'congruence_K.npy',allow_pickle=False),2)==t['K_binary_sha256']
                and t['root_or_minimum_claimed'] is False,'synthetic exact K archive roundtrip failed')
    for bad_data in (I.astype(bool),np.array([[np.nan,0.],[0.,1.]],dtype='<f8'),I[:1],I.astype('>f8')):
        require(reject(lambda:array_hash(bad_data,2)),'wrong type/shape/nonfinite matrix accepted')
    row=p['ordered_targets'][0]['identity']
    for field in ('field_sha256','cache_packet_sha256','R_binary_sha256','drive_exact','center_binary_sha256'):
        require(reject(lambda:require_identity(row,{**row,field:'wrong'})),'mixed same-target identity accepted')
    dummy_h='b'*64
    binding=args_binding(p,ph,dummy_h)
    pf=preflight(ph,{'fixture':True})
    valid={'status':GATE_STATUS,**binding,'checks':{'fixture':True},'permitted_case_ids':p['ordered_case_ids'],
           'root_or_minimum_claimed':False}
    check_phase(valid,pf,p,binding)
    tests=({**valid,'root_or_minimum_claimed':True},{**valid,'checks':{'truthy':1}},
           {**valid,'permitted_case_ids':p['ordered_case_ids'][1:]},{**valid,'preflight_sha256':'a'*64})
    require(all(reject(lambda:check_phase(test,pf,p,binding)) for test in tests)
            and reject(lambda:check_phase(valid,{**pf,'new_production_matrix_jobs':False},p,binding)),
            'mixed external SPD cold stage flags or preflight self-link accepted')
    return {'exact_binary_sympart_distinguishes_rounded_float_proposal':True,
            'nonsymmetric_R_with_positive_exact_symmetric_part_certificates':True,
            'indefinite_part_rejected_despite_positive_R_eigenvalues':True,
            'singular_K_and_exact_defect_boundary_one_rejected':True,
            'exact_integer_R_row_norm_matches_synthetic_oracle':True,
            'frozen_K_binary_and_file_hash_roundtrip':True,
            'wrong_type_shape_nonfinite_or_mixed_target_rejected':True,
            'positive_external_gate_roundtrip_and_literal_negative_fixtures':True,
            'all_32_same_immutable_exact_metadata_inputs_preflighted':True}


def prepare():
    require(not PREP.exists() and not ACTUAL.exists(),'cold artifacts and actual attempts are write-once')
    with patch.object(adapter,'open_context',side_effect=AssertionError('cold phase cannot decode any actual R or u')):
        rows,links=authorities()
        p=protocol(rows,links)
        ph=sha256((json.dumps(p,indent=2,allow_nan=False)+'\n').encode('utf-8')).hexdigest()
        checks=tiny_fixtures(p,ph)
    PREP.mkdir(exist_ok=False)
    require(save(PREP/'protocol.json',p)==ph,'serialized protocol hash changed')
    pf=save(PREP/'preflight.json',preflight(ph,checks))
    with (PREP/'producer_frozen.py').open('xb') as f:
        f.write(SOURCE.read_bytes())
    print(json.dumps({'source_sha256':sha(SOURCE),'protocol_sha256':ph,'preflight_sha256':pf,
                      'cold_targets':len(rows),'new_production_matrix_jobs':0}),flush=True)


def gate(args):
    p=read(PREP/'protocol.json',args.protocol_sha256)
    pf=read(PREP/'preflight.json',args.preflight_sha256)
    rpath=args.review.resolve()
    require(rpath==PREP/'independent_prelaunch/review.json','canonical separate matrix source review required')
    reviewer=read(rpath,args.review_sha256)
    binding=args_binding(p,args.protocol_sha256,args.preflight_sha256)
    check_phase(reviewer,pf,p,binding)
    require(sha(PREP/'producer_frozen.py')==binding['source_sha256'],'frozen matrix source drift')
    rows,links=authorities()
    require(p==protocol(rows,links),'current exact model, source or numeric backend changed')
    return rows,links,{**binding,'independent_prelaunch_review_sha256':args.review_sha256}


def run(args):
    rows,links,binding=gate(args)
    ACTUAL.mkdir(exist_ok=False)
    start=perf_counter()
    mark={'pid':os.getpid(),'bindings':binding,'ordered_case_ids':[r['pilot_id'] for r in rows]}
    save(ACTUAL/'started.json',{'status':'RUNNING',**mark})
    terminal={'status':'FAILED','exit_code':1,**mark}
    results=[]
    try:
        with ExitStack() as guards:
            guards.enter_context(adapter.X.no_integrals())
            for obj,k in ((core,'variation_rows'),(core,'bulk_local'),(core,'surface_local'),(maps,'terms_from_arrays')):
                guards.enter_context(patch.object(obj,k,side_effect=AssertionError('no physical force or derivative jobs in matrix stage')))
            for row in rows:
                pid=row['pilot_id']
                folder=ACTUAL/pid
                folder.mkdir(exist_ok=False)
                R=context_matrix(row)
                result={**matrix_candidate(R,row,folder),'case_id':pid,'pid':mark['pid'],'bindings':binding}
                path=folder/'result.json'
                rh=save(path,result)
                results.append({'case_id':pid,'file':rel(path),'sha256':rh,'status':result['status'],
                                'R_binary_sha256':row['identity']['R_binary_sha256'],
                                'K_file_sha256':result.get('K_file_sha256')})
                print(json.dumps({'case_id':pid,'status':result['status'],'result_sha256':rh}),flush=True)
        for path,h in links.items():
            checked(canon(path),h)
        candidates=sum(r['status']==CANDIDATE for r in results)
        require(len(results)==32 and [r['case_id'] for r in results]==mark['ordered_case_ids'],'complete explicit matrix stage required')
        result={'status':'COMPLETE_REMAINING32_SYMMETRIC_PART_ATTEMPTS_AWAITING_INDEPENDENT_REPLAY',**mark,
                'cases':results,'candidate_spd_count':candidates,'unresolved_count':len(results)-candidates,
                'externally_observed_exit_required':True,'external_independent_actual_review_required':True,
                'root_or_minimum_claimed':False,'scope':SCOPE,'seconds':perf_counter()-start}
        mh=save(ACTUAL/'manifest.json',result)
        terminal.update(status='FINISHED',exit_code=0,manifest_sha256=mh)
        print(json.dumps({'manifest_sha256':mh,'candidates':candidates,'unresolved':len(results)-candidates}),flush=True)
    except BaseException as error:
        save(ACTUAL/'failure.json',{'error':str(error),'traceback':traceback.format_exc(),
                                    'completed_case_results':results,**mark})
        raise
    finally:
        save(ACTUAL/'terminal.json',terminal)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest='mode',required=True)
    modes.add_parser('prepare')
    actual=modes.add_parser('run')
    actual.add_argument('--review',type=Path,required=True)
    for k in ('review-sha256','protocol-sha256','preflight-sha256'):
        actual.add_argument('--'+k,required=True)
    args=parser.parse_args()
    if args.mode=='prepare':
        prepare()
    else:
        run(args)


if __name__=='__main__':
    main()
