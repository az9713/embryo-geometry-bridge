"""Cold review of the first-wave independent signed endpoint replay."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import combinations_with_replacement
import json
from pathlib import Path
from random import Random
from types import SimpleNamespace
from unittest.mock import patch


HERE=Path(__file__).resolve().parent
PREP=HERE.parent
SOURCE=PREP/'replay_global.py'
PINS={'source_sha256':'4682165bf59bd633092c042a7df619c9bc1ad272dc80d1755404c2999a65287f',
      'protocol_sha256':'9ca811a3e7f9e35346dd28d040ff3f370599e014ce3df942c48d347f6bd2f24a',
      'preflight_sha256':'90242d2cde1a81c075cf342f716cf8a229f45b5ce3f73778fdd9db6ade249d91'}
STATUS='PASS_INDEPENDENT_COLD_FIRST_WAVE_GLOBAL_ENDPOINT_SOURCE_REVIEW'


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def finite_images(module):
    rng=Random(4016)
    def x():
        return Q(rng.randint(-12,12),rng.choice([1,2,3,5,7,11]))
    def interval():
        return tuple(sorted((x(),x())))
    def interval_product_dot(row,values):
        products=[tuple(sorted((a*l,a*u))) for a,(l,u) in zip(row,values,strict=True)]
        return sum(t[0] for t in products),sum(t[1] for t in products)
    cases=70
    for sample in range(cases):
        n=sample%4+1
        R=[[x() for _ in range(n)] for _ in range(n)]
        F=[interval() for _ in range(n)]
        H={(i,j):interval() for i,j in combinations_with_replacement(range(n),2) if (i+j+sample)%3!=0}
        exact_f=[interval_product_dot(row,F) for row in R]
        exact_ht=[[interval_product_dot(row,[H.get((min(k,j),max(k,j)),(0,0)) for k in range(n)])
                   for j in range(n)] for row in R]
        hnorm=[sum(max(abs(x),abs(y)) for x,y in row) for row in exact_ht]
        result_f,result_h=module.bound(R,F,H)
        require(all(a<=c<=d<=b and max(c-a,b-d)<Q(1,2**700)
                    for (a,b),(c,d) in zip(result_f,exact_f,strict=True))
                and all(b<=a and a-b<Q(1,2**700) for a,b in zip(result_h,hnorm,strict=True)),
                'directed positive/negative split misses an exact-rational endpoint or is unexpectedly wide')
    return cases


def scatter_and_logic(m):
    zero=(Q(0),Q(0)); force=[zero]*3; H={}
    # Two actual-order maps use both (1,2) and (2,1). Each local upper
    # off-diagonal contributes exactly once to the canonical global pair.
    m.scatter_add(force,H,[2,1],[(3,4),(-5,-3)],[(1,2),(-2,3),(7,8)])
    m.scatter_add(force,H,[1,2],[(0,2),(-2,1)],[(2,3),(-4,-2),(-2,0)])
    require(force==[zero,(-5,-1),(1,5)] and H=={(2,2):(-1,2),(1,2):(-6,1),(1,1):(9,11)},
            'global diagonal or one-time symmetry mirror/scatter misrouting')
    require(reject(lambda:m.scatter_add(force,H,[False],[(0,0)],[(0,0)]))
            and reject(lambda:m.scatter_add(force,H,[2,2],[(0,0)]*2,[(0,0)]*3)),
            'typed distinct global free projection required')
    f=[(-Q(1,8),-Q(1,8))]; anchor=[(Q(1,4),Q(1,4))]
    s=m.classify(f,{(0,0):(Q(1,8),Q(1,8))},anchor,[[2]],Q(1,4),m.S.RADIUS)
    require(s['signed_R_F_I_rows_exact']==[['0','0']] and s['z_upper_exact']=='1/2'
            and s['strict_selfmap'] is True and s['strict_cube_exclusion'] is False,
            'already preconditioned cached anchor must cancel the signed defect once')
    target={'archived_center_shift_exact':'1/64','original_addition_gap_exact':'1/128','observer_goal_exact':'1/64'}
    with_spd=m.observer_minimum(s,target|{'observer_goal_exact':'1/32'},True,Q(2))
    no_spd=m.observer_minimum(s,target|{'observer_goal_exact':'1/32'},False,Q(2))
    require(with_spd['strict_minimum_proved'] is True and with_spd['normalized_H_euclidean_lower_exact']=='1/4'
            and with_spd['physical_nodal_error_um_upper_exact']=='201/256'
            and no_spd['strict_minimum_proved'] is False and no_spd['normalized_H_euclidean_lower_exact'] is None,
            'shift plus posterior plus original xi, or same-R SPD dependency, lost')
    rr=m.S.RADIUS
    tied=m.classify([(rr,rr)],{},[(0,0)],[[1]],0,rr)
    excluded=m.classify([(2*rr,2*rr)],{},[(0,0)],[[1]],0,rr)
    require(tied['strict_selfmap'] is tied['strict_cube_exclusion'] is False
            and excluded['strict_selfmap'] is False and excluded['strict_cube_exclusion'] is True,
            'strict self-map or only-cube reverse triangle criterion mislabeled')
    src={'index':[0,1],'independent_exact':['-3','1'],'producer_exact':['0','2'],'hull_exact':['-3','2'],
         'producer_contains_independent':False,'disjoint':False}
    def signed(row):
        return {'signed_channels':{'H':{'rows':[row],'producer_containment_misses':[[0,1]],'disjoint_indices':[]}}}
    require(m.validated_channel(signed(src),'H',[[0,1]],[['0','2']])==([(Q(-3),Q(2))],1)
            and all(reject(lambda:m.validated_channel(signed({**src,**v}),'H',[[0,1]],[['0','2']])) for v in
                    ({'hull_exact':['0','1']},{'index':[False,1]},{'independent_exact':['-3','-2'],'disjoint':True})),
            'nonnested certified union cannot be replaced by overlap or mismatched local index')


def root_gate_roundtrip(m):
    class StopBeforeActualProducerAdmission(Exception):
        pass

    root_h='a'*64; external_h='b'*64; source_review_sha='c'*64
    args=SimpleNamespace(protocol_sha256=PINS['protocol_sha256'],preflight_sha256=PINS['preflight_sha256'],
                         review=m.SOURCE_REVIEW,review_sha256=source_review_sha,
                         producer_observed_sha256=root_h,producer_admission_sha256=external_h)
    valid={'status':m.REVIEW,**PINS,'producer_root_terminal_admission_sha256':root_h,
           'producer_reduction_input_admission_sha256':external_h,
           'checks':{'synthetic_phase_roundtrip':True},'new_spatial_integrations_authorized':False}
    choice=[valid]; real_read=m.read
    def read(path,h,*a,**k):
        if path==m.SOURCE_REVIEW:
            require(h==source_review_sha,'wrong source authorizer')
            return choice[0]
        return real_read(path,h,*a,**k)
    with patch.object(m,'read',side_effect=read), \
         patch.object(m.D,'gate',side_effect=StopBeforeActualProducerAdmission) as stopped, \
         patch.object(m,'observed_producer',side_effect=AssertionError('cold gate read actual output')):
        try:
            m.run(args)
        except StopBeforeActualProducerAdmission:
            pass
        else:
            raise AssertionError('frozen source gate did not reach the separate producer admission seam')
        require(stopped.call_count==1,'frozen source gate did not bind real preflight and source first')
        for bad in ({**valid,'preflight_sha256':'d'*64},{**valid,'checks':{'truthy':1}},
                    {**valid,'producer_root_terminal_admission_sha256':'e'*64},
                    {**valid,'producer_reduction_input_admission_sha256':'f'*64},
                    {**valid,'new_spatial_integrations_authorized':True}):
            choice[0]=bad
            require(reject(lambda:m.run(args)) and stopped.call_count==1,
                    'global replay crossed its own real source gate with synthetic malformed authority')


def native_os_fixture(m):
    p=m.D.ACTUAL; order=list(m.W.FIRST); b={'pilot':'test'}
    hashes={m.rel(p/(k+'.json')):'0'*64 for k in ('started','terminal','results')}
    hashes.update({m.rel(p/c/(k+'.json')):'0'*64 for c in order for k in ('assembled_operands','results')})
    obs={'status':m.OS_STATUS,'pid':123,'unified_exec_session':78,'observed_exit_code':0,
         'producer_reduction_input_admission_sha256':'f'*64,'sha256':hashes}
    comm={'pid':123,'bindings':b,'case_order':order}
    terminal={**comm,'status':'FINISHED','exit_code':0,'results_sha256':'0'*64}
    combined={**comm,'status':'COMPLETE_FIRST_WAVE_SIGNED_REDUCTION_AWAITS_INDEPENDENT_GLOBAL_REVIEW',
              'root_asserted':False,'full_five_gate_complete':False,'cases':[{'case_id':c} for c in order]}
    documents={m.PRODUCER_OS:obs,p/'started.json':{**comm,'status':'RUNNING'},
               p/'terminal.json':terminal,p/'results.json':combined}
    def read(path,h):
        require(path in documents,'premature local or numeric payload open')
        return documents[path]
    with patch.object(m,'read',side_effect=read):
        a,top,_=m.observed_producer('a'*64,'f'*64,b)
        require(a==obs and top==combined,'11-pin synthetic admission fixture not executable')
        for bad in ({**obs,'observed_exit_code':False},{**obs,'unified_exec_session':78.0},
                    {**obs,'sha256':dict(list(hashes.items())[:-1])}):
            documents[m.PRODUCER_OS]=bad
            require(reject(lambda:m.observed_producer('a'*64,'f'*64,b)), 'foreign native process or incomplete file set accepted')
        documents[m.PRODUCER_OS]=obs
        documents[p/'terminal.json']={**terminal,'pid':123.0}
        require(reject(lambda:m.observed_producer('a'*64,'f'*64,b)),'parent float PID alias accepted')


def exact_spd_context(m,a):
    required=m.PINS[m.D.SPD/'actual_replay/os_exit_observation.json']
    replay=m.D.optional_spd({'observed_os_exit_sha256':required,'replay_phase':{
        'source_sha256':'28c511fc6087d130f353f67c3cb655aa729031b557ae0bd84c9714b632435a13',
        'protocol_sha256':'a740f6a3305421de51a8fa7ce1e5da94bf308ae35ab7226c0ef3336e1beed8c7',
        'preflight_sha256':'c8150950f934719f6212b579bcac066573c94b09ac5943cd7f4c9b5d489339d7'}})
    authority=m.read(m.V.A.PREP/'manifest.json','dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242')
    checked=[]
    for item in a['plans']:
        case=item['case_id']; t=item['target']
        entry=next(x for x in authority['targets'] if x['case_id']==case)
        c=m.read(m.V.A.canonical(entry['model_contract']),entry['model_contract_sha256'])
        p=m.read(m.S.PLAN/case/'plan.json',item['plan_sha256'])
        require(t['model_contract_sha256']==entry['model_contract_sha256']
                and t['trial_box_sha256']==p['trial_box_sha256'] and Q(p['whole_box_radius_exact'])==m.S.RADIUS
                and len(p['signed_cached_RF_center_rows_exact'])==611
                and all(len(m.pair(x))==2 for x in p['signed_cached_RF_center_rows_exact'])
                and 0<=Q(p['sizing']['ZQ_exact'])<1 and Q(p['sizing']['R_inf_exact'])>0,
                'not the same exact planned r34 target and already preconditioned cache vector')
        ctx=SimpleNamespace(row={'frozen_binary_cache_sha256':c['frozen_binary_cache_sha256'],
                                 'parameter_cell_sha256':c['parameter_cell_sha256'],'point_cache_file':c['point_cache']},
                            dose=Q(c['exact_requested_drive']['optical_dose_exact']),
                            base={'common_gain':Q(c['exact_requested_drive']['common_gain_exact'])})
        norm=Q(p['sizing']['R_inf_exact'])
        q=m.own_spd(replay,t,ctx,p,norm)
        require(q['exact_Sym_R_native_receipt_sha256']=='de31a30f3270829c7fd89ebc75db58e45d2d9a911cd02dfc6e564493748d8615'
                and m.own_spd(None,t,ctx,p,norm) is None
                and reject(lambda:m.own_spd(replay,{**t,'R_binary_sha256':'0'*64},ctx,p,norm))
                and reject(lambda:m.own_spd(replay,{**t,'archived_center_shift_exact':'1'},ctx,p,norm))
                and reject(lambda:m.own_spd(replay,t,ctx,p,norm+1)),
                'same binary R, shift or exact norm necessary even with valid SPD of another target')
        checked.append({'case_id':case,'plan_sha256':item['plan_sha256'],'trial_box_sha256':t['trial_box_sha256'],**q})
    return checked


def main():
    require(sha(SOURCE)==sha(PREP/'replay_frozen.py')==PINS['source_sha256']
            and sha(PREP/'protocol.json')==PINS['protocol_sha256'] and sha(PREP/'preflight.json')==PINS['preflight_sha256'],
            'exact frozen local, producer, and review phases required')
    spec=importlib.util.spec_from_file_location('independent_global_cold_subject',SOURCE)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    require(not module.OUT.exists() and not module.D.ACTUAL.exists(),'this review admits no live actual result')
    protocol=module.read(PREP/'protocol.json',PINS['protocol_sha256'])
    require(protocol==module.static() and all(v is True for v in module.fixtures().values()),
            'authentic frozen upstream phases and prior tiny exercises required')
    raw,runtime=module.D.static_phase()
    require(runtime==protocol['worker_binding'] and tuple(x['case_id'] for x in raw['plans'])==module.W.FIRST,
            'same four authorized selector targets and scope required')
    trials=finite_images(module)
    scatter_and_logic(module)
    root_gate_roundtrip(module)
    native_os_fixture(module)
    same=exact_spd_context(module,raw)
    require(not module.OUT.exists() and not module.D.ACTUAL.exists(),'no cold result payload created')
    result={'status':STATUS,**PINS,'independent_review_source_sha256':sha(Path(__file__)),
            'case_order':list(module.W.FIRST),'true_os_and_matrix_actual_launch_authorized':False,
            'future_required_root_source_status':module.REVIEW,'native_producer11_pin_shape_tested':True,
            'separate_true_producer_terminal_and_scope_required':True,'actual_hull_leaves_opened':0,
            'actual_global_outputs_opened':0,'actual_global_matrix_jobs':0,'new_spatial_integration_or_AD':0,
            'same_target_plan_and_SPD_links':same,'synthetic_exact_fraction_image_cases':trials,
            'checks':{'full_same_four_hull_truth_and_point_force_whole_H_role_source_read':True,
                      'two_independent_valid_interval_union_not_intersection':True,
                      'exact_F_H_scatter_and_one_time_symmetry_mirror':True,
                      'signed_split_Arb768_exact_fraction_oracle_and_cached_RF_once':True,
                      'strict_Banach_reverse_cube_exclusion_observer_sigma_plus_e_plus_xi':True,
                      'optional_native_same_R_congruence_and_no_minimum_without_SPD':True,
                      'frozen_actual_phase_gate_roundtrip_real_PF_schema':True,
                      'synthetic_complete_11_pin_native_OS_chain_and_rejections':True},
            'scope':'Cold source and synthetic arithmetic review only. Validates no actual hull leaf or global numerical result. Reuses typed upstream source/metadata readers and closed independently reviewed per-leaf mathematical enclosure authority.'}
    receipt=HERE/'receipt.json'
    with receipt.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2); f.write('\n')
    print(json.dumps({'status':STATUS,'receipt_sha256':sha(receipt),'synthetic_image_cases':trials,'matched_targets':len(same)}))


if __name__=='__main__':
    main()
