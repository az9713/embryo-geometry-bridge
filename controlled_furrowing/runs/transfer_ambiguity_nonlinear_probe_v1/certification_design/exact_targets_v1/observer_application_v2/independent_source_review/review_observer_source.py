"""Cold review for the exact-spatial four-observer scalar transfer."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


HERE=Path(__file__).resolve().parent
APP=HERE.parent
EXACT=APP.parent
DESIGN=EXACT.parent
RUN=DESIGN.parent
ROOT=RUN.parents[1]
SOURCE=ROOT/'scripts/transfer_ambiguity_exact_four_observers.py'
OLD=DESIGN/'cached_pilots_v1/observer_application_v1'
FOUR=EXACT/'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json'
PHASE={'source_sha256':'4779f97e37fdb7b4bd6cb66dc5a7934ebd6adf12f5c007ea88e0c26244845cb1',
       'protocol_sha256':'080729ec73de6d92d1f2e4ceeaa0ac748bf252a3c56295b19eeccb8338bb3a14',
       'preflight_sha256':'43f518b1c2fcefee7812214e02c8bcd1a7995189603783cd14582e44d5feb3d9'}
STATUS='PASS_INDEPENDENT_EXACT_FOUR_OBSERVER_SOURCE'
CASES=tuple((side,role) for side in ('max_contrast_change','max_negative_contrast_change')
            for role in ('small_17','large_03'))
VARIANTS=('current_disk__spatial_max','current_disk__regional_mean',
          'current_radial_line__spatial_max','current_radial_line__regional_mean')
OBS=('maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um')


def check(condition,msg):
    if not condition:
        raise ValueError(msg)


def sha(path):
    h=sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def load(path,h):
    raw=path.read_bytes()
    check(sha256(raw).hexdigest()==h,f'changed input {path}')
    return json.loads(raw)


def one(rows,predicate):
    row,=(row for row in rows if predicate(row))
    return row


def executable_gate(source,protocol,pref):
    """Run the actual CLI function through the admission predicates before shadowing any actual transfer."""
    class StopBeforeApplication(Exception):
        pass

    gate_file=HERE/'receipt.json'
    args=SimpleNamespace(**PHASE,review=gate_file,review_sha256='1'*64)
    valid={'status':STATUS,**PHASE,'checks':{'fixture':True}}
    reads={APP/'protocol.json':protocol,APP/'preflight.json':pref,gate_file:valid}
    with patch.object(source,'read',side_effect=lambda path,h,links:reads[path]), \
         patch.object(source,'plan',return_value=protocol), \
         patch.object(source,'authorities',side_effect=StopBeforeApplication), \
         patch.object(source,'save',side_effect=AssertionError('cold run attempted to persist observer results')):
        try:
            source.run(args)
        except StopBeforeApplication:
            pass
        else:
            raise AssertionError('positive fixture did not reach first actual instruction')
        for bad in ({**valid,'status':'PASS'}, {**valid,'preflight_sha256':'2'*64}, {**valid,'checks':{'x':1}}):
            reads[gate_file]=bad
            try:
                source.run(args)
            except ValueError:
                continue
            raise AssertionError('malformed gate accepted')
        reads[gate_file]=valid
        reads[APP/'preflight.json']={**pref,'actual_observer_applications':False}
        try:
            source.run(args)
        except ValueError:
            return
    raise AssertionError('Boolean masquerading as zero accepted')


def math_fixtures(source):
    # On each half the negative endpoint shape has minimum -1/8 at t=1/4 or 3/4.
    for t in (Q(1,4),Q(3,4)):
        shape=(2*(t-Q(1,2))*(t-1),4*t*(1-t),2*t*(t-Q(1,2)))
        check(sum(shape)==1 and sum(map(abs,shape))==Q(5,4),'P2 shape 5/4 bound')
    for t in (0,1):
        check(sum(map(abs,(4*t-3,4-8*t,4*t-1)))==8,'endpoint derivative-sum bound')
    check(source.inflate(('1/3','1/2'),Q(1,5))==['2/15','7/10'],'signed outward inflation mismatch')
    expected=[('IN',['0','2']),('UNRESOLVED',['2','3']),('OUT',['5/2','4'])]
    fixtures=[(('-1','2'),Q(0)),(('3','4'),Q(1)),(('-4','-5/2'),Q(0))]
    check([source.classify(v,c,Q(2)) for v,c in fixtures]==expected,'absolute residual boundary misuse')
    # A new nodal error can invalidate a previously valid monotonicity guard.
    previous,larger,b=Q(1,40),Q(1,8),Q(1)
    check(b-8*previous>0 and b-8*larger==0,'strict monotonicity boundary fixture incorrect')
    e1,e2=Q(1,100),Q(1,10)
    check(Q(5,4)*(1+Q(3))*e1==Q(1,20) and Q(5,4)*(1+Q(3))*e2==Q(1,2),'linear Eulerian graph bound fixture')


def main():
    check(not (APP/'results.json').exists(),'source review must precede actual application')
    check(sha(SOURCE)==PHASE['source_sha256'],'wrong frozen producer source')
    protocol=load(APP/'protocol.json',PHASE['protocol_sha256'])
    pf=load(APP/'preflight.json',PHASE['preflight_sha256'])
    check(protocol['status']=='FROZEN_FOUR_EXACT_SPATIAL_OBSERVER_TRANSFER'
          and protocol['source_sha256']==pf['source_sha256']==PHASE['source_sha256']
          and pf['protocol_sha256']==PHASE['protocol_sha256'] and pf['status']=='PASS'
          and type(pf['actual_observer_applications']) is int and pf['actual_observer_applications']==0,
          'unbound preparation phase')
    links=protocol['inputs_sha256']
    def opened(p):
        return load(p,links[rel(p)])
    for name,h in links.items():
        path=(ROOT/name).resolve()
        check(path.is_relative_to(ROOT) and rel(path)==name and sha(path)==h,f'changed closure {name}')
    old=opened(OLD/'results.json')
    review=opened(OLD/'independent_readback/receipt.json')
    four=opened(FOUR)
    theorem=opened(ROOT/'evidence/nodal_observer_error_independent_review.json')
    manifest=opened(EXACT/'preparation/target_manifest.json')
    case_ids=tuple('_'.join(x) for x in CASES)
    check(review['status']=='PASS_INDEPENDENT_CACHED_ROOT_OBSERVER_APPLICATION_READBACK'
          and review['producer_application_functions_called'] is False and len(review['all_32_comparisons'])==32
          and review['input_sha256'][rel(OLD/'results.json')]==links[rel(OLD/'results.json')]
          and Q(review['conditional_delta_um_exact'])==Q(protocol['allowance_um_exact'])==Q(old['frozen_model_discrepancy_allowance_um_exact'])
          and four['status']=='PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA'
          and four['independent_exact_spatial_root_count']==4 and four['full_five_gate_goal_achieved'] is False
          and len(four['current_file_sha256'])==four['current_file_binding_count']==3006
          and all(links[name]==h for name,h in four['current_file_sha256'].items())
          and tuple(c['pilot_id'] for c in four['cases'])==tuple(r['pilot_id'] for r in old['four_pilots'])
          ==tuple(manifest['pilot_order'])==case_ids and all(v is True for v in theorem['checks'].values()),
          'four current certifications and independently reviewed graph theorem required')
    cases=[]
    for root,prior in zip(four['cases'],old['four_pilots'],strict=True):
        pid=root['pilot_id']; ident=root['target_identity']; sm=root['selfmap']
        model=one(manifest['targets'],lambda c:c['pilot_id']==pid)
        contract=opened(ROOT/model['model_contract'])
        cache=opened(ROOT/contract['point_cache'])
        geo=opened(EXACT/'independent_actual/replay/receipt.json')
        whole=one(geo['geometry_records'],lambda c:c['pilot_id']==pid and Q(c['radius_exact'])==Q(1,2**34))
        native=opened(ROOT/root['independent_proof_file'])
        os=opened(ROOT/root['independent_os_observation_file'])
        check(root['independent_proof_sha256']==links[root['independent_proof_file']]
              and root['independent_os_observation_sha256']==links[root['independent_os_observation_file']]
              and native['status']=='PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM'
              and os['status']=='PASS_OBSERVED_INDEPENDENT_GLOBAL_REDUCTION_TRUE_OS_EXIT'
              and type(os['observed_exit_code']) is int and os['observed_exit_code']==0 and type(os['pid']) is int
              and os['pid']==native['pid'] and os['receipt_sha256']==root['independent_proof_sha256']
              and native['target_identity']==ident and native['conservative_compared_bounds']==sm
              and native['same_geometry_admissibility_sha256']==whole['sha256']
              and Q(root['radius_exact'])==Q(native['radius_exact'])==Q(1,2**34)
              and ident['field_sha256']==contract['source_field_sha256']==cache['input_state_sha256']
              and ident['cache_packet_sha256']==contract['point_cache_sha256']==links[contract['point_cache']]
              and ident['embedded_cache_sha256']==contract['frozen_binary_cache_sha256']==cache['frozen_binary_cache']['sha256']
              and ident['R_binary_sha256']==contract['binary_preconditioner_sha256']
              and ident['model_contract_sha256']==model['model_contract_sha256']==links[model['model_contract']]
              and ident['parameter_cell_sha256']==contract['parameter_cell_sha256']
              and Q(ident['dose_exact'])==Q(prior['exact_dose'])==Q(cache['optical_dose_exact'])
              ==Q(contract['exact_requested_drive']['optical_dose_exact'])
              and Q(ident['gain_exact'])==Q(prior['exact_gain'])==Q(cache['common_gain_exact'])
              ==Q(contract['exact_requested_drive']['common_gain_exact'])
              and Q(ident['drive_exact'])==Q(prior['exact_drive'])
              ==Q(contract['exact_requested_drive']['activation_strength_exact'])
              ==Q(cache['activation_strength_exact'])==Q(ident['dose_exact'])*Q(ident['gain_exact']),
              f'different field, proof, cache, R, radius or drive: {pid}')
        if prior['role']=='small_record':
            pairs={(r['variant'],r['observable']) for r in prior['root_observer_comparisons']}
            check(len(prior['root_observer_comparisons'])==16 and pairs==set((v,o) for v in VARIANTS for o in OBS)
                  and prior['record_id']=='optoGEF2_small_col_L','expected all sixteen distinct small-L comparisons')
        else:
            check(prior['role']=='large_final' and prior['record_id']=='fig4f_col_I','unexpected large observation')
        cases.append({'pilot_id':pid,'proof_sha256':root['independent_proof_sha256'],
                      'contract_sha256':ident['model_contract_sha256'],'field_sha256':ident['field_sha256']})
    spec=importlib.util.spec_from_file_location('observer_cold_peer_review',SOURCE)
    candidate=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(candidate)
    executable_gate(candidate,protocol,pf)
    math_fixtures(candidate)
    check(not (APP/'results.json').exists(),'cold review wrote actual application data')
    result={'status':STATUS,**PHASE,'review_source_sha256':sha(Path(__file__)),'checks':{
        'pinned_complete_spatial_closure_and_native_zero_exits':True,
        'aligned_four_target_field_cache_R_and_exact_gain_times_dose':True,
        'same_reviewed_saved_observers_and_original_xi_with_enlarged_radius_guard':True,
        'P2_shape_and_derivative_bounds_graph_coefficient_and_strict_margin':True,
        'signed_large_gap_and_all_32_L_scope_no_other_state_admission':True,
        'actual_gate_entry_fixture_and_malformed_receipt_rejections':True,
        'exact_outward_and_closed_IN_strict_OUT_fixtures':True},
        'rehash_input_count':len(links),'proof_cases':cases,'actual_observer_applications':0,
        'scope':'Cold source and exact authority review. Independent numeric arithmetic on stored fields follows the admitted producer result.'}
    path=HERE/'receipt.json'
    with path.open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':STATUS,'review_sha256':sha(path),'inputs':len(links)}))


if __name__=='__main__':
    main()
