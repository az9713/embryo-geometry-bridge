"""Cold review of the frozen two-proof v3 caller and its metadata-only phase gate."""

from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch


HERE=Path(__file__).resolve().parent
PREP=HERE.parent
REDUCE=PREP.parent
EXACT=REDUCE.parent
ROOT=EXACT.parents[3]
SOURCE=ROOT/'scripts/transfer_ambiguity_exact_bridge_reduce_driver_v3.py'
MAIN='445fef269a49f5af4e4a32cf491995f87172bcda214520ed4574d355fb8405d8'
PROTOCOL='1819e5d09d74029c4dbb8ef222c0a8753e00b2768a16a69c359dbe8c18ff7442'
PREFLIGHT='1aa2e62af2703546df1b526e2009fefee4bdaa3274b933647ce17d85464df7b5'
RAW_OS='02940fd7afeaf820f217be711ad286902e8b18970759a21140c7bcd44c3031c6'
V3_OS='802ddb8839f07995dcd5b50b99235df0694f8e84c69c18c6d2fae1aae46c593e'
RAW_TOP='fdaee6daea7ebc2ac9df554281076c3648f6d3858b5cac6c62f254d70a95bf63'
V3_TOP='579d5b81dcad739ad363e708e2d354609f68d412f26f5eb898a2b0a644cab29b'
PURE='d31ca58afe2e1b935d723513fb9fa18f2f978885f8f6ec32e48766f8b15bf0e8'
PILOT='max_contrast_change_small_17'


def require(ok,msg):
    if not ok:
        raise ValueError(msg)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path,digest):
    raw=path.read_bytes()
    require(sha256(raw).hexdigest()==digest,f'changed file {path}')
    return json.loads(raw)


def rejected(callback):
    try:
        callback()
    except ValueError:
        return True
    return False


def gate_roundtrip(v):
    authoritative={'source_sha256':MAIN,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,
                   'old_caller_sha256':v.SOURCE_PINS[v.base.SOURCE],
                   'v3_cold_peer_review_sha256':v.SOURCE_PINS[v.COMPANION/'independent_review_v3/receipt.json'],
                   'v3_source_sha256':v.PHASE['source_sha256'],'v3_protocol_sha256':v.PHASE['protocol_sha256'],
                   'v3_preflight_sha256':v.PHASE['preflight_sha256']}
    probe=HERE/'synthetic_in_memory_review.json'
    mark='7'*64
    cold={'status':v.REVIEW_STATUS,'checks':{'fixture':True},**authoritative}
    original=v.pinned
    def reader(path,h):
        if path==probe:
            require(h==mark,'wrong synthetic approval digest')
            return cold
        return original(path,h)
    stop=AssertionError('cold gate attempted actual domain or packet admission')
    with patch.object(v,'pinned',side_effect=reader),patch.object(v,'admitted',side_effect=stop),\
         patch.object(v.base,'context',side_effect=stop):
        p,links=v.gate(probe,mark)
        require(p==read(PREP/'protocol.json',PROTOCOL)
                and links=={**authoritative,'prelaunch_review_sha256':mark},
                'actual gate does not round trip through frozen protocol and preflight')
        for field in ('protocol_sha256','preflight_sha256','v3_source_sha256','v3_cold_peer_review_sha256'):
            old=cold[field]
            cold[field]='0'*64
            require(rejected(lambda:v.gate(probe,mark)),'wrong reviewed input accepted')
            cold[field]=old
    return authoritative


def phase_metadata(v):
    root=v.base.BRIDGE/'actual'/PILOT/'root_terminal_admission.json'
    args=SimpleNamespace(pilot_id=PILOT,admission=root,admission_sha256=sha(root),
                        raw_observation=v.RAW/PILOT/'os_exit_observation.json',raw_observation_sha256=RAW_OS,
                        companion_observation=v.V3/PILOT/'os_exit_observation.json',companion_observation_sha256=V3_OS)
    a=read(root,args.admission_sha256)
    require(a['status']=='PASS_ROOT_OBSERVED_SINGLETON_BRIDGE_COMPONENT_TERMINALS'
            and a['pilot_id']==PILOT and type(a['observed_exit_code']) is int and a['observed_exit_code']==0,
            'literal observed original producer authority')
    unopened=[]
    original=v.pinned
    def no_leaves(path,h):
        require(path.name in {'os_exit_observation.json','receipt.json','started.json','terminal.json'},
                'metadata-only phase visited a payload')
        unopened.append(path.name)
        return original(path,h)
    with patch.object(v,'pinned',side_effect=no_leaves):
        old,new,raw,fresh,_,_=v.phases(args)
    require(len(unopened)==8 and raw['receipt_sha256']==RAW_TOP and fresh['receipt_sha256']==V3_TOP
            and (raw['unified_exec_session'],raw['observed_exit_code'],raw['native_documented_exit_code'])==(90695,1,2)
            and (fresh['unified_exec_session'],fresh['observed_exit_code'])==(31584,0)
            and old['status']=='HOLD_INDEPENDENT_COMPONENT_CONTAINMENT'
            and new['status']==v.base.NATIVE_STATUS
            and len(old['components'])==len(new['components'])==168
            and new['fresh_surface_H_jobs']==24 and new['remaining_noncontained']==[],
            'actual +L metadata differs from the approved two-proof scope')
    for folder,ob,rec in ((v.RAW/PILOT,raw,old),(v.V3/PILOT,fresh,new)):
        term=read(folder/'terminal.json',ob['terminal_sha256'])
        start=read(folder/'started.json',sha(folder/'started.json'))
        require(all(type(d['pid']) is int and d['pid']==ob['pid'] for d in (term,start,rec,ob)),
                'literal-int top replay identifiers required in the observed chain')
    return {'pilot_id':PILOT,'producer_admission_sha256':args.admission_sha256,
            'raw_os_sha256':RAW_OS,'raw_receipt_sha256':RAW_TOP,'raw_observed_host_exit':1,'raw_native_declared_exit':2,
            'fresh_os_sha256':V3_OS,'fresh_receipt_sha256':V3_TOP,'fresh_observed_exit':0,
            'actual_component_packets_opened':0}


def catalogue_fixture(v):
    h='1'*64
    b=v.base
    native={'components':[],'checked_packet_sha256':{},'plan_packet_sha256':{}}
    new={'components':[]}
    plan={'components':[]}
    manifest={'files':[]}
    for kind,i in b.ORDER:
        fn=f'{kind}_{i:03d}.json'
        pf=b.rel(b.BRIDGE/'plans'/PILOT/fn)
        orig=b.rel(b.BRIDGE/'actual'/PILOT/kind/fn)
        common={'kind':kind,'slot':i,'global_indices':[min(i,600)],'point_force_scope':b.POINT,
                'whole_box_H_scope':b.BOX,'packet_file':orig,'packet_sha256':h,'plan_file_sha256':h,'replay_sha256':h}
        native['components'].append({**common,'replay_file':b.rel(v.RAW/PILOT/fn)})
        new['components'].append({**common,'replay_file':b.rel(v.V3/PILOT/fn)})
        native['checked_packet_sha256'][orig]=native['plan_packet_sha256'][pf]=h
        plan['components'].append({'kind':kind,'slot':i,'file':fn,'sha256':h})
        manifest['files'].append({'kind':kind,'slot':i,'file':orig,'sha256':h})
    def load(fresh):
        return v.adopt_catalogue(PILOT,native,fresh,plan,manifest,native['checked_packet_sha256'],native['plan_packet_sha256'])
    with patch.object(v,'sha',return_value=h):
        entries=load(new)
        require(len(entries)==168 and sum(e[0]['kind']=='surface' for e in entries)==24,'complete synthetic dual catalogue')
        for key,val in (('replay_file',native['components'][0]['replay_file']),('global_indices',[1]),
                        ('slot',True),('whole_box_H_scope',b.POINT)):
            wrong=deepcopy(new)
            wrong['components'][0][key]=val
            require(rejected(lambda:load(wrong)),f'mixed replay coordinate, scope or path accepted: {key}')
        require(rejected(lambda:load({'components':new['components'][:-1]})),'partial companion coverage accepted')
    wrongfile=entries[-1][-1]
    with patch.object(v,'sha',side_effect=lambda path: '2'*64 if path==wrongfile else h):
        require(rejected(lambda:load(new)),'changed last dual leaf was not rejected before payload parse')


def leaf_fixtures(v):
    b=v.base
    for kind in ('bulk','surface'):
        row={'kind':kind,'slot':0,'global_indices':[11],'point_force_scope':b.POINT,
             'whole_box_H_scope':b.BOX,'plan_file_sha256':'a'*64,'packet_file':'product',
             'packet_sha256':'b'*64,'replay_file':'prior','replay_sha256':'c'*64}
        direct={'direct_point_cache':[],'point_force':[0],'whole_box_H':[] if kind=='bulk' else [[0,0]]}
        raw={**row,'checked_force_rows':1,'checked_cache_sites':64 if kind=='bulk' else 12,
             'checked_upper_H_entries':1,'containment_failures':direct,
             'independent_point_F_exact':[['-3','1']],
             'independent_whole_box_H_exact':[[0,0,'-1','1']] if kind=='bulk' else [[0,0,'-3','2']]}
        t=Q(1,2) if kind=='surface' else Q(1)
        shifted=v.refined_endpoints((Q(-3),Q(1)),Q(3,2),Q(1),t,Q(1),Q(1),kind)
        # Exact scalar Minkowski witness for the validated tail erosion.
        require(shifted==(Q(-3)+Q(3,2)-t,Q(1)-Q(3,2)+t),'exact erode/reinflate formula')
        fresh={**raw,'containment_failures':dict.fromkeys(b.FAILURE_KEYS,[]),
               'raw_component_replay_file':'prior','raw_component_replay_sha256':'c'*64,
               'independent_point_F_exact':[list(map(str,shifted))],
               'independent_whole_box_H_exact':[[0,0,'-1','1']],
               'point_original_excess_rounded_tails':[[0,'3/2','1',str(t),True]],
               'fresh_surface_H_tails':[] if kind=='bulk' else [[0,0,'1/4']]}
        prod={'point_force':{v.FKEY:[list(map(str,shifted))],'force_tail_upper_exact':['1']},
              'whole_box_hessian':{v.HKEY:[[0,0,'-2','2']]}}
        def verify(old=raw,new=fresh,outer=prod):
            return v.leaf_refinement(row,row,old,new,outer,(Q(1),Q(1)))
        require(verify()==(1,0 if kind=='bulk' else 1),'valid raw and refined synthetic intervals rejected')
        altered={**fresh,'independent_point_F_exact':[['-1','-1']]}
        require(rejected(lambda:verify(new=altered)),'arbitrary force narrowing accepted')
        wrong={**fresh,'fresh_surface_H_tails':[[0,0,'1']] if kind=='bulk' else []}
        require(rejected(lambda:verify(new=wrong)),'bulk H substitution or missing fresh surface H accepted')
        overlap={**fresh,'independent_whole_box_H_exact':[[0,0,'-3','1']]}
        require(rejected(lambda:verify(new=overlap)),'mere H overlap accepted')


def main():
    require(sha(SOURCE)==sha(PREP/'transfer_ambiguity_exact_bridge_reduce_driver_v3_frozen.py')==MAIN,
            'source and frozen bytes differ')
    p=read(PREP/'protocol.json',PROTOCOL)
    pf=read(PREP/'preflight.json',PREFLIGHT)
    require(p['source_sha256']==pf['source_sha256']==MAIN and pf['protocol_sha256']==PROTOCOL
            and pf['status']=='PASS' and all(v is True for v in pf['checks'].values())
            and pf['actual_packet_payloads_opened']==pf['actual_matrix_work']==0,'changed cold preflight')
    sys.path.insert(0,str(ROOT/'scripts'))
    import transfer_ambiguity_exact_bridge_reduce_driver_v3 as v
    require(Path(v.__file__).resolve()==SOURCE and sha(v.base.KERNEL_SOURCE)==PURE,'unchanged reviewed pure kernel required')
    checks=v.fixtures()
    require(all(value is True for value in checks.values()),'caller cold fixtures failed')
    bindings=gate_roundtrip(v)
    metadata=phase_metadata(v)
    catalogue_fixture(v)
    leaf_fixtures(v)
    result={'status':v.REVIEW_STATUS,**bindings,
            'checks':{'frozen_caller_and_all_pure_and_companion_source_pins_verified':True,
                      'real_frozen_preflight_gate_roundtrip_accepts_review_and_rejects_mixed_phase':True,
                      'actual_plus_L_raw_HOLD2_host1_and_v3_PASS_exit0_metadata_verified':True,
                      'dual_full_168_canonical_scope_route_and_digest_fixtures_pass':True,
                      'strict_point_F_narrowing_and_surface_H_containment_rechecked_in_synthetic_leaves':True,
                      'old_direct_cache_and_bulk_H_remain_required_with_fresh_surface_H':True,
                      'same_binary_R_cache_cell_radius_and_ZQ_via_reviewed_factory':True,
                      'signed_cache_anchor_already_preconditioned_and_added_once_by_pinned_kernel':True},
            'source_of_independent_review_sha256':sha(Path(__file__)),
            'metadata_only_actual_phase':metadata,'pure_kernel_sha256':PURE,
            'actual_packet_payloads_opened':0,'actual_611_matrix_work':0,
            'scope':'Cold prelaunch approval. The companion is a separate actual proof for the refined rows. '
                    'It supplies the actual inserted tail radii and fresh surface H used by the caller leaf checks. '
                    'The original raw HOLD remains unchanged. No root, minimum, or actual numerical reduction reviewed here.',
            'nested_pid_typing_note':'Observed process PIDs and exits have literal integer guards. Original and companion '
                                      'start/end/native PID equality also has exact-byte binding. All actual +L '
                                      'nested PID types were verified as literal integers in this metadata readback.'}
    path=HERE/'review.json'
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':result['status'],'receipt_sha256':sha(path),'actual_leaf_packets_opened':0}),flush=True)


if __name__=='__main__':
    main()
