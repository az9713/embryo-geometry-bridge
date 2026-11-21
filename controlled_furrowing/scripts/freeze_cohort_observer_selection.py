"""Translate the same immutable all17 provider list; never select on outcomes."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_static_observers_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    output=RUN/'provider_selection.json'
    if output.exists():raise RuntimeError('preserve frozen provider selection')
    mf=ROOT/'runs/cohort_uniform_root_review_v1/provider_manifest.json';source=read(mf)
    cp=ROOT/'runs/cohort_uniform_roots_v2/protocol.json';p=read(cp)
    terminal=ROOT/'runs/cohort_uniform_roots_v2/results.json';t=read(terminal)
    job=read(ROOT/'runs/cohort_uniform_roots_v2/job.json')
    supplement=ROOT/'runs/cohort_radius_supplement_v1/results.json'
    if job['status']!='FINISHED' or job['results_sha256']!=sha(terminal):raise RuntimeError('true terminal cohort required')
    if t['certified']!=15 or t['unresolved']!=2 or len(source['records'])!=17:raise RuntimeError('expected retained15/2 state required')
    rows=[];checks={}
    for declared,entry,completed in zip(p['records'],source['records'],t['records']):
        rid=declared['record_id'];f=ROOT/entry['certificate_file'];r=read(f)
        if not rid==entry['record_id']==completed['record_id']==r['record_id']:raise RuntimeError('record order/identity changed')
        if sha(f)!=entry['certificate_sha256'] or sha(f)!=completed['sha256']:raise RuntimeError('provider bytes changed')
        if entry['provider_protocol_file']!=cp.relative_to(ROOT).as_posix() or r['protocol_sha256']!=sha(cp):raise RuntimeError('provider protocol changed')
        if Q(r['dose_exact'])!=Q(declared['dose_exact']) or r['parameter_cell_sha256']!=p['shared_parameter_cell_sha256'] or r['model_observation_contract_sha256']!=p['common_model_observation_contract_sha256']:raise RuntimeError('numeric target/cell/contract changed')
        rows.append({'record_id':rid,'root_result_file':entry['certificate_file'],'root_result_sha256':sha(f),
            'root_protocol_file':entry['provider_protocol_file'],'root_protocol_sha256':sha(cp)})
        checks[rid+'_identity_hash_target']=True
    packet={'version':1,'parameter_cell_sha256':p['shared_parameter_cell_sha256'],'common_contract_sha256':p['common_model_observation_contract_sha256'],
        'source_provider_manifest_file':mf.relative_to(ROOT).as_posix(),'source_provider_manifest_sha256':sha(mf),
        'terminal_cohort_file':terminal.relative_to(ROOT).as_posix(),'terminal_cohort_sha256':sha(terminal),
        'failed_radius_supplement_file':supplement.relative_to(ROOT).as_posix(),'failed_radius_supplement_sha256':sha(supplement),
        'records':rows,'selection_rule':'Exactly every v2 terminal provider, in source order. F/K remain UNRESOLVED; failed radius-only supplement is retained separately. No subset, retuning or inferred root.',
        'independent_root_replay_status':'PENDING_AT_SELECTION_FREEZE','source_sha256':sha(__file__)}
    output.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
    checks.update(all17_preserved=len(rows)==17,no_provider_replaced_or_dropped=True,
        no_all17_claim=t['all17_certified'] is False,exact15_certified_2_unresolved=True,
        failed_radius_history_bound=read(supplement)['certified']==0 and read(supplement)['unresolved']==2)
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'selection_sha256':sha(output),
        'adapter_source_sha256':sha(__file__),'source_manifest_sha256':sha(mf),
        'scope':'Exact schema translation and immutable selected-input provenance, not an independent root or observer proof.'}
    (RUN/'selection_verification.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'selection_sha256':sha(output)}))
    if receipt['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
