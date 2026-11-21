"""Review the independent assembly protocol, not its still-unexecuted results."""
from pathlib import Path
from fractions import Fraction as Q
import ast,hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_uniform_root_review_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    out=RUN/'prelaunch_review.json'
    if out.exists():raise RuntimeError('preserve concrete replay review')
    module=importlib.import_module('review_cohort_uniform_roots')
    p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');manifest=read(ROOT/p['manifest'])
    source=Path(module.__file__);tree=ast.parse(source.read_text(encoding='utf-8'))
    calls=[n.func.attr for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)]
    guards=module.guard_tests()
    checks={
        'review_source_protocol_and_preflight_bound':p['review_source_sha256']==pre['review_source_sha256']==sha(source) and pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'preflight_nonempty_boolean_PASS':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'all_frozen_files_current':all(sha(ROOT/k)==h for k,h in p['file_hashes'].items()),
        'one_expected_parameter_cell':p['shared_parameter_cell_sha256']=='64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf',
        'one_expected_common_contract':p['model_observation_contract_sha256']=='2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82',
        'all17_providers_retained':len(p['records'])==17 and p['records']==manifest['records'],
        'separate_assembly_no_producer_point_or_cell_calls':'build_point' not in calls and 'evaluate_cell' not in calls,
        'fresh_192bit_precision_declared':p['precision_bits']==192,
        'lightweight_acceptance_normalization_symbolic_guards_reexecuted':bool(guards) and all(v is True for v in guards.values()),
        'guards_serializable_without_numpy_booleans':bool(json.dumps(guards)),
        'strict_full_radius_policy':p['acceptance'].startswith('Strict independent Z<1 and eta+Z*r<r at EXACT producer FULL corrector radius.'),
        'no_claim_of_fresh_M_assembly':p['actual_M_reuse_review_sha256']=='4b03db70583be266a8d4ae581f19c8db0099bc028ddd6af82a5c990a91d8e600',
        'no_replay_result_written_yet':not (ROOT/'evidence/cohort_uniform_roots_independent_review.json').exists() and not (RUN/'records').exists()}
    statuses={};record_checks=[]
    for selected in p['records']:
        f=ROOT/selected['certificate_file'];r=read(f);provider=ROOT/selected['provider_protocol_file']
        statuses[selected['record_id']]=r['status']
        record_checks.append(sha(f)==selected['certificate_sha256'] and r['protocol_sha256']==sha(provider) and r['parameter_cell_sha256']==p['shared_parameter_cell_sha256'] and r['model_observation_contract_sha256']==p['model_observation_contract_sha256'])
    checks['all17_current_provider_certificate_protocol_bindings']=all(record_checks)
    checks['15_certified_and_two_unresolved_not_silently_deleted']=sum(s=='CERTIFIED_UNIFORM_LOCAL_ROOT' for s in statuses.values())==15 and {k for k,v in statuses.items() if v=='UNRESOLVED'}=={'optoGEF2_small_col_F','optoGEF2_small_col_K'}
    mindex=next(i for i,r in enumerate(p['records']) if r['record_id']=='optoGEF2_small_col_M')
    original=read(ROOT/'runs/cohort_uniform_roots_v1/protocol.json')['records'][mindex]
    checks['M_original_independent_cache_predictor_review_reuse_passes']=all(v is True for v in module.reuse_M(original,read(ROOT/p['records'][mindex]['certificate_file']))['checks'].values())
    text=source.read_text(encoding='utf-8')
    checks['failed_producer_status_cannot_become_uniform_audit_PASS']="checks['producer_declares_certified_uniform_root']=result['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'" in text and "'uniform_quantifier_verified':alltrue(checks)" in text
    module.bindings(p);checks['frozen_environment_and_input_guards_reexecuted']=True
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'review_source_sha256':sha(source),'protocol_sha256':sha(RUN/'protocol.json'),
        'preflight_sha256':sha(RUN/'preflight.json'),'parent_review_source_sha256':sha(__file__),
        'mathematical_review':'Analytic bulk and surface coefficient derivatives agree with the declared affine five-parameter force law. The frozen binary predictor defect H*S+B is explicitly retained; Taylor remainder is bounded on the parameter/displacement hull and the preconditioned absolute-row majorant bounds the infinity operator norm. Strict independent self-map/contraction at the full radius proves only the accepted nearby root families. Exact upper endpoints determine all norm maxima. Shared reviewed local algebra is reused; assembly/columns/majorant implementation are separate from producer.',
        'scope':'Prelaunch review of source/protocol/math/lineage; expected completed result is partial15/17 with F/K unresolved. No actual root replay, source-matched fit, continuum bound or goal completion is claimed.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':len(checks),'guard_checks':len(guards),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
