"""Review an explicit new-center to original-graph bound adapter."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs/cohort_recentered_observers_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    out=RUN/'prelaunch_independent_review.json'
    if out.exists():raise RuntimeError('preserve concrete adapter review')
    m=importlib.import_module('recentered_static_observers');p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');tests=read(RUN/'synthetic_checks.json')
    rootp=read(ROOT/p['root_protocol_file']);text=Path(m.__file__).read_text(encoding='utf-8')
    checks={
        'source_protocol_preflight_and_synthetics_bound':pre['source_sha256']==tests['source_sha256']==sha(m.__file__) and pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'current_preflight_all_boolean_PASS':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'current_synthetic_checks_all_boolean_PASS':tests['status']=='PASS' and bool(tests['checks']) and all(v is True for v in tests['checks'].values()),
        'all_frozen_sources_and_inputs_current':all(sha(ROOT/k)==h for k,h in {**p['source_hashes'],**p['input_hashes']}.items()),
        'same_original_contract_and_full_parameter_cell':p['common_contract_sha256']=='2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82' and p['parameter_cell_sha256']=='64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf',
        'new_explicit_adapter_kind':p['adapter_kind']=='RECENTERED_ROOT_TO_ORIGINAL_SEED_GRAPH',
        'normal_corrector_separate_from_center_translation':p['normal_corrector_radius_is_separate'] is True,
        'original_strict_consumer_not_weakened':sha(ROOT/'scripts/cohort_static_observers.py')=='981d7442aefc793a986b9ecb36da175b631d2cd15fa5f3e93754e818b6aefc14' and tests['checks']['reject_invalid_3'] is True,
        'finite_shape_dtype_and_R_array_changes_rejected':all(tests['checks'][k] is True for k in ('reject_invalid_0','reject_invalid_1','reject_invalid_2')),
        'no_new_composition_already_written':not (RUN/'results.json').exists() and not (RUN/'records').exists(),
        'new_root_and_application_independent_review_required':p['independent_root_and_application_review_required'] is True,
        'no_shared_time_or_stability_claim':p['no_common_path_or_stability_claim'] is True}
    # Reexecute pure, no-assembly adversarial tests; immutable original output.
    replay=m.tests();checks['synthetic_checks_reexecuted']=replay==tests
    for row in rootp['records']:
        rid=row['record_id'];root=read(ROOT/row['output_folder']/'result.json')
        seed=read(ROOT/row['original_seed_state_file']);prior=read(ROOT/row['original_point_cache_file']);new=read(ROOT/root['point_cache_file'])
        a=m.old.decode_cache(prior);b=m.old.decode_cache(new);fp=m.unchanged_cache(a,b)
        center=read(ROOT/root['center_refinement_file']);physical=m.full_potential_fingerprint(a,prior['parameters'],rootp,row)
        checks[rid+'_every_physical_cache_array_and_R_unchanged']=fp==center['cache_arrays_except_u_fingerprint']==center['new_cache_arrays_except_u_fingerprint']
        checks[rid+'_full_target_including_parameters_and_dose_unchanged']=physical==center['full_physical_potential_fingerprint']==center['new_full_physical_potential_fingerprint'] and prior['parameters']==new['parameters']==seed['parameters']
        checks[rid+'_separate_source_and_aux_center_identity']=root['original_seed_state_sha256']==row['original_seed_state_sha256'] and new['input_state_sha256']==root['proof_center_state_sha256'] and root['proof_center_state_sha256']!=row['original_seed_state_sha256']
        checks[rid+'_unchanged_binary_S']=sha(ROOT/root['predictor_file'])==sha(ROOT/row['predictor_file'])
        checks[rid+'_normal_radius_small_positive_and_separate']=m.old.accepted_radius(root['uniform_cell_result'])==Q(1,10**11) and root['normal_corrector_radius_is_not_center_translation'] is True
    checks['actual_translation_present_in_coordinate_error_formula']='abs(new_exact[k//2][k%2]-centers[k//2][k%2])+max(map(abs,hulls[i]))+radius' in text
    checks['reviewed_geometry_and_graph_lemmas_reused']='old.certify_box(' in text and 'old.curve_bound(' in text and 'errors=dict(zip(old.OBSERVABLES,(ea+eb,ea+eb,ea,eb)))' in text
    packet={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'observer_source_sha256':sha(m.__file__),'observer_protocol_sha256':sha(RUN/'protocol.json'),
        'observer_preflight_sha256':sha(RUN/'preflight.json'),'review_source_sha256':sha(__file__),
        'mathematical_review':'Fresh root is centered at exactX+u1+S*delta with normal radius r. Relative to original stored graph nodes the error is bounded by the exact center-to-original serialization gap plus max|S*delta| plus r, coordinate by coordinate. Physical cache arrays and full parameter/input fingerprints are unchanged. The reviewed P2 whole-body and graph-transfer lemmas therefore apply to this enlarged nodal box. Translation is never inserted into the fresh Banach self-map radius. Baseline graph observables use the original seed; no source data or material parameters are changed.',
        'scope':'Explicit offset-adapter prelaunch review only; independent fresh F/K root replay and actual geometry/output application review are still required.'}
    out.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':packet['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if packet['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
