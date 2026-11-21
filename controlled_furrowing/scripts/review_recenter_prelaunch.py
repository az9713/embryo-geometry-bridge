"""Independent prelaunch review: new numerical center, identical potential/C."""
from pathlib import Path
from fractions import Fraction as Q
import ast,hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs/cohort_recentered_roots_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    out=RUN/'prelaunch_independent_review.json'
    if out.exists():raise RuntimeError('preserve concrete center review')
    m=importlib.import_module('recenter_uniform_roots');p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');old=read(ROOT/p['original_protocol_file'])
    source=Path(m.__file__);text=source.read_text(encoding='utf-8');tree=ast.parse(text)
    guards=m.guard_tests();m.helpers.assert_frozen(p)
    checks={
        'current_source_protocol_preflight_bindings':sha(source)==p['source_hashes']['recenter_uniform_roots.py']==pre['source_sha256'] and pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'preflight_boolean_PASS':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'all_sources_and_inputs_current':all(sha(ROOT/'scripts'/k)==h for k,h in p['source_hashes'].items()) and all(sha(ROOT/k)==h for k,h in p['input_hashes'].items()),
        'all_single_rounding_fingerprint_guards_reexecuted':bool(guards) and all(v is True for v in guards.values()),
        'same_shared_cell_and_parameter_center':p['shared_parameter_cell']==old['shared_parameter_cell'] and p['base_parameters_exact']==old['base_parameters_exact'],
        'same_material_geometry_quadrature_and_observation_contract':all(p[k]==old[k] for k in ('fixed_full_parameters','integration_policy','common_model_observation_contract_sha256')),
        'only_F_K_same_source_doses':{r['record_id'] for r in p['records']}=={'optoGEF2_small_col_F','optoGEF2_small_col_K'} and all(any(r['record_id']==s['record_id'] and Q(r['dose_exact'])==Q(s['dose_exact']) and r['original_seed_state_sha256']==s['saved_state_sha256'] for s in old['records']) for r in p['records']),
        'one_correction_not_iteration':p['center_recipe']['number_of_corrections']==1,
        'same_full_normal_corrector_test_policy':p['attempt_policy']==old['attempt_policy'],
        'translation_not_added_to_Banach_radius':'DO NOT add center translation to Banach corrector r' in p['center_recipe']['translation_policy'],
        'fresh_assembly_and_root_test_required':'newpoint=kernel.build_point' in text and 'proof=kernel.evaluate_cell(newpoint' in text,
        'original_R_and_S_preserved_by_actual_code':"fingerprint(values,{'u'})!=fingerprint(newvalues,{'u'})" in text and "helpers.binary(work/'affine_predictor.npz',(ROOT/row['predictor_file']).read_bytes())" in text,
        'numerical_auxiliary_state_not_mislabelled_as_experiment_or_converged_solver':"'converged':False" in text and "'solver_status':'NO_NONLINEAR_SOLVER_RUN'" in text,
        'original_seed_identity_and_new_proof_center_separated':all(k in text for k in ("'original_seed_state_file'","'proof_center_state_file'","'observer_requires_explicit_versioned_center_offset_adapter':True")),
        'physical_potential_array_and_parameter_fingerprints_preserved':all(k in text for k in ("'full_physical_potential_fingerprint':physical(values)","'new_full_physical_potential_fingerprint':physical(newvalues)","'theta0_exact':p['base_parameters_exact']")),
        'no_model_or_center_work_before_review':not (RUN/'records').exists() and pre['actual_RF_or_center_calculations']==0,
        'launch_review_and_exclusive_job_guard_present':"with (OUT/'active.lock').open('x'" in text and "helpers.binary(OUT/'launch_review.json',raw)" in text}
    for row in p['records']:
        oldcache=read(ROOT/row['original_point_cache_file']);values=m.arrays(oldcache['frozen_binary_cache']);state=read(ROOT/row['original_seed_state_file'])
        checks[row['record_id']+'_seed_cache_dimensions_and_original_failed_context']=values['u'].shape==(611,) and values['preconditioner'].shape==(611,611) and oldcache['input_state_sha256']==row['original_seed_state_sha256'] and state['parameters']==p['fixed_full_parameters'] and read(ROOT/row['original_result_file'])['status']=='UNRESOLVED'
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(source),
        'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'review_source_sha256':sha(__file__),
        'mathematical_review':'The correction is solely a choice of auxiliary binary proof center, not an assumed exact Newton iterate/root. The unchanged frozen potential and same parameter cell must be proved by all-array-except-u and full physical target fingerprints. A fresh Banach theorem around u1+S*delta uses its own normal corrector radius. Any future observer based on original stored curves must separately include center translation/serialization in its nodal hull; original same-center consumers must not silently accept the new cache.',
        'scope':'Prelaunch recipe/lineage/target/guard review only. No corrected centers, root existence, all17 tuple, observer result, biological fit or continuum result is claimed.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':len(checks),'guards':len(guards),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
