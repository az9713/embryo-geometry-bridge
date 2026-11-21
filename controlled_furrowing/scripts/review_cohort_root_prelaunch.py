"""Independent source/data/launch review of the shared-cell cohort wrapper.

No Hessian, equilibrium optimization or interval assembly is performed here.
The completed numerical proof still requires its own independent review.
"""
from pathlib import Path
from fractions import Fraction as Q
import csv,hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_uniform_roots_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
cell=lambda c:{k:list(map(str,map(Q,c[k]))) for k in sorted(c)}
digest=lambda c:hashlib.sha256(json.dumps(cell(c),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def main():
    out=RUN/'prelaunch_independent_review.json'
    if out.exists():raise RuntimeError('preserve launch review; do not overwrite')
    p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');guards=read(RUN/'guard_tests.json')
    contract=read(ROOT/p['common_model_observation_contract_file'])
    direct=read(ROOT/'runs/material_cohort_direct_v1/protocol.json')
    results=read(ROOT/'runs/material_cohort_direct_v1/results.json')
    old=read(ROOT/'runs/parametric_root_cells_v1/actual_M/root_certificate.json')
    audit=read(ROOT/'evidence/parametric_actual_M_independent_review.json')
    with (ROOT/'data/derived/calibration_small_18um.csv').open(encoding='utf-8-sig',newline='') as f:data=list(csv.DictReader(f))
    checks={}
    checks['frozen_wrapper_source_and_preflight_bound']=pre['source_sha256']==p['source_hashes']['cohort_uniform_roots.py']==sha(ROOT/'scripts/cohort_uniform_roots.py') and pre['protocol_sha256']==sha(RUN/'protocol.json')
    checks['producer_preflight_all_boolean_pass']=pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values())
    checks['guard_tests_all_boolean_pass_bound']=guards['status']=='PASS' and all(v is True for v in guards['checks'].values()) and guards['source_sha256']==sha(ROOT/'scripts/cohort_uniform_roots.py') and sha(RUN/'guard_tests.json')==pre['guard_tests_sha256']
    checks['all_direct_source_kernel_hashes_current']=all(sha(ROOT/'scripts'/n)==h for n,h in p['source_hashes'].items())
    checks['all_contract_kernel_hashes_current']=all(sha(ROOT/'scripts'/n)==h for n,h in p['common_contract_source_hashes'].items())
    checks['all_frozen_inputs_current']=all(sha(ROOT/n)==h for n,h in p['input_hashes'].items())
    checks['no_generated_record_or_live_job_before_launch']=not (RUN/'results.json').exists() and not (RUN/'job.json').exists() and not (RUN/'records').exists()
    checks['exact_training_records_and_order']=len(data)==17 and [r['record_id'] for r in p['records']]==[r['record_id'] for r in data] and len({r['record_id'] for r in data})==17
    checks['same_cell_as_independently_reviewed_M']=cell(p['shared_parameter_cell'])==cell(old['parameter_cell']) and digest(p['shared_parameter_cell'])==p['shared_parameter_cell_sha256']==audit['root_audit']['parameter_cell_sha256']
    checks['same_cell_in_common_contract']=cell(contract['parameter_cell'])==cell(p['shared_parameter_cell']) and contract['parameter_cell_sha256']==p['shared_parameter_cell_sha256']
    checks['common_contract_hash_bound']=sha(ROOT/p['common_model_observation_contract_file'])==p['common_model_observation_contract_sha256']
    checks['contract_primary_static_operator_semantics']=contract['primary_variant']=='current_disk__spatial_max' and len(contract['endpoint_definitions'])==4 and contract['source_temporal_semantics_validated'] is False and contract['common_path_or_energetic_stability_certified'] is False
    checks['fixed_coefficient_geometry_and_normalization']=p['fixed_full_parameters']==direct['new_parameters'] and all(direct['new_parameters'][k]==v for k,v in contract['fixed_full_parameters'].items()) and p['fixed_full_parameters']['bulk']==5000 and p['fixed_full_parameters']['mu']==1
    checks['one_actual_M_reuse_only']=sum(r['mode']=='reuse_immutable_actual_M' and r['record_id']=='optoGEF2_small_col_M' for r in p['records'])==1 and sum(r['mode']=='prepare_and_one_affine_test' for r in p['records'])==16
    checks['reuse_root_audit_acceptance_and_binding']=audit['status']=='PASS' and audit['root_audit']['status']=='PASS' and bool(audit['root_audit']['checks']) and all(v is True for v in audit['root_audit']['checks'].values()) and audit['root_audit']['certificate_sha256']==sha(ROOT/'runs/parametric_root_cells_v1/actual_M/root_certificate.json')
    per_record={}
    for row,datum in zip(p['records'],data):
        r=next(x for x in results['records'] if x['record_id']==row['record_id'])
        state=read(ROOT/row['saved_state']);gain=Q(direct['new_gain_exact']);dose=Q(datum['activation_magnitude'])
        c={'state_hash_and_direct_record':sha(ROOT/row['saved_state'])==row['saved_state_sha256']==r['sha256'],
           'exact_declared_dose_gain_product':Q(row['dose_exact'])==dose and Q(row['base_gain_exact'])==gain and Q(row['exact_requested_strength'])==gain*dose==Q(r['requested_strength_exact']),
           'binary_solver_input_separately_bound':Q.from_float(float(state['activation_strength']))==Q(row['solver_strength_binary_exact'])==Q(r['solver_input_strength_exact']),
           'correct_material_dimensions_and_policy':state['parameters']==p['fixed_full_parameters'] and len(state['free_displacement'])==611 and len(state['triangles'])==144 and state['converged'] is True,
           'record_output_owned':(ROOT/row['output_folder']).resolve().is_relative_to(RUN.resolve()),
           'source_state_input_frozen':p['input_hashes'][row['saved_state']]==row['saved_state_sha256']}
        per_record[row['record_id']]=c
    checks['all17_source_state_and_parameter_associations']=all(all(v is True for v in c.values()) for c in per_record.values())
    policy=p['attempt_policy']
    checks['one_affine_attempt_no_cell_shrink']=policy['method']=='affine' and policy['exactly_one_shared_cell_per_record'] is True and policy['radius_multipliers']==[1]
    checks['failed_bound_tests_unresolved_not_OUT']='UNRESOLVED, never OUT' in policy['failed_sufficient_test']
    checks['stationary_tuple_scope_explicit']=p['branch_law_semantics']=='independent_stationary_root_per_record' and p['no_new_equilibrium_optimization'] is True and p['heldout_response_files_read'] is False and p['response_values_used_for_root_acceptance'] is False
    checks['validated_kernel_unchanged']=sha(ROOT/'scripts/parametric_root_cells.py')=='362a19548aa60c4234017316e22deef4a80e93b73c725f47132573a1a9103b60' and sha(ROOT/'scripts/algebraic_root_feasibility.py')=='7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf'
    module=importlib.import_module('cohort_uniform_roots')
    checks['selfcontained_guard_tests_independently_reexecuted']=module.guard_tests()==guards['checks']
    module.assert_frozen(p)
    checks['producer_integrity_and_environment_guards_reexecuted']=True
    text=(ROOT/'scripts/cohort_uniform_roots.py').read_text(encoding='utf-8')
    checks['prepared_metadata_guard_in_actual_execution_path']='validate_prepared_metadata(p,state,model.reference,model.triangles,u,asdict(model.integration_policy),verification,local_terms)' in text
    checks['review_bytes_archived_at_actual_launch']='reviewbytes=reviewpath.read_bytes()' in text and "binary(launchdir/f'review_{j:03}.json',reviewbytes)" in text
    checks['exclusive_lease_and_creation_identity_guard_present']="with lock.open('x'" in text and "identity['creation_ticks']==prior['creation_ticks']" in text
    packet={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'per_record_checks':per_record,
        'source_sha256':sha(ROOT/'scripts/cohort_uniform_roots.py'),'protocol_sha256':sha(RUN/'protocol.json'),
        'preflight_sha256':sha(RUN/'preflight.json'),'guard_tests_sha256':sha(RUN/'guard_tests.json'),
        'review_source_sha256':sha(__file__),'shared_cell_sha256':p['shared_parameter_cell_sha256'],
        'scope':'Independent wrapper/source/data/launch review only. No fresh interval proof performed here; completed cohort bounds require independent reconstruction and output/geometry composition review.'}
    out.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':packet['status'],'checks':len(checks),'record_checks':sum(map(len,per_record.values())),'receipt_sha256':sha(out)}))
    if packet['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
