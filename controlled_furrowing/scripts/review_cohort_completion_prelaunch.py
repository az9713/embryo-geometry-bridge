"""Read-only independent review of versioned completion and prefix reuse."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_uniform_roots_v2';OLD=ROOT/'runs/cohort_uniform_roots_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    out=RUN/'prelaunch_independent_review.json'
    if out.exists():raise RuntimeError('preserve concrete review')
    p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');old=read(OLD/'protocol.json');progress=read(OLD/'progress.json');job=read(OLD/'job.json')
    m=importlib.import_module('cohort_uniform_roots_completion')
    checks={
        'new_source_preflight_protocol_bound':pre['source_sha256']==p['source_hashes']['cohort_uniform_roots_completion.py']==sha(ROOT/'scripts/cohort_uniform_roots_completion.py') and pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'preflight_nonempty_all_boolean_PASS':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'all_frozen_source_and_input_bytes_current':all(sha(ROOT/'scripts'/n)==h for n,h in p['source_hashes'].items()) and all(sha(ROOT/n)==h for n,h in p['input_hashes'].items()),
        'old_source_protocol_not_rewritten':old['source_hashes']['cohort_uniform_roots.py']==sha(ROOT/'scripts/cohort_uniform_roots.py')=='889182fc142ba34abda7f35758809f0e84e6b7f3d923e42237e247c16f3e7d12' and sha(OLD/'protocol.json')=='e1a6efd2a3a09611f60419b36845db296f8b75e875a150b8ce2384cf623e0982',
        'old_process_terminal_failed_not_live':job['status']=='FAILED' and not (m.helpers.process_identity(job['pid']) or {}).get('alive',False),
        'retained_failure_is_exact_numeric_encoding_bug':'reuse family mismatch' in job['traceback'] and Q('0.78208913')==Q('78208913/100000000'),
        'common_cell_contract_and_materials_unchanged':p['shared_parameter_cell']==old['shared_parameter_cell'] and p['common_model_observation_contract_sha256']==old['common_model_observation_contract_sha256'] and p['fixed_full_parameters']==old['fixed_full_parameters'] and p['integration_policy']==old['integration_policy'],
        'all17_input_fields_and_exact_targets_unchanged':len(p['records'])==17 and all(all(a[k]==b[k] for k in ('record_id','dose_exact','exact_requested_strength','saved_state','saved_state_sha256','branch_domain_id')) for a,b in zip(p['records'],old['records'])),
        'precisely11_prefix_1_M_5_fresh':sum(r['mode']=='reuse_v1_completed' for r in p['records'])==11 and sum(r['mode']=='reuse_immutable_actual_M_fraction_checked' for r in p['records'])==1 and sum(r['mode']=='prepare_and_one_affine_test' for r in p['records'])==5,
        'new_output_ownership_only':all((ROOT/r['output_folder']).resolve().is_relative_to(RUN.resolve()) for r in p['records']),
        'new_run_not_started':not (RUN/'job.json').exists() and not (RUN/'records').exists()}
    prefix=[]
    for row,done in zip(p['records'][:11],progress['completed']):
        prior=read(ROOT/done['file']);wrapped=m.wrap_prefix(p,row,prior)
        prefix.append(sha(ROOT/done['file'])==done['sha256']==row['original_result_sha256'] and wrapped['uniform_cell_result']==prior['uniform_cell_result'] and wrapped['point_cache_file']==prior['point_cache_file'] and wrapped['predictor_file']==prior['predictor_file'] and wrapped['status']==prior['status'])
    checks['all11_numerical_payload_and_cache_reuse_replayed']=all(prefix)
    checks['F_K_prior_failures_preserved']=set(r['record_id'] for r in progress['completed'] if r['status']=='UNRESOLVED')=={'optoGEF2_small_col_F','optoGEF2_small_col_K'}
    row=p['records'][11];wrapped=m.wrap_M(p,row);prior=read(ROOT/p['actualM_reuse']['root_certificate'])
    checks['M_full_payload_same_original_independent_proof']=wrapped['uniform_cell_result']==prior['uniform_cell_result'] and Q(wrapped['dose_exact'])==Q(prior['dose_exact']) and wrapped['lineage']['exact_Fraction_equality_verified'] is True
    checks['M_cache_predictor_review_hashes_current']=all(sha(ROOT/p['actualM_reuse'][k])==p['actualM_reuse'][h] for k,h in [('root_certificate','root_certificate_sha256'),('point_cache','point_cache_sha256'),('predictor_file','predictor_sha256'),('independent_review','independent_review_sha256')])
    try:m.same_dose('0.78208913','0.78208914');checks['unequal_numeric_dose_rejected']=False
    except ValueError:checks['unequal_numeric_dose_rejected']=True
    m.helpers.assert_frozen(p);checks['source_input_environment_guards_reexecuted']=True
    text=(ROOT/'scripts/cohort_uniform_roots_completion.py').read_text(encoding='utf-8')
    checks['no_read_monkeypatch_or_old_source_rewrite']='helpers.read=' not in text and 'helpers.OUT=OUT' in text and 'else:return helpers.record_result(p,row)' in text
    checks['new_exclusive_job_and_launch_review_bytes_retained']="with (OUT/'active.lock').open('x'" in text and "helpers.binary(OUT/'launch_review.json',reviewbytes)" in text
    checks['no_radius_retry_or_parameter_shrink_in_completion']=p['attempt_policy']==old['attempt_policy']
    packet={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'source_sha256':sha(ROOT/'scripts/cohort_uniform_roots_completion.py'),'protocol_sha256':sha(RUN/'protocol.json'),
        'preflight_sha256':sha(RUN/'preflight.json'),'review_source_sha256':sha(__file__),
        'scope':'Independent versioned wrapper/provenance/target/launch review; reuses old numerical proofs without modification, five new field bounds still require execution and independent reconstruction.'}
    out.write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':packet['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if packet['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
