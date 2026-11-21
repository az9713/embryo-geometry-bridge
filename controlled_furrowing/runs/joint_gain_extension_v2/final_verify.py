"""Receipt-link and exact-arithmetic final bookkeeping; no new solves."""
from pathlib import Path
from fractions import Fraction as Q
import csv,hashlib,json
import numpy as np
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vector_sha(v):return hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest()
job=read(out/'direct_job.json')
if job.get('status')!='FINISHED':raise RuntimeError('Terminal full direct job required')
result=read(out/'direct_results.json');protocol=read(out/'protocol.json');dp=read(out/'direct_protocol.json')
progress=read(out/'direct_progress.json');cand=read(out/'surrogate_candidates.json')
review=read(root/'evidence/piecewise_gain_minimax_v2_independent_review.json')
checks={}
checks['terminal_job_68']=job['targets_verified']==68
checks['complete_result_status']=result['status']=='COMPLETE_ALL_FOUR_CONDITIONAL_CANDIDATE_REPLAYS'
checks['complete_progress']=progress['status']=='COMPLETE' and len(progress['records'])==68
checks['four_unique_variants']=len(result['variants'])==4 and len({v['variant'] for v in result['variants']})==4
checks['all_frozen_sources_unchanged']=all(sha(root/'scripts'/n)==h for n,h in result['source_hashes'].items())
checks['source_sets_match']=result['source_hashes']==job['source_hashes']==dp['source_hashes']==protocol['source_hashes']
checks['protocol_hashes_match']=result['protocol_sha256']==sha(out/'protocol.json')==dp['protocol_sha256']
checks['direct_protocol_bound']=result['direct_protocol_sha256']==sha(out/'direct_protocol.json')
checks['candidate_table_bound']=dp['candidate_file_sha256']==sha(out/'surrogate_candidates.json') and dp['response_table_sha256']==sha(out/'response_table.json')
checks['independent_surrogate_review_passes']=review['status']=='PASS' and all(review['checks'].values()) and dp['independent_review_sha256']==sha(root/'evidence/piecewise_gain_minimax_v2_independent_review.json')
checks['parent_artifacts_preserved']=all(sha(root/'runs/joint_gain_calibration'/n)==h for n,h in protocol['parent_artifact_hashes'].items())
checks['seed_fields_preserved']=all(sha(root/p)==h for p,h in protocol['seed_field_hashes'].items())
checks['data_hash_unchanged']=sha(root/'data/derived/calibration_small_18um.csv')==protocol['calibration_csv_sha256']
with (root/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig') as s:rows=list(csv.DictReader(s))
byid={r['record_id']:r for r in rows};keys=('maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um')
record_checks=[];jlo=[];jhi=[];q8=[];q12=[];triangles=0;accepted=rejected=0
previous=read(root/'runs/joint_gain_calibration/zero/final.json')
for i,(r,request) in enumerate(zip(result['records'],dp['requests'])):
    field=root/r['file'];s=read(field);vpath=field.with_name('verification.json');v=read(vpath)
    local={}
    local['field_and_verification_hash']=sha(field)==r['sha256'] and sha(vpath)==r['verification_sha256']
    local['immutable_parameters']=s['parameters']==protocol['physical_parameters']
    local['target_binding']=r['index']==i and r['requested_strength_exact']==request['strength_exact'] and r['solver_input_strength_exact']==request['actual_solver_strength_exact'] and r['uses']==request['uses']
    local['binary64_target_rounding']=Q(r['solver_input_strength_exact'])==Q.from_float(s['activation_strength']) and Q(r['strength_rounding_difference_exact'])==Q(r['solver_input_strength_exact'])-Q(r['requested_strength_exact'])
    local['convergence_force_gates']=s['converged'] and s['reached_target'] and r['reached'] and v['Q8_force_linf']<1e-7 and v['Q12_cross_evaluation']['scaled_force_linf']<1e-7
    local['stored_force_agreement']=s['scaled_free_force_linf']==r['force_Q8']==v['Q8_force_linf'] and r['force_Q12']==v['Q12_cross_evaluation']['scaled_force_linf']
    g=v['geometry'];local['exact_J_global_embedding']=g['all_exact_J_positive'] and g['global_embedding_certified'] and Q(g['J_lower_exact'])>0 and all(t['positive_J_certified'] for t in g['triangle_certificates'])
    triangles+=len(g['triangle_certificates']);jlo.append(Q(g['J_lower_exact']));jhi.append(Q(g['J_upper_exact']));q8.append(r['force_Q8']);q12.append(r['force_Q12'])
    at=s['robust_continuation']['attempts'];last=previous;chain=True
    for a in at:
        chain=chain and a['initial_verified_vector_sha256']==vector_sha(last['free_displacement']) and a['from_strength']==last['activation_strength']
        ap=Path(a['result_path']);ap=ap if ap.is_absolute() else root.parent/ap
        trial=read(ap);chain=chain and a['final_vector_sha256']==vector_sha(trial['free_displacement'])
        if a['accepted']:last=trial;accepted+=1
        else:rejected+=1
    local['last_verified_only_chain']=chain and not s['robust_continuation']['failed_state_reuse'] and vector_sha(last['free_displacement'])==vector_sha(s['free_displacement']) and r['previous_verified_strength']==previous['activation_strength']
    local['dose_common_gain_binding']=all(Q(r['requested_strength_exact'])==Q(next(c['surrogate_fit']['gain_exact'] for c in cand['variants'] if c['variant']==u['variant']))*Q(byid[u['record_id']]['activation_magnitude']) for u in r['uses'])
    local['attempt_count_binding']=len(at)==r['attempts'] and sum(not a['accepted'] for a in at)==r['rejected_attempts']
    record_checks.append({'index':i,'checks':local});previous=s
checks['all68_records_present']=len(result['records'])==len(dp['requests'])==68
checks['every_record_receipt_passes']=all(all(r['checks'].values()) for r in record_checks)
checks['68_unique_ascending_strengths']=all(Q(a['requested_strength_exact'])<Q(b['requested_strength_exact']) for a,b in zip(result['records'],result['records'][1:]))
variant_summaries=[]
for variant in result['variants']:
    name=variant['variant'];rs=variant['all68_residuals']['constraint_rows'];errors=variant['all68_surrogate_direct_errors']
    expected={(r,k) for r in byid for k in keys}
    checks[name+'_all68_coverage']=len(rs)==68 and {(r['record_id'],r['observable']) for r in rs}==expected and len(errors)==68 and {(r['record_id'],r['observable']) for r in errors}==expected
    lowers=[];uppers=[];maxwidth=Q(0);channel={}
    for r in rs:
        y=Q(byid[r['record_id']][r['observable']]);a,b=map(Q,r['prediction_interval_exact']);lo=max(a-y,y-b,Q(0));hi=max(abs(a-y),abs(b-y))
        assert a<=b and Q(r['observed_exact'])==y and Q(r['minimum_interval_residual_exact'])==lo and Q(r['maximum_interval_residual_exact'])==hi
        lowers.append(lo);uppers.append(hi);maxwidth=max(maxwidth,b-a)
        channel.setdefault(r['observable'],[]).append((lo,hi))
    discrepancy=[max(lowers),max(uppers)]
    checks[name+'_summary_residuals']=list(map(Q,variant['direct_candidate_discrepancy_requirement_um_exact']))==discrepancy
    checks[name+'_surrogate_gain_unchanged']=variant['selected_exact_SURROGATE_gain']==next(v['surrogate_fit']['gain_exact'] for v in cand['variants'] if v['variant']==name)
    checks[name+'_all_error_intervals_consistent']=all((lambda a,b,y,e: Q(e[0])==max(a-y,y-b,Q(0)) and Q(e[1])==max(abs(a-y),abs(b-y)))(*map(Q,e['direct_interval_exact']),Q(e['surrogate_midpoint_prediction_exact']),e['absolute_surrogate_to_direct_error_interval_exact']) for e in errors)
    variant_summaries.append({'variant':name,'gain':float(Q(variant['selected_exact_SURROGATE_gain'])),'surrogate_loss_um':float(Q(variant['surrogate_loss_exact'])),'direct_discrepancy_interval_um':list(map(float,discrepancy)),'max_observer_interval_width_um':float(maxwidth),'maximum_sampled_surrogate_direct_error_um':variant['maximum_surrogate_to_direct_error_upper_um'],'channel_discrepancy_intervals_um':{k:[float(max(a for a,b in values)),float(max(b for a,b in values))] for k,values in channel.items()}})
checks['no_original_goal_overclaim']=result['original_five_gate_goal_achieved'] is False
checks['evidence_copy_byte_identical']=sha(out/'direct_results.json')==sha(root/'evidence/joint_gain_extension_direct.json')
receipt={'status':'PASS' if all(checks.values()) else 'FAIL','scope':'Final receipt-link and exact-arithmetic bookkeeping on all68 direct states. No new optimization; independent physics/observer replay supplied separately by data_audit. Stored discrete certificates and force cross-checks do not certify continuum error or physical observation assumptions.','checks':checks,'check_count':len(checks),'per_record_checks':record_checks,'counts':{'direct_targets':68,'retained_variants':4,'residuals_per_variant':68,'total_residual_entries':272,'triangle_certificates':triangles,'accepted_attempts':accepted,'rejected_attempts':rejected},'max_Q8_force_linf':max(q8),'max_Q12_force_linf':max(q12),'all_stored_geometry_J_enclosure_exact':[str(min(jlo)),str(max(jhi))],'all_stored_geometry_J_enclosure_display':[float(min(jlo)),float(max(jhi))],'variants':variant_summaries,'hashes':{p.name:sha(p) for p in [Path(__file__),out/'direct_results.json',out/'protocol.json',out/'direct_protocol.json',out/'response_table.json',out/'surrogate_candidates.json']}}
(root/'evidence/joint_gain_extension_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ('per_record_checks','checks')},indent=2))
if not all(checks.values()):raise RuntimeError('Final bookkeeping verification failed')
