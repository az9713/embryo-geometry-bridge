"""Terminal-only producer comparison, provenance and visual receipts."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
COHORT=ROOT/'runs/material_cohort_descent_v1'
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
vsha=lambda v:hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest()
job=read(OUT/'job.json')
if job['status']!='FINISHED':raise RuntimeError('Complete direct job required; partial results cannot generate final comparison')
results=read(OUT/'results.json');p=read(OUT/'protocol.json');cohort=read(COHORT/'results.json')
baseline_rows={}
for record in cohort['records']:
    baseline_path=ROOT/record['file']
    if sha(baseline_path)!=record['sha256']:raise RuntimeError('Baseline observer artifact hash changed')
    baseline_rows[record['record_id']]=read(baseline_path)['baseline_endpoint_intervals']
exact=read(COHORT/'exact_extraction_v2/exact_candidates.json')
primary=next(c for c in exact['candidates'] if c['variant']==p['primary_variant'])
step=list(map(Q,primary['exact_affine_result']['step_exact']))
checks={'complete17':len(results['records'])==job['records_verified']==17,
    'four_variants_same_new_primary_fields':len(results['variants'])==4 and len({v['variant'] for v in results['variants']})==4,
    'protocol_bound':results['protocol_sha256']==sha(OUT/'protocol.json'),
    'all21_source_hashes_current':all(sha(ROOT/'scripts'/n)==h for n,h in results['source_hashes'].items()),
    'no_original_goal_overclaim':results['original_five_gate_goal_achieved'] is False,
    'v2_candidate_bound':sha(COHORT/'exact_extraction_v2/exact_candidates.json')==p['candidate_file_sha256'],
    'cohort_result_bound':sha(COHORT/'results.json')==p['cohort_results_sha256']}
launch=read(OUT/'launch_job.json');launch_review=read(OUT/'launch_review_exact.json');current_review=read(ROOT/'evidence/material_cohort_direct_prelaunch_review.json')
checks['exact_launch_review_snapshot_bound']=sha(OUT/'launch_review_exact.json')==launch['prelaunch_review_sha256']
checks['final_documentary_review_bound']=results['prelaunch_review_sha256']==sha(ROOT/'evidence/material_cohort_direct_prelaunch_review.json')
checks['review_change_only_archive_metadata']=set(k for k in set(launch_review)|set(current_review) if launch_review.get(k)!=current_review.get(k))=={'prelaunch_revision_archive_hashes'}
checks['same_launch_driver_protocol_bindings']=launch_review['driver_sha256']==sha(ROOT/'scripts/material_cohort_direct.py') and launch_review['protocol_sha256']==sha(OUT/'protocol.json')
zero=read(OUT/'zero/final.json');zv=read(OUT/'zero/verification.json')
checks['new_material_zero_verified']=zero['parameters']==p['new_parameters'] and zero['activation_strength']==0 and zero['converged'] and zv['Q8_force_linf']<1e-7 and zv['Q12_cross_evaluation']['scaled_force_linf']<1e-7 and zv['geometry']['global_embedding_certified'] and zv['geometry']['all_exact_J_positive']
previous=zero;records=[];accepted=rejected=0;triangles=len(zv['geometry']['triangle_certificates']);q8=[];q12=[];jlo=[];jhi=[]
for r,request in zip(results['records'],p['records']):
    f=ROOT/r['file'];s=read(f);vp=f.with_name('verification.json');v=read(vp);g=v['geometry']
    local={'field_hash':sha(f)==r['sha256'],'verification_hash':sha(vp)==r['verification_sha256'],
        'new_material_fixed':s['parameters']==p['new_parameters'],
        'record_dose_target_binding':r['record_id']==request['record_id'] and r['dose_exact']==request['dose_exact'] and Q(r['requested_strength_exact'])==Q(p['new_gain_exact'])*Q(request['dose_exact']),
        'actual_binary64_strength_bound':Q.from_float(s['activation_strength'])==Q(r['solver_input_strength_exact'])==Q(request['solver_input_strength_exact']),
        'requested_actual_rounding_retained':Q(r['strength_rounding_difference_exact'])==Q(r['solver_input_strength_exact'])-Q(r['requested_strength_exact']),
        'Q8_Q12_gates':s['converged'] and s['reached_target'] and v['Q8_force_linf']<1e-7 and v['Q12_cross_evaluation']['scaled_force_linf']<1e-7,
        'exact_geometry_gates':g['all_exact_J_positive'] and g['global_embedding_certified'] and Q(g['J_lower_exact'])>0 and all(t['positive_J_certified'] for t in g['triangle_certificates'])}
    last=previous;chain=True
    for a in s['robust_continuation']['attempts']:
        chain=chain and a['from_strength']==last['activation_strength'] and a['initial_verified_vector_sha256']==vsha(last['free_displacement'])
        if 'result_path' in a:
            path=Path(a['result_path']);trial=read(path if path.is_absolute() else ROOT.parent/path)
            chain=chain and a['final_vector_sha256']==vsha(trial['free_displacement'])
            if a['accepted']:last=trial
        if a['accepted']:accepted+=1
        else:rejected+=1
    local['last_verified_only_warm_chain']=chain and not s['robust_continuation']['failed_state_reuse'] and vsha(last['free_displacement'])==vsha(s['free_displacement']) and r['previous_verified_strength']==previous['activation_strength']
    records.append({'record_id':r['record_id'],'checks':local});previous=s
    triangles+=len(g['triangle_certificates']);q8.append(v['Q8_force_linf']);q12.append(v['Q12_cross_evaluation']['scaled_force_linf']);jlo.append(Q(g['J_lower_exact']));jhi.append(Q(g['J_upper_exact']))
checks['all17_record_receipts_pass']=all(all(r['checks'].values()) for r in records)
checks['all17_source_ids_unique']=len({r['record_id'] for r in records})==17
summaries=[]
for variant in results['variants']:
    name=variant['variant'];aff=next(c['affine_LP'] for c in cohort['candidates'] if c['variant']==name)
    exact_affine=[Q(r)+sum(Q(a)*x for a,x in zip(row,step)) for r,row in zip(aff['affine_residual_intercepts_exact'],aff['affine_derivative_coefficients_exact'])]
    affine_labels=next(c['all68_constraint_order'] for c in cohort['candidates'] if c['variant']==name)
    amap={(r['record_id'],r['observable']):a for r,a in zip(affine_labels,exact_affine)}
    lo=hi=Q(0);baseline_lo=baseline_hi=Q(0);errors=[]
    for r in variant['all68_residuals']['constraint_rows']:
        a,b=map(Q,r['prediction_interval_exact']);y=Q(r['observed_exact']);minus=max(a-y,y-b,Q(0));plus=max(abs(a-y),abs(b-y))
        assert minus==Q(r['minimum_interval_residual_exact']) and plus==Q(r['maximum_interval_residual_exact'])
        lo=max(lo,minus);hi=max(hi,plus);prediction=y+amap[r['record_id'],r['observable']]
        old_a,old_b=map(Q,baseline_rows[r['record_id']][name][r['observable']])
        baseline_lo=max(baseline_lo,old_a-y,y-old_b,Q(0))
        baseline_hi=max(baseline_hi,abs(old_a-y),abs(old_b-y))
        error=[max(a-prediction,prediction-b,Q(0)),max(abs(a-prediction),abs(b-prediction))]
        errors.append({'record_id':r['record_id'],'observable':r['observable'],
            'primary_step_affine_prediction_exact':str(prediction),
            'actual_observer_interval_exact':r['prediction_interval_exact'],
            'absolute_affine_direct_error_interval_exact':list(map(str,error))})
    assert [lo,hi]==list(map(Q,variant['candidate_discrepancy_interval_exact']))
    affloss=max(map(abs,exact_affine))
    summaries.append({'variant':name,'same_baseline_loss_um':aff['baseline_loss_um'],
        'primary_step_affine_loss_exact':str(affloss),'primary_step_affine_loss_um':float(affloss),
        'same_baseline_discrepancy_interval_exact':[str(baseline_lo),str(baseline_hi)],
        'same_baseline_discrepancy_interval_um':[float(baseline_lo),float(baseline_hi)],
        'direct_discrepancy_interval_exact':[str(lo),str(hi)],'direct_discrepancy_interval_um':[float(lo),float(hi)],
        'baseline_minus_direct_improvement_interval_exact':[str(baseline_lo-hi),str(baseline_hi-lo)],
        'baseline_minus_direct_improvement_interval_um':[float(baseline_lo-hi),float(baseline_hi-lo)],
        'all68_primary_affine_direct_errors':errors,
        'maximum_sampled_primary_affine_direct_error_um':float(max(Q(e['absolute_affine_direct_error_interval_exact'][1]) for e in errors)),
        'scope':'Same ONE primary material proposal across all variants, not all4 variant-specific LP steps. Improvement includes possible loading-history effects.'})
checks['all272_residual_and_error_entries']=len(summaries)==4 and all(len(s['all68_primary_affine_direct_errors'])==68 for s in summaries)
receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'per_record_checks':records,
    'counts':{'new_nonzero_fields':17,'new_material_zero_fields':1,'accepted_nonzero_attempts':accepted,'rejected_nonzero_attempts':rejected,'exact_triangle_certificates_including_zero':triangles,'source_observations':68,'retained_operator_variants':4},
    'max_Q8_force_linf':max(q8),'max_Q12_force_linf':max(q12),
    'nonzero_field_exact_J_enclosure_display':[float(min(jlo)),float(max(jhi))],
    'nonzero_field_exact_J_enclosure_exact':[str(min(jlo)),str(max(jhi))],
    'variants':summaries,'hashes':{x.name:sha(x) for x in [Path(__file__),OUT/'results.json',OUT/'protocol.json',OUT/'launch_review_exact.json',ROOT/'evidence/material_cohort_direct_prelaunch_review.json']},
    'scope':'Producer receipt bookkeeping and sampled affine/direct comparison. Independent field replay is separate. No uniform nonlinear/continuum error, physical identified set, causal material attribution or novelty is established.'}
(OUT/'verification_and_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n')
if not all(checks.values()):raise RuntimeError('Final direct producer verification failed')
plt.style.use('dark_background');fig,ax=plt.subplots(figsize=(12,6),layout='constrained');x=np.arange(4)
base=[v['same_baseline_loss_um'] for v in summaries];aff=[v['primary_step_affine_loss_um'] for v in summaries];direct=[v['direct_discrepancy_interval_um'][1] for v in summaries]
for offset,values,color,label in [(-.25,base,'#5cc8ff','verified old common baseline'),(0,aff,'#ffb86b','AFFINE prediction of primary step'),(.25,direct,'#50fa7b','new directly equilibrated primary step')]:
    ax.bar(x+offset,values,.25,color=color,label=label)
    for i,value in enumerate(values):ax.text(i+offset,value+.035,f'{value:.4f}',ha='center',fontsize=9)
ax.set_xticks(x,['Disk / spatial max','Disk / regional mean','Line / spatial max','Line / regional mean'],fontsize=10)
ax.set_ylabel('uniform discrepancy requirement (µm)');ax.grid(axis='y',alpha=.15);ax.set_ylim(0,max(base+aff+direct)+1)
ax.legend(fontsize=9);ax.set_title('One primary material proposal: all 17 new nonlinear fields\nSame 68 observations, all 4 operators; different continuation paths prevent isolated causal attribution',fontsize=13)
figure=OUT/'material-primary-step-direct-comparison.png';fig.savefig(figure,dpi=170);plt.close(fig)
(OUT/'visual_verification.json').write_text(json.dumps({'generator_sha256':sha(__file__),'results_sha256':sha(OUT/'results.json'),'comparison_sha256':sha(OUT/'verification_and_comparison.json'),'figure_sha256':sha(figure),'scope':'Completed primary-step comparison; no invented experimental or physical uncertainty bars'},indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'counts':receipt['counts'],'variants':[{k:v for k,v in s.items() if k!='all68_primary_affine_direct_errors'} for s in summaries]},indent=2))
