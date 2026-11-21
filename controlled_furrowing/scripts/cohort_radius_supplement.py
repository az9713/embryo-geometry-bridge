"""Retain failed first tests; retry corrector radii on the identical F/K cell.

Reuses frozen saved-state caches, predictors and reviewed interval algebra.
No optimization, new Hessian preparation, changed parameter cell or OUT claim.
"""
from pathlib import Path
import argparse,copy,hashlib,json,os,time,traceback
import cohort_uniform_roots as helpers

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'runs/cohort_uniform_roots_v1'
OUT=ROOT/'runs/cohort_radius_supplement_v1'
helpers.OUT=OUT
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()

def freeze():
    if (OUT/'protocol.json').exists():raise RuntimeError('preserve frozen supplemental protocol')
    old=read(OLD/'protocol.json');job=read(OLD/'job.json')
    identity=helpers.process_identity(job['pid'])
    if job['status']!='FAILED' or (identity and identity['alive'] and identity['creation_ticks']==job['creation_ticks']):raise RuntimeError('original attempt must be terminal')
    p=copy.deepcopy(old);p.update(version=1,created_utc=helpers.stamp(),scope='Same-cell radius-only supplement for prior F/K UNRESOLVED sufficient tests')
    p['source_hashes']['cohort_radius_supplement.py']=sha(__file__)
    p['input_hashes'].update({rel(x):sha(x) for x in [OLD/'protocol.json',OLD/'job.json',OLD/'prelaunch_independent_review.json']})
    p['original_protocol_file']=rel(OLD/'protocol.json');p['original_protocol_sha256']=sha(OLD/'protocol.json')
    p['attempt_policy']['radius_multipliers']=[2,4,8]
    p['attempt_policy']['radius_rule']='Keep the previous exact shared parameter cell and affine predictor; test r=max(floor,2*validated eta upper) times2,4,8, stopping at the first strict success.'
    p['records']=[]
    for row in old['records']:
        if row['record_id'] not in ('optoGEF2_small_col_F','optoGEF2_small_col_K'):continue
        priorpath=ROOT/row['output_folder']/'result.json';prior=read(priorpath)
        if prior['status']!='UNRESOLVED' or 'uniform_cell_result' not in prior or prior['uniform_cell_result']['status']!='UNRESOLVED':raise RuntimeError('only retained failed sufficient tests allowed')
        if not prior['uniform_cell_result']['radius_attempts'] or any(a.get('strict_selfmap_contraction') for a in prior['uniform_cell_result']['radius_attempts']):raise RuntimeError('prior failed-radius lineage required')
        r=copy.deepcopy(row);r.update(mode='same_cell_radius_only',original_result_file=rel(priorpath),original_result_sha256=sha(priorpath),point_cache_file=prior['point_cache_file'],predictor_file=prior['predictor_file'],output_folder=rel(OUT/'records'/row['record_id']))
        p['records'].append(r);p['input_hashes'][rel(priorpath)]=sha(priorpath);p['input_hashes'].update(prior['artifact_hashes'])
    if len(p['records'])!=2:raise RuntimeError('exactly F and K required')
    helpers.immutable(OUT/'protocol.json',p);helpers.assert_frozen(p)
    checks={'same_parameter_cell':p['shared_parameter_cell']==old['shared_parameter_cell'],
        'same_common_contract':p['common_model_observation_contract_sha256']==old['common_model_observation_contract_sha256'],
        'same_model_parameters_and_geometry':p['fixed_full_parameters']==old['fixed_full_parameters'],
        'same_affine_method_and_exact_floor':p['attempt_policy']['method']=='affine' and p['attempt_policy']['error_radius_floor_exact']==old['attempt_policy']['error_radius_floor_exact'],
        'only_radius_policy_changes':{k:v for k,v in p['attempt_policy'].items() if k not in ('radius_multipliers','radius_rule')}=={k:v for k,v in old['attempt_policy'].items() if k not in ('radius_multipliers','radius_rule')},
        'exact_frozen_multipliers':p['attempt_policy']['radius_multipliers']==[2,4,8],
        'only_F_K_and_prior_failures_retained':{r['record_id'] for r in p['records']}=={'optoGEF2_small_col_F','optoGEF2_small_col_K'},
        'all_cache_predictor_and_source_hashes_current':True,
        'no_parameter_cell_shrink_or_equilibrium_optimization':p['no_new_equilibrium_optimization'] is True,
        'no_output_work_before_review':not (OUT/'records').exists()}
    helpers.immutable(OUT/'preflight.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'scope':'Metadata/cache lineage only; no interval assembly or radius test performed.'})
    print(json.dumps({'status':'PREFLIGHT_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}))
    if not all(checks.values()):raise RuntimeError('preflight failure')

def run():
    p=read(OUT/'protocol.json');helpers.assert_frozen(p)
    reviewpath=OUT/'prelaunch_independent_review.json';raw=reviewpath.read_bytes();review=json.loads(raw)
    if review.get('status')!='PASS' or not review.get('checks') or any(v is not True for v in review['checks'].values()) or review.get('source_sha256')!=sha(__file__) or review.get('protocol_sha256')!=sha(OUT/'protocol.json') or review.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise RuntimeError('concrete independent review required')
    if (OUT/'job.json').exists() or (OUT/'results.json').exists():raise RuntimeError('preserve prior job; no duplicate')
    if any(os.environ.get(k)!=v for k,v in p['required_threads'].items() if k!='torch_num_threads'):raise RuntimeError('thread policy mismatch')
    identity=helpers.process_identity(os.getpid())
    if not identity or not identity['alive']:raise RuntimeError('actual worker identity required')
    job={'status':'RUNNING',**identity,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':hashlib.sha256(raw).hexdigest(),'started_utc':helpers.stamp()}
    with (OUT/'active.lock').open('x',encoding='utf-8') as f:json.dump(job,f)
    helpers.binary(OUT/'launch_review.json',raw);helpers.atomic(OUT/'job.json',job)
    import parametric_root_cells as kernel
    from flint import ctx
    ctx.prec=p['precision_bits'];records=[];start=time.time()
    try:
        for row in p['records']:
            helpers.assert_frozen(p);prior=read(ROOT/row['original_result_file'])
            kernel.OUT=(ROOT/row['predictor_file']).parent;kernel.A.cache_clear()
            local={'point_certificate_file':row['point_cache_file'],'predictor_file':Path(row['predictor_file']).name,
                'dose_exact':row['dose_exact'],'base_parameters_exact':p['base_parameters_exact'],
                'calibration_record_ids_covered':[row['record_id']],'input_kind':'source_record_dose',
                'strength_relation':'common_gain_times_source_dose','branch_domain_id':row['branch_domain_id'],**p['attempt_policy']}
            point=kernel.build_point(local)
            if point['n']!=611 or len(point['cache'])!=9504:raise RuntimeError('cache dimensions changed')
            result=kernel.evaluate_cell(point,local,'unchanged_shared_cell_radius_supplement',p['shared_parameter_cell'],'affine')
            answer=copy.deepcopy(prior)
            answer.update(protocol_sha256=sha(OUT/'protocol.json'),source_hashes=p['source_hashes'],mode='same_cell_radius_only',status=result['status'],uniform_cell_result=result,
                point_defect_Z0=kernel.pack(point['Z0']),artifact_hashes={**prior['artifact_hashes'],row['original_result_file']:row['original_result_sha256']},
                lineage={'kind':'SAME_CELL_RADIUS_ONLY_SUPPLEMENT','original_result_file':row['original_result_file'],'original_result_sha256':row['original_result_sha256'],'original_protocol_sha256':p['original_protocol_sha256'],'parameter_cell_and_predictor_unchanged':True,'no_new_Hessian_preparation_or_equilibrium_solve':True})
            output=ROOT/row['output_folder']/'result.json';helpers.immutable(output,answer)
            records.append({'record_id':row['record_id'],'status':answer['status'],'file':rel(output),'sha256':sha(output),
                'provider_protocol_file':rel(OUT/'protocol.json'),'provider_protocol_sha256':sha(OUT/'protocol.json')})
            print(json.dumps({'record':row['record_id'],'status':answer['status'],'radius_attempts':len(result['radius_attempts']),'seconds':time.time()-start}),flush=True)
        helpers.assert_frozen(p);summary={'status':'COMPLETE_RADIUS_SUPPLEMENT','records':records,'certified':sum(r['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' for r in records),
            'unresolved':sum(r['status']=='UNRESOLVED' for r in records),'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),
            'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],'common_model_observation_contract_sha256':p['common_model_observation_contract_sha256'],
            'no_OUT_classifications':True,'no_cell_shrink':True,'no_new_equilibrium_optimization':True,'elapsed_seconds':time.time()-start}
        helpers.immutable(OUT/'results.json',summary);job.update(status='FINISHED',results_sha256=sha(OUT/'results.json'),finished_utc=helpers.stamp());helpers.atomic(OUT/'job.json',job)
    except Exception:
        job.update(status='FAILED',traceback=traceback.format_exc());helpers.atomic(OUT/'job.json',job);raise

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','run']);args=ap.parse_args()
    freeze() if args.stage=='freeze' else run()
