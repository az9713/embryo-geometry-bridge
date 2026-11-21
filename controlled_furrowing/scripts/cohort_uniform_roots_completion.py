"""Versioned completion after the v1 dose-string guard stopped at record M.

No monkeypatch, old source edit, old run mutation, repeated prefix Hessian,
parameter-cell change, physical optimization or unresolved-radius retry.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,copy,hashlib,json,os,time,traceback
import cohort_uniform_roots as helpers

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'runs/cohort_uniform_roots_v1'
OUT=ROOT/'runs/cohort_uniform_roots_v2'
helpers.OUT=OUT  # only output ownership/context of existing explicit helpers


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def relative(p):return Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
def payload_hash(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def common(p,row):
    return {'record_id':row['record_id'],'protocol_sha256':sha(OUT/'protocol.json'),'parameter_cell':p['shared_parameter_cell'],'parameter_cell_sha256':p['shared_parameter_cell_sha256'],
        'model_observation_contract_sha256':p['common_model_observation_contract_sha256'],'quantifier':p['quantifier'],'input_kind':'source_record_dose',
        'dose_exact':row['dose_exact'],'strength_relation':'common_gain_times_source_dose','base_strength_exact':row['exact_requested_strength'],
        'saved_state':row['saved_state'],'saved_state_sha256':row['saved_state_sha256'],'branch_scope':'selected_unique_root_in_declared_domain',
        'branch_domain_id':row['branch_domain_id'],'branch_law_semantics':p['branch_law_semantics'],'common_path_or_energetic_stability_certified':False,
        'source_hashes':p['source_hashes'],'no_new_equilibrium_optimization':True,'mode':row['mode']}


def same_dose(a,b):
    if Q(a)!=Q(b):raise ValueError('source doses differ as exact rationals')
    return True


def wrap_prefix(p,row,old):
    same_dose(old['dose_exact'],row['dose_exact'])
    if old['parameter_cell_sha256']!=p['shared_parameter_cell_sha256'] or helpers.canonical_cell(old['parameter_cell'])!=p['shared_parameter_cell'] or old['saved_state_sha256']!=row['saved_state_sha256']:raise RuntimeError('prefix family changed')
    if old['status'] not in ('CERTIFIED_UNIFORM_LOCAL_ROOT','UNRESOLVED') or old['source_hashes']!=p['original_source_hashes']:raise RuntimeError('invalid prefix result')
    answer={**common(p,row),'status':old['status'],'uniform_cell_result':copy.deepcopy(old['uniform_cell_result']),
        'point_cache_file':old['point_cache_file'],'predictor_file':old['predictor_file'],'free_dofs':old['free_dofs'],'local_terms':old['local_terms'],
        'artifact_hashes':{**old['artifact_hashes'],row['original_result_file']:row['original_result_sha256']},
        'lineage':{'kind':'IMMUTABLE_V1_COMPLETED_RESULT','original_result_file':row['original_result_file'],'original_result_sha256':row['original_result_sha256'],
            'original_protocol_sha256':p['original_protocol_sha256'],'original_source_sha256':p['original_source_hashes']['cohort_uniform_roots.py'],
            'numerical_payload_sha256':payload_hash(old['uniform_cell_result']),'no_Hessian_or_root_replay':True}}
    if 'point_defect_Z0' in old:answer['point_defect_Z0']=copy.deepcopy(old['point_defect_Z0'])
    if payload_hash(answer['uniform_cell_result'])!=answer['lineage']['numerical_payload_sha256']:raise RuntimeError('prefix numerical payload changed')
    return answer


def wrap_M(p,row):
    binding=p['actualM_reuse'];old=read(ROOT/binding['root_certificate']);same_dose(old['dose_exact'],row['dose_exact'])
    for key,hkey in [('root_certificate','root_certificate_sha256'),('point_cache','point_cache_sha256'),('predictor_file','predictor_sha256'),('independent_review','independent_review_sha256')]:
        if sha(ROOT/binding[key])!=binding[hkey]:raise RuntimeError('immutable actualM cache/predictor/review lineage changed')
    if sha(ROOT/binding['root_certificate'])!=binding['root_certificate_sha256'] or helpers.cell_hash(old['parameter_cell'])!=p['shared_parameter_cell_sha256']:raise RuntimeError('immutable actualM family changed')
    answer={**common(p,row),'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','uniform_cell_result':copy.deepcopy(old['uniform_cell_result']),
        'point_cache_file':binding['point_cache'],'predictor_file':binding['predictor_file'],'free_dofs':611,'local_terms':9504,
        'artifact_hashes':{binding[k]:binding[h] for k,h in [('root_certificate','root_certificate_sha256'),('point_cache','point_cache_sha256'),('predictor_file','predictor_sha256'),('independent_review','independent_review_sha256')]},
        'reuse_bindings':binding,'lineage':{'kind':'IMMUTABLE_APPROVED_ACTUAL_M','original_dose_spelling':old['dose_exact'],'cohort_dose_spelling':row['dose_exact'],
            'exact_Fraction_equality_verified':True,'numerical_payload_sha256':payload_hash(old['uniform_cell_result']),'original_certificate_model_contract_sha256':old['model_observation_contract_sha256'],
            'common_cohort_model_contract_sha256':p['common_model_observation_contract_sha256'],'same_mechanical_target_new_common_contract_wrapper':True,'no_Hessian_or_root_replay':True}}
    return answer


def freeze():
    if (OUT/'protocol.json').exists():raise RuntimeError('v2 protocol exists; preserve it')
    oldp=read(OLD/'protocol.json');job=read(OLD/'job.json');progress=read(OLD/'progress.json');identity=helpers.process_identity(job['pid'])
    if job['status']!='FAILED' or (identity and identity['alive'] and identity['creation_ticks']==job['creation_ticks']):raise RuntimeError('original job must be failed and dead')
    if len(progress['completed'])!=11:raise RuntimeError('exact eleven-record immutable prefix required')
    if any(sha(ROOT/'scripts'/k)!=v for group in ('source_hashes','common_contract_source_hashes') for k,v in oldp[group].items()) or any(sha(ROOT/k)!=v for k,v in oldp['input_hashes'].items()) or helpers.versions()!=oldp['environment_versions']:raise RuntimeError('old source/input/environment changed')
    p=copy.deepcopy(oldp);p['version']=2;p['created_utc']=helpers.stamp();p['scope']='Complete all17 shared-cell records with immutable11-record v1 prefix reuse, Fraction-checked M reuse and only5 remaining fresh calculations'
    p['original_protocol_file']=relative(OLD/'protocol.json');p['original_protocol_sha256']=sha(OLD/'protocol.json');p['original_source_hashes']=copy.deepcopy(oldp['source_hashes'])
    p['source_hashes']['cohort_uniform_roots_completion.py']=sha(__file__)
    p['input_hashes'].update({relative(x):sha(x) for x in [OLD/'protocol.json',OLD/'job.json',OLD/'progress.json',OLD/'prelaunch_independent_review.json']})
    p['recovery_reason']='v1 compared decimal and reduced-fraction dose spellings literally at M; version2 uses exact Fraction equality. No numerical result or physical target changes.'
    for i,row in enumerate(p['records']):
        row['output_folder']=relative(OUT/'records'/f'{i:02d}_{row["record_id"]}')
        if i<11:
            record=progress['completed'][i]
            if record['record_id']!=row['record_id'] or sha(ROOT/record['file'])!=record['sha256']:raise RuntimeError('prefix order/hash changed')
            old=read(ROOT/record['file']);row.update(mode='reuse_v1_completed',original_result_file=record['file'],original_result_sha256=record['sha256'])
            p['input_hashes'][record['file']]=record['sha256'];p['input_hashes'].update(old['artifact_hashes'])
        elif i==11:row['mode']='reuse_immutable_actual_M_fraction_checked'
        else:row['mode']='prepare_and_one_affine_test'
    helpers.immutable(OUT/'protocol.json',p);preflight(p)
    print(json.dumps({'status':'FROZEN_V2_PREFLIGHT_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)


def preflight(p):
    helpers.assert_frozen(p);checks={'unchanged_shared_cell_and_common_contract':p['shared_parameter_cell_sha256']==read(OLD/'protocol.json')['shared_parameter_cell_sha256'] and p['common_model_observation_contract_sha256']==read(OLD/'protocol.json')['common_model_observation_contract_sha256'],
        'exact11_prefix_1_M_5_remaining':sum(r['mode']=='reuse_v1_completed' for r in p['records'])==11 and sum(r['mode']=='reuse_immutable_actual_M_fraction_checked' for r in p['records'])==1 and sum(r['mode']=='prepare_and_one_affine_test' for r in p['records'])==5,
        'all_original_fields_and_exact_targets_unchanged':all(a['saved_state_sha256']==b['saved_state_sha256'] and a['exact_requested_strength']==b['exact_requested_strength'] and Q(a['dose_exact'])==Q(b['dose_exact']) for a,b in zip(p['records'],read(OLD/'protocol.json')['records']))}
    same_dose('0.78208913','78208913/100000000');checks['decimal_fraction_equal_dose_accepted']=True
    try:same_dose('0.78208913','78208914/100000000');checks['unequal_dose_rejected']=False
    except ValueError:checks['unequal_dose_rejected']=True
    payloads=True;unresolved=[]
    for row in p['records']:
        if row['mode']=='reuse_v1_completed':
            old=read(ROOT/row['original_result_file']);new=wrap_prefix(p,row,old);payloads &= new['uniform_cell_result']==old['uniform_cell_result']
            if new['status']=='UNRESOLVED':unresolved.append(row['record_id'])
    checks['all11_numerical_payloads_preserved_exactly']=payloads;checks['F_K_UNRESOLVED_preserved_without_retest']=set(unresolved)=={'optoGEF2_small_col_F','optoGEF2_small_col_K'}
    m=next(r for r in p['records'] if r['mode']=='reuse_immutable_actual_M_fraction_checked');mr=wrap_M(p,m);checks['M_approved_payload_preserved_with_numeric_dose_check']=mr['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' and mr['lineage']['exact_Fraction_equality_verified']
    tests=helpers.guard_tests();checks['old_metadata_and_path_guards_still_pass']=all(tests.values());identity=helpers.process_identity(os.getpid());checks['actual_PID_creation_identity_readable']=bool(identity and identity['alive'])
    helpers.immutable(OUT/'guard_tests.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':{**tests,**checks},'source_sha256':sha(__file__),'no_model_or_hessian_work':True})
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'guard_tests_sha256':sha(OUT/'guard_tests.json'),
        'new_model_or_Hessian_calculations':0,'original_v1_files_mutated':False,'launch_requires':'new prelaunch_independent_review.json binding exact v2 source/protocol/preflight'}
    helpers.immutable(OUT/'preflight.json',receipt)
    if not all(checks.values()):raise RuntimeError('v2 preflight failed')


def result_record(p,row):
    path=ROOT/row['output_folder']/'result.json'
    if path.exists():
        answer=read(path)
        if answer['protocol_sha256']!=sha(OUT/'protocol.json') or answer['parameter_cell_sha256']!=p['shared_parameter_cell_sha256'] or any(sha(ROOT/k)!=v for k,v in answer.get('artifact_hashes',{}).items()):raise RuntimeError('existing v2 record changed')
        return answer
    if row['mode']=='reuse_v1_completed':answer=wrap_prefix(p,row,read(ROOT/row['original_result_file']))
    elif row['mode']=='reuse_immutable_actual_M_fraction_checked':answer=wrap_M(p,row)
    else:return helpers.record_result(p,row)
    helpers.immutable(path,answer);return answer


def run():
    p=read(OUT/'protocol.json');helpers.assert_frozen(p);reviewbytes=(OUT/'prelaunch_independent_review.json').read_bytes();review=json.loads(reviewbytes)
    if review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()) or review.get('source_sha256')!=sha(__file__) or review.get('protocol_sha256')!=sha(OUT/'protocol.json') or review.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise RuntimeError('bound independent v2 launch review required')
    if (OUT/'job.json').exists() or (OUT/'results.json').exists():raise RuntimeError('v2 job exists; preserve/no duplicate')
    if any(os.environ.get(k)!=v for k,v in p['required_threads'].items() if k!='torch_num_threads'):raise RuntimeError('singlethread policy required')
    identity=helpers.process_identity(os.getpid());job={'status':'RUNNING',**identity,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':hashlib.sha256(reviewbytes).hexdigest(),'started_utc':helpers.stamp()}
    with (OUT/'active.lock').open('x',encoding='utf8') as f:json.dump(job,f)
    helpers.binary(OUT/'launch_review.json',reviewbytes);helpers.atomic(OUT/'job.json',job);started=time.time();completed=[]
    try:
        for row in p['records']:
            helpers.assert_frozen(p);answer=result_record(p,row);path=ROOT/row['output_folder']/'result.json';completed.append({'record_id':row['record_id'],'status':answer['status'],'file':relative(path),'sha256':sha(path),'mode':row['mode'],'provider_protocol_file':relative(OUT/'protocol.json'),'provider_protocol_sha256':sha(OUT/'protocol.json')})
            helpers.atomic(OUT/'progress.json',{'status':'RUNNING','pid':os.getpid(),'completed':completed,'total':17,'elapsed_seconds':time.time()-started})
            print(json.dumps({'record_id':row['record_id'],'status':answer['status'],'mode':row['mode'],'completed':len(completed),'total':17,'elapsed_seconds':time.time()-started}),flush=True)
        helpers.assert_frozen(p);certified=sum(x['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' for x in completed)
        summary={'status':'CERTIFIED_SHARED_CELL_STATIONARY_ROOT_TUPLE' if certified==17 else 'COMPLETE_WITH_UNRESOLVED_ROOT_RECORDS','protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),
            'records':completed,'shared_parameter_cell':p['shared_parameter_cell'],'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],
            'common_model_observation_contract_sha256':p['common_model_observation_contract_sha256'],'all17_certified':certified==17,'certified':certified,'unresolved':17-certified,
            'branch_law_semantics':p['branch_law_semantics'],'quantified_conclusion':'All17 shared-cell roots certified.' if certified==17 else 'No all17 root tuple certificate; unresolved records retained.',
            'original_v1_immutable_protocol_sha256':p['original_protocol_sha256'],'no_equilibrium_optimization':True,'no_OUT_classifications':True,
            'elapsed_seconds':time.time()-started,'completed_utc':helpers.stamp(),'not_claimed':p['not_claimed']}
        helpers.immutable(OUT/'results.json',summary);job.update(status='FINISHED',finished_utc=helpers.stamp(),results_sha256=sha(OUT/'results.json'));helpers.atomic(OUT/'job.json',job)
        print(json.dumps({'status':summary['status'],'certified':certified,'unresolved':17-certified,'results_sha256':sha(OUT/'results.json')}),flush=True)
    except Exception:
        job.update(status='FAILED',traceback=traceback.format_exc(),failed_utc=helpers.stamp());helpers.atomic(OUT/'job.json',job);raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','run']);args=parser.parse_args()
    freeze() if args.stage=='freeze' else run()
