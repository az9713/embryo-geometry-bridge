"""Uniform roots on ONE shared parameter cell for all17 saved source doses.

Freeze/preflight precede expensive work and require a separate bound launch
review. No equilibrium optimization, physical-source edits, heldout reads,
response fitting, silent per-record cell shrink or OUT classification.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,csv,datetime,hashlib,importlib.metadata,io,json,os,platform,sys,time,traceback

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/cohort_uniform_roots_v1'
DIRECT=ROOT/'runs/material_cohort_direct_v1'
ACTUAL_M=ROOT/'runs/parametric_root_cells_v1/actual_M'
REVIEW=ROOT/'evidence/parametric_actual_M_independent_review.json'
SOURCE_NAMES=('cohort_uniform_roots.py','parametric_root_cells.py','algebraic_root_feasibility.py','material_tangent_inference.py',
 'sparse_quadratic.py','quadratic_quadrature_study.py','robust_quadratic.py','axisymmetric_quadratic.py',
 'axisymmetric_mixed.py','axisymmetric_prestress.py','axisymmetric_mechanics.py')
QUANTIFIER='for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def relative(p):return Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
def safe_path(p,write=False):
    p=Path(p);p=p if p.is_absolute() else ROOT/p;p=p.resolve()
    if not p.is_relative_to(ROOT.resolve()) or (write and not p.is_relative_to(OUT.resolve())):raise ValueError('path outside declared workspace/output ownership')
    return p
def canonical_cell(c):return {k:list(map(str,map(Q,c[k]))) for k in sorted(c)}
def cell_hash(c):return hashlib.sha256(json.dumps(canonical_cell(c),sort_keys=True,separators=(',',':')).encode()).hexdigest()


def atomic(p,value):
    p=safe_path(p,write=True);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.partial')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8');tmp.replace(p)


def immutable(p,value):
    p=Path(p)
    if p.exists():raise RuntimeError('preserve immutable artifact: '+str(p))
    atomic(p,value)


def binary(p,data):
    p=safe_path(p,write=True)
    if p.exists():raise RuntimeError('preserve immutable binary artifact')
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.partial');tmp.write_bytes(data);tmp.replace(p)


def versions():
    return {'python':platform.python_version(),**{k:importlib.metadata.version(k) for k in ('numpy','scipy','torch','python-flint')}}


def freeze(contract_file):
    if (OUT/'protocol.json').exists():raise RuntimeError('protocol already frozen')
    contract_file=safe_path(contract_file)
    contract=read(contract_file);directp=read(DIRECT/'protocol.json');directr=read(DIRECT/'results.json');audit=read(REVIEW)
    oldroot=read(ACTUAL_M/'root_certificate.json');oldp=read(ACTUAL_M/'protocol.json')
    if directr['status']!='COMPLETE_PRIMARY_MATERIAL_STEP_DIRECT_REPLAY' or audit['status']!='PASS':raise RuntimeError('reviewed direct cohort and actualM proof required')
    if audit['root_audit']['status']!='PASS' or audit['root_audit']['certificate_sha256']!=sha(ACTUAL_M/'root_certificate.json') or not all(v is True for v in audit['root_audit']['checks'].values()):raise RuntimeError('actualM root audit not bound')
    with (ROOT/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    if len(rows)!=17 or len({r['record_id'] for r in rows})!=17:raise RuntimeError('exact17 source-column records required')
    records=[];gain=Q(directp['new_gain_exact'])
    for i,row in enumerate(rows):
        rec=next(r for r in directr['records'] if r['record_id']==row['record_id']);statepath=safe_path(rec['file']);state=read(statepath);dose=Q(row['activation_magnitude'])
        if sha(statepath)!=rec['sha256'] or state['parameters']!=directp['new_parameters'] or dose!=Q(rec['dose_exact']):raise RuntimeError('new-material source record mismatch')
        if Q(rec['requested_strength_exact'])!=gain*dose or Q.from_float(float(state['activation_strength']))!=Q(rec['solver_input_strength_exact']):raise RuntimeError('exact source-strength law mismatch')
        if not state['converged'] or len(state['free_displacement'])!=611 or len(state['triangles'])!=144:raise RuntimeError('expected archived field dimensions/convergence required')
        records.append({'index':i,'record_id':row['record_id'],'dose_exact':str(dose),'base_gain_exact':str(gain),
            'exact_requested_strength':str(gain*dose),'solver_strength_binary_exact':rec['solver_input_strength_exact'],
            'saved_state':rec['file'],'saved_state_sha256':rec['sha256'],
            'mode':'reuse_immutable_actual_M' if row['record_id']=='optoGEF2_small_col_M' else 'prepare_and_one_affine_test',
            'output_folder':relative(OUT/'records'/f'{i:02d}_{row["record_id"]}'),
            'branch_domain_id':oldroot['branch_domain_id'] if row['record_id']=='optoGEF2_small_col_M' else 'cohort_shared_cell_v1:'+row['record_id']})
    C=canonical_cell(oldroot['parameter_cell'])
    if cell_hash(C)!=audit['root_audit']['parameter_cell_sha256']:raise RuntimeError('same shared actualM cell required')
    if canonical_cell(contract['parameter_cell'])!=C or contract['parameter_cell_sha256']!=cell_hash(C) or contract['branch_law_semantics']!='independent_stationary_root_per_record':raise RuntimeError('common contract does not describe this shared family')
    if any(directp['new_parameters'][k]!=v for k,v in contract['fixed_full_parameters'].items()):raise RuntimeError('common fixed geometry/material differs')
    if any(sha(ROOT/'scripts'/n)!=v for n,v in contract['source_kernel_hashes'].items()):raise RuntimeError('common contract dependency changed')
    if sha(ROOT/contract['calibration_source_file'])!=contract['calibration_source_sha256']:raise RuntimeError('common contract data changed')
    if contract['primary_variant']!='current_disk__spatial_max' or set(contract['endpoint_definitions'])!={'maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um'} or contract['source_temporal_semantics_validated'] is not False or contract['common_path_or_energetic_stability_certified'] is not False:raise RuntimeError('common observer semantics differ from approved conditional primary contract')
    inputs=[DIRECT/'protocol.json',DIRECT/'results.json',ROOT/'evidence/material_cohort_direct_independent_review.json',
        ROOT/'data/derived/calibration_small_18um.csv',REVIEW,ACTUAL_M/'root_certificate.json',ACTUAL_M/'point_cache.json',
        ACTUAL_M/'protocol.json',ACTUAL_M/'affine_predictor.npz',contract_file]
    inputs.extend(ROOT/r['saved_state'] for r in records)
    p={'version':1,'created_utc':stamp(),'scope':'Uniform stationary roots for all17 saved new-material source inputs over ONE unchanged shared five-parameter cell',
       'source_hashes':{n:sha(ROOT/'scripts'/n) for n in SOURCE_NAMES},'input_hashes':{relative(x):sha(x) for x in inputs},
       'shared_parameter_cell':C,'shared_parameter_cell_sha256':cell_hash(C),'common_model_observation_contract_file':relative(contract_file),
       'common_model_observation_contract_sha256':sha(contract_file),'base_parameters_exact':oldp['base_parameters_exact'],
       'common_contract_source_hashes':contract['source_kernel_hashes'],
       'fixed_full_parameters':directp['new_parameters'],'record_order':'canonical training CSV order B:R; no nonlinear warm-chain operations','records':records,
       'integration_policy':directp['integration_policy'],'precision_bits':160,'expected_free_dofs':611,'expected_local_terms':9504,
       'attempt_policy':{'method':'affine','exactly_one_shared_cell_per_record':True,'radius_multipliers':[1],
            'error_radius_floor_exact':'1/100000000000','maximum_test_radius_exact':'1/100',
            'radius_rule':'max(floor,2*validated eta upper); ONE radius test; never shrink/change the parameter cell',
            'failed_sufficient_test':'UNRESOLVED, never OUT; retained terminal record; continue other records',
            'integrity_failure':'stop job FAILED; preserve completed records; never consume stale source/cache',
            'restart_policy':'skip completed hash-bound record results; reuse completed prepared manifests; new preparation-attempt folder after interrupted preparation'},
       'branch_law_semantics':'independent_stationary_root_per_record','quantifier':QUANTIFIER,
       'actualM_reuse':{'root_certificate':relative(ACTUAL_M/'root_certificate.json'),'root_certificate_sha256':sha(ACTUAL_M/'root_certificate.json'),
            'point_cache':relative(ACTUAL_M/'point_cache.json'),'point_cache_sha256':sha(ACTUAL_M/'point_cache.json'),
            'predictor_file':relative(ACTUAL_M/'affine_predictor.npz'),'predictor_sha256':sha(ACTUAL_M/'affine_predictor.npz'),
            'independent_review':relative(REVIEW),'independent_review_sha256':sha(REVIEW)},
       'environment_versions':versions(),'required_threads':{'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','torch_num_threads':1},
       'no_new_equilibrium_optimization':True,'heldout_response_files_read':False,'response_values_used_for_root_acceptance':False,
       'not_claimed':['all68 fit or biological validation','stability/commontrajectory','other stationary branches','ideal quadrature/continuum PDE','full identified set/unbounded domain','novel concept']}
    immutable(OUT/'protocol.json',p);preflight(p)
    print(json.dumps({'status':'FROZEN_PREFLIGHT_ONLY','protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'records':17,'new_cache_records':16,'reused_records':1}),flush=True)


def assert_frozen(p):
    if any(sha(ROOT/'scripts'/n)!=v for n,v in p['source_hashes'].items()):raise RuntimeError('frozen source changed')
    if any(sha(ROOT/'scripts'/n)!=v for n,v in p['common_contract_source_hashes'].items()):raise RuntimeError('common contract source changed')
    if any(sha(ROOT/n)!=v for n,v in p['input_hashes'].items()):raise RuntimeError('frozen input changed')
    if versions()!=p['environment_versions']:raise RuntimeError('numerical environment versions changed')
    for n in p['input_hashes']:safe_path(n)
    for row in p['records']:
        safe_path(row['saved_state']);folder=safe_path(row['output_folder'],write=True)
        if folder.parent!=OUT.resolve()/'records':raise RuntimeError('record output must be under owned records directory')


def validate_prepared_metadata(p,state,reference,triangles,u,integration_policy,verification,local_terms):
    import numpy as np
    checks=verification.get('checks')
    if not isinstance(checks,dict) or not checks or not all(v is True for v in checks.values()):raise ValueError('nonempty all-bool passing tangent verification required')
    if not np.array_equal(np.asarray(u),np.asarray(state['free_displacement'])):raise ValueError('prepared free vector differs from exact saved vector')
    if not np.array_equal(np.asarray(reference),np.asarray(state['reference_nodes'])) or not np.array_equal(np.asarray(triangles),np.asarray(state['triangles'])):raise ValueError('prepared reference/mesh differs from saved field')
    if integration_policy!=p['integration_policy']:raise ValueError('prepared integration rule differs from frozen actual policy')
    if len(u)!=p['expected_free_dofs'] or local_terms!=p['expected_local_terms']:raise ValueError('unexpected field dimensions/quadrature')


def guard_tests():
    fakep={'integration_policy':{'volume':8,'surface':12},'expected_free_dofs':2,'expected_local_terms':3}
    state={'free_displacement':[1.,2.],'reference_nodes':[[0.,0.],[1.,1.]],'triangles':[[0,1]]}
    def check(verification=None,u=None,policy=None,reference=None,triangles=None,terms=3):
        validate_prepared_metadata(fakep,state,reference if reference is not None else state['reference_nodes'],triangles if triangles is not None else state['triangles'],
            u if u is not None else state['free_displacement'],policy if policy is not None else fakep['integration_policy'],verification if verification is not None else {'checks':{'good':True}},terms)
    check();tests={'passing_metadata_accepted':True}
    mutations={'empty_checks':{'verification':{'checks':{}}},'failed_check':{'verification':{'checks':{'bad':False}}},'truthy_nonbool_check':{'verification':{'checks':{'bad':1}}},
        'missing_checks':{'verification':{}},'different_vector':{'u':[1.,3.]},'different_policy':{'policy':{'volume':4,'surface':8}},'different_reference':{'reference':[[0,0],[2,2]]},'different_mesh':{'triangles':[[1,0]]},'different_term_count':{'terms':4}}
    for name,kwargs in mutations.items():
        try:check(**kwargs);tests[name+'_rejected']=False
        except ValueError:tests[name+'_rejected']=True
    for name,path,write in [('outside_root',ROOT.parent/'not_owned.json',False),('outside_output',ROOT/'scripts/not_owned.json',True)]:
        try:safe_path(path,write);tests[name+'_rejected']=False
        except ValueError:tests[name+'_rejected']=True
    tests['owned_output_accepted']=safe_path(OUT/'records/test.json',True).is_relative_to(OUT.resolve())
    if not all(tests.values()):raise RuntimeError('prelaunch guard tests failed')
    return tests


def preflight(p):
    assert_frozen(p);C=p['shared_parameter_cell'];checks={
        'all17_distinct_source_columns':len(p['records'])==17 and len({r['record_id'] for r in p['records']})==17,
        'one_identical_shared_cell':cell_hash(C)==p['shared_parameter_cell_sha256'],
        'same_exact_gain_for_all_doses':len({r['base_gain_exact'] for r in p['records']})==1 and all(Q(r['exact_requested_strength'])==Q(r['base_gain_exact'])*Q(r['dose_exact']) for r in p['records']),
        'all5_coordinates_positive_width':len(C)==5 and all(Q(a)<Q(b) for a,b in C.values()),
        'apical_zero_face_in_shared_cell':Q(C['apical_surface_tension'][0])==0,
        'precisely_one_bound_actualM_reuse':sum(r['mode']=='reuse_immutable_actual_M' for r in p['records'])==1,
        'no_parameter_cell_shrink_or_nonlinear_optimization_policy':p['attempt_policy']['exactly_one_shared_cell_per_record'] and p['no_new_equilibrium_optimization'] and p['attempt_policy']['radius_multipliers']==[1],
        'fixed611dof9504terms_required':p['expected_free_dofs']==611 and p['expected_local_terms']==9504,
        'sources_data_state_and_common_contract_bound':True,'no_heldout_reads':p['heldout_response_files_read'] is False,
        'all_unresolved_failures_retained_never_OUT':'UNRESOLVED, never OUT' in p['attempt_policy']['failed_sufficient_test']}
    identity=process_identity(os.getpid());checks['actual_windows_PID_creation_identity_readable_without_signals']=bool(identity and identity['alive'] and identity['creation_ticks']>0)
    tests=guard_tests();checks['all_selfcontained_guard_tests_pass']=all(tests.values())
    immutable(OUT/'guard_tests.json',{'status':'PASS','checks':tests,'scope':'Synthetic metadata/path guards only; no model/Hessian calculation','source_sha256':sha(__file__)})
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'protocol_sha256':sha(OUT/'protocol.json'),
        'source_sha256':sha(__file__),'guard_tests_sha256':sha(OUT/'guard_tests.json'),'preflight_only':True,'new_model_or_hessian_calculations':0,
        'launch_requires':'prelaunch_independent_review.json PASS with nonempty bool checks and exact source/protocol/preflight hashes'}
    immutable(OUT/'preflight.json',result)
    if result['status']!='PASS':raise RuntimeError('preflight failed')


def prepare_record(p,row):
    import numpy as np,torch
    from dataclasses import asdict
    from scipy.sparse.linalg import splu
    import algebraic_root_feasibility as core
    from material_tangent_inference import mechanical_tangent
    folder=safe_path(row['output_folder'],write=True);manifest=folder/'prepared.json'
    if manifest.exists():
        done=read(manifest)
        if done['protocol_sha256']!=sha(OUT/'protocol.json') or any(sha(safe_path(k))!=v for k,v in done['artifact_hashes'].items()):raise RuntimeError('prepared record altered')
        work=safe_path(done['work_folder'],write=True)
        if work.parent!=folder/'preparation_attempts' or safe_path(done['point_cache_file'])!=work/'point_cache.json' or safe_path(done['predictor_file'])!=work/'affine_predictor.npz':raise RuntimeError('prepared paths differ from owned cache layout')
        return done
    attempts=folder/'preparation_attempts';attempts.mkdir(parents=True,exist_ok=True)
    index=1+max([int(d.name) for d in attempts.iterdir() if d.is_dir() and d.name.isdigit()]+[-1]);work=attempts/f'{index:03d}';work.mkdir()
    immutable(work/'started.json',{'protocol_sha256':sha(OUT/'protocol.json'),'record_id':row['record_id'],'started_utc':stamp(),'no_equilibrium_optimization':True})
    torch.set_num_threads(1);state=read(ROOT/row['saved_state']);model,u,H,B,S,verification=mechanical_tangent(state)
    local_terms=len(model.triangles)*len(model.shape)+sum(len(s[0])*len(s[1]) for s in model.surfaces)
    validate_prepared_metadata(p,state,model.reference,model.triangles,u,asdict(model.integration_policy),verification,local_terms)
    R=splu(H.tocsc()).solve(np.eye(len(u)))
    buf=io.BytesIO();np.savez_compressed(buf,free_vector=u,free_indices=model.free,H_data=H.data,H_indices=H.indices,H_indptr=H.indptr,H_shape=np.array(H.shape),parameter_force_derivatives=B,displacement_tangents=S);binary(work/'tangents.npz',buf.getvalue())
    predictor=S[:,[1,2,3,4,5]].copy();base=p['base_parameters_exact'];predictor/=np.array([float(Q(base[k])) for k in ('apical_area_modulus','basal_area_modulus')]+[1.,float(Q(base['basal_surface_tension'])),float(Q(base['common_gain']))])[None,:]
    buf=io.BytesIO();np.savez_compressed(buf,S=predictor);binary(work/'affine_predictor.npz',buf.getvalue())
    immutable(work/'tangent_verification.json',verification)
    point={'status':'FROZEN_SAVED_STATE_CACHE_NO_ROOT_CLAIM','parameters':state['parameters'],'activation_strength_binary_exact':row['solver_strength_binary_exact'],
        'input_state':row['saved_state'],'input_state_sha256':row['saved_state_sha256'],'matrix_archive':relative(work/'tangents.npz'),'matrix_archive_sha256':sha(work/'tangents.npz'),
        'frozen_binary_cache':core.freeze_cache(model,u,R),'protocol_sha256':sha(OUT/'protocol.json'),'record_id':row['record_id'],
        'dose_exact':row['dose_exact'],'exact_requested_strength':row['exact_requested_strength'],'no_new_equilibrium_optimization':True,
        'integration_policy':asdict(model.integration_policy),'uniform_target_warning':'Use exact shared common_gain*dose, not the archived rounded solver strength'}
    immutable(work/'point_cache.json',point)
    info={'status':'PREPARED_NO_ROOT_CLAIM','record_id':row['record_id'],'protocol_sha256':sha(OUT/'protocol.json'),
        'point_cache_file':relative(work/'point_cache.json'),'predictor_file':relative(work/'affine_predictor.npz'),'work_folder':relative(work),
        'free_dofs':len(u),'local_terms':p['expected_local_terms'],'artifact_hashes':{relative(work/n):sha(work/n) for n in ('tangents.npz','affine_predictor.npz','tangent_verification.json','point_cache.json','started.json')}}
    immutable(manifest,info);return info


def record_result(p,row):
    folder=safe_path(row['output_folder'],write=True);path=folder/'result.json'
    if path.exists():
        done=read(path)
        if done['protocol_sha256']!=sha(OUT/'protocol.json') or done['saved_state_sha256']!=row['saved_state_sha256'] or done['parameter_cell_sha256']!=p['shared_parameter_cell_sha256']:raise RuntimeError('completed record context changed')
        for k,v in done.get('artifact_hashes',{}).items():
            if sha(ROOT/k)!=v:raise RuntimeError('completed record artifact changed')
        return done
    started=time.time();common={'record_id':row['record_id'],'protocol_sha256':sha(OUT/'protocol.json'),'parameter_cell':p['shared_parameter_cell'],'parameter_cell_sha256':p['shared_parameter_cell_sha256'],
        'model_observation_contract_sha256':p['common_model_observation_contract_sha256'],'quantifier':p['quantifier'],'input_kind':'source_record_dose',
        'dose_exact':row['dose_exact'],'strength_relation':'common_gain_times_source_dose','base_strength_exact':row['exact_requested_strength'],
        'saved_state':row['saved_state'],'saved_state_sha256':row['saved_state_sha256'],'branch_scope':'selected_unique_root_in_declared_domain',
        'branch_domain_id':row['branch_domain_id'],'branch_law_semantics':p['branch_law_semantics'],'common_path_or_energetic_stability_certified':False,
        'source_hashes':p['source_hashes'],'no_new_equilibrium_optimization':True,'mode':row['mode']}
    if row['mode']=='reuse_immutable_actual_M':
        old=read(ROOT/p['actualM_reuse']['root_certificate']);oldcell=old['uniform_cell_result']
        if cell_hash(oldcell['cell_bounds_exact'])!=p['shared_parameter_cell_sha256'] or old['dose_exact']!=row['dose_exact']:raise RuntimeError('reuse family mismatch')
        answer={**common,'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','uniform_cell_result':oldcell,'reuse_bindings':p['actualM_reuse'],
            'point_cache_file':p['actualM_reuse']['point_cache'],'predictor_file':p['actualM_reuse']['predictor_file'],
            'free_dofs':611,'local_terms':9504,'artifact_hashes':{p['actualM_reuse'][key]:p['actualM_reuse'][hashkey] for key,hashkey in (('root_certificate','root_certificate_sha256'),('point_cache','point_cache_sha256'),('predictor_file','predictor_sha256'),('independent_review','independent_review_sha256'))},
            'reuse_scope':'Mechanically identical actualM root family; common cohort contract wrapper supplied here, original immutable certificate/review retained.'}
    else:
        try:
            prepared=prepare_record(p,row)
            import parametric_root_cells as kernel
            from flint import ctx
            ctx.prec=p['precision_bits'];kernel.A.cache_clear();kernel.OUT=ROOT/prepared['work_folder']
            local={'point_certificate_file':prepared['point_cache_file'],'predictor_file':'affine_predictor.npz','dose_exact':row['dose_exact'],
                'base_parameters_exact':p['base_parameters_exact'],'calibration_record_ids_covered':[row['record_id']],'input_kind':'source_record_dose',
                'strength_relation':'common_gain_times_source_dose','branch_domain_id':row['branch_domain_id'],**p['attempt_policy']}
            point=kernel.build_point(local)
            if point['n']!=611 or len(point['cache'])!=9504:raise RuntimeError('unexpected assembled point size')
            result=kernel.evaluate_cell(point,local,'shared_actual_M_scale_3',p['shared_parameter_cell'],'affine')
            answer={**common,'status':result['status'],'uniform_cell_result':result,'point_cache_file':prepared['point_cache_file'],
                'predictor_file':prepared['predictor_file'],'prepared_manifest':relative(folder/'prepared.json'),'prepared_manifest_sha256':sha(folder/'prepared.json'),
                'free_dofs':point['n'],'local_terms':len(point['cache']),'point_defect_Z0':kernel.pack(point['Z0']),'artifact_hashes':prepared['artifact_hashes']}
        except (ValueError,AssertionError) as error:
            answer={**common,'status':'UNRESOLVED','failure_category':'preparation_or_sufficient_bound_failed','error':repr(error),'traceback':traceback.format_exc(),
                'parameter_cell_unchanged':True,'no_OUT_classification':True}
    answer['elapsed_seconds']=time.time()-started;immutable(path,answer);return answer


def process_identity(pid):
    """Read-only Windows PID+creation-time guard; never sends a signal."""
    if os.name!='nt':raise RuntimeError('this launch controller requires Windows process identity API')
    import ctypes
    from ctypes import wintypes
    lib=ctypes.WinDLL('kernel32',use_last_error=True);lib.OpenProcess.restype=wintypes.HANDLE
    lib.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD]
    lib.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)]
    lib.GetProcessTimes.argtypes=[wintypes.HANDLE,*([ctypes.POINTER(wintypes.FILETIME)]*4)]
    lib.CloseHandle.argtypes=[wintypes.HANDLE]
    handle=lib.OpenProcess(0x1000,False,int(pid))
    if not handle:return None
    try:
        code=wintypes.DWORD();creation=wintypes.FILETIME();exit_=wintypes.FILETIME();kernel=wintypes.FILETIME();user=wintypes.FILETIME()
        if not lib.GetExitCodeProcess(handle,ctypes.byref(code)) or not lib.GetProcessTimes(handle,ctypes.byref(creation),ctypes.byref(exit_),ctypes.byref(kernel),ctypes.byref(user)):return None
        return {'pid':int(pid),'alive':code.value==259,'creation_ticks':(creation.dwHighDateTime<<32)+creation.dwLowDateTime}
    finally:lib.CloseHandle(handle)


def run(resume=False):
    p=read(OUT/'protocol.json');assert_frozen(p);pre=read(OUT/'preflight.json');reviewpath=OUT/'prelaunch_independent_review.json';reviewbytes=reviewpath.read_bytes();review=json.loads(reviewbytes)
    if pre['status']!='PASS' or review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()):raise RuntimeError('bound independent prelaunch review required')
    if review.get('source_sha256')!=sha(__file__) or review.get('protocol_sha256')!=sha(OUT/'protocol.json') or review.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise RuntimeError('prelaunch review source/protocol/preflight mismatch')
    if any(os.environ.get(k)!=v for k,v in p['required_threads'].items() if k!='torch_num_threads'):raise RuntimeError('declared BLAS/OMP thread policy required')
    if (OUT/'results.json').exists():raise RuntimeError('cohort terminal receipt already exists; preserve it')
    lock=OUT/'active.lock';jobpath=OUT/'job.json'
    if lock.exists():
        prior=read(lock);identity=process_identity(prior['pid'])
        if identity and identity['alive'] and identity['creation_ticks']==prior['creation_ticks']:raise RuntimeError('same cohort job is still live; never duplicate')
        if not resume:raise RuntimeError('stale lease retained; explicit --resume required')
        old=OUT/'launch_history';old.mkdir(exist_ok=True);j=1+max([int(x.stem.split('_')[-1]) for x in old.glob('lease_*.json')]+[-1]);immutable(old/f'lease_{j:03}.json',prior)
        if jobpath.exists():immutable(old/f'job_{j:03}.json',read(jobpath))
        lock.unlink()  # one fixed file in this explicitly owned workspace; archived above
    identity=process_identity(os.getpid())
    if not identity or not identity['alive']:raise RuntimeError('actual worker identity unavailable')
    with lock.open('x',encoding='utf8') as f:json.dump({**identity,'protocol_sha256':sha(OUT/'protocol.json'),'started_utc':stamp()},f)
    launchdir=OUT/'launch_history';launchdir.mkdir(exist_ok=True);j=1+max([int(x.stem.split('_')[-1]) for x in launchdir.glob('review_*.json')]+[-1]);binary(launchdir/f'review_{j:03}.json',reviewbytes)
    job={'status':'RUNNING',**identity,'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),'prelaunch_review_sha256':hashlib.sha256(reviewbytes).hexdigest(),'started_utc':stamp()};atomic(jobpath,job)
    completed=[];started=time.time()
    try:
        for row in p['records']:
            assert_frozen(p);answer=record_result(p,row);completed.append({'record_id':row['record_id'],'status':answer['status'],'file':relative(ROOT/row['output_folder']/'result.json'),'sha256':sha(ROOT/row['output_folder']/'result.json')})
            atomic(OUT/'progress.json',{'status':'RUNNING','pid':os.getpid(),'completed':completed,'total':17,'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],'elapsed_seconds':time.time()-started})
            print(json.dumps({'record_id':row['record_id'],'status':answer['status'],'completed':len(completed),'total':17,'elapsed_seconds':time.time()-started}),flush=True)
        assert_frozen(p);allgood=all(x['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' for x in completed)
        summary={'status':'CERTIFIED_SHARED_CELL_STATIONARY_ROOT_TUPLE' if allgood else 'COMPLETE_WITH_UNRESOLVED_ROOT_RECORDS','protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),
            'records':completed,'shared_parameter_cell':p['shared_parameter_cell'],'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],
            'common_model_observation_contract_sha256':p['common_model_observation_contract_sha256'],'all17_certified':allgood,'certified':sum(x['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' for x in completed),'unresolved':sum(x['status']=='UNRESOLVED' for x in completed),
            'branch_law_semantics':p['branch_law_semantics'],'quantified_conclusion':'For every parameter in the shared cell, a stationary root tuple exists uniquely within the product of the17 declared nearby domains.' if allgood else 'No all17 root tuple certificate; unresolved records retained.',
            'no_equilibrium_optimization':True,'no_OUT_classifications':True,'elapsed_seconds':time.time()-started,'completed_utc':stamp(),
            'not_claimed':p['not_claimed']}
        immutable(OUT/'results.json',summary);job.update(status='FINISHED',finished_utc=stamp(),results_sha256=sha(OUT/'results.json'));atomic(jobpath,job)
        print(json.dumps({'status':summary['status'],'certified':summary['certified'],'unresolved':summary['unresolved'],'results_sha256':sha(OUT/'results.json')}),flush=True)
    except Exception:
        job.update(status='FAILED',failed_utc=stamp(),traceback=traceback.format_exc());atomic(jobpath,job);raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','run']);parser.add_argument('--contract');parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    if args.stage=='freeze':
        if not args.contract:raise SystemExit('--contract path required before protocol freeze')
        freeze(args.contract)
    else:run(args.resume)
