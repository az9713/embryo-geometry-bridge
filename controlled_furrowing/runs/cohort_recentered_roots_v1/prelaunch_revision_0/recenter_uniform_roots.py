"""ONE validated-Newton proof-center correction for F/K on the SAME cell.

Not an optimizer, data fit, physical change or silent seed-state replacement.
All original caches/R/S/states stay immutable. Fresh interval proofs required.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,copy,hashlib,io,json,os,time,traceback
import numpy as np
import cohort_uniform_roots as helpers

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'runs/cohort_uniform_roots_v1'
OUT=ROOT/'runs/cohort_recentered_roots_v1'
helpers.OUT=OUT
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()


def arrays(packet):
    raw=base64.b64decode(packet['data'])
    if hashlib.sha256(raw).hexdigest()!=packet['sha256']:raise RuntimeError('frozen cache checksum mismatch')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}


def fingerprint(values,exclude):
    h=hashlib.sha256()
    for k in sorted(values):
        if k in exclude:continue
        a=np.ascontiguousarray(values[k]);h.update(k.encode()+b'\0'+a.dtype.str.encode()+b'\0'+json.dumps(list(a.shape)).encode()+b'\0'+a.tobytes())
    return h.hexdigest()


def packet(values):
    buf=io.BytesIO();np.savez_compressed(buf,**values);raw=buf.getvalue()
    return {'encoding':'base64 np.savez_compressed exact binary arrays; ONLYu changed from seed cache','sha256':hashlib.sha256(raw).hexdigest(),
        'bytes':len(raw),'data':base64.b64encode(raw).decode(),'array_shapes':{k:list(v.shape) for k,v in values.items()}}


def one_round_center(old,midpoints):
    if len(old)!=len(midpoints):raise ValueError('one midpoint per free coordinate required')
    return np.array([float(Q.from_float(float(u))-Q(m)) for u,m in zip(old,midpoints)],dtype=np.float64)


def guard_tests():
    tests={};u=np.array([.5,-.25]);m=[Q(1,10),Q(-1,20)];v=one_round_center(u,m)
    tests['single_final_rounding_from_exact_binary_minus_rational']=all(Q.from_float(float(x))==Q.from_float(float(float(Q.from_float(float(y))-z))) for x,y,z in zip(v,u,m))
    tests['nonzero_correction_changes_center']=not np.array_equal(u,v)
    old={'u':u,'preconditioner':np.eye(2),'reference':np.ones((2,2))};new={**old,'u':v}
    tests['all_cache_arrays_except_u_fingerprint_invariant']=fingerprint(old,{'u'})==fingerprint(new,{'u'})
    tests['changed_R_detected']=fingerprint(old,{'u'})!=fingerprint({**new,'preconditioner':2*np.eye(2)},{'u'})
    tests['changed_reference_detected']=fingerprint(old,{'u'})!=fingerprint({**new,'reference':np.zeros((2,2))},{'u'})
    tests['metadata_guards_from_reviewed_driver_pass']=all(helpers.guard_tests().values())
    if not all(tests.values()):raise RuntimeError('synthetic recenter guards failed')
    return tests


def freeze():
    if (OUT/'protocol.json').exists():raise RuntimeError('recenter protocol exists; retain')
    old=read(OLD/'protocol.json');p=copy.deepcopy(old);p.update(version=1,created_utc=helpers.stamp(),scope='One proof-center refinement per F/K, unchanged shared physical model and parameter cell; fresh interval existence test')
    p['source_hashes']['recenter_uniform_roots.py']=sha(__file__);p['input_hashes'][rel(OLD/'protocol.json')]=sha(OLD/'protocol.json')
    diag=ROOT/'runs/cohort_bound_diagnostics_v1/results.json';p['input_hashes'][rel(diag)]=sha(diag)
    p['original_protocol_file']=rel(OLD/'protocol.json');p['original_protocol_sha256']=sha(OLD/'protocol.json');p['records']=[]
    for row in old['records']:
        if row['record_id'] not in ('optoGEF2_small_col_F','optoGEF2_small_col_K'):continue
        priorpath=ROOT/row['output_folder']/'result.json';prior=read(priorpath)
        if prior['status']!='UNRESOLVED':raise RuntimeError('retained unresolved seed input required')
        r=copy.deepcopy(row);r.update(mode='single_proof_center_refinement',original_result_file=rel(priorpath),original_result_sha256=sha(priorpath),
            original_point_cache_file=prior['point_cache_file'],predictor_file=prior['predictor_file'],
            original_seed_state_file=row['saved_state'],original_seed_state_sha256=row['saved_state_sha256'],
            branch_domain_id='cohort_recentered_v1:'+row['record_id'],output_folder=rel(OUT/'records'/row['record_id']))
        p['records'].append(r);p['input_hashes'][rel(priorpath)]=sha(priorpath);p['input_hashes'].update(prior['artifact_hashes'])
        cache=read(ROOT/r['original_point_cache_file']);a=arrays(cache['frozen_binary_cache'])
        if a['u'].shape!=(611,) or a['preconditioner'].shape!=(611,611):raise RuntimeError('cached seed dimensions unexpected')
    p['center_recipe']={'number_of_corrections':1,'formula':'u1_i=float64(Fraction(binary_u0_i)-exact_dyadic_midpoint((R*F_exact(u0,theta0))_i))',
        'rounding':'ONE nearest-even binary64 rounding of the exact subtraction; a separately rounded correction is not used in subtraction',
        'midpoint_precision_bits':p['precision_bits'],'R_policy':'preserve original exact binary preconditioner','S_policy':'preserve original exact binary affine predictor file',
        'physical_target':'all cached arrays exceptu, modelparams, theta0, exactgain*dose, geometry/reference/quadrature and sharedC remain unchanged',
        'acceptance':'Fresh Arb center and affine-cell theorem around u1. Correction itself is not a root certificate. One normal radius test using original policy.',
        'translation_policy':'Keep exact u1-u0 separate. DO NOT add center translation to Banach corrector r or its eta/Z calculation.',
        'auxiliary_state':'Numerical proof center, NOT raw observation or converged nonlinear-solver state; original seed identity remains separate.'}
    helpers.immutable(OUT/'protocol.json',p);helpers.assert_frozen(p);tests=guard_tests()
    helpers.immutable(OUT/'guard_tests.json',{'status':'PASS','checks':tests,'source_sha256':sha(__file__),'synthetic_only':True})
    checks={'only_F_K_unresolved_seed_records':len(p['records'])==2,'same_exact_sharedC':p['shared_parameter_cell']==old['shared_parameter_cell'],
        'same_theta0_model_geometry_and_integration':all(p[k]==old[k] for k in ('base_parameters_exact','fixed_full_parameters','integration_policy')),
        'same_normal_affine_cell_radius_policy':p['attempt_policy']==old['attempt_policy'],'exactly_one_correction_no_iteration':p['center_recipe']['number_of_corrections']==1,
        'all_seed_cache_R_S_state_and_source_hashes_current':True,'all_synthetic_guards_pass':all(tests.values()),'no_real_center_or_assembly_before_review':not (OUT/'records').exists()}
    helpers.immutable(OUT/'preflight.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),
        'protocol_sha256':sha(OUT/'protocol.json'),'guard_tests_sha256':sha(OUT/'guard_tests.json'),'actual_RF_or_center_calculations':0})
    print(json.dumps({'status':'FROZEN_PREFLIGHT_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)


def local_protocol(p,row,cachefile,predictor_name):
    return {'point_certificate_file':cachefile,'predictor_file':predictor_name,'dose_exact':row['dose_exact'],'base_parameters_exact':p['base_parameters_exact'],
        'calibration_record_ids_covered':[row['record_id']],'input_kind':'source_record_dose','strength_relation':'common_gain_times_source_dose',
        'branch_domain_id':row['branch_domain_id'],**p['attempt_policy']}


def calculate(p,row):
    import parametric_root_cells as kernel
    from flint import ctx
    ctx.prec=p['precision_bits'];kernel.A.cache_clear();kernel.OUT=(ROOT/row['predictor_file']).parent
    oldpoint=kernel.build_point(local_protocol(p,row,row['original_point_cache_file'],Path(row['predictor_file']).name))
    oldcache=read(ROOT/row['original_point_cache_file']);values=arrays(oldcache['frozen_binary_cache']);u0=values['u'].copy()
    mids=[str(oldpoint['gR'][i,0].mid().fmpq()) for i in range(611)];u1=one_round_center(u0,mids)
    translations=[Q.from_float(float(b))-Q.from_float(float(a)) for a,b in zip(u0,u1)];newvalues={k:v.copy() for k,v in values.items()};newvalues['u']=u1.copy()
    if fingerprint(values,{'u'})!=fingerprint(newvalues,{'u'}) or not np.array_equal(newvalues['preconditioner'],values['preconditioner']):raise RuntimeError('physical cache or R changed')
    work=ROOT/row['output_folder'];work.mkdir(parents=True,exist_ok=True)
    buf=io.BytesIO();np.savez_compressed(buf,old_u=u0,new_u=u1,RF_midpoint_float64=np.array([float(Q(m)) for m in mids]),actual_center_translation_float64=np.array([float(q) for q in translations]));helpers.binary(work/'center_arrays.npz',buf.getvalue())
    seed=read(ROOT/row['original_seed_state_file']);current=values['reference'].copy().ravel();current[values['free']]+=u1
    aux={'kind':'AUXILIARY_NUMERICAL_PROOF_CENTER_NOT_EXPERIMENT_OR_SOLVER_RESULT','converged':False,'solver_status':'NO_NONLINEAR_SOLVER_RUN',
        'parameters':copy.deepcopy(oldcache['parameters']),'reference_nodes':values['reference'].tolist(),'triangles':values['triangles'].tolist(),
        'free_displacement':u1.tolist(),'current_nodes':current.reshape(-1,2).tolist(),
        'original_seed_state_file':row['original_seed_state_file'],'original_seed_state_sha256':row['original_seed_state_sha256'],
        'proof_activation_strength_exact':row['exact_requested_strength'],'original_solver_strength_binary_exact':row['solver_strength_binary_exact']}
    helpers.immutable(work/'proof_center_state.json',aux)
    cache=copy.deepcopy(oldcache);cache.update(status='FROZEN_AUXILIARY_PROOF_CENTER_NO_ROOT_CLAIM',input_state=rel(work/'proof_center_state.json'),input_state_sha256=sha(work/'proof_center_state.json'),
        frozen_binary_cache=packet(newvalues),original_seed_state_file=row['original_seed_state_file'],original_seed_state_sha256=row['original_seed_state_sha256'],
        original_point_cache_file=row['original_point_cache_file'],original_point_cache_sha256=sha(ROOT/row['original_point_cache_file']),
        matrix_archive_role='Original seed floatH/tangent archive, NOT re-evaluated at u1. R/S remain fixed; freshArbH at u1 validates their defect.')
    helpers.immutable(work/'point_cache.json',cache)
    center={'status':'ONE_NUMERICAL_PROOF_CENTER_CORRECTION','center_arrays_file':rel(work/'center_arrays.npz'),'center_arrays_sha256':sha(work/'center_arrays.npz'),
        'original_seed_state_file':row['original_seed_state_file'],'original_seed_state_sha256':row['original_seed_state_sha256'],
        'new_proof_center_state_file':rel(work/'proof_center_state.json'),'new_proof_center_state_sha256':sha(work/'proof_center_state.json'),
        'RF_exact_enclosures':[kernel.pack(oldpoint['gR'][i,0]) for i in range(611)],'RF_midpoints_exact':mids,'actual_center_translation_exact':list(map(str,translations)),
        'translation_linf_exact':str(max(map(abs,translations))),'cache_arrays_except_u_fingerprint':fingerprint(values,{'u'}),
        'new_cache_arrays_except_u_fingerprint':fingerprint(newvalues,{'u'}),'potential_arrays_excluding_u_R_fingerprint':fingerprint(values,{'u','preconditioner'}),
        'new_potential_arrays_excluding_u_R_fingerprint':fingerprint(newvalues,{'u','preconditioner'}),
        'R_binary_unchanged':True,'S_file':row['predictor_file'],'S_sha256':sha(ROOT/row['predictor_file']),
        'all_parameters_reference_quadrature_and_cell_unchanged':True,'root_existence_not_inferred_from_correction':True}
    helpers.immutable(work/'center_refinement.json',center)
    # Predictor is copied bit-for-bit to this new owned proof work directory.
    helpers.binary(work/'affine_predictor.npz',(ROOT/row['predictor_file']).read_bytes());kernel.OUT=work;kernel.A.cache_clear()
    newpoint=kernel.build_point(local_protocol(p,row,rel(work/'point_cache.json'),'affine_predictor.npz'))
    if newpoint['n']!=611 or len(newpoint['cache'])!=9504:raise RuntimeError('fresh proof dimensions changed')
    proof=kernel.evaluate_cell(newpoint,local_protocol(p,row,rel(work/'point_cache.json'),'affine_predictor.npz'),'same_shared_cell_recentered_once',p['shared_parameter_cell'],'affine')
    answer={'status':proof['status'],'record_id':row['record_id'],'protocol_sha256':sha(OUT/'protocol.json'),'source_hashes':p['source_hashes'],
        'parameter_cell':p['shared_parameter_cell'],'parameter_cell_sha256':p['shared_parameter_cell_sha256'],'model_observation_contract_sha256':p['common_model_observation_contract_sha256'],
        'quantifier':p['quantifier'],'branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':row['branch_domain_id'],'branch_law_semantics':p['branch_law_semantics'],
        'dose_exact':row['dose_exact'],'input_kind':'source_record_dose','strength_relation':'common_gain_times_source_dose','base_strength_exact':row['exact_requested_strength'],
        'original_seed_state_file':row['original_seed_state_file'],'original_seed_state_sha256':row['original_seed_state_sha256'],
        'proof_center_state_file':rel(work/'proof_center_state.json'),'proof_center_state_sha256':sha(work/'proof_center_state.json'),
        'point_cache_file':rel(work/'point_cache.json'),'predictor_file':rel(work/'affine_predictor.npz'),'center_refinement_file':rel(work/'center_refinement.json'),'center_refinement_sha256':sha(work/'center_refinement.json'),
        'uniform_cell_result':proof,'normal_corrector_radius_is_not_center_translation':True,'point_defect_Z0':kernel.pack(newpoint['Z0']),
        'free_dofs':611,'local_terms':9504,'common_path_or_energetic_stability_certified':False,'no_parameter_change_or_equilibrium_optimizer':True,
        'artifact_hashes':{rel(work/n):sha(work/n) for n in ('center_arrays.npz','proof_center_state.json','point_cache.json','center_refinement.json','affine_predictor.npz')},
        'lineage':{'kind':'ONE_VALIDATED_NEWTON_PROOF_CENTER_REFINEMENT','original_result_file':row['original_result_file'],'original_result_sha256':row['original_result_sha256'],
            'original_point_cache_file':row['original_point_cache_file'],'original_point_cache_sha256':sha(ROOT/row['original_point_cache_file']),
            'original_seed_state_is_not_new_proof_center':True,'observer_requires_explicit_versioned_center_offset_adapter':True},
        'scope':'Numerical proofcenter changed once; physical model/sharedcell/gain*dose/R/S unchanged. Fresh proof determines root existence; oldsame-centerobserver cannotconsume this as originalstate.'}
    helpers.immutable(work/'result.json',answer);return answer


def run():
    p=read(OUT/'protocol.json');helpers.assert_frozen(p);raw=(OUT/'prelaunch_independent_review.json').read_bytes();audit=json.loads(raw)
    if audit.get('status')!='PASS' or not audit.get('checks') or not all(v is True for v in audit['checks'].values()) or audit.get('source_sha256')!=sha(__file__) or audit.get('protocol_sha256')!=sha(OUT/'protocol.json') or audit.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise RuntimeError('independent concrete recenter prelaunch review required')
    if (OUT/'job.json').exists():raise RuntimeError('recenter job exists; preserve/no duplicate')
    if any(os.environ.get(k)!=v for k,v in p['required_threads'].items() if k!='torch_num_threads'):raise RuntimeError('thread policy mismatch')
    identity=helpers.process_identity(os.getpid());job={'status':'RUNNING',**identity,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':hashlib.sha256(raw).hexdigest(),'started_utc':helpers.stamp()}
    with (OUT/'active.lock').open('x',encoding='utf8') as f:json.dump(job,f)
    helpers.binary(OUT/'launch_review.json',raw);helpers.atomic(OUT/'job.json',job);results=[];started=time.time()
    try:
        for row in p['records']:
            helpers.assert_frozen(p);answer=calculate(p,row);path=ROOT/row['output_folder']/'result.json';results.append({'record_id':row['record_id'],'status':answer['status'],'file':rel(path),'sha256':sha(path),'provider_protocol_file':rel(OUT/'protocol.json'),'provider_protocol_sha256':sha(OUT/'protocol.json')})
            print(json.dumps({'record':row['record_id'],'status':answer['status'],'seconds':time.time()-started}),flush=True)
        helpers.assert_frozen(p);result={'status':'COMPLETE_SINGLE_CENTER_REFINEMENT_TESTS','records':results,'certified':sum(r['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT' for r in results),'unresolved':sum(r['status']=='UNRESOLVED' for r in results),
            'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],'common_model_observation_contract_sha256':p['common_model_observation_contract_sha256'],
            'no_cell_or_physics_change':True,'one_center_correction_each_no_iteration':True,'no_OUT_classification':True,'elapsed_seconds':time.time()-started}
        helpers.immutable(OUT/'results.json',result);job.update(status='FINISHED',results_sha256=sha(OUT/'results.json'));helpers.atomic(OUT/'job.json',job)
    except Exception:
        job.update(status='FAILED',traceback=traceback.format_exc());helpers.atomic(OUT/'job.json',job);raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','run']);args=parser.parse_args();freeze() if args.stage=='freeze' else run()
