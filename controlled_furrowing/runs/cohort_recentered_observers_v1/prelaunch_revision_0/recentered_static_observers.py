"""Explicit corrected-proof-center -> ORIGINAL seed-graph observer adapter.

No root calculation or optimization. All old receipts and strict same-center
guards stay unchanged; fresh root receipts must supply the new corrector cube.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json
import numpy as np
import cohort_static_observers as old

ROOT=old.ROOT
OUT=ROOT/'runs/cohort_recentered_observers_v1'
ROOTRUN=ROOT/'runs/cohort_recentered_roots_v1'
IDS=('optoGEF2_small_col_F','optoGEF2_small_col_K')

def fingerprint(values,exclude):
    digest=hashlib.sha256()
    for key in sorted(values):
        if key in exclude:continue
        a=np.ascontiguousarray(values[key]);digest.update(key.encode()+b'\0'+a.dtype.str.encode()+b'\0'+json.dumps(list(a.shape)).encode()+b'\0'+a.tobytes())
    return digest.hexdigest()

def unchanged_cache(old_arrays,new_arrays):
    if set(old_arrays)!=set(new_arrays):raise ValueError('cache array keys changed')
    for key in old_arrays:
        a,b=old_arrays[key],new_arrays[key]
        if a.dtype!=b.dtype or a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('cache dtype/shape/finite guard: '+key)
        if key!='u' and a.tobytes()!=b.tobytes():raise ValueError('cached physical array or preconditioner changed: '+key)
    return fingerprint(old_arrays,{'u'})

def translated_box(seed,old_arrays,new_arrays,S,cell,base,radius):
    if radius<=0:raise ValueError('positive NORMAL corrector radius required')
    unchanged_cache(old_arrays,new_arrays)
    # Deliberately validate the ORIGINAL seed against the ORIGINAL cache.
    centers,unused,gaps,hulls,free=old.nodal_box(seed,old_arrays,S,cell,base,Q(0))
    radii=[row[:] for row in gaps];new_exact=[[old.q(v) for v in row] for row in seed['reference_nodes']]
    translations=[]
    for i,k in enumerate(free):
        new_exact[k//2][k%2]+=old.q(new_arrays['u'][i])
        translation=old.q(new_arrays['u'][i])-old.q(old_arrays['u'][i]);translations.append(translation)
        radii[k//2][k%2]=abs(new_exact[k//2][k%2]-centers[k//2][k%2])+max(map(abs,hulls[i]))+radius
    center_gaps=[[abs(new_exact[i][j]-centers[i][j]) for j in range(2)] for i in range(len(centers))]
    return centers,radii,center_gaps,hulls,free,translations

def full_potential_fingerprint(arrays,parameters,rootp,row):
    payload={'potential_arrays':fingerprint(arrays,{'u','preconditioner'}),'model_parameters':parameters,
        'theta0_exact':rootp['base_parameters_exact'],'dose_exact':row['dose_exact'],'shared_parameter_cell':rootp['shared_parameter_cell']}
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def tests():
    cell={k:['0','0'] for k in old.ORDER};base={k:'0' for k in old.ORDER};S=np.zeros((1,5))
    seed={'reference_nodes':[[0.,0.]],'current_nodes':[[0.,.1]],'triangles':[], 'free_displacement':[.1]}
    a={'reference':np.array([[0.,0.]]),'triangles':np.empty((0,)),'free':np.array([1]),'u':np.array([.1]),'preconditioner':np.eye(1)}
    b={k:v.copy() for k,v in a.items()};b['u']=np.array([.11]);r=Q(1,1000)
    centers,radii,gaps,hulls,free,translation=translated_box(seed,a,b,S,cell,base,r)
    checks={'only_u_may_change':unchanged_cache(a,b)==fingerprint(a,{'u'}),
        'normal_corrector_separate_from_translation':radii[0][1]==abs(old.q(.11)-old.q(.1))+r,
        'exact_translation_retained':translation==[old.q(.11)-old.q(.1)],'fixed_coordinate_error_stays_zero':radii[0][0]==0,
        'new_center_gap_is_not_old_serialization_gap':gaps[0][1]>0}
    bad=[]
    for key in ('preconditioner','reference','free'):
        changed={k:v.copy() for k,v in b.items()};changed[key]=changed[key]+1;bad.append(lambda changed=changed:unchanged_cache(a,changed))
    bad.append(lambda:old.nodal_box(seed,b,S,cell,base,r))
    for i,call in enumerate(bad):
        try:call();checks['reject_invalid_'+str(i)]=False
        except ValueError:checks['reject_invalid_'+str(i)]=True
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':old.sha(__file__),'synthetic_only':True}
    old.write(OUT/'synthetic_checks.json',result)
    if not all(checks.values()):raise ValueError('center-offset synthetic guard failed')
    return result

def freeze():
    rootp=old.read(ROOTRUN/'protocol.json');contract=old.read(old.CONTRACT)
    sources=old.dependencies();sources[old.rel(__file__)]=old.sha(__file__)
    sources.update({'scripts/'+k:v for k,v in rootp['source_hashes'].items()})
    inputs=dict(rootp['input_hashes'])
    for file in (ROOTRUN/'protocol.json',ROOTRUN/'preflight.json',ROOTRUN/'prelaunch_independent_review.json',old.CONTRACT,
                 ROOT/'runs/cohort_static_observers_v1/protocol.json',ROOT/'runs/cohort_static_observers_v1/composition_final/results.json'):
        inputs[old.rel(file)]=old.sha(file)
    # Bind terminal provider bytes only; mathematical payloads are consumed after
    # this adapter's independent prelaunch review, never during freeze.
    inputs[old.rel(ROOTRUN/'results.json')]=old.sha(ROOTRUN/'results.json')
    for row in rootp['records']:
        file=old.path(row['output_folder'])/'result.json';inputs[old.rel(file)]=old.sha(file)
    old.assert_hashes({**sources,**inputs})
    checks={'exactly_F_K':tuple(row['record_id'] for row in rootp['records'])==IDS,
        'same_shared_cell':rootp['shared_parameter_cell_sha256']==contract['parameter_cell_sha256'],
        'same_common_contract':rootp['common_model_observation_contract_sha256']==old.sha(old.CONTRACT),
        'one_numerical_center_correction_declared':rootp['center_recipe']['number_of_corrections']==1,
        'no_new_root_receipts_consumed':True}
    for row in rootp['records']:
        seed=old.read(old.path(row['original_seed_state_file']));packet=old.read(old.path(row['original_point_cache_file']));arrays=old.decode_cache(packet)
        old.validate_point_target(packet,seed,row)
        with np.load(old.path(row['predictor_file']),allow_pickle=False) as z:S=z['S'].copy()
        old.nodal_box(seed,arrays,S,rootp['shared_parameter_cell'],rootp['base_parameters_exact'],Q(0))
        checks[row['record_id']+'_original_seed_cache_exact']=old.sha(old.path(row['original_seed_state_file']))==row['original_seed_state_sha256']
    if not all(checks.values()):raise ValueError('recentered observer preflight failed')
    protocol={'version':1,'created_utc':old.stamp(),'scope':'F/K only, new proof center mapped to ORIGINAL seed graphs; no root optimization',
        'source_hashes':sources,'input_hashes':inputs,'root_protocol_file':old.rel(ROOTRUN/'protocol.json'),'root_protocol_sha256':old.sha(ROOTRUN/'protocol.json'),
        'common_contract_sha256':old.sha(old.CONTRACT),'parameter_cell_sha256':contract['parameter_cell_sha256'],
        'adapter_kind':'RECENTERED_ROOT_TO_ORIGINAL_SEED_GRAPH','normal_corrector_radius_is_separate':True,
        'coordinate_error_formula':'abs(exact_binary_reference+new_cache_u-original_serialized_seed_node)+max_abs(S*parameter_delta)+FULL_normal_corrector_r',
        'fixed_coordinates':'exact gap only; no corrector assigned to constrained coordinates',
        'physical_target_guard':'all old/new cache arrays exceptu byte-identical INCLUDING preconditioner, same parameters/theta0/dose/C/S and full-potential fingerprints',
        'auxiliary_state_rule':'kind explicitly numerical, converged false, no nonlinear solver result or raw-data identity',
        'maximum_reporting_policy':'request1e-7um/8192cells; valid wider intervals retained with actual reporting certificate',
        'branch_law_semantics':'independent_stationary_root_per_record','no_common_path_or_stability_claim':True,
        'prelaunch_required_fields':['status=PASS','checks all literal true','observer_source_sha256','observer_protocol_sha256'],
        'unresolved_policy':'No endpoint certificate or interval fabricated if fresh root proof or geometry/ROI validity is unresolved',
        'independent_root_and_application_review_required':True}
    old.write(OUT/'protocol.json',protocol)
    old.write(OUT/'preflight.json',{'status':'PASS','checks':checks,'source_sha256':old.sha(__file__),'protocol_sha256':old.sha(OUT/'protocol.json')})
    return protocol

def compose(p,rootp,row,root):
    expected={'record_id':row['record_id'],'protocol_sha256':p['root_protocol_sha256'],'parameter_cell_sha256':p['parameter_cell_sha256'],
        'model_observation_contract_sha256':p['common_contract_sha256'],'quantifier':rootp['quantifier'],'branch_scope':'selected_unique_root_in_declared_domain',
        'branch_domain_id':row['branch_domain_id'],'branch_law_semantics':'independent_stationary_root_per_record','dose_exact':row['dose_exact'],
        'base_strength_exact':row['exact_requested_strength'],'input_kind':'source_record_dose','strength_relation':'common_gain_times_source_dose',
        'common_path_or_energetic_stability_certified':False,'normal_corrector_radius_is_not_center_translation':True,
        'original_seed_state_file':row['original_seed_state_file'],'original_seed_state_sha256':row['original_seed_state_sha256']}
    old.validate_root_context(root,expected,rootp['shared_parameter_cell']);old.assert_hashes(root['artifact_hashes'])
    old.assert_hashes({'scripts/'+k:v for k,v in root['source_hashes'].items()})
    if root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT':return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Fresh corrected-center root not certified','no_endpoint_certificate':True}
    proof=root['uniform_cell_result']
    if old.canonical(proof['cell_bounds_exact'])!=old.canonical(rootp['shared_parameter_cell']) or old.q(proof['dose_exact'])!=old.q(row['dose_exact']):raise ValueError('inner corrected root family mismatch')
    radius=old.accepted_radius(proof);seed=old.read(old.path(row['original_seed_state_file']))
    originalpacket=old.read(old.path(row['original_point_cache_file']));newpacket=old.read(old.path(root['point_cache_file']))
    old.validate_point_target(originalpacket,seed,row)
    lineage=root['lineage']
    if lineage['original_result_file']!=row['original_result_file'] or lineage['original_result_sha256']!=row['original_result_sha256'] or lineage['original_point_cache_file']!=row['original_point_cache_file'] or lineage['original_point_cache_sha256']!=old.sha(old.path(row['original_point_cache_file'])):raise ValueError('original proof-cache/result lineage differs')
    a=old.decode_cache(originalpacket);b=old.decode_cache(newpacket);unchanged=fingerprint(a,{'u'});unchanged_cache(a,b)
    if newpacket['parameters']!=originalpacket['parameters'] or newpacket['parameters']!=seed['parameters']:raise ValueError('material/geometry metadata changed')
    if newpacket['original_point_cache_file']!=row['original_point_cache_file'] or newpacket['original_point_cache_sha256']!=old.sha(old.path(row['original_point_cache_file'])):raise ValueError('old-cache lineage differs')
    if newpacket['input_state']!=root['proof_center_state_file'] or newpacket['input_state_sha256']!=root['proof_center_state_sha256']:raise ValueError('new numerical center identity differs')
    if old.q(newpacket['exact_base_gain_times_dose'])!=old.q(row['exact_requested_strength']) or old.q(newpacket['activation_strength_binary_exact'])!=old.q(row['solver_strength_binary_exact']):raise ValueError('exact drive or archived strength changed')
    for key in ('matrix_archive','matrix_archive_sha256'):
        if newpacket[key]!=originalpacket[key]:raise ValueError('original matrix lineage changed')
    old.assert_hashes({newpacket['matrix_archive']:newpacket['matrix_archive_sha256']})
    aux=old.read(old.path(root['proof_center_state_file']));center=old.read(old.path(root['center_refinement_file']))
    old.assert_hashes({root['proof_center_state_file']:root['proof_center_state_sha256'],root['center_refinement_file']:root['center_refinement_sha256'],center['center_arrays_file']:center['center_arrays_sha256']})
    if aux.get('kind')!='AUXILIARY_NUMERICAL_PROOF_CENTER_NOT_EXPERIMENT_OR_SOLVER_RESULT' or aux.get('converged') is not False or aux.get('solver_status')!='NO_NONLINEAR_SOLVER_RUN':raise ValueError('auxiliary center misrepresented as solver/data')
    if aux['original_seed_state_file']!=row['original_seed_state_file'] or aux['original_seed_state_sha256']!=row['original_seed_state_sha256']:raise ValueError('auxiliary seed identity changed')
    if aux['parameters']!=seed['parameters'] or not np.array_equal(np.asarray(aux['free_displacement']),b['u']) or not np.array_equal(np.asarray(aux['reference_nodes']),b['reference']) or not np.array_equal(np.asarray(aux['triangles']),b['triangles']):raise ValueError('auxiliary numeric field differs from new cache')
    serialized=b['reference'].copy().ravel();serialized[b['free']]+=b['u']
    if not np.array_equal(np.asarray(aux['current_nodes']).ravel(),serialized):raise ValueError('auxiliary serialization differs from declared numerical center')
    with np.load(old.path(root['predictor_file']),allow_pickle=False) as z:S=z['S'].copy()
    if old.sha(old.path(root['predictor_file']))!=old.sha(old.path(row['predictor_file'])):raise ValueError('binary predictor changed')
    centers,radii,gaps,hulls,free,translation=translated_box(seed,a,b,S,rootp['shared_parameter_cell'],rootp['base_parameters_exact'],radius)
    if list(map(Q,center['actual_center_translation_exact']))!=translation:raise ValueError('exact translation metadata differs')
    mids=list(map(Q,center['RF_midpoints_exact']))
    if len(mids)!=len(a['u']) or np.asarray([float(old.q(v)-m) for v,m in zip(a['u'],mids)],dtype=np.float64).tobytes()!=b['u'].tobytes():raise ValueError('declared single-round proof-center recipe differs')
    physical=full_potential_fingerprint(a,originalpacket['parameters'],rootp,row)
    if any(center[k]!=unchanged for k in ('cache_arrays_except_u_fingerprint','new_cache_arrays_except_u_fingerprint')) or any(center[k]!=physical for k in ('full_physical_potential_fingerprint','new_full_physical_potential_fingerprint')):raise ValueError('full target fingerprint differs')
    with np.load(old.path(center['center_arrays_file']),allow_pickle=False) as z:
        if not np.array_equal(z['old_u'],a['u']) or not np.array_equal(z['new_u'],b['u']):raise ValueError('center array archive differs')
    ids,endpoints=old.fixed_endpoints(seed,radii,gaps,free)
    body=old.certify_box(seed['reference_nodes'],centers,seed['triangles'],radii)
    if body['status']!='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING':return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Entire translated root box geometry not certified','body_geometry':body,'no_endpoint_certificate':True}
    contract=old.read(old.CONTRACT);L=Q(contract['physical_length_um_exact']);roi=list(map(Q,contract['roi_um_exact']));controls,H=old._state_controls(seed,L)
    transfer={}
    for name in ('top','bottom'):
        epsilon=L*max(radii[int(node)][j] for edge in ids[name] for node in edge for j in range(2))
        transfer[name]=old.curve_bound(controls[name],epsilon,fixed_radial_endpoints=True,roi_um=roi)
    if any(t['status']!='CONDITIONAL_EXACT_ERROR_TRANSFER' for t in transfer.values()):return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Translated graph/ROI bound unresolved','transfer':transfer,'no_endpoint_certificate':True}
    ea=Q(transfer['top']['Eulerian_vertical_error_um_exact']);eb=Q(transfer['bottom']['Eulerian_vertical_error_um_exact']);errors=dict(zip(old.OBSERVABLES,(ea+eb,ea+eb,ea,eb)))
    baseline=old.observe_state(seed,roi_um=roi,length_unit_um=L,reference_apical_um=L,reference_basal_um=0,weighting='current_disk',component_aggregation='spatial_max',tolerance_um='1/10000000')
    report=old.accept_observer_report(baseline)
    values=(baseline['four_endpoints']['maximum_Eulerian_height_shortening'],baseline['mean_alternatives']['current_disk']['height_regional_mean_shortening'],baseline['spatial_max_components']['apical_downward'],baseline['spatial_max_components']['basal_upward'])
    predictions={k:[v.lo-errors[k],v.hi+errors[k]] for k,v in zip(old.OBSERVABLES,values)}
    return {'status':'COMPOSED_UNIFORM_STATIC_ENDPOINTS_PRODUCER','adapter_kind':p['adapter_kind'],'record_id':row['record_id'],
        'root_certificate_sha256':old.sha(old.path(row['output_folder'])/'result.json'),'branch_domain_id':root['branch_domain_id'],
        'parameter_cell_sha256':p['parameter_cell_sha256'],'common_contract_sha256':p['common_contract_sha256'],
        'original_seed_state_file':row['original_seed_state_file'],'original_seed_state_sha256':row['original_seed_state_sha256'],
        'proof_center_state_file':root['proof_center_state_file'],'proof_center_state_sha256':root['proof_center_state_sha256'],
        'full_physical_potential_fingerprint':physical,'all_cache_arrays_except_u_fingerprint':unchanged,
        'full_NORMAL_corrector_radius_exact':radius,'exact_center_translation':translation,'exact_new_center_to_original_seed_gaps':gaps,
        'predictor_hulls_exact':hulls,'coordinate_radii_normalized_exact':radii,'fixed_radial_endpoints':endpoints,
        'body_geometry':body,'surface_error_transfer':transfer,'endpoint_error_um_exact':errors,'base_maximum_reporting':report,
        'base_endpoints_exact':{k:[v.lo,v.hi] for k,v in zip(old.OBSERVABLES,values)},'prediction_intervals_exact':predictions,
        'branch_law_semantics':'independent_stationary_root_per_record','center_translation_not_added_to_Banach_radius':True,
        'independent_root_and_application_review_required':True,'lineage':root['lineage']}

def run(review_file):
    p=old.read(OUT/'protocol.json');old.assert_hashes({**p['source_hashes'],**p['input_hashes']});review=old.read(old.path(review_file))
    if review.get('status')!='PASS' or not review.get('checks') or any(v is not True for v in review['checks'].values()) or review.get('observer_source_sha256')!=old.sha(__file__) or review.get('observer_protocol_sha256')!=old.sha(OUT/'protocol.json'):raise ValueError('matching corrected-center adapter prelaunch review required')
    rootp=old.read(ROOTRUN/'protocol.json');terminal=old.read(ROOTRUN/'results.json')
    if terminal['protocol_sha256']!=p['root_protocol_sha256'] or len(terminal['records'])!=2 or {r['record_id'] for r in terminal['records']}!=set(IDS):raise ValueError('terminal root tuple binding differs')
    locks={old.rel(ROOTRUN/'results.json'):old.sha(ROOTRUN/'results.json'),review_file:old.sha(old.path(review_file))}
    for row in rootp['records']:
        file=row['output_folder']+'/result.json';digest=old.sha(old.path(file));reported=next(r for r in terminal['records'] if r['record_id']==row['record_id'])
        if reported['file']!=file or reported['sha256']!=digest:raise ValueError('terminal root record hash differs')
        locks[file]=digest
    old.write(OUT/'composition_inputs.json',{'created_utc':old.stamp(),'observer_protocol_sha256':old.sha(OUT/'protocol.json'),'root_result_hashes':locks,'producer_results_pending_independent_review':True})
    records=[];predictions={}
    for row in rootp['records']:
        root=old.read(old.path(row['output_folder'])/'result.json')
        try:answer=compose(p,rootp,row,root)
        except (ValueError,KeyError,AssertionError) as error:answer={'status':'UNRESOLVED','record_id':row['record_id'],'reason':repr(error),'no_endpoint_certificate':True}
        output=OUT/'records'/(row['record_id']+'.json');old.write(output,answer)
        records.append({'record_id':row['record_id'],'status':answer['status'],'file':old.rel(output),'sha256':old.sha(output)})
        if 'prediction_intervals_exact' in answer:predictions[row['record_id']]=answer['prediction_intervals_exact']
        print(json.dumps({'record_id':row['record_id'],'status':answer['status']}),flush=True)
    old.assert_hashes(locks);old.assert_hashes({**p['source_hashes'],**p['input_hashes']})
    observations={k:v for k,v in old.load_observations().items() if k in IDS};aggregate=old.aggregate(observations,predictions)
    result={'status':'PRODUCER_TWO_RECORDS_COMPOSED' if aggregate['covered']==8 else 'PRODUCER_PARTIAL','records':records,'aggregate_F_K_only':aggregate,
        'protocol_sha256':old.sha(OUT/'protocol.json'),'composition_inputs_sha256':old.sha(OUT/'composition_inputs.json'),
        'scope':'Only F/K new centers to original seed graphs; existing60 outputs unchanged. No all68 or fullidentifiedset certificate from this adapter alone.',
        'independent_root_and_application_review_required':True}
    old.write(OUT/'results.json',result);return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','freeze','run']);parser.add_argument('--review');args=parser.parse_args()
    result=tests() if args.stage=='test' else freeze() if args.stage=='freeze' else run(args.review) if args.review else parser.error('--review required')
    print(json.dumps({'stage':args.stage,'status':result.get('status','FROZEN')}))
