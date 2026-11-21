"""Independent F/K center-offset application audit, not a root recomputation.

Boxes/trace constants/residuals rebuilt independently. Whole-body and baseline
static endpoints are DISCLOSED replays of reviewed components, not reimplemented.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,copy,csv,hashlib,io,json,time
import numpy as np
from root_body_geometry import certify_box
from quadratic_endpoint_observer import observe_state
from joint_calibration_contract import OBSERVABLES

ROOT=Path(__file__).resolve().parents[1]
PROD=ROOT/'runs/cohort_recentered_observers_v1'
ROOTRUN=ROOT/'runs/cohort_recentered_roots_v1'
ROOTREVIEW=ROOT/'evidence/cohort_recentered_roots_independent_review.json'
ROOTREVIEWRUN=ROOT/'runs/cohort_recentered_root_review_v1'
RUN=ROOT/'runs/recentered_observer_review_v1'
OUT=ROOT/'evidence/recentered_observers_independent_review.json'
ORDER=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
q=lambda x:Q.from_float(float(x)) if isinstance(x,(float,np.floating)) else Q(x)


def save(p,v):
    if p.exists():raise RuntimeError('retain review artifact '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')


def bindings(p):
    if any(sha(ROOT/k)!=v for k,v in p['file_hashes'].items()):raise RuntimeError('frozen review input/source changed')


def decode(c):
    raw=base64.b64decode(c['frozen_binary_cache']['data'])
    if hashlib.sha256(raw).hexdigest()!=c['frozen_binary_cache']['sha256']:raise RuntimeError('cache checksum mismatch')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}


def fingerprints(a,exclude):
    h=hashlib.sha256()
    for k in sorted(a):
        if k in exclude:continue
        v=np.ascontiguousarray(a[k]);h.update(k.encode()+b'\0'+v.dtype.str.encode()+b'\0'+json.dumps(list(v.shape)).encode()+b'\0'+v.tobytes())
    return h.hexdigest()


def boundary_edges(seed):
    edges={};ref=seed['reference_nodes']
    for a,b,c,d,e,f in seed['triangles']:
        for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):
            key=tuple(sorted((i,j)))
            if key in edges and edges[key]!=m:raise RuntimeError('nonconforming shared midpoint')
            edges[key]=m
    result={}
    for name,z in [('top',q(seed['parameters']['height'])),('bottom',Q(0))]:
        selected=[]
        for (i,j),m in edges.items():
            if q(ref[i][1])==q(ref[j][1])==z:
                if q(ref[i][0])>q(ref[j][0]):i,j=j,i
                selected.append((q(ref[i][0]),[i,m,j]))
        result[name]=[edge for _,edge in sorted(selected)]
    return result


def prepare():
    p=read(PROD/'protocol.json');rp=read(ROOTRUN/'protocol.json');rr=read(ROOTREVIEW);rrp=read(ROOTREVIEWRUN/'protocol.json')
    files={**p['source_hashes'],**p['input_hashes'],**rrp['file_hashes']}
    for path in [Path(__file__),PROD/'results.json',PROD/'results_scoped_v2.json',PROD/'scope_correction.json',PROD/'composition_inputs.json',PROD/'protocol.json',ROOTRUN/'results.json',ROOTREVIEW,ROOTREVIEWRUN/'protocol.json',ROOTREVIEWRUN/'preflight.json',ROOTREVIEWRUN/'prelaunch_review.json']:
        files[rel(path)]=sha(path)
    for entry in read(PROD/'results.json')['records']:
        path=ROOT/entry['file'];o=read(path);files[entry['file']]=sha(path)
        rrow=next(x for x in rp['records'] if x['record_id']==entry['record_id']);rpath=ROOT/rrow['output_folder']/'result.json';r=read(rpath);files[rel(rpath)]=sha(rpath);files.update(r['artifact_hashes'])
    protocol={'source_sha256':sha(__file__),'file_hashes':files,'scope':'F/K8 outputs only. Independent box/trace/residual arithmetic plus disclosed body/baseline component replay. No root assembly.',
        'root_review_file':rel(ROOTREVIEW),'root_review_sha256':sha(ROOTREVIEW),'root_review_protocol_sha256':sha(ROOTREVIEWRUN/'protocol.json'),
        'record_ids':['optoGEF2_small_col_F','optoGEF2_small_col_K'],'common_contract_sha256':p['common_contract_sha256'],'parameter_cell_sha256':p['parameter_cell_sha256'],
        'reporting_tolerance_um':'1/10000000','conditional_until_independent_new_roots_bound':True,'scope_metadata_issue':'Raw producer aggregate incorrectly inherited all68=true for8of8. Scopedv2 wrapper must reverseonlydeclaredmetadata changes and forbidglobal68claim.'}
    bindings(protocol)
    if rr['status']!='PASS' or not rr['checks'] or not all(v is True for v in rr['checks'].values()):raise RuntimeError('fresh independently reviewed new roots required')
    if rr['protocol_sha256']!=sha(ROOTREVIEWRUN/'protocol.json') or rr['review_source_sha256']!=sha(ROOT/'scripts/review_recentered_roots.py'):raise RuntimeError('new root review source/protocol mismatch')
    save(RUN/'protocol.json',protocol);print(json.dumps({'status':'REVIEW_PROTOCOL_FROZEN','source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),'root_review_sha256':sha(ROOTREVIEW)}),flush=True)


def arithmetic():
    start=time.time();p=read(RUN/'protocol.json');bindings(p);rp=read(ROOTRUN/'protocol.json');raw=read(PROD/'results.json');scoped=read(PROD/'results_scoped_v2.json');fix=read(PROD/'scope_correction.json')
    contract=read(ROOT/'runs/cohort_static_observers_v1/model_observation_contract.json');files=dict(p['file_hashes']);audits=[];residuals=[]
    with (ROOT/contract['calibration_source_file']).open(newline='',encoding='utf-8-sig') as f:observed={r['record_id']:r for r in csv.DictReader(f)}
    for entry in raw['records']:
        rid=entry['record_id'];row=next(x for x in rp['records'] if x['record_id']==rid);o=read(ROOT/entry['file']);rfile=ROOT/row['output_folder']/'result.json';root=read(rfile)
        seed=read(ROOT/row['original_seed_state_file']);oldpc=read(ROOT/row['original_point_cache_file']);newpc=read(ROOT/root['point_cache_file']);a=decode(oldpc);b=decode(newpc)
        aux=read(ROOT/root['proof_center_state_file']);center=read(ROOT/root['center_refinement_file'])
        with np.load(ROOT/root['predictor_file'],allow_pickle=False) as z:S=z['S'].copy()
        checks={'observer_and_new_root_bytes_bound':sha(ROOT/entry['file'])==entry['sha256'] and o['root_certificate_sha256']==sha(rfile),
            'new_root_same_domain_cell_contract_and_source_dose':root['record_id']==rid and o['branch_domain_id']==root['branch_domain_id']==row['branch_domain_id'] and root['parameter_cell_sha256']==o['parameter_cell_sha256']==p['parameter_cell_sha256'] and root['model_observation_contract_sha256']==o['common_contract_sha256']==p['common_contract_sha256'] and Q(root['dose_exact'])==Q(row['dose_exact']),
            'all_cache_arrays_except_u_byteidentical_including_R':set(a)==set(b) and all(a[k].dtype==b[k].dtype and a[k].shape==b[k].shape and a[k].tobytes()==b[k].tobytes() for k in a if k!='u'),
            'all_old_new_cache_arrays_finite':all(np.isfinite(v).all() for x in (a,b) for v in x.values()),
            'original_seed_matches_original_cache':np.asarray(seed['free_displacement'],dtype=np.float64).tobytes()==a['u'].tobytes() and np.asarray(seed['reference_nodes'],dtype=np.float64).tobytes()==a['reference'].tobytes(),
            'S_bytes_unchanged':sha(ROOT/root['predictor_file'])==sha(ROOT/row['predictor_file']),
            'new_auxiliary_identity_not_seed_or_solver_result':newpc['input_state']==root['proof_center_state_file'] and newpc['input_state_sha256']==sha(ROOT/root['proof_center_state_file']) and root['proof_center_state_file']!=row['original_seed_state_file'] and aux['converged'] is False and aux['kind']=='AUXILIARY_NUMERICAL_PROOF_CENTER_NOT_EXPERIMENT_OR_SOLVER_RESULT',
            'new_auxiliary_u_matches_new_cache':np.asarray(aux['free_displacement'],dtype=np.float64).tobytes()==b['u'].tobytes(),
            'physical_parameters_reference_and_mesh_unchanged':oldpc['parameters']==newpc['parameters']==seed['parameters']==aux['parameters'] and aux['reference_nodes']==seed['reference_nodes'] and aux['triangles']==seed['triangles']}
        physical=lambda x:hashlib.sha256(json.dumps({'potential_arrays':fingerprints(x,{'u','preconditioner'}),'model_parameters':oldpc['parameters'],'theta0_exact':rp['base_parameters_exact'],'dose_exact':row['dose_exact'],'shared_parameter_cell':rp['shared_parameter_cell']},sort_keys=True,separators=(',',':')).encode()).hexdigest()
        checks['full_physical_family_fingerprints_independently_reproduced']=physical(a)==physical(b)==o['full_physical_potential_fingerprint']==center['full_physical_potential_fingerprint']==center['new_full_physical_potential_fingerprint'] and fingerprints(a,{'u'})==o['all_cache_arrays_except_u_fingerprint']
        normal=Q(root['uniform_cell_result']['selected_radius_exact']);attempt=next(x for x in root['uniform_cell_result']['radius_attempts'] if Q(x['radius_exact'])==normal and x['strict_selfmap_contraction'])
        checks['normal_full_corrector_not_translation_or_posterior']=normal==Q(o['full_NORMAL_corrector_radius_exact']) and Q(attempt['Z_uniform']['upper_exact'])<1 and Q(root['uniform_cell_result']['eta_uniform']['upper_exact'])+Q(attempt['Z_uniform']['upper_exact'])*normal<normal and o['center_translation_not_added_to_Banach_radius'] is True
        delta=[(Q(root['parameter_cell'][k][0])-Q(rp['base_parameters_exact'][k]),Q(root['parameter_cell'][k][1])-Q(rp['base_parameters_exact'][k])) for k in ORDER];hulls=[]
        for coef in S:
            terms=[(q(s)*lo,q(s)*hi) for s,(lo,hi) in zip(coef,delta)];hulls.append((sum(min(v) for v in terms),sum(max(v) for v in terms)))
        free=list(map(int,b['free']));lookup={k:i for i,k in enumerate(free)};nominal=[[q(v) for v in pair] for pair in seed['current_nodes']];refs=[[q(v) for v in pair] for pair in b['reference']]
        radii=[];gaps=[];translations=[q(v)-q(u) for u,v in zip(a['u'],b['u'])]
        for k,v in enumerate(b['reference'].ravel()):
            actual=q(v)+(q(b['u'][lookup[k]]) if k in lookup else Q(0));gap=abs(actual-nominal[k//2][k%2]);gaps.append(gap)
            radii.append(gap+(max(map(abs,hulls[lookup[k]]))+normal if k in lookup else Q(0)))
        checks['all611_exact_predictor_hulls_rebuilt']=list(map(lambda x:tuple(map(Q,x)),o['predictor_hulls_exact']))==hulls
        checks['all650_exact_new_center_to_original_seed_gaps_rebuilt']=[Q(v) for pair in o['exact_new_center_to_original_seed_gaps'] for v in pair]==gaps
        checks['all650_total_nodal_radii_rebuilt']=[Q(v) for pair in o['coordinate_radii_normalized_exact'] for v in pair]==radii
        checks['exact611_translations_and_single_round_recipe_verified']=list(map(Q,o['exact_center_translation']))==list(map(Q,center['actual_center_translation_exact']))==translations and np.asarray([float(q(u)-Q(m)) for u,m in zip(a['u'],center['RF_midpoints_exact'])],dtype=np.float64).tobytes()==b['u'].tobytes()
        body=certify_box(refs,nominal,seed['triangles'],[radii[i:i+2] for i in range(0,len(radii),2)])
        checks['disclosed_reviewed_body_component_replay']=body==o['body_geometry'] and body['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING'
        edges=boundary_edges(seed);L=Q(contract['physical_length_um_exact']);roi=list(map(Q,contract['roi_um_exact']));errors={}
        for name,triples in edges.items():
            eps=L*max(radii[2*n+j] for nodes in triples for n in nodes for j in (0,1));minimum=None;lip=Q(0)
            for nodes in triples:
                points=[tuple(L*q(v) for v in seed['current_nodes'][n]) for n in nodes];linear=[-3*points[0][j]+4*points[1][j]-points[2][j] for j in (0,1)];quad=[2*points[0][j]-4*points[1][j]+2*points[2][j] for j in (0,1)]
                mr=min(linear[0],linear[0]+2*quad[0]);minimum=mr if minimum is None else min(minimum,mr);lip=max(lip,max(abs(linear[1]),abs(linear[1]+2*quad[1]))/mr)
            checks[name+'_fixed_radial_domain_endpoints']=all(2*n not in lookup and radii[2*n]==0 for n in (triples[0][0],triples[-1][-1]))
            E=Q(5,4)*eps*(1+lip);errors[name]=E;t=o['surface_error_transfer'][name]
            checks[name+'_independent_exact_Eulerian_transfer']=minimum>8*eps and eps==Q(t['nodal_coordinate_error_um_exact']) and lip==Q(t['base_graph_Lipschitz_upper_exact']) and E==Q(t['Eulerian_vertical_error_um_exact']) and roi==list(map(Q,t['ROI_um_exact'])) and roi[0]<roi[1]
        ea,eb=errors['top'],errors['bottom'];E=dict(zip(OBSERVABLES,(ea+eb,ea+eb,ea,eb)))
        checks['all_four_endpoint_transport_assignments_correct']=all(Q(o['endpoint_error_um_exact'][k])==v for k,v in E.items())
        baseline=observe_state(seed,roi_um=roi,length_unit_um=L,reference_apical_um=L,reference_basal_um=0,weighting='current_disk',component_aggregation='spatial_max',tolerance_um=p['reporting_tolerance_um'])
        vals=(baseline['four_endpoints']['maximum_Eulerian_height_shortening'],baseline['mean_alternatives']['current_disk']['height_regional_mean_shortening'],baseline['spatial_max_components']['apical_downward'],baseline['spatial_max_components']['basal_upward'])
        checks['disclosed_reviewed_all_four_base_endpoint_replay']=all(tuple(map(Q,o['base_endpoints_exact'][k]))==(v.lo,v.hi) for k,v in zip(OBSERVABLES,vals))
        checks['actual_maximum_reporting_width_preserved']=Q(o['base_maximum_reporting']['actual_interval_width_um_exact'])==baseline['maximum_shortening_certificate']['interval'].hi-baseline['maximum_shortening_certificate']['interval'].lo
        for key,v in zip(OBSERVABLES,vals):
            lo,hi=map(Q,o['prediction_intervals_exact'][key]);checks[key+'_uniform_interval_exact']=lo==v.lo-E[key] and hi==v.hi+E[key];y=Q(observed[rid][key])
            residuals.append({'record_id':rid,'observable':key,'observed_exact':str(y),'prediction_interval_exact':[str(lo),str(hi)],'minimum_residual_exact':str(max(lo-y,y-hi,Q(0))),'maximum_residual_exact':str(max(abs(lo-y),abs(hi-y)))})
        audits.append({'record_id':rid,'status':'PASS_ARITHMETIC_CONDITIONAL_ON_NEW_ROOT_REVIEW' if all(checks.values()) else 'FAIL','checks':checks,'root_certificate_sha256':sha(rfile),'observer_certificate_sha256':sha(ROOT/entry['file']),'branch_domain_id':root['branch_domain_id']})
        print(json.dumps({'record':rid,'checks':len(checks),'status':audits[-1]['status']}),flush=True)
    lower=max(Q(v['minimum_residual_exact']) for v in residuals);upper=max(Q(v['maximum_residual_exact']) for v in residuals)
    keys=('record_id','observable','observed_exact','prediction_interval_exact','minimum_residual_exact','maximum_residual_exact')
    reported=[{k:x[k] for k in keys} for x in raw['aggregate_F_K_only']['rows']]
    undo=copy.deepcopy(scoped);undo.pop('metadata_scope_correction');ag=undo['aggregate_F_K_only'];ag['all68_uniform_static_tuple_enclosed']=True
    for k in ('all_required_FK_channels_uniform_static_tuple_enclosed','required_F_K_observables','global_required_observables','global_all68_IN_certified'):ag.pop(k)
    ag['common_delta_family'].pop('IN_static_stationary_root_tuple_scope')
    scope=scoped['aggregate_F_K_only'];checks={'all_two_records_application_checks_pass':len(audits)==2 and all(all(v is True for v in a['checks'].values()) for a in audits),
        'all8_exact_residuals_reproduced':reported==residuals and len(residuals)==8,
        'exact_FK_only_residual_aggregate_reproduced':lower==Q(scope['maximum_residual_lower_exact']) and upper==Q(scope['maximum_residual_upper_exact']),
        'scope_issue_raw_legacy_flag_identified':raw['aggregate_F_K_only']['all68_uniform_static_tuple_enclosed'] is True,
        'metadata_wrapper_reverses_exactly_to_original_packet':undo==raw and scoped['metadata_scope_correction']['original_result_sha256']==sha(PROD/'results.json'),
        'authoritative_wrapper_is_only8_not_all68':scope['all68_uniform_static_tuple_enclosed'] is False and scope['all_required_FK_channels_uniform_static_tuple_enclosed'] is True and scope['required_F_K_observables']==scope['covered']==scope['required']==8 and scope['global_required_observables']==68 and scope['global_all68_IN_certified'] is False,
        'scope_correction_receipt_bound_and_allchecks_true':fix['scoped_result_sha256']==sha(PROD/'results_scoped_v2.json') and fix['original_result_sha256']==sha(PROD/'results.json') and all(v is True for v in fix['checks'].values()),
        'all_review_inputs_unchanged':all(sha(ROOT/k)==v for k,v in files.items())}
    save(RUN/'arithmetic.json',{'status':'PASS_ARITHMETIC_PENDING_INDEPENDENT_ROOT_BINDING' if all(checks.values()) else 'FAIL','checks':checks,'record_audits':audits,'exact_residuals':residuals,
        'FK_only_lower_exact':str(lower),'FK_only_upper_exact':str(upper),'file_hashes':files,'elapsed_seconds':time.time()-start,
        'disclosure':'Body and base static-observer evaluations are replays of reviewed components, not independent reimplementations. Node hulls, translated gaps, trace derivatives/transport and residual arithmetic rebuilt independently.'})
    print(json.dumps({'status':read(RUN/'arithmetic.json')['status'],'checks':sum(len(v['checks']) for v in audits)+len(checks),'elapsed':time.time()-start}),flush=True)
    if not all(checks.values()):raise RuntimeError('application arithmetic finding retained')


def finalize():
    p=read(RUN/'protocol.json');bindings(p);a=read(RUN/'arithmetic.json');rootreview=read(ROOTREVIEW);rp=read(ROOTREVIEWRUN/'protocol.json')
    checks={'arithmetic_gate_passed':a['status']=='PASS_ARITHMETIC_PENDING_INDEPENDENT_ROOT_BINDING','new_independent_root_batch_PASS':rootreview['status']=='PASS' and bool(rootreview['checks']) and all(v is True for v in rootreview['checks'].values()),
        'new_independent_root_receipt_source_protocol_bound':sha(ROOTREVIEW)==p['root_review_sha256'] and rootreview['protocol_sha256']==p['root_review_protocol_sha256']==sha(ROOTREVIEWRUN/'protocol.json') and rootreview['review_source_sha256']==sha(ROOT/'scripts/review_recentered_roots.py'),
        'new_rootreview_same_cell_and_common_contract':rootreview['parameter_cell_sha256']==p['parameter_cell_sha256'] and rootreview['model_observation_contract_sha256']==p['common_contract_sha256']}
    capsule={v['record_id']:v for v in rootreview['record_audits']};out=[]
    for record in a['record_audits']:
        rid=record['record_id'];entry=capsule[rid];r=entry['audit']
        checks[rid+'_fresh_root_capsule_PASS_and_all_checks_literal_true']=r['status']=='PASS' and bool(r['checks']) and all(v is True for v in r['checks'].values()) and r['uniform_quantifier_verified'] is True
        checks[rid+'_fresh_root_capsule_exact_certificate_domain_cell_contract']=r['certificate_sha256']==record['root_certificate_sha256'] and r['branch_domain_id']==record['branch_domain_id'] and r['parameter_cell_sha256']==p['parameter_cell_sha256'] and r['model_observation_contract_sha256']==p['common_contract_sha256'] and r['branch_law_semantics']=='independent_stationary_root_per_record'
        checks[rid+'_individual_audit_file_and_protocol_bound']=sha(ROOT/entry['review_file'])==entry['review_sha256'] and r['review_protocol_sha256']==p['root_review_protocol_sha256']
        out.append({**record,'status':'PASS_CONDITIONAL_UNIFORM_OBSERVER_APPLICATION','independent_new_root_review_sha256':sha(ROOTREVIEW),'parameter_cell_sha256':p['parameter_cell_sha256'],
            'independent_root_review_file':rel(ROOTREVIEW),'independent_root_review_sha256':sha(ROOTREVIEW),
            'common_model_observation_contract_sha256':p['common_contract_sha256'],'common_contract_sha256':p['common_contract_sha256'],'uniform_quantifier_verified':True,'forward_family_admissible':True,'branch_law_semantics':'independent_stationary_root_per_record'})
    checks['all_frozen_arithmetic_and_root_inputs_current']=all(sha(ROOT/k)==v for k,v in a['file_hashes'].items()) and all(sha(ROOT/k)==v for k,v in rp['file_hashes'].items())
    receipt={'status':'PASS_F_K_EIGHT_OUTPUT_APPLICATION' if all(checks.values()) else 'FAIL','checks':checks,'record_audits':out,'arithmetic_sha256':sha(RUN/'arithmetic.json'),
        'independent_new_root_review_file':rel(ROOTREVIEW),'independent_new_root_review_sha256':sha(ROOTREVIEW),'root_review_protocol_sha256':p['root_review_protocol_sha256'],
        'independent_root_review_file':rel(ROOTREVIEW),'independent_root_review_sha256':sha(ROOTREVIEW),
        'independent_root_review_bindings':{rel(path):sha(path) for path in (ROOTREVIEW,ROOTREVIEWRUN/'protocol.json',ROOTREVIEWRUN/'preflight.json',ROOTREVIEWRUN/'prelaunch_review.json',ROOT/'scripts/review_recentered_roots.py')},
        'authoritative_scoped_results_sha256':sha(PROD/'results_scoped_v2.json'),'scope_correction_sha256':sha(PROD/'scope_correction.json'),'raw_results_sha256':sha(PROD/'results.json'),
        'parameter_cell_sha256':p['parameter_cell_sha256'],'common_model_observation_contract_sha256':p['common_contract_sha256'],'covered_observables':8,'required_scope_observables':8,
        'global_required_observables':68,'global_all68_fit_or_IN_certified':False,'FK_only_necessary_common_discrepancy_lower_exact':a['FK_only_lower_exact'],'FK_only_sufficient_common_discrepancy_upper_exact':a['FK_only_upper_exact'],
        'review_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),
        'scope':'Verified only F/K8 conditional static outputs on independently reviewed new root domains. Root owns combining60+8. This receipt alone is not all68/globalfit, physicalidentifiedset, trajectory, stability, PDE or novelty.',
        'disclosure':a['disclosure']}
    save(OUT,receipt);print(json.dumps({'status':receipt['status'],'final_checks':len(checks),'covered':8,'global68_IN':False,'receipt_sha256':sha(OUT)}),flush=True)
    if not all(checks.values()):raise RuntimeError('final new-root binding finding retained')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','arithmetic','finalize']);args=parser.parse_args();{'prepare':prepare,'arithmetic':arithmetic,'finalize':finalize}[args.stage]()
