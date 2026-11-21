"""Uniform static endpoints for an independently parameterized stationary-root tuple.

No equilibrium computation is performed here. Producer root certificates are
composed conditionally and must subsequently receive independent review.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,base64,datetime,hashlib,io,json
import numpy as np
from quadratic_endpoint_observer import _state_controls,observe_state
from quadratic_mean_derivatives import _boundary_node_ids
from nodal_observer_error_transfer import curve_bound
from root_body_geometry import certify_box
from joint_calibration_contract import OBSERVABLES,load_observations,interval_residual

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/cohort_static_observers_v1'
ROOTRUN=ROOT/'runs/cohort_uniform_roots_v1'
CONTRACT=OUT/'model_observation_contract.json'
ORDER=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def rel(p):return Path(p).resolve().relative_to(ROOT.resolve()).as_posix()
def path(p):
    target=(ROOT/p).resolve()
    target.relative_to(ROOT.resolve())
    return target
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def serial(v):
    if isinstance(v,Q):return str(v)
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [serial(x) for x in v]
    if hasattr(v,'lo') and hasattr(v,'hi'):return [str(v.lo),str(v.hi)]
    if isinstance(v,np.integer):return int(v)
    return v
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    data=json.dumps(serial(v),indent=2,allow_nan=False)+'\n'
    if p.exists():
        if p.read_text(encoding='utf8')!=data:raise ValueError('immutable artifact already exists: '+str(p))
    else:p.write_text(data,encoding='utf8')
def q(v):
    if isinstance(v,bool):raise ValueError('boolean is not numeric')
    return Q.from_float(float(v)) if isinstance(v,(float,np.floating)) else Q(v)
def canonical(c):return {k:[str(q(v)) for v in c[k]] for k in sorted(c)}
def cell_hash(c):return hashlib.sha256(json.dumps(canonical(c),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def assert_hashes(bindings):
    for name,digest in bindings.items():
        if sha(path(name))!=digest:raise ValueError('changed bound artifact: '+name)

def dependencies():
    """Bind the recursive local runtime import closure, without importing solvers."""
    found={};pending=[Path(__file__)]
    while pending:
        p=pending.pop()
        if rel(p) in found:continue
        found[rel(p)]=sha(p)
        for node in ast.walk(ast.parse(p.read_text(encoding='utf8'))):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                if name:
                    child=ROOT/'scripts'/(name.split('.')[0]+'.py')
                    if child.exists():pending.append(child)
    return found

def predictor_hulls(S,cell,base):
    if np.asarray(S).ndim!=2 or np.asarray(S).shape[1]!=5 or not np.isfinite(S).all():raise ValueError('finite n by5 predictor required')
    bounds=[(q(cell[k][0])-q(base[k]),q(cell[k][1])-q(base[k])) for k in ORDER]
    if any(lo>hi for lo,hi in bounds):raise ValueError('ordered parameter box required')
    answer=[]
    for row in S:
        terms=[(q(s)*lo,q(s)*hi) for s,(lo,hi) in zip(row,bounds)]
        answer.append((sum((min(t) for t in terms),Q(0)),sum((max(t) for t in terms),Q(0))))
    return answer

def accepted_radius(result):
    if result['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT' or result['method']!='affine':raise ValueError('certified affine family required')
    radius=q(result['selected_radius_exact']);eta=q(result['eta_uniform']['upper_exact'])
    matches=[a for a in result['radius_attempts'] if q(a['radius_exact'])==radius and a.get('strict_selfmap_contraction') is True]
    if radius<=0 or len(matches)!=1:raise ValueError('one accepted matching full corrector radius required')
    Z=q(matches[0]['Z_uniform']['upper_exact'])
    if not 0<=Z<1 or eta<0 or not eta+Z*radius<radius:raise ValueError('strict rational Banach inequalities failed')
    if q(matches[0]['eta_plus_Zr']['upper_exact'])>=radius:raise ValueError('reported selfmap bound failed')
    return radius

def decode_cache(packet):
    frozen=packet['frozen_binary_cache'];raw=base64.b64decode(frozen['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=frozen['sha256']:raise ValueError('cache bytes changed')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
    if any(not np.isfinite(v).all() for v in arrays.values()):raise ValueError('nonfinite cache')
    return arrays

def nodal_box(state,arrays,S,cell,base,radius):
    ref=np.asarray(state['reference_nodes']);tri=np.asarray(state['triangles']);u=np.asarray(state['free_displacement'])
    if not np.array_equal(ref,arrays['reference']) or not np.array_equal(tri,arrays['triangles']) or not np.array_equal(u,arrays['u']):raise ValueError('cache versus saved field mismatch')
    free=[int(x) for x in arrays['free']]
    if len(free)!=len(u) or len(set(free))!=len(free) or min(free)<0 or max(free)>=ref.size or S.shape!=(len(free),5):raise ValueError('free-coordinate map mismatch')
    centers=[[q(v) for v in row] for row in state['current_nodes']]
    exact=[[q(v) for v in row] for row in state['reference_nodes']]
    for i,k in enumerate(free):exact[k//2][k%2]+=q(u[i])
    hulls=predictor_hulls(S,cell,base)
    radii=[[abs(exact[i][j]-centers[i][j]) for j in range(2)] for i in range(len(ref))]
    gaps=[r[:] for r in radii]
    for i,k in enumerate(free):radii[k//2][k%2]+=max(map(abs,hulls[i]))+radius
    return centers,radii,gaps,hulls,free

def fixed_endpoints(state,radii,gaps,free):
    ids=_boundary_node_ids(state);result={}
    for name in ('top','bottom'):
        endpoints=[ids[name][0][0],ids[name][-1][-1]]
        for node in endpoints:
            if 2*node in free or radii[node][0]!=0 or gaps[node][0]!=0:raise ValueError('radial endpoint not exactly fixed')
        result[name]=endpoints
    return ids,result

def aggregate(observations,predictions):
    if set(predictions)-set(observations):raise ValueError('unexpected source record')
    rows=[];missing=[]
    for rid,obs in observations.items():
        if set(predictions.get(rid,{}))-set(OBSERVABLES):raise ValueError('unexpected endpoint')
        for key in OBSERVABLES:
            if key not in predictions.get(rid,{}):missing.append([rid,key]);continue
            lo,hi=map(q,predictions[rid][key]);lower,upper=interval_residual(obs[key],[lo,hi])
            rows.append({'record_id':rid,'observable':key,'observed_exact':obs[key],'prediction_interval_exact':[lo,hi],
                         'minimum_residual_exact':lower,'maximum_residual_exact':upper})
    lower=max((r['minimum_residual_exact'] for r in rows),default=Q(0))
    upper=max((r['maximum_residual_exact'] for r in rows),default=Q(0)) if not missing else None
    return {'rows':rows,'missing_rows':missing,'covered':len(rows),'required':4*len(observations),
        'maximum_residual_lower_exact':lower,'maximum_residual_upper_exact':upper,
        'lower_scope':'Every hypothetical complete tuple extending the certified selected record branches has D>=lower.',
        'common_delta_family':{'OUT_selected_root_tuple_if':'delta < lower',
            'IN_static_stationary_root_tuple_if':'delta >= upper' if upper is not None else None,
            'otherwise':'UNRESOLVED','biological_tolerance_assumed':False},
        'all68_uniform_static_tuple_enclosed':not missing}

def tests():
    checks={};C={k:['-1','2'] for k in ORDER};base={k:'0' for k in ORDER}
    S=np.array([[1.,-2.,0.,.5,-.25]])
    hull=predictor_hulls(S,C,base)[0]
    import itertools
    vals=[sum(q(s)*v for s,v in zip(S[0],vtx)) for vtx in itertools.product((Q(-1),Q(2)),repeat=5)]
    checks['affine_hull_matches_all32_vertices']=hull==(min(vals),max(vals))
    obs={'B':{k:Q(2) for k in OBSERVABLES}}
    full=aggregate(obs,{'B':{k:['3','4'] for k in OBSERVABLES}})
    checks['all_four_residuals_exact']=full['maximum_residual_lower_exact']==1 and full['maximum_residual_upper_exact']==2
    partial=aggregate(obs,{'B':{OBSERVABLES[0]:['3','4']}})
    checks['partial_lower_retained']=partial['maximum_residual_lower_exact']==1
    checks['partial_upper_forbidden']=partial['maximum_residual_upper_exact'] is None and partial['common_delta_family']['IN_static_stationary_root_tuple_if'] is None
    for i,(y,I,expected) in enumerate([(2,[1,3],(0,1)),(4,[1,3],(1,3)),(0,[1,3],(1,3))]):checks['absolute_interval_'+str(i)]=interval_residual(y,I)==expected
    trial={'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','method':'affine','selected_radius_exact':'1/100','eta_uniform':{'upper_exact':'1/1000'},
        'radius_attempts':[{'radius_exact':'1/100','strict_selfmap_contraction':True,'Z_uniform':{'upper_exact':'1/10'},'eta_plus_Zr':{'upper_exact':'1/500'}}]}
    checks['full_corrector_not_posterior_radius']=accepted_radius(trial)==Q(1,100)
    import copy
    bad=[]
    t=copy.deepcopy(trial);t['radius_attempts'][0]['Z_uniform']['upper_exact']='1';bad.append(lambda t=t:accepted_radius(t))
    t=copy.deepcopy(trial);t['selected_radius_exact']='1/10000';bad.append(lambda t=t:accepted_radius(t))
    bad.extend([lambda:predictor_hulls(np.ones((2,4)),C,base),lambda:aggregate(obs,{'X':{}}),lambda:q(True)])
    for i,call in enumerate(bad):
        try:call();checks['invalid_'+str(i)]=False
        except (ValueError,KeyError):checks['invalid_'+str(i)]=True
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),'no_root_result_consumption':True}
    write(OUT/'synthetic_checks.json',result)
    if result['status']!='PASS':raise ValueError('synthetic failure')
    return result

def freeze():
    rootp=read(ROOTRUN/'protocol.json');contract=read(CONTRACT)
    if rootp['common_model_observation_contract_sha256']!=sha(CONTRACT) or contract['parameter_cell_sha256']!=cell_hash(rootp['shared_parameter_cell']):raise ValueError('shared contract/cell mismatch')
    sources=dependencies();sources.update({'scripts/'+k:v for k,v in rootp['source_hashes'].items()})
    inputs={rel(ROOTRUN/'protocol.json'):sha(ROOTRUN/'protocol.json'),rel(CONTRACT):sha(CONTRACT),**rootp['input_hashes']}
    for name in ('evidence/cohort_uniform_roots_prelaunch_review.json',):
        if path(name).exists():inputs[name]=sha(path(name))
    assert_hashes({**sources,**inputs})
    checks={'all17_distinct':len(rootp['records'])==17 and len({r['record_id'] for r in rootp['records']})==17,
        'all_source_records':{r['record_id'] for r in rootp['records']}==set(load_observations()),'no_root_results_read':True}
    for row in rootp['records']:
        state=read(path(row['saved_state']));controls,H=_state_controls(state,Q(contract['physical_length_um_exact']))
        ids=_boundary_node_ids(state)
        checks[row['record_id']+'_saved_controls']=H==Q(contract['physical_length_um_exact']) and all(len(ids[k])==12 for k in ('top','bottom'))
    if not all(checks.values()):raise ValueError('preflight failure')
    protocol={'version':1,'created_utc':stamp(),'source_hashes':sources,'input_hashes':inputs,
        'root_protocol_file':rel(ROOTRUN/'protocol.json'),'root_protocol_sha256':sha(ROOTRUN/'protocol.json'),
        'common_contract_file':rel(CONTRACT),'common_contract_sha256':sha(CONTRACT),'parameter_cell_sha256':contract['parameter_cell_sha256'],
        'static_endpoint_tolerance_um_exact':'1/10000000','full_corrector_cube_required':True,'producer_results_require_independent_review':True,
        'residual_lower_formula':'max(l-y,y-h,0)','residual_upper_formula':'max(abs(l-y),abs(h-y))',
        'aggregation':'D(theta)=max_{17 records,4 endpoints}|prediction(theta)-observation|',
        'branch_law_semantics':'independent_stationary_root_per_record','missing_row_policy':'retain lower; no sufficient upper or IN if any of68 missing',
        'prelaunch_required_fields':['status=PASS','checks all boolean true','observer_source_sha256','observer_protocol_sha256'],
        'no_new_physical_solve':True,'no_biological_tolerance_assumed':True}
    write(OUT/'protocol.json',protocol)
    write(OUT/'preflight.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    return protocol

def compose_record(protocol,rootp,contract,row,root):
    expected={'record_id':row['record_id'],'protocol_sha256':protocol['root_protocol_sha256'],'parameter_cell_sha256':protocol['parameter_cell_sha256'],
        'model_observation_contract_sha256':protocol['common_contract_sha256'],'dose_exact':row['dose_exact'],'base_strength_exact':row['exact_requested_strength'],
        'saved_state_sha256':row['saved_state_sha256'],'branch_domain_id':row['branch_domain_id'],'input_kind':'source_record_dose',
        'strength_relation':'common_gain_times_source_dose','branch_law_semantics':'independent_stationary_root_per_record','common_path_or_energetic_stability_certified':False}
    if any(root.get(k)!=v for k,v in expected.items()):raise ValueError('root context mismatch')
    if canonical(root['parameter_cell'])!=canonical(rootp['shared_parameter_cell']):raise ValueError('cell mismatch')
    assert_hashes(root.get('artifact_hashes',{}));assert_hashes({'scripts/'+k:v for k,v in root['source_hashes'].items()})
    if root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT':return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Root producer unresolved','root_status':root['status']}
    result=root['uniform_cell_result']
    if canonical(result['cell_bounds_exact'])!=canonical(rootp['shared_parameter_cell']) or result['dose_exact']!=row['dose_exact']:raise ValueError('inner uniform cell mismatch')
    radius=accepted_radius(result);state=read(path(row['saved_state']));packet=read(path(root['point_cache_file']));arrays=decode_cache(packet)
    with np.load(path(root['predictor_file']),allow_pickle=False) as z:S=z['S'].copy()
    centers,radii,gaps,hulls,free=nodal_box(state,arrays,S,rootp['shared_parameter_cell'],rootp['base_parameters_exact'],radius)
    ids,endpoints=fixed_endpoints(state,radii,gaps,free)
    body=certify_box(state['reference_nodes'],centers,state['triangles'],radii)
    if body['status']!='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING':return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Whole-body sufficient geometry bound failed','body':body}
    L=Q(contract['physical_length_um_exact']);roi=list(map(Q,contract['roi_um_exact']));controls,H=_state_controls(state,L)
    transfer={}
    for name in ('top','bottom'):
        epsilon=L*max(radii[int(n)][j] for edge in ids[name] for n in edge for j in range(2))
        transfer[name]=curve_bound(controls[name],epsilon,fixed_radial_endpoints=True,roi_um=roi)
    if any(t['status']!='CONDITIONAL_EXACT_ERROR_TRANSFER' for t in transfer.values()):return {'status':'UNRESOLVED','record_id':row['record_id'],'reason':'Radial/ROI error transfer unresolved','transfer':transfer}
    ea=Q(transfer['top']['Eulerian_vertical_error_um_exact']);eb=Q(transfer['bottom']['Eulerian_vertical_error_um_exact'])
    errors=dict(zip(OBSERVABLES,(ea+eb,ea+eb,ea,eb)))
    observed=observe_state(state,roi_um=roi,length_unit_um=L,reference_apical_um=L,reference_basal_um=0,
        weighting='current_disk',component_aggregation='spatial_max',tolerance_um=protocol['static_endpoint_tolerance_um_exact'])
    if not observed['maximum_shortening_certificate']['requested_width_achieved']:raise ValueError('base max observer requested width not achieved')
    values=(observed['four_endpoints']['maximum_Eulerian_height_shortening'],observed['mean_alternatives']['current_disk']['height_regional_mean_shortening'],
            observed['spatial_max_components']['apical_downward'],observed['spatial_max_components']['basal_upward'])
    predictions={k:[v.lo-errors[k],v.hi+errors[k]] for k,v in zip(OBSERVABLES,values)}
    return {'status':'COMPOSED_UNIFORM_STATIC_ENDPOINTS_PRODUCER','record_id':row['record_id'],'dose_exact':row['dose_exact'],'branch_domain_id':row['branch_domain_id'],
        'root_result_sha256':sha(path(row['output_folder'])/'result.json'),'parameter_cell_sha256':protocol['parameter_cell_sha256'],
        'common_contract_sha256':protocol['common_contract_sha256'],'full_corrector_radius_exact':radius,
        'predictor_hulls_exact':hulls,'serialization_gaps_exact':gaps,'coordinate_radii_normalized_exact':radii,'fixed_radial_endpoints':endpoints,
        'body_geometry':body,'surface_error_transfer':transfer,'endpoint_error_um_exact':errors,'base_endpoints_exact':{k:[v.lo,v.hi] for k,v in zip(OBSERVABLES,values)},
        'prediction_intervals_exact':predictions,'point_cache_sha256':sha(path(root['point_cache_file'])),'predictor_sha256':sha(path(root['predictor_file'])),
        'no_common_physical_path_or_stability_claim':True,'independent_review_required':True}

def run(review_file):
    protocol=read(OUT/'protocol.json');assert_hashes({**protocol['source_hashes'],**protocol['input_hashes']})
    review=read(path(review_file))
    if review.get('status')!='PASS' or not review.get('checks') or any(v is not True for v in review['checks'].values()) or review.get('observer_source_sha256')!=sha(__file__) or review.get('observer_protocol_sha256')!=sha(OUT/'protocol.json'):raise ValueError('matching concrete prelaunch review required')
    rootp=read(path(protocol['root_protocol_file']));contract=read(CONTRACT);terminal=read(ROOTRUN/'results.json')
    if terminal.get('protocol_sha256')!=protocol['root_protocol_sha256']:raise ValueError('root terminal protocol mismatch')
    locks={rel(ROOTRUN/'results.json'):sha(ROOTRUN/'results.json'),review_file:sha(path(review_file))}
    for row in rootp['records']:locks[row['output_folder']+'/result.json']=sha(path(row['output_folder'])/'result.json')
    write(OUT/'composition_inputs.json',{'created_utc':stamp(),'observer_protocol_sha256':sha(OUT/'protocol.json'),'root_result_hashes':locks,'scope':'Producer root results; independent review remains required.'})
    predictions={};records=[]
    for i,row in enumerate(rootp['records']):
        output=OUT/'records'/f'{i:02d}_{row["record_id"].split("_")[-1]}.json'
        root=read(path(row['output_folder'])/'result.json')
        try:composed=compose_record(protocol,rootp,contract,row,root)
        except (ValueError,KeyError,AssertionError) as error:composed={'status':'UNRESOLVED','record_id':row['record_id'],'reason':repr(error)}
        write(output,composed);records.append({'record_id':row['record_id'],'status':composed['status'],'file':rel(output),'sha256':sha(output)})
        if 'prediction_intervals_exact' in composed:predictions[row['record_id']]=composed['prediction_intervals_exact']
        print(json.dumps({'finished':i+1,'record_id':row['record_id'],'status':composed['status']}),flush=True)
    assert_hashes(locks);assert_hashes({**protocol['source_hashes'],**protocol['input_hashes']})
    summary=aggregate(load_observations(),predictions)
    answer={'status':'PRODUCER_COMPLETE' if summary['all68_uniform_static_tuple_enclosed'] else 'PRODUCER_PARTIAL','protocol_sha256':sha(OUT/'protocol.json'),
        'composition_inputs_sha256':sha(OUT/'composition_inputs.json'),'common_contract_sha256':sha(CONTRACT),'parameter_cell_sha256':protocol['parameter_cell_sha256'],
        'records':records,'aggregate':summary,'quantifier':'forall theta in C, independently for each covered record there exists its selected stationary root with these static outputs',
        'branch_law_semantics':'independent_stationary_root_per_record','source_time_and_roi_and_pairing_conditional':True,
        'PDE_or_all_branch_or_full_identified_set_claim':False,'independent_review_required':True}
    write(OUT/'results.json',answer);return answer

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','freeze','run']);parser.add_argument('--review');args=parser.parse_args()
    result=tests() if args.stage=='test' else freeze() if args.stage=='freeze' else run(args.review) if args.review else parser.error('--review required')
    print(json.dumps({'stage':args.stage,'status':result.get('status','FROZEN')}))
