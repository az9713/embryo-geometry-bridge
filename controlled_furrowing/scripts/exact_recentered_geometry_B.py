"""Preparation-first common-domain geometry for the translated B candidate.
No roots are assumed in these geometric boxes. Preparation DOES NOT call the
actual whole-body certificate. Run is separately gated by an exact review.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,datetime,hashlib,json
import numpy as np
import cohort_static_observers as cohort
import root_body_geometry as body
from nodal_observer_error_transfer import curve_bound
from quadratic_endpoint_observer import _state_controls,Curve,positive_thickness
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry'
CANDIDATE=ROOT/'runs/exact_integration_bridge_v1/retube_B_preparation'
REVIEW='evidence/exact_retube_B_preparation_independent_review.json'
REVIEW_SHA='71a96ef4acd3525eccf252058fd245889b6449f3dea47b7e99bca2c282db4981'
CONTRACT='runs/cohort_static_observers_v1/model_observation_contract.json'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def q(v):return cohort.q(v)
def canonical_path(p):
    value=str(p).replace(chr(92),'/');target=(ROOT/value).resolve();target.relative_to(ROOT.resolve());return target.relative_to(ROOT.resolve()).as_posix()
def bind(bindings,name,digest):
    name=canonical_path(name)
    if name in bindings and bindings[name]!=digest:raise ValueError('conflicting immutable binding: '+name)
    if sha(ROOT/name)!=digest:raise ValueError('changed frozen input: '+name)
    bindings[name]=digest

def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.write_text(json.dumps(cohort.serial(value),indent=2,allow_nan=False)+'\n',encoding='utf8')
def pass_checks(checks):return isinstance(checks,dict)and bool(checks)and all(v is True for v in checks.values())
def validate_prelaunch(review,source,protocol,preflight):
    if not isinstance(review,dict)or review.get('status')!='PASS'or not pass_checks(review.get('checks'))or review.get('source_sha256')!=source or review.get('protocol_sha256')!=protocol or review.get('preflight_sha256')!=preflight:raise ValueError('nonempty literal-True exact-bound independent prelaunch review required')
def assert_bindings(protocol):
    for name,h in protocol['source_and_input_hashes'].items():
        if sha(ROOT/name)!=h:raise ValueError('frozen source/input changed: '+name)
def contains(outer,inner):return outer[0]<=inner[0]<=inner[1]<=outer[1]
def interval_add(a,b):return a[0]+b[0],a[1]+b[1]
def fingerprint(array):return {'dtype':np.asarray(array).dtype.str,'shape':list(np.asarray(array).shape),'C_bytes_sha256':hashlib.sha256(np.asarray(array).tobytes(order='C')).hexdigest()}

def dependencies():
    found={};pending=[Path(__file__)]
    while pending:
        file=pending.pop();name=file.relative_to(ROOT).as_posix()
        if name in found:continue
        found[name]=sha(file)
        for node in ast.walk(ast.parse(file.read_text(encoding='utf8'))):
            names=([node.module]if isinstance(node,ast.ImportFrom)else[a.name for a in node.names]if isinstance(node,ast.Import)else[])
            for module in names:
                if module:
                    child=ROOT/'scripts'/(module.split('.')[0]+'.py')
                    if child.exists():pending.append(child)
    return found

def tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    # Exact scalar boxes, including signed/asymmetric predictor hulls.
    for i,(center,hull,d,r)in enumerate([(Q(2),(-Q(1,3),Q(2,5)),Q(7,10),Q(1,9)),(Q(-1),(-Q(3,5),Q(1,7)),Q(-2,3),Q(1,10)),(Q(0),(Q(0),Q(0)),Q(0),Q(1,11))]):
        common=r+abs(d);box=(center+hull[0]-common,center+hull[1]+common)
        candidate=(center-d+hull[0]-r,center-d+hull[1]+r);path=(center+hull[0]+min(Q(0),-d),center+hull[1]+max(Q(0),-d));swept=(path[0]-r,path[1]+r)
        ck('toy_candidate_cube_containment_'+str(i),contains(box,candidate));ck('toy_translation_path_containment_'+str(i),contains(box,path));ck('toy_swept_cube_containment_'+str(i),contains(box,swept))
    ck('old_box_escape_not_geometry_failure',not contains((-Q(1),Q(1)),(Q(0),Q(2)))and contains((-Q(2),Q(2)),(Q(0),Q(2))))
    valid={'status':'PASS','checks':{'one':True},'source_sha256':'s','protocol_sha256':'p','preflight_sha256':'f'};validate_prelaunch(valid,'s','p','f');ck('valid_review_accepted',True)
    for name,c in [('empty',{}),('list',[True]),('false',{'x':False}),('truthy_integer',{'x':1}),('truthy_string',{'x':'yes'}),('missing',None)]:
        try:validate_prelaunch({**valid,'checks':c},'s','p','f');ck('bad_review_rejected_'+name,False)
        except ValueError:checks['bad_review_rejected_'+name]=True
    try:validate_prelaunch(valid,'bad','p','f');ck('review_source_mismatch_rejected',False)
    except ValueError:checks['review_source_mismatch_rejected']=True
    # The body kernel call must exist ONLY inside the gated run function.
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'));body_call_functions=[]
    for node in tree.body:
        if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)):
            if any(isinstance(call,ast.Call)and isinstance(call.func,ast.Attribute)and isinstance(call.func.value,ast.Name)and call.func.value.id=='body'and call.func.attr=='certify_box'for call in ast.walk(node)):body_call_functions.append(node.name)
    ck('actual_body_certificate_call_only_in_gated_run',body_call_functions==['run'])
    return checks

def prepare():
    if(OUT/'protocol.json').exists():raise ValueError('write-once preparation exists')
    checks=tests();bindings={};bind(bindings,REVIEW,REVIEW_SHA);review=read(ROOT/REVIEW)
    if review['status']!='PASS_CONDITIONAL_RETUBE_PREPARATION_ONLY'or not pass_checks(review['checks'])or len(review['checks'])!=23:raise ValueError('current independent retube preparation review required')
    candidate_receipt=read(CANDIDATE/'preparation_receipt.json');candidate_protocol=read(CANDIDATE/'protocol.json');row=candidate_protocol['selected_record']
    for name,key in[('prepare.py','candidate_source_sha256'),('protocol.json','candidate_protocol_sha256'),('preparation_receipt.json','candidate_receipt_sha256')]:
        if sha(CANDIDATE/name)!=review[key]:raise ValueError('candidate/review binding differs: '+name)
    for group in(review['input_hashes'],candidate_protocol['source_and_input_hashes']):
        for name,h in group.items():bind(bindings,name,h)
    for name,h in candidate_receipt['artifact_hashes'].items():bind(bindings,(CANDIDATE/name).relative_to(ROOT),h)
    for name,h in dependencies().items():bind(bindings,name,h)
    bind(bindings,(CANDIDATE/'preparation_receipt.json').relative_to(ROOT),sha(CANDIDATE/'preparation_receipt.json'));bind(bindings,CONTRACT,sha(ROOT/CONTRACT))
    state=read(ROOT/candidate_protocol['selected_record']['proof_center_state_file']);packet=read(ROOT/row['point_cache_file']);arrays=cohort.decode_cache(packet)
    if sha(ROOT/packet['input_state'])!=packet['input_state_sha256']or packet['input_state']!=row['proof_center_state_file']:raise ValueError('original B seed/cache association differs')
    with np.load(ROOT/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    with np.load(CANDIDATE/'candidate_vectors.npz',allow_pickle=False)as z:u0=z['u0'].copy();u1=z['u1'].copy()
    if not np.array_equal(u0,arrays['u'])or u1.shape!=u0.shape or not np.isfinite(u1).all():raise ValueError('original/candidate binary centers differ')
    ledger=read(CANDIDATE/'rounding_ledger.json');applied=[q(a)-q(b)for a,b in zip(u0,u1)]
    if applied!=list(map(Q,ledger['applied_decrement_exact'])):raise ValueError('exact applied shift changed')
    rold=Q(candidate_protocol['old_radius_exact']);dnorm=max(map(abs,applied));common=rold+dnorm
    if common!=Q(candidate_receipt['proposed_common_radius_exact']):raise ValueError('common radius must use actual binary-applied shift')
    C=row['parameter_cell'];base=row['base_parameters_exact'];contract=read(ROOT/CONTRACT)
    if cohort.cell_hash(C)!=row['parameter_cell_sha256']or cohort.canonical(C)!=cohort.canonical(contract['parameter_cell']):raise ValueError('same exact parameter cell required')
    invariant=read(CANDIDATE/'unchanged_target.json')
    if invariant['unchanged_R']!=fingerprint(arrays['preconditioner'])or invariant['unchanged_S']!=fingerprint(S)or invariant['full_parameters']!=packet['parameters']:raise ValueError('R/S/physical model changed')
    # PURE factory; radius is a NEW explicitly declared geometry-domain radius.
    centers,radii,gaps,hulls,free=cohort.nodal_box(state,arrays,S,C,base,common)
    ids,endpoints=cohort.fixed_endpoints(state,radii,gaps,free)
    index={int(k):i for i,k in enumerate(free)};rows=[]
    for node in range(len(centers)):
        for coord in(0,1):
            k=2*node+coord;i=index.get(k);outer=(centers[node][coord]-radii[node][coord],centers[node][coord]+radii[node][coord]);X=q(arrays['reference'][node,coord]);isfree=i is not None
            if isfree:
                oldcenter=X+q(u0[i]);newcenter=X+q(u1[i]);lo,hi=hulls[i];d=applied[i]
                oldcube=(oldcenter+lo-rold,oldcenter+hi+rold);candidate=(newcenter+lo-rold,newcenter+hi+rold);path=(oldcenter+lo+min(Q(0),-d),oldcenter+hi+max(Q(0),-d));swept=(path[0]-rold,path[1]+rold)
            else:oldcube=candidate=path=swept=(X,X)
            included={name:contains(outer,value)for name,value in[('old_cube',oldcube),('candidate_cube',candidate),('translation_paths',path),('swept_translation_cubes',swept)]}
            if not all(included.values()):raise ValueError('new common factory box misses required domain at DOF'+str(k))
            rows.append({'global_dof':k,'node':node,'component':'R'if coord==0 else'Z','free':isfree,'common_center_exact':centers[node][coord],'common_radius_exact':radii[node][coord],'common_box_exact':outer,'old_cube_exact':oldcube,'candidate_cube_exact':candidate,'translation_path_hull_exact':path,'swept_translation_cube_hull_exact':swept,'containment':included})
    checks['all650_new_common_boxes_frozen']=len(rows)==650;checks['all2600_exact_domain_inclusions']=all(all(x['containment'].values())for x in rows);checks['all611_applied_shifts_recomputed']=len(applied)==611
    checks['same_C_R_S_physical_model']=True;checks['axis_and_boundary_endpoint_radial_coordinates_fixed']=all(radii[node][0]==0 and centers[node][0]==0 for node,rowref in enumerate(arrays['reference'])if q(rowref[0])==0)and bool(endpoints)
    L=Q(contract['physical_length_um_exact']);roi=list(map(Q,contract['roi_um_exact']));controls,height=_state_controls(state,L)
    boundary={'physical_length_um_exact':L,'same_current_radius_ROI_um_exact':roi,'original_baseline_state_file':packet['input_state'],'original_baseline_state_sha256':packet['input_state_sha256'],'edge_node_ids':ids,'fixed_radial_endpoints':endpoints,'baseline_controls_um':controls,'reference_height_um_exact':height,'boundary_specific_nodal_error_um_exact':{name:L*max(radii[node][c]for edge in ids[name]for node in edge for c in(0,1))for name in('top','bottom')},'candidate_is_NOT_marked_converged':True,'future_uniform_surface_G_upper_not_established':True}
    box={'domain_id':'B_exact_integral_retube_common_geometry_v1','domain_definition':'p0(theta)+B_infinity(r_common), theta in SAME C; contains old normal cubes, translated-candidate same-radius cubes, translation paths and swept cubes.','root_membership_assumed':False,'parameter_cell_sha256':row['parameter_cell_sha256'],'reference_nodes_exact':[[q(v)for v in rr]for rr in arrays['reference']],'triangles':arrays['triangles'].tolist(),'current_centers_exact':centers,'coordinate_radii_exact':radii,'predictor_hulls_exact':hulls,'rows':rows}
    protocol={'version':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_and_input_hashes':bindings,'record_id':row['record_id'],'geometry_domain_id':box['domain_id'],'parameter_cell':C,'parameter_cell_sha256':row['parameter_cell_sha256'],'base_parameters_exact':base,'old_radius_exact':rold,'actual_applied_decrement_norm_exact':dnorm,'new_common_enclosing_radius_exact':common,'domain_centering':'ORIGINAL B exact predictor p0(theta), represented around original serialized seed with exact gaps; candidate p1=p0-d is a contained proposal, not a root assertion.','physical_model_unchanged':True,'exact_geometry_target':'Affine material corners and exact P2 current coordinates, including all fixed binary offsets; same canonical field as reviewed exact-FE target.','body_kernel':'root_body_geometry.certify_box; independently validates supplied nodal boxes without requiring any root.','body_acceptance':'Fresh positive-J/all-element embedding+axis+strong-monotonicity, plus explicit Jlower>1e-6. Derive rho ONLY from returned fresh global J interval.','boundary_acceptance':'Fresh fixed endpoints, curve_bound derivative margins/common ROI, and baseline positive thickness minus apical/basal Eulerian errors. No unproved-root endpoint values.','surface_domain_corollary':'If uniform meridional symmetric-gradient lower eigenvalue m>0 and axis radius fixed, r/R>=m and dr/dR>=m; hence Js>=m² and G>=m⁴ on canonical boundary traces. Does not supply upperG/rhoG or cache-strain-line domain.','separate_missing_checks':['Uniform surface upperG/expansion bounds need a NEW interval-coefficient trace wrapper; point-only positive_G_bounds cannot be promoted.','Uniform HQ/HI/cache convex-strain-line bounds, parameter-discrepancy variation and root selfmap/contraction remain separate.','No PDE, biological observation, source temporal or trajectory validation follows.'],'prelaunch_review_required':{'file':'prelaunch_independent_review.json','status':'PASS','nonempty_literal_true_checks':True,'hashes':['source_sha256','protocol_sha256','preflight_sha256']},'body_certificate_called_during_prepare':False,'actual_root_certified':False}
    write(OUT/'coordinate_boxes.json',box);write(OUT/'boundary_preparation.json',boundary);protocol['prepared_box_sha256']=sha(OUT/'coordinate_boxes.json');protocol['boundary_preparation_sha256']=sha(OUT/'boundary_preparation.json');write(OUT/'protocol.json',protocol)
    for name,h in bindings.items():
        if sha(ROOT/name)!=h:raise ValueError('source/input changed during preparation '+name)
    if not all(checks.values()):raise AssertionError(checks)
    write(OUT/'synthetic_and_containment_tests.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'scope':'Synthetic scalar/formula/review guards plus exact new-domain650coordinate containment. No actual full-body or boundary graph certificate run.'})
    preflight={'status':'PASS_GEOMETRY_PREPARATION_ONLY','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'coordinate_boxes_sha256':sha(OUT/'coordinate_boxes.json'),'boundary_preparation_sha256':sha(OUT/'boundary_preparation.json'),'tests_sha256':sha(OUT/'synthetic_and_containment_tests.json'),'all_frozen_bindings_current':True,'actual_whole_body_certify_box_called':False,'actual_boundary_curve_bound_called':False,'new_root_claim':False};write(OUT/'preflight.json',preflight)
    print(json.dumps({'status':preflight['status'],'checks':len(checks),'coordinates':len(rows),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)

def run():
    p=read(OUT/'protocol.json');validate_prelaunch(read(OUT/'prelaunch_independent_review.json'),sha(__file__),sha(OUT/'protocol.json'),sha(OUT/'preflight.json'));assert_bindings(p)
    if sha(OUT/'coordinate_boxes.json')!=p['prepared_box_sha256']or sha(OUT/'boundary_preparation.json')!=p['boundary_preparation_sha256']:raise ValueError('prepared geometry inputs changed')
    box=read(OUT/'coordinate_boxes.json');boundary=read(OUT/'boundary_preparation.json');folder=OUT/'certification_run';write(folder/'started.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':sha(OUT/'protocol.json'),'review_sha256':sha(OUT/'prelaunch_independent_review.json')})
    geometry=body.certify_box(box['reference_nodes_exact'],box['current_centers_exact'],box['triangles'],box['coordinate_radii_exact']);write(folder/'whole_body_geometry.json',geometry)
    success=geometry['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING';J=geometry.get('J_bound');analytic=success and J is not None and Q(J['lower_exact'])>Q(1,1000000);rho=max(abs(Q(J['lower_exact'])-1),abs(Q(J['upper_exact'])-1))if J else None
    transfer={name:curve_bound(boundary['baseline_controls_um'][name],Q(boundary['boundary_specific_nodal_error_um_exact'][name]),fixed_radial_endpoints=True,roi_um=boundary['same_current_radius_ROI_um_exact'])for name in('top','bottom')};graphs=all(v['status']=='CONDITIONAL_EXACT_ERROR_TRANSFER'for v in transfer.values());separation=None;baseline_thickness=None
    if graphs:
        baseline_thickness=positive_thickness(Curve(boundary['baseline_controls_um']['top']),Curve(boundary['baseline_controls_um']['bottom']),list(map(Q,boundary['same_current_radius_ROI_um_exact'])),Q(1,10000000),8192);separation=baseline_thickness-Q(transfer['top']['Eulerian_vertical_error_um_exact'])-Q(transfer['bottom']['Eulerian_vertical_error_um_exact'])
    m=Q(geometry['strong_monotonicity_lower_exact'])if success else None
    output={'status':'CERTIFIED_COMMON_NODAL_GEOMETRY_UNREVIEWED'if analytic and graphs and separation is not None and separation>0 else'UNRESOLVED_COMMON_GEOMETRY_COMPONENT','protocol_sha256':sha(OUT/'protocol.json'),'body_geometry_file':'whole_body_geometry.json','body_geometry_sha256':sha(folder/'whole_body_geometry.json'),'new_J_interval':J,'fresh_rho_about1_exact':rho,'analytic_no_penalty_J_gt1e6_inverse':analytic,'body_embedding_certified':success,'boundary_graph_transfer':transfer,'baseline_thickness_lower_um_exact':baseline_thickness,'fresh_uniform_Eulerian_thickness_lower_um_exact':separation,'uniform_surface_positive_lower_corollary':{'meridional_m_exact':m,'hoop_lower_exact':m,'canonical_radial_derivative_lower_exact':m,'Js_lower_exact':m*m,'G_lower_exact':m**4}if m is not None else None,'uniform_surface_G_upper_or_binomial_rho_established':False,'root_membership_or_new_integral_root_certified':False,'scope':'Fresh geometry of every canonical field in a newly declared containing box ONLY. No exact-integral/cache Hessian bound, parameter discrepancy variation, root, endpoint-data fit, PDE or biological trajectory claim.'}
    assert_bindings(p);write(folder/'results.json',output);print(json.dumps({'status':output['status'],'root_certified':False}),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run','test']);args=ap.parse_args()
    if args.command=='prepare':prepare()
    elif args.command=='run':run()
    else:print(json.dumps({'status':'PASS','checks':tests()}))
