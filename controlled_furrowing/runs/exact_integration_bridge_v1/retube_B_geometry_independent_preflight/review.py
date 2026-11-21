"""Independent exact650-box/2600-containment audit; no geometry evaluation."""
from pathlib import Path
from fractions import Fraction as F
import ast,base64,hashlib,io,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];OUT=CF/'runs/exact_integration_bridge_v1/retube_B_geometry';CAND=CF/'runs/exact_integration_bridge_v1/retube_B_preparation'
SOURCE=CF/'scripts/exact_recentered_geometry_B.py';EXPECTED='73c7765c53781e2f0fc5337dc252c61481175c84e21fbd6faba9821a2ea4ed80';PROTOCOL='c220487fc55b1f057ef6bfe6bc1a8b01216ee8556c6d7db7506d95c0ce29bb6e';PREFLIGHT='55b33469c6c8fc05bc9cfde65cd22b9d2e1d4fc63a1fce068b94e236b334ce4b'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def q(x):return F.from_float(float(x))if isinstance(x,(float,np.floating))else F(x)
def inside(a,b):return a[0]<=b[0]<=b[1]<=a[1]
def main():
    if sha(SOURCE)!=EXPECTED or sha(OUT/'protocol.json')!=PROTOCOL or sha(OUT/'preflight.json')!=PREFLIGHT:raise ValueError('frozen geometry preparation changed')
    p=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');box=read(OUT/'coordinate_boxes.json');boundary=read(OUT/'boundary_preparation.json');cp=read(CAND/'protocol.json');row=cp['selected_record'];checks={};inputs={}
    for k,d in p['source_and_input_hashes'].items():
        if sha(CF/k)!=d:raise ValueError('frozen source/input differs '+k)
        inputs[k]=d
    checks['box_and_boundary_packets_hash_bound']=sha(OUT/'coordinate_boxes.json')==p['prepared_box_sha256']==pre['coordinate_boxes_sha256']and sha(OUT/'boundary_preparation.json')==p['boundary_preparation_sha256']==pre['boundary_preparation_sha256']
    checks['preparation_checks_literal_true_and_no_geometry_run']=pre['status']=='PASS_GEOMETRY_PREPARATION_ONLY'and bool(pre['checks'])and all(v is True for v in pre['checks'].values())and pre['actual_whole_body_certify_box_called']is False and pre['actual_boundary_curve_bound_called']is False and pre['new_root_claim']is False
    cache=read(CF/row['point_cache_file']);raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('cache corrupt')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:a={k:z[k].copy()for k in z.files}
    with np.load(CF/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    with np.load(CAND/'candidate_vectors.npz',allow_pickle=False)as z:u0=z['u0'].copy();u1=z['u1'].copy()
    state=read(CF/cache['input_state']);physical=read(CAND/'unchanged_target.json');contract=read(CF/'runs/cohort_static_observers_v1/model_observation_contract.json');rold=F(cp['old_radius_exact']);applied=[q(x)-q(y)for x,y in zip(u0,u1)];dnorm=max(map(abs,applied));common=rold+dnorm
    checks['same_recordB_seed_u0_free_map']=row['record_id']==p['record_id']=='optoGEF2_small_col_B'and cache['input_state']==row['proof_center_state_file']and sha(CF/cache['input_state'])==cache['input_state_sha256']and u0.tobytes()==a['u'].tobytes()and a['free'].shape==(611,)
    checks['exact_common_radius_uses_applied_not_requested_shift']=common==F(p['new_common_enclosing_radius_exact'])and dnorm==F(p['actual_applied_decrement_norm_exact'])and rold==F(p['old_radius_exact'])
    C=row['parameter_cell'];base=row['base_parameters_exact'];order=['apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain'];deltas=[(F(C[k][0])-F(base[k]),F(C[k][1])-F(base[k]))for k in order];hulls=[]
    for ss in S:
        terms=[(q(s)*lo,q(s)*hi)for s,(lo,hi)in zip(ss,deltas)];hulls.append((sum(min(x)for x in terms),sum(max(x)for x in terms)))
    checks['same_C_base_and_S_hulls']=p['parameter_cell']==C==physical['parameter_cell']and p['base_parameters_exact']==base and hulls==[tuple(map(F,x))for x in box['predictor_hulls_exact']]
    checks['reference_mesh_unchanged']=np.array_equal(np.asarray([[q(v)for v in r]for r in a['reference']],dtype=object),np.asarray([[F(v)for v in r]for r in box['reference_nodes_exact']],dtype=object))and np.array_equal(a['triangles'],np.asarray(box['triangles']))
    lookup={int(k):i for i,k in enumerate(a['free'])};rows={x['global_dof']:x for x in box['rows']};centers=[];radii=[];count=0;axisok=True;boundariesok=True
    for node in range(325):
        cc=[];rr=[]
        for component in range(2):
            k=2*node+component;ctr=q(state['current_nodes'][node][component]);ref=q(a['reference'][node,component]);exact=ref;rad=abs(ref-ctr);i=lookup.get(k)
            if i is not None:exact+=q(u0[i]);lo,hi=hulls[i];rad=abs(exact-ctr)+max(abs(lo),abs(hi))+common;oldcube=(exact+lo-rold,exact+hi+rold);candidate=(oldcube[0]-applied[i],oldcube[1]-applied[i]);path=(exact+lo+min(F(0),-applied[i]),exact+hi+max(F(0),-applied[i]));swept=(path[0]-rold,path[1]+rold)
            else:oldcube=candidate=path=swept=(ref,ref)
            outer=(ctr-rad,ctr+rad);families={'old_cube_exact':oldcube,'candidate_cube_exact':candidate,'translation_path_hull_exact':path,'swept_translation_cube_hull_exact':swept};entry=rows[k]
            valid=entry['node']==node and bool(entry['free'])==(i is not None)and tuple(map(F,entry['common_box_exact']))==outer and all(tuple(map(F,entry[name]))==v and inside(outer,v)for name,v in families.items())and F(box['current_centers_exact'][node][component])==ctr and F(box['coordinate_radii_exact'][node][component])==rad
            checks[f'coordinate{k:03d}_four_exact_family_containments_and_factory_values']=valid;count+=4;cc.append(ctr);rr.append(rad)
            if component==0 and ref==0:axisok &= rad==0 and ctr==0
        centers.append(cc);radii.append(rr)
    checks['exact650_unique_rows_and2600_containments']=len(rows)==len(box['rows'])==650 and count==2600
    checks['axis_radial_coordinates_exactly_fixed']=axisok
    L=F(boundary['physical_length_um_exact']);checks['same_length_and_current_radius_ROI']=L==F(contract['physical_length_um_exact'])and list(map(F,boundary['same_current_radius_ROI_um_exact']))==list(map(F,contract['roi_um_exact']))
    for name,face in[('bottom',0),('top',1)]:
        edges=a[f'surface_{face}_edges'].tolist();controls=[[[q(v)*L for v in state['current_nodes'][node]]for node in edge]for edge in edges];eps=L*max(radii[node][c]for edge in edges for node in edge for c in range(2));fixed=[edges[0][0],edges[-1][-1]]
        checks[name+'_boundary_controls_nodal_error_and_fixed_radial_endpoints']=edges==boundary['edge_node_ids'][name]and controls==[[list(map(F,v))for v in e]for e in boundary['baseline_controls_um'][name]]and eps==F(boundary['boundary_specific_nodal_error_um_exact'][name])and fixed==boundary['fixed_radial_endpoints'][name]and all(2*node not in lookup and radii[node][0]==0 for node in fixed)
    # Source closure and the absence of actual certificate calls in preparation.
    tree=ast.parse(SOURCE.read_text());calls={}
    for n in tree.body:
        if isinstance(n,ast.FunctionDef):calls[n.name]=[(x.func.value.id,x.func.attr)for x in ast.walk(n)if isinstance(x,ast.Call)and isinstance(x.func,ast.Attribute)and isinstance(x.func.value,ast.Name)]
    checks['whole_body_certificate_only_inside_gated_run']=[name for name,x in calls.items()if('body','certify_box')in x]==['run']
    checks['prepare_uses_original_pure_nodal_box_factory']=('cohort','nodal_box')in calls['prepare']and('body','certify_box')not in calls['prepare']
    # Verify local runtime-module source closure by an independent queue walk.
    todo=[SOURCE];visited=set()
    while todo:
        f=todo.pop();key=f.relative_to(CF).as_posix()
        if key in visited:continue
        visited.add(key)
        for n in ast.walk(ast.parse(f.read_text())):
            names=[n.module]if isinstance(n,ast.ImportFrom)else[x.name for x in n.names]if isinstance(n,ast.Import)else[]
            for module in names:
                if module:
                    child=CF/'scripts'/(module.split('.')[0]+'.py')
                    if child.exists():todo.append(child)
    checks['all_runtime_local_dependency_sources_locked']=all(k in inputs and sha(CF/k)==inputs[k]for k in visited)
    checks['future_run_requires_fresh_J_and_separate_boundary_checks']='J[\'lower_exact\']'in SOURCE.read_text()and 'rho=max('in SOURCE.read_text()and 'curve_bound('in SOURCE.read_text()and 'positive_thickness('in SOURCE.read_text()
    checks['no_surface_G_upper_or_root_assertion_in_preparation']=p['actual_root_certified']is False and boundary['future_uniform_surface_G_upper_not_established']is True and p['body_certificate_called_during_prepare']is False
    # Gate validator only: import definitions, never execute prepare/run.
    sys.path.insert(0,str(CF/'scripts'));import exact_recentered_geometry_B as G
    good={'status':'PASS','checks':{'one':True},'source_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT};G.validate_prelaunch(good,EXPECTED,PROTOCOL,PREFLIGHT)
    checks['strict_literal_review_gate']=False
    try:G.validate_prelaunch({**good,'checks':{'one':1}},EXPECTED,PROTOCOL,PREFLIGHT)
    except ValueError:checks['strict_literal_review_gate']=True
    checks['all_frozen_bindings_still_current']=all(sha(CF/k)==v for k,v in inputs.items())and sha(SOURCE)==EXPECTED
    out={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'source_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'preflight_sha256':PREFLIGHT,'review_source_sha256':sha(Path(__file__)),'input_hashes':inputs,'coordinate_boxes_checked':650,'individual_family_containments_checked':2600,'runtime_local_source_count':len(visited),'actual_body_certificate_or_curve_bound_called':False,'mathematical_review':['Common radius r+||d_applied|| contains oldcube,translatedcandidatecube,translationpathsandsweptcubes atsameC.','Serializedseedcentering includes exactX+u0gap, predictorhulls andnewcommonradius onlyonfreecoordinates.','Freshsym(F)>=mIplusfixedaxis impliesRcur/Rref>=m andradialboundaryderivative>=m, henceJs>=m² andG>=m⁴; noGupper/rhoG follows.','FutureJ/rho mustcome fromfreshreturnedbodybounds; boundarygraph/thickness transfer remainsseparate anddoesnotassume roots.'],'scope':'PRELAUNCHgeometrypreparation only. No actualwholebodycertificate, graphcertification, uniformcache/H/integration/root/PDEclaim.'}
    HERE.mkdir(parents=True,exist_ok=True);(HERE/'review.json').write_text(json.dumps(out,indent=2)+'\n');(OUT/'prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(OUT/'prelaunch_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
