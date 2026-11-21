"""Read-only exact arithmetic inventory. No model/Hessian assembly or solves."""
from pathlib import Path
from fractions import Fraction as Q
import base64, hashlib, io, json
import numpy as np

HERE=Path(__file__).resolve().parent
CF=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def q(x):return Q.from_float(float(x)) if isinstance(x,(float,np.floating)) else Q(x)
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')

def main():
    bindings={}
    def load(name,digest=None):
        p=CF/name; h=sha(p)
        if digest is not None and h!=digest:raise ValueError('hash mismatch: '+name)
        if name in bindings and bindings[name]!=h:raise ValueError('conflicting snapshot')
        bindings[name]=h;return read(p)
    join=load('evidence/shared_cell_cohort_inference.json')
    records=[]; reference_fingerprint=None; cache_geometry=None
    checks={'17_distinct_source_records':len(join['proof_provenance'])==17 and len({x['record_id'] for x in join['proof_provenance']})==17}
    for item in join['proof_provenance']:
        rid=item['record_id'];root=load(item['root_certificate_file'],item['root_certificate_sha256'])
        obs=load(item['observer_certificate_file'],item['observer_certificate_sha256'])
        review=load(item['root_review_file'],item['root_review_sha256'])
        cap=next(x['audit'] for x in review['record_audits'] if x['record_id']==rid)
        if 'bounds' in cap:b=cap['bounds']; eta=q(b['eta']['upper_exact'])
        else:
            reused=load(cap['reuse_review_file'],cap['reuse_review_sha256'])
            eta=q(reused['independent_bounds']['eta_uniform']['upper_exact'])
        s=load('runs/cohort_energy_stability_review_v1/records/'+rid+'.json')
        r=q(s['normal_root_radius_exact']);Z=q(s['independent_fulltube_Z_upper_exact']);N=q(s['exact_R_infinity_norm']);m=q(s['uniform_Hessian_lower_bound_exact'])
        slack=r-eta-Z*r
        if not (s['status']=='PASS' and all(v is True for v in s['checks'].values()) and 0<=Z<1 and slack>0 and m==(1-Z)/N):raise ValueError('existing proof mismatch')
        cache=load(root['point_cache_file']);raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
        if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('cache checksum')
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
        fp=hashlib.sha256(b''.join(np.ascontiguousarray(a[k]).tobytes() for k in ['reference','triangles','shape','shape_grad','qr','volume_weights','active'])).hexdigest()
        if reference_fingerprint is None:reference_fingerprint=fp;cache_geometry=(a,cache['parameters'],cache['integration_policy'])
        if fp!=reference_fingerprint:raise ValueError('different geometry/rule across records')
        J=obs['body_geometry']['J_bound'];low=q(J['lower_exact']);hi=q(J['upper_exact'])
        row={'record_id':rid,'dose_exact':root['dose_exact'],'radius_exact':str(r),'eta_exact':str(eta),'Z_exact':str(Z),'R_infinity_exact':str(N),'m_exact':str(m),'remaining_selfmap_slack_exact':str(slack),'sufficient_separate_unpreconditioned_force_budget_exact':str(slack/(4*N)),'sufficient_preconditioned_Hessian_defect_budget_exact':str(slack/(4*r)),'sufficient_separate_unpreconditioned_Hessian_budget_exact':str(min(slack/(4*r*N),m/2)),'wholebody_J_lower_exact':str(low),'wholebody_J_upper_exact':str(hi),'strong_monotonicity_lower_exact':obs['body_geometry']['strong_monotonicity_lower_exact']}
        row['display']={k:float(q(v)) for k,v in row.items() if k.endswith('_exact') and k!='dose_exact'}
        records.append(row)
    a,p,policy=cache_geometry; X=a['reference'];tri=a['triangles'];gap=Q(0);nonzero=0
    for t in tri:
        for i,j,k in [(0,1,3),(1,2,4),(2,0,5)]:
            for d in range(2):
                v=abs(q(X[t[k],d])-(q(X[t[i],d])+q(X[t[j],d]))/2);gap=max(gap,v);nonzero+=v!=0
    shape_partition=max(abs(sum(map(q,row),Q(0))-1) for row in a['shape'])
    grad_partition=max(abs(sum((q(v) for v in a['shape_grad'][t,iq,:,d]),Q(0))) for t in range(len(tri)) for iq in range(a['shape'].shape[0]) for d in range(2))
    counts={'axis_edge_triangles':0,'axis_vertex_only_triangles':0,'off_axis_triangles':0,'activation_interfaces_cutting_element_interior':0,'active_triangles':int(a['active'].sum())}
    pr=q(p['patch_radius']);zlo=q(p['activation_z_min']);zhi=q(p['activation_z_max'])
    for t,flag in zip(tri,a['active'].reshape(-1)):
        pts=[[q(v) for v in row] for row in X[t[:3]]];zeros=sum(v[0]==0 for v in pts)
        counts['axis_edge_triangles' if zeros==2 else 'axis_vertex_only_triangles' if zeros==1 else 'off_axis_triangles']+=1
        for axis,c in [(0,pr),(1,zlo),(1,zhi)]:
            if min(v[axis] for v in pts)<c<max(v[axis] for v in pts):counts['activation_interfaces_cutting_element_interior']+=1
        expected=max(v[0] for v in pts)<=pr and min(v[1] for v in pts)>=zlo and max(v[1] for v in pts)<=zhi
        if bool(flag)!=expected:raise ValueError('cache activation vs aligned rectangle differs')
    checks.update({'same_geometry_and_cache_rule_all17':True,'all17_independent_selfmap_slacks_positive':True,'all17_existing_m_reproduced_exactly':True,'activation_mask_matches_exact_closed_rectangle_elements':counts['activation_interfaces_cutting_element_interior']==0,'no_Hessian_assembly_or_solve':True})
    geometry={'cached_integration_policy':policy,'cache_shapes':{k:list(v.shape) for k,v in a.items()},'triangle_counts':counts,'binary_midpoint_max_gap_exact':str(gap),'binary_midpoint_nonzero_element_coordinate_occurrences':int(nonzero),'cached_shape_partition_unity_max_defect_exact':str(shape_partition),'cached_gradient_partition_zero_max_defect_exact':str(grad_partition),'fingerprint':reference_fingerprint,'qualification':'Cache partition/midpoint defects are exact rational arithmetic facts, not yet propagated force/Hessian bounds.'}
    extrema={}
    for key in ['remaining_selfmap_slack_exact','sufficient_separate_unpreconditioned_force_budget_exact','sufficient_preconditioned_Hessian_defect_budget_exact','m_exact','wholebody_J_lower_exact','strong_monotonicity_lower_exact']:
        row=min(records,key=lambda z:q(z[key]));extrema['minimum_'+key]={'record_id':row['record_id'],'exact':row[key],'display':float(q(row[key]))}
    for key in ['Z_exact','R_infinity_exact','wholebody_J_upper_exact']:
        row=max(records,key=lambda z:q(z[key]));extrema['maximum_'+key]={'record_id':row['record_id'],'exact':row[key],'display':float(q(row[key]))}
    out={'status':'PASS_READ_ONLY_SCALE_INVENTORY' if all(checks.values()) else 'FAIL','checks':checks,'bindings':bindings,'records':records,'geometry':geometry,'extrema':extrema,'source_sha256':sha(Path(__file__)),'exact_integration_bridge_certified':False,'no_error_bound_computed':True}
    save(HERE/'scales.json',out); print(json.dumps({'status':out['status'],'checks':len(checks),'extrema':extrema,'geometry_summary':{k:v for k,v in geometry.items()if k!='cache_shapes'}},indent=2))
if __name__=='__main__':main()
