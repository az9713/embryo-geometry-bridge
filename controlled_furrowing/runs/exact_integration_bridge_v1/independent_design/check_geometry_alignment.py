"""Check canonical-field inclusion in all17 existing whole-body boxes; no solve."""
from pathlib import Path
from fractions import Fraction as Q
import base64, hashlib, io, json
import numpy as np
H=Path(__file__).resolve().parent;C=H.parents[2]
def q(x):return Q.from_float(float(x)) if isinstance(x,(float,np.floating)) else Q(x)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    joined=json.loads((C/'evidence/shared_cell_cohort_inference.json').read_text()); bindings={}; rows=[]
    def load(name,digest=None):
        p=C/name; d=sha(p)
        if digest is not None and d!=digest:raise ValueError('hash differs')
        bindings[name]=d;return json.loads(p.read_text())
    for p in joined['proof_provenance']:
        root=load(p['root_certificate_file'],p['root_certificate_sha256']);obs=load(p['observer_certificate_file'],p['observer_certificate_sha256'])
        cache=load(root['point_cache_file']);blob=cache['frozen_binary_cache'];raw=base64.b64decode(blob['data'],validate=True)
        if hashlib.sha256(raw).hexdigest()!=blob['sha256']:raise ValueError('cache corrupt')
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
        statepath=root.get('original_seed_state_file',root.get('saved_state'))
        seed=load(statepath,root.get('original_seed_state_sha256',root.get('saved_state_sha256')))
        radius=q(root['uniform_cell_result']['selected_radius_exact']);hulls=[[q(x) for x in v]for v in obs['predictor_hulls_exact']]
        lookup={int(k):i for i,k in enumerate(a['free'])};count=0
        for node in range(len(a['reference'])):
            for d in range(2):
                k=2*node+d;exact=q(a['reference'][node,d]);extent=Q(0)
                if k in lookup:
                    i=lookup[k];exact+=q(a['u'][i]);extent=max(map(abs,hulls[i]))+radius
                need=abs(exact-q(seed['current_nodes'][node][d]))+extent
                if need>q(obs['coordinate_radii_normalized_exact'][node][d]):raise ValueError('canonical field not enclosed')
                count+=1
        body=obs['body_geometry'];jl=q(body['J_bound']['lower_exact']);ju=q(body['J_bound']['upper_exact'])
        if not body['all_element_J_positive'] or not body['axis_preserved_for_entire_box'] or not 0<jl<=ju:raise ValueError('geometry does not PASS')
        rho=max(abs(jl-1),abs(ju-1));rows.append({'record_id':p['record_id'],'canonical_coordinate_inclusions':count,'same_affine_reference_and_P2_current_interpretation':True,'J_lower_exact':str(jl),'J_upper_exact':str(ju),'rho_exact':str(rho)})
    rho=max(q(x['rho_exact']) for x in rows);upper=Q(21,500)
    if not rho<upper<1:raise ValueError('rho budget fails')
    remainders=[]
    for n in [8,10,12,16,20]:
        rem=upper**(n+1)/((n+1)*(1-upper));dr=upper**n/(1-upper);d2=n*upper**(n-1)/(1-upper)+upper**n/(1-upper)**2
        remainders.append({'N':n,'log_remainder_exact':str(rem),'first_J_derivative_remainder_exact':str(dr),'second_J_derivative_remainder_exact':str(d2),'display':[float(rem),float(dr),float(d2)]})
    sourcefiles=['scripts/root_body_geometry.py','scripts/cohort_static_observers.py','scripts/recentered_static_observers.py','scripts/mixed_geometry_certificate_v2.py']
    for f in sourcefiles:bindings[f]=sha(C/f)
    result={'status':'PASS_CANONICAL_P2_GEOMETRY_ALIGNMENT','checks':{'all17canonical_fields_enclosed_by_existing_geometry_boxes':len(rows)==17,'all11050coordinate_inclusions_checked':sum(x['canonical_coordinate_inclusions'] for x in rows)==11050,'all_root_observer_cache_hashes_verified':True,'rho_strictly_below_21_over_500':rho<upper,'no_new_geometry_or_Hessian_assembly_or_solve':True},'records':rows,'rho_exact':str(rho),'safe_rho_upper_exact':str(upper),'safe_scalar_log_series_remainders':remainders,'source_interpretation':'root_body_geometry.triangle_box constructs exact affine map from first3reference vertices and exactP2current from6nodal intervals; J denominator is affine reference radius elevated todegree4. Neither quadratic reference-midpoint interpolation nor cached shape_grad defines this geometry.','axis_interpretation':'Bounds hold for positive reference radius and extend as inequalities in axis-adapted coordinates; no unique physical-corner limit or global polynomial divisibility is assumed for axis-vertex-only triangles.','coefficient_preservation':'Stored binary nodal X retained; not silently replaced by affine midpoint values. Fixed coordinates remain fixed binary values.','remaining_scope':'Scalar log and J bounds do not bound complete FE gradient/Hessian integral defects, establish an exact-integration root, or supply PDE errors.','bindings':bindings,'script_sha256':sha(Path(__file__))}
    (H/'geometry_alignment.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(result['checks']),'rho_display':float(rho),'safe_upper':str(upper),'receipt_sha256':sha(H/'geometry_alignment.json'),'scalar_remainders_display':[{k:v for k,v in x.items()if k in ['N','display']}for x in remainders]},indent=2))
if __name__=='__main__':main()
