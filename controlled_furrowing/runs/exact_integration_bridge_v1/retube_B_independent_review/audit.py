"""Exact arithmetic review of a candidate center translation; no field calls."""
from pathlib import Path
from fractions import Fraction as F
import base64,hashlib,io,json
import numpy as np
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];PROD=CF/'runs/exact_integration_bridge_v1/retube_B_preparation'
EXPECTED_SOURCE='e118f5f8c6e9e919847bfd596b439eaf99fbe945a224da8771a8bffed47e8131'
EXPECTED_PROTOCOL='b71a90b45bae25a792f602d22915791601141dcf90000fd09afefa403a09f394'
EXPECTED_RECEIPT='5620374b95deaf5362c0e4818804226e217a9ffbb4d6522f3eebd3d83a974ef0'
ORDER=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def q(x):return F.from_float(float(x))if isinstance(x,(float,np.floating))else F(x)
def rounding_proved(exact,binary):
    y=np.float64(binary)
    if not np.isfinite(y):return False
    a,b=np.nextafter(y,np.float64(-np.inf)),np.nextafter(y,np.float64(np.inf))
    if not np.isfinite(a)or not np.isfinite(b):return False
    lo,hi=(q(a)+q(y))/2,(q(y)+q(b))/2
    even=(int(y.view(np.uint64))&1)==0
    return lo<=exact<=hi and((exact!=lo and exact!=hi)or even)
def fingerprint(a):return{'dtype':np.asarray(a).dtype.str,'shape':list(np.asarray(a).shape),'C_bytes_sha256':hashlib.sha256(np.asarray(a).tobytes(order='C')).hexdigest()}
def contains(a,b):return a[0]<=b[0]and b[1]<=a[1]
def main():
    if sha(PROD/'prepare.py')!=EXPECTED_SOURCE or sha(PROD/'protocol.json')!=EXPECTED_PROTOCOL or sha(PROD/'preparation_receipt.json')!=EXPECTED_RECEIPT:raise ValueError('frozen preparation changed')
    p=read(PROD/'protocol.json');r=read(PROD/'preparation_receipt.json');checks={};bindings={}
    for k,v in p['source_and_input_hashes'].items():
        if sha(CF/k)!=v:raise ValueError('changed frozen input '+k)
        bindings[k]=v
    for k,v in r['artifact_hashes'].items():
        if sha(PROD/k)!=v:raise ValueError('changed prepared artifact '+k)
    row=p['selected_record'];cache=read(CF/row['point_cache_file']);blob=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(blob).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('cache corrupt')
    with np.load(io.BytesIO(blob),allow_pickle=False)as z:old={k:z[k].copy()for k in z.files}
    with np.load(PROD/'candidate_vectors.npz',allow_pickle=False)as z:new={k:z[k].copy()for k in z.files}
    with np.load(CF/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    point=read(CF/p['exact_point_discrepancy_file']);observer=read(CF/'runs/cohort_static_observers_v1/composition_final/records/00_B.json');state=read(CF/cache['input_state']);candidate=read(PROD/'candidate_state.json');nodes=read(PROD/'candidate_exact_nodes.json');ledger=read(PROD/'rounding_ledger.json');domains=read(PROD/'domain_obligations.json');physical=read(PROD/'unchanged_target.json');lemma=read(PROD/'conditional_reuse_lemma.json')
    checks['preparation_receipt_true_checks_only']=r['status']=='PREPARED_TRANSLATED_CENTER_FRESH_PROOFS_REQUIRED'and bool(r['checks'])and all(v is True for v in r['checks'].values())
    checks['all611_original_binary_vectors_preserved']=new['u0'].shape==new['u1'].shape==new['requested_binary_decrement'].shape==(611,)and new['u0'].tobytes()==old['u'].tobytes()
    D=[(F(v['lower_exact']),F(v['upper_exact']))for v in point['preconditioned_force_defect']];applied=[q(a)-q(b)for a,b in zip(new['u0'],new['u1'])]
    checks['all611_midpoint_float_roundings_proved']=all(rounding_proved((lo+hi)/2,v)for(lo,hi),v in zip(D,new['requested_binary_decrement']))
    checks['all611_binary_subtractions_proved_exactly']=all(rounding_proved(q(a)-q(d),b)for a,d,b in zip(new['u0'],new['requested_binary_decrement'],new['u1']))
    checks['round_to_nearest_tie_even_fixture']=rounding_proved(F(1)+F(1,2**53),np.float64(1))and not rounding_proved(F(1)+F(1,2**53),np.nextafter(np.float64(1),np.float64(2)))
    checks['all611_applied_intercept_identities_exact']=all(q(b)==q(a)-d for a,b,d in zip(new['u0'],new['u1'],applied))and list(map(F,ledger['applied_decrement_exact']))==applied and list(map(F,ledger['new_predictor_intercept_delta_exact']))==[-d for d in applied]
    dnorm=max(map(abs,applied));mismatch=max(max(abs(lo-d),abs(hi-d))for(lo,hi),d in zip(D,applied));rad=F(p['old_radius_exact']);common=rad+dnorm
    checks['D_minus_actual_applied_shift_bound_recomputed']=mismatch==F(r['Dpoint_minus_applied_d_upper_exact'])==F(ledger['complete_point_D_minus_applied_d_norm_upper_exact'])==F(lemma['Dpoint_minus_applied_d_upper_exact'])
    checks['applied_shift_norm_and_common_radius_recomputed']=dnorm==F(r['applied_decrement_norm_exact'])==F(ledger['applied_decrement_norm_exact'])and common==F(r['proposed_common_radius_exact'])==F(domains['proposed_common_enclosing_radius_exact'])
    checks['C_and_base_and_dose_and_strength_unchanged']=physical['parameter_cell']==row['parameter_cell']and physical['parameter_cell_sha256']==row['parameter_cell_sha256']and physical['base_parameters_exact']==row['base_parameters_exact']and physical['source_dose_exact']==row['dose_exact']and F(physical['exact_integral_base_strength_exact'])==F(row['base_parameters_exact']['common_gain'])*F(row['dose_exact'])and physical['full_parameters']==cache['parameters']==candidate['parameters']and F(physical['energy_scale_exact'])==q(old['energy_scale'][0])
    checks['same_R_and_S_bit_fingerprints']=physical['unchanged_R']==fingerprint(old['preconditioner'])and physical['unchanged_S']==fingerprint(S)
    checks['every_reference_auxiliary_array_unchanged_fingerprint']=set(physical['unchanged_geometry_and_auxiliary_arrays'])==set(old)-{'u','preconditioner'}and all(physical['unchanged_geometry_and_auxiliary_arrays'][k]==fingerprint(v)for k,v in old.items()if k not in('u','preconditioner'))
    delta=[(F(row['parameter_cell'][k][0])-F(row['base_parameters_exact'][k]),F(row['parameter_cell'][k][1])-F(row['base_parameters_exact'][k]))for k in ORDER]
    hulls=[]
    for ss in S:
        ll=hh=F(0)
        for s,(a,b)in zip(ss,delta):aa,bb=q(s)*a,q(s)*b;ll+=min(aa,bb);hh+=max(aa,bb)
        hulls.append((ll,hh))
    checks['all611_predictor_hulls_unchanged']=hulls==[tuple(map(F,x))for x in observer['predictor_hulls_exact']]
    free={int(k):i for i,k in enumerate(old['free'])};escaped=[];contained=0;exact_nodes=[];fresh_hulls=[];fixed_ok=True;old_ok=True;node_round=True
    for node in range(325):
        pair=[];coordh=[]
        for c in range(2):
            key=2*node+c;ref=q(old['reference'][node,c]);oldcenter=ref;newcenter=ref
            if key in free:
                i=free[key];oldcenter+=q(new['u0'][i]);newcenter+=q(new['u1'][i]);lo,hi=hulls[i];oldtrue=(oldcenter+lo-rad,oldcenter+hi+rad);shifted=(oldtrue[0]-applied[i],oldtrue[1]-applied[i])
            else:oldtrue=(ref,ref);shifted=(ref,ref);fixed_ok &= newcenter==ref
            pair.append(str(newcenter));coordh.append(list(map(str,shifted)));ctr=q(state['current_nodes'][node][c]);rr=F(observer['coordinate_radii_normalized_exact'][node][c]);oldbox=(ctr-rr,ctr+rr);old_ok &= contains(oldbox,oldtrue)
            if contains(oldbox,shifted):contained+=1
            else:escaped.append((node,'R'if c==0 else'Z'))
            node_round &= rounding_proved(newcenter,candidate['current_nodes'][node][c])and F(ledger['new_exact_node_vs_serialization_gaps'][node][c])==abs(newcenter-q(candidate['current_nodes'][node][c]))
        exact_nodes.append(pair);fresh_hulls.append(coordh)
    checks['all650_authoritative_exact_nodes_match']=exact_nodes==nodes['nodes_exact']
    checks['all650_serialized_node_roundings_and_gaps_accounted']=node_round
    checks['old_body_box_contains_actual_old_tube']=old_ok
    checks['fixed_coordinate_values_unchanged']=fixed_ok
    checks['exact396_escaped_254_contained_coordinate_hulls']=len(escaped)==396 and contained==254 and domains['failed_coordinate_count']==396 and domains['contained_coordinate_count']==254 and set(escaped)=={(x['node'],x['component'])for x in domains['failed_coordinates']}
    checks['old_geometry_never_silently_reused']=domains['old_whole_coordinate_box_reused']is False and domains['old_body_boxes_enclose_new_predictor_plus_old_radius']is False and r['new_cube_contained_in_old_body_coordinate_boxes']is False
    checks['common_radius_encloses_entire_translation_and_new_cube']=common==rad+dnorm and all(abs(d)+rad<=common for d in applied)
    checks['old_scalar_domain_insufficient_for_translation']=dnorm>rad and domains['old_Q_scalar_proof_cube_contains_translation_path']is False and domains['old_Q_scalar_proof_cube_contains_new_same_radius_cube']is False
    checks['fresh_ZQ_b_V_eta_and_domain_NOT_claimed']=all(lemma[k]is None for k in['ZQ_fresh','b_fresh','V_parameter_predictor','eta_new_certified'])and lemma['common_domain_proved']is False and lemma['root_assertion']is False and r['new_integral_root_certified']is False
    checks['candidate_explicitly_not_solved']=candidate['converged']is False and candidate['no_gradient_Hessian_or_equilibrium_computed']is True and p['no_solve_or_residual_evaluation']is True
    checks['all_frozen_bindings_still_match']=all(sha(CF/k)==v for k,v in bindings.items())and sha(PROD/'prepare.py')==EXPECTED_SOURCE
    details={'applied_shift_norm_exact':str(dnorm),'D_minus_applied_shift_norm_upper_exact':str(mismatch),'old_radius_exact':str(rad),'common_comparison_radius_exact':str(common),'escaped_coordinate_hulls':escaped,'contained_count':contained,'new_hulls_exact':fresh_hulls}
    HERE.mkdir(parents=True,exist_ok=True);(HERE/'exact_arithmetic.json').write_text(json.dumps(details,indent=2)+'\n')
    out={'status':'PASS_CONDITIONAL_RETUBE_PREPARATION_ONLY'if all(checks.values())else'FAIL_PREPARATION_REVIEW','checks':checks,'candidate_source_sha256':EXPECTED_SOURCE,'candidate_protocol_sha256':EXPECTED_PROTOCOL,'candidate_receipt_sha256':EXPECTED_RECEIPT,'review_source_sha256':sha(Path(__file__)),'input_hashes':bindings,'exact_arithmetic_sha256':sha(HERE/'exact_arithmetic.json'),'arithmetic_summary':{k:v for k,v in details.items()if k.endswith('_exact')or k=='contained_count'},'conditional_lemma':'Forp0(theta), p1=p0-d: ||R FI(p1,theta)|| <= eta_old+||D_base-d||+(ZQ_fresh+b_fresh)||d||+V, onlyiffreshcommontranslationdomainandparameter/predictorvariationboundsareproved.','scope':'611binarycandidatecenter arithmetic and650coordinatehulls only. No fieldgradient/Hessian/integration/geometry/root/PDE assertion.'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');(CF/'evidence/exact_retube_B_preparation_independent_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'shift':float(dnorm),'mismatch':float(mismatch),'common_radius':float(common),'escaped':len(escaped),'receipt_sha256':sha(CF/'evidence/exact_retube_B_preparation_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
