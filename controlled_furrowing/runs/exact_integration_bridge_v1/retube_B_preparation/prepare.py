"""One explicit proof-center translation; no field evaluation or solve."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target as T
POINT='runs/exact_integration_bridge_v1/complete_point_B/results.json'
REVIEW='evidence/exact_complete_point_B_independent_review.json'
POINT_SHA='57bb157ba4ea40d790cb01c3775e03b6dbb7a811c8961f78d32167532af3b98d'
REVIEW_SHA='bc22a23e7af96d0d1a99836401112378600c14be61d5cba957e6f15a12d699c5'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,obj):
    if p.exists():raise ValueError('write-once file exists: '+str(p))
    p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')
def q(v):return T.rational(v)
def canonical(f):return str(f).replace('\\','/')
def bind(d,f,h):
    f=canonical(f);T.bind(d,f,h)
    if sha(ROOT/f)!=h:raise ValueError('changed immutable input: '+f)
def fingerprint(array):return {'dtype':np.asarray(array).dtype.str,'shape':list(np.asarray(array).shape),'C_bytes_sha256':hashlib.sha256(np.asarray(array).tobytes(order='C')).hexdigest()}
def contains(outer,inner):return outer[0]<=inner[0]<=inner[1]<=outer[1]
def conditional_eta(eta,mismatch,ZQ,b,dnorm,V,domain_verified):
    vals=list(map(F,(eta,mismatch,ZQ,b,dnorm,V)))
    if domain_verified is not True or any(v<0 for v in vals):raise ValueError('fresh common domains and nonnegative certified inputs required')
    eta,mismatch,ZQ,b,dnorm,V=vals;return eta+mismatch+(ZQ+b)*dnorm+V

def tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    # Exact scalar affine families test every term in the reuse estimate.
    for i,(a,beta,gamma,epsilon,d)in enumerate([(F(2),F(1,5),F(1,7),F(1,100),F(1,200)),(F(3),F(-1,9),F(2,11),F(-1,50),F(-1,170))]):
        theta0=F(1,3);delta=F(1,100);u0=theta0/a;D0=(epsilon+beta*u0+gamma*theta0)/a;V=abs(beta/a**2+gamma/a)*delta;bound=conditional_eta(0,abs(D0-d),0,abs(beta/a),abs(d),V,True)
        for sign in(-1,0,1):
            th=theta0+sign*delta;old=th/a;new=old-d;res=abs((a*new-th+epsilon+beta*new+gamma*th)/a);ck(f'affine_residual_lemma_{i}_{sign}',res<=bound)
    ck('parameter_variation_is_not_optional',F(1,4)*F(1,10)>conditional_eta(0,0,0,0,F(1,100),0,True))
    try:conditional_eta(0,0,0,0,1,0,False);ck('old_domain_not_silently_reused',False)
    except ValueError:checks['old_domain_not_silently_reused']=True
    ck('translated_cube_needs_fresh_domain',not contains((-F(1),F(1)),(F(-3),F(-1))))
    ck('common_radius_covers_translation',contains((-F(3),F(3)),(-F(3),F(-1))))
    d0=np.float64(2.**-60);u0=np.float64(1.);u1=np.subtract(u0,d0,dtype=np.float64);applied=q(u0)-q(u1)
    ck('rounded_requested_step_can_differ_from_applied',q(d0)!=applied and applied==0)
    ck('applied_intercept_identity_exact',q(u1)==q(u0)-applied)
    return checks

def main():
    if(OUT/'protocol.json').exists():raise ValueError('write-once preparation exists')
    checks=tests();bindings={}
    for f,h in[(POINT,POINT_SHA),(REVIEW,REVIEW_SHA)]:bind(bindings,f,h)
    point=read(ROOT/POINT);review=read(ROOT/REVIEW)
    if review['status']!='PASS_COMPLETE_POINT_DISCREPANCY_AND_NEGATIVE_SELF_MAP_CERTIFICATE'or not isinstance(review['checks'],dict)or len(review['checks'])!=195 or any(v is not True for v in review['checks'].values())or review['producer_point_results_sha256']!=POINT_SHA:raise ValueError('matching independently reviewed complete point discrepancy required')
    for f,h in review['source_and_input_hashes'].items():bind(bindings,f,h)
    prepfile='runs/exact_integration_bridge_v1/target_preparation/protocol.json';prep=read(ROOT/prepfile);bind(bindings,prepfile,sha(ROOT/prepfile));row=next(x for x in prep['records']if x['record_id']=='optoGEF2_small_col_B')
    for group in(prep['source_hashes'],prep['input_hashes']):
        for f,h in group.items():bind(bindings,f,h)
    bind(bindings,str(Path(__file__).relative_to(ROOT)),sha(__file__))
    root=read(ROOT/row['root_certificate_file']);packet=read(ROOT/row['point_cache_file']);target=T.Target(packet);a=target.a
    if root['artifact_hashes'][row['point_cache_file']]!=row['point_cache_sha256']or root['artifact_hashes'][row['predictor_file']]!=row['predictor_sha256']:raise ValueError('authoritative unchanged cache/predictor required')
    with np.load(ROOT/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    observerfile='runs/cohort_static_observers_v1/composition_final/records/00_B.json';observer=read(ROOT/observerfile);bind(bindings,observerfile,sha(ROOT/observerfile));oldstate=read(ROOT/packet['input_state']);bind(bindings,packet['input_state'],packet['input_state_sha256'])
    oldreview=read(ROOT/row['independent_root_capsule_file']);capsules=[v['audit']for v in oldreview['record_audits']if v['record_id']==row['record_id']]
    if len(capsules)!=1 or capsules[0]['status']!='PASS'or capsules[0]['uniform_quantifier_verified']is not True or capsules[0]['certificate_sha256']!=row['root_certificate_sha256']:raise ValueError('selected old uniform root audit required')
    capsule=capsules[0];rold=F(point['old_radius_exact']);eta_old=max(F(point['old_uniform_eta_upper_exact']),F(capsule['bounds']['eta']['upper_exact']))
    if rold!=F(row['selected_corrector_radius_exact'])or rold!=F(observer['full_corrector_radius_exact']):raise ValueError('same full old normal radius required')
    if any(F(point['base_parameters_exact'][k])!=F(row['base_parameters_exact'][k])for k in T.ORDER):raise ValueError('same physical base parameters required')
    if len(point['preconditioned_force_defect'])!=len(a['u'])or len(a['u'])!=611:raise ValueError('complete611-coordinate discrepancy required')
    protocol={'version':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_and_input_hashes':bindings,'selected_record':row,'exact_point_discrepancy_file':POINT,'exact_point_discrepancy_sha256':POINT_SHA,'independent_point_review_file':REVIEW,'independent_point_review_sha256':REVIEW_SHA,'recipe':'For each exact interval D_j=R(FI-FQ)(u0,theta0): d_requested=float64((Dlo+Dhi)/2); u1=float64(u0-d_requested); d_applied=Q(u0)-Q(u1) exactly. NEWpredictor=u1+S(theta-theta0), unchanged S/R/C.','no_solve_or_residual_evaluation':True,'new_root_assertion':False,'old_radius_exact':str(rold),'eta_old_uniform_exact':str(eta_old),'domain_policy':'Old whole-coordinate geometry boxes and old root proof domains are checked, never enlarged by implication. New uniform HI/HQ/parameter variation and geometry proofs are explicit future obligations.'}
    write(OUT/'protocol.json',protocol)
    D=[(F(v['lower_exact']),F(v['upper_exact']))for v in point['preconditioned_force_defect']]
    if any(l>h for l,h in D):raise ValueError('ordered exact discrepancy intervals required')
    requested=np.array([float((lo+hi)/2)for lo,hi in D],dtype=np.float64);u0=a['u'].copy();u1=np.subtract(u0,requested,dtype=np.float64)
    if not np.isfinite(requested).all()or not np.isfinite(u1).all():raise ValueError('finite binary center required')
    applied=[q(x)-q(y)for x,y in zip(u0,u1)];dnorm=max(map(abs,applied));mismatch=max(max(abs(lo-d),abs(hi-d))for(lo,hi),d in zip(D,applied));intercept=[-v for v in applied]
    if any(float(q(x)-q(d))!=float(y)for x,d,y in zip(u0,requested,u1)):raise ValueError('binary subtraction differs from independent exact-rational rounding')
    checks['all611_midpoints_and_binary_subtractions_accounted']=True;checks['exact_applied_intercept_translation']=all(q(y)==q(x)+t for x,y,t in zip(u0,u1,intercept))
    cell=row['parameter_cell'];base=row['base_parameters_exact'];delta=[(F(cell[k][0])-F(base[k]),F(cell[k][1])-F(base[k]))for k in T.ORDER];hulls=[]
    for ss in S:
        terms=[(q(v)*lo,q(v)*hi)for v,(lo,hi)in zip(ss,delta)];hulls.append((sum(min(t)for t in terms),sum(max(t)for t in terms)))
    if hulls!=[tuple(map(F,p))for p in observer['predictor_hulls_exact']]:raise ValueError('old exact predictor hulls changed')
    exact0=target.nodes();exact1=target.nodes(u1);oldcenters=[[q(v)for v in row]for row in oldstate['current_nodes']];oldradii=[[F(v)for v in row]for row in observer['coordinate_radii_normalized_exact']]
    oldbox=[[(c-r,c+r)for c,r in zip(cc,rr)]for cc,rr in zip(oldcenters,oldradii)];newbox=[[(v,v)for v in row]for row in exact1];oldactual=[[(v,v)for v in row]for row in exact0]
    for i,k in enumerate(a['free']):
        node,coord=divmod(int(k),2);lo,hi=hulls[i];newbox[node][coord]=(exact1[node][coord]+lo-rold,exact1[node][coord]+hi+rold);oldactual[node][coord]=(exact0[node][coord]+lo-rold,exact0[node][coord]+hi+rold)
    checks['old_body_box_encloses_its_actual_old_root_tube']=all(contains(oldbox[n][c],oldactual[n][c])for n in range(325)for c in(0,1))
    failed=[];contained=0
    for n in range(325):
        for c in(0,1):
            if contains(oldbox[n][c],newbox[n][c]):contained+=1
            else:failed.append({'node':n,'component':'R'if c==0 else'Z','global_dof':2*n+c,'old_box_exact':list(map(str,oldbox[n][c])),'new_cube_hull_exact':list(map(str,newbox[n][c])),'lower_excess_exact':str(max(F(0),oldbox[n][c][0]-newbox[n][c][0])),'upper_excess_exact':str(max(F(0),newbox[n][c][1]-oldbox[n][c][1]))})
    common=rold+dnorm;checks['common_radius_encloses_shifted_same_radius_cube']=all(abs(d)+rold<=common for d in applied);checks['fixed_coordinates_unchanged']=all(exact0[n][c]==exact1[n][c]for n in range(325)for c in(0,1)if 2*n+c not in target.index)
    arrfile=OUT/'candidate_vectors.npz';np.savez_compressed(arrfile,u0=u0,u1=u1,requested_binary_decrement=requested)
    coordinate_packet={'definition':'Exact rational current nodes Xbinary+P*u1; fixed coordinates unchanged. These are proof-center candidates, not a solved or experimentally observed field.','nodes_exact':[[str(v)for v in row]for row in exact1]};write(OUT/'candidate_exact_nodes.json',coordinate_packet)
    current=np.array(a['reference'],copy=True);current.ravel()[a['free']]+=u1
    state={'status':'NUMERICAL_PROOF_CENTER_CANDIDATE_NOT_SOLVED','converged':False,'original_seed_state_file':packet['input_state'],'original_seed_state_sha256':packet['input_state_sha256'],'parameters':packet['parameters'],'reference_nodes':a['reference'].tolist(),'triangles':a['triangles'].tolist(),'free_displacement':u1.tolist(),'current_nodes':current.tolist(),'current_nodes_are_rounded_serialization_only':True,'authoritative_exact_nodes_file':'candidate_exact_nodes.json','no_gradient_Hessian_or_equilibrium_computed':True};write(OUT/'candidate_state.json',state)
    roundledger={'requested_binary_decrement_exact':list(map(lambda v:str(q(v)),requested)),'applied_decrement_exact':list(map(str,applied)),'new_predictor_intercept_delta_exact':list(map(str,intercept)),'applied_decrement_norm_exact':str(dnorm),'complete_point_D_minus_applied_d_norm_upper_exact':str(mismatch),'requested_minus_applied_norm_exact':str(max(abs(q(r)-d)for r,d in zip(requested,applied))),'nonzero_applied_coordinates':sum(d!=0 for d in applied),'new_exact_node_vs_serialization_gaps':[[str(abs(v-q(w)))for v,w in zip(row,rr)]for row,rr in zip(exact1,current)]};write(OUT/'rounding_ledger.json',roundledger)
    physical={'parameter_cell':cell,'parameter_cell_sha256':row['parameter_cell_sha256'],'base_parameters_exact':base,'source_dose_exact':row['dose_exact'],'exact_integral_base_strength_exact':str(F(base['common_gain'])*F(row['dose_exact'])),'original_binary_solver_strength_exact':packet['activation_strength_binary_exact'],'energy_scale_exact':str(q(a['energy_scale'][0])),'full_parameters':packet['parameters'],'unchanged_R':fingerprint(a['preconditioner']),'unchanged_S':fingerprint(S),'unchanged_geometry_and_auxiliary_arrays':{k:fingerprint(v)for k,v in a.items()if k not in('u','preconditioner')},'same_physical_energy_target':True,'changed_only_proof_center_free_intercept':True};write(OUT/'unchanged_target.json',physical)
    domains={'old_whole_coordinate_box_reused':not failed,'old_body_boxes_enclose_new_predictor_plus_old_radius':not failed,'contained_coordinate_count':contained,'failed_coordinate_count':len(failed),'failed_coordinates':failed,'old_Q_scalar_proof_cube_contains_translation_path':dnorm<=rold,'old_Q_scalar_proof_cube_contains_new_same_radius_cube':dnorm+rold<=rold,'old_normal_radius_exact':str(rold),'proposed_common_enclosing_radius_exact':str(common),'new_candidate_root_radius_not_changed_or_certified':str(rold),'fresh_obligations':['Uniform HQ preconditioned defect on a common domain enclosing old predictor, translation path and new same-radius cube; old ZQ bound is not automatically valid outside old radius.','Uniform HI-HQ Hessian defect b including bulk AND surface, all quadrature/integral remainders and finite-cache corrections on the fresh domain.','Parameter/predictor variation V of R(FI-FQ)(p0(theta),theta) relative to its BASE point discrepancy; point H bounds alone do not prove it.','Whole canonical P2 positive-J/no-penalty and global embedding on fresh coordinate boxes; preserve exact axis/fixed constraints.','Canonical surface positive hoop/tangent and G/Js analytic domains on all fresh boundary tubes.','Fresh observer graph monotonicity, top/bottom separation and same-current-radius ROI transport/endpoint bounds if the new root is to support measurements.','Any exact-integral energy minimum needs new uniform HI contraction/positivity proof; fixed-R SPD witness alone is insufficient.']};write(OUT/'domain_obligations.json',domains)
    theorem={'eta_new_inequality':'eta_new <= eta_old + norm(Dpoint-d_applied) + (ZQ_fresh+b_fresh)*norm(d_applied) + V_parameter_predictor','eta_old_uniform_exact':str(eta_old),'Dpoint_minus_applied_d_upper_exact':str(mismatch),'applied_d_norm_exact':str(dnorm),'available_old_ZQ_upper_exact_NOT_fresh':capsule['bounds']['Z']['upper_exact'],'ZQ_fresh':None,'b_fresh':None,'V_parameter_predictor':None,'eta_new_certified':None,'common_domain_proved':False,'root_assertion':False};write(OUT/'conditional_reuse_lemma.json',theorem)
    for f,h in bindings.items():
        if sha(ROOT/f)!=h:raise ValueError('frozen source/input changed during arithmetic preparation '+f)
    checks['R_S_C_physics_dose_normalization_preserved']=True;checks['no_old_geometry_certificate_silently_widened']=True;checks['no_new_root_residual_or_Hessian_evaluation']=True
    write(OUT/'synthetic_tests.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    result={'status':'PREPARED_TRANSLATED_CENTER_FRESH_PROOFS_REQUIRED','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'checks':checks,'artifact_hashes':{f.name:sha(f)for f in OUT.iterdir()if f.is_file()and f.name not in('preparation_receipt.json',)},'applied_decrement_norm_exact':str(dnorm),'Dpoint_minus_applied_d_upper_exact':str(mismatch),'proposed_common_radius_exact':str(common),'new_cube_contained_in_old_body_coordinate_boxes':not failed,'old_body_coordinate_box_failures':len(failed),'all_frozen_bindings_current':True,'new_integral_root_certified':False,'scope':'One explicit proof-center/predictor-intercept candidate only. No new solve, residual, Hessian, integration, uniform transfer, observer theorem, full611 job or scientific parameter change.'};write(OUT/'preparation_receipt.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'receipt_sha256':sha(OUT/'preparation_receipt.json'),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'shift_norm_approx':float(dnorm),'rounding_mismatch_upper_approx':float(mismatch),'common_radius_approx':float(common),'old_body_coordinate_box_failures':len(failed)}),flush=True)
if __name__=='__main__':main()
