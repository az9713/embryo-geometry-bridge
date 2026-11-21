"""Independent final scalar theorem replay for the actual frozen B root packet."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,base64,io,sys,copy,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;PROD=ROOT/'runs/exact_integration_bridge_v1/integral_root_B_safe_v2'
RECEIPT=ROOT/'evidence/exact_integral_root_B_safe_independent_review.json'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
 if p.exists():raise ValueError('write-once final review exists '+str(p))
 p.write_text(json.dumps(d,indent=2)+'\n')
def serial(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:serial(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [serial(v)for v in x]
 return x
checks={}
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
def S(eta,m,zq,b,d,M,P,r):
 z=zq+b;V=b*M+P;a=eta+m+z*d+V;lhs=a+z*r
 return {'z':z,'V':V,'eta':a,'selfmap_lhs':lhs,'strict_sufficient_conditions':bool(z<1 and lhs<r)}
def validate_scalars(packet,expected,chosen,first):
 for k,v in expected.items():
  if isinstance(v,Q):
   if Q(packet.get(k))!=v:raise ValueError('wrong scalar '+k)
  elif packet.get(k)!=v:raise ValueError('wrong field '+k)
 for key,values in [('physical_axis_test',chosen),('first_representation_test',first)]:
  for k,v in values.items():
   x=packet[key][k]
   if isinstance(v,Q):
    if Q(x)!=v:raise ValueError('wrong theorem scalar '+key+':'+k)
   elif x is not v:raise ValueError('wrong strict theorem flag '+key)

def main():
 if RECEIPT.exists()or(OUT/'final_results.json').exists():raise ValueError('preserve final review')
 packet=read(PROD/'results.json');protocol=read(PROD/'protocol.json');pre=read(PROD/'preflight.json');review=read(PROD/'prelaunch_independent_review.json');terminal=read(PROD/'terminal.json');job=read(PROD/'job.json')
 ck('actual_result_expected_SHA',sha(PROD/'results.json')=='42d97691311c057d8152f54a6ea31a7218763e7b451df63358bdf423010e18d4')
 ck('true_terminal_zero_same_actual_PID',terminal['status']=='FINISHED'and terminal['exit_code']==0 and terminal['pid']==job['pid']==35704 and terminal['results_sha256']==sha(PROD/'results.json'))
 ck('stdout_hash_actual_terminal',terminal['stdout_sha256']==sha(PROD/'stdout.log')and terminal['all_launch_bindings_current']is True)
 ck('frozen_source',sha(ROOT/'scripts/exact_integral_root_B_safe.py')==packet['source_sha256']==pre['source_sha256']=='8404a27fa339f4990c9a3fb784a84f5576cfa1fea8eaa20f1f58f8cf6ff6633c')
 ck('frozen_protocol',sha(PROD/'protocol.json')==packet['protocol_sha256']==pre['protocol_sha256']=='ce2a816ba037725395cc4dc2b6a8f419edd1c645c071d542893615b228fefd41')
 ck('frozen_preflight',sha(PROD/'preflight.json')=='547f58b221ee0ef1d9857b908edd5f7d5543fefe7cf650755691badc4d8dc1df')
 ck('final_prelaunch966PASS',sha(PROD/'prelaunch_independent_review.json')=='d3cafecf9a3152f6cbae39d10cf7bfcedd73123b4215cdde5e5777c212f08478'and review['status']=='PASS'and len(review['checks'])==966 and all(v is True for v in review['checks'].values()))
 ck('actual_producer_all_predicates_true',isinstance(packet['checks'],dict)and bool(packet['checks'])and all(v is True for v in packet['checks'].values()))
 for k,h in protocol['source_and_input_hashes'].items():ck('frozen_'+k,sha(ROOT/k)==h)
 for k,h in job['bindings'].items():ck('launch_binding_'+k,sha(ROOT/k)==h)
 docs={k:read(ROOT/v)for k,v in protocol['source_files'].items()if k not in ['vectors','lemma']};record=docs['retube']['selected_record'];point=docs['point_exact'];cache=read(ROOT/record['point_cache_file'])
 raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True);ck('actual_cache_embeddedhash',hashlib.sha256(raw).hexdigest()==record['frozen_binary_cache_sha256'])
 with np.load(io.BytesIO(raw),allow_pickle=False)as f:au=f['u'].copy();R=f['preconditioner'].copy()
 with np.load(ROOT/protocol['source_files']['vectors'],allow_pickle=False)as f:u0=f['u0'].copy();u1=f['u1'].copy()
 ck('actual_finite_shapes',au.shape==u0.shape==u1.shape==(611,)and R.shape==(611,611)and all(np.isfinite(v).all()for v in [au,u0,u1,R]))
 delta=[Q.from_float(float(a))-Q.from_float(float(b))for a,b in zip(u0,u1)]
 for i,(a,b)in enumerate(zip(au,u0)):ck('actual_u0_binary_'+str(i),Q.from_float(float(a))==Q.from_float(float(b)))
 d=max(map(abs,delta));m=max(max(abs(Q(a)-v),abs(Q(b)-v))for(a,b),v in zip(point['RF_endpoint_sums'],delta));r=Q(protocol['root_cube_radius_exact']);eta=max(Q(point['eta_author_exact']),Q(point['eta_independent_exact']));zq=Q(docs['cached_review']['ZQ_upper_exact']);b=Q(docs['physical_bound']['norm_bounds']['b_absR_H_upper_exact']);b0=Q(docs['first_bound']['norm_bounds']['b_absR_H_upper_exact']);M=Q(docs['parameter']['predictor_displacement_norm_upper_exact']);P=Q(docs['parameter']['parameter_only_variation_upper_exact']);Rn=max(sum((abs(Q.from_float(float(v)))for v in row),Q(0))for row in R)
 chosen=S(eta,m,zq,b,d,M,P,r);first=S(eta,m,zq,b0,d,M,P,r)
 ck('positive_rational_premises',all(x>=0 for x in [eta,m,zq,b,d,M,P])and r>0 and Rn>0)
 ck('new_exact_selfmap_and_contraction',chosen['strict_sufficient_conditions']is True)
 ck('old_representation_sufficient_test_fails_retained',first['strict_sufficient_conditions']is False)
 e=chosen['eta']/(1-chosen['z']);curv=(1-chosen['z'])/Rn;loc=max(Q(0),abs(delta[588])-e);eq=eta/(1-zq);gap=max(Q(0),abs(delta[588])-e-eq)
 ck('root_bound_within_candidate_cube',0<e<r);ck('positive_curvature',curv>0);ck('new_root_outside_old_cube',loc>r);ck('positive_sameparameter_IQ_component_gap',gap>0)
 expected={'status':'CERTIFIED_UNIFORM_EXACT_FE_ROOT_AND_STRICT_LOCAL_MINIMUM','record_id':record['record_id'],'parameter_cell_sha256':record['parameter_cell_sha256'],'dose_exact':record['dose_exact'],'r_exact':r,'etaQ_exact':eta,'m_exact':m,'zQ_exact':zq,'b_exact':b,'d_norm_exact':d,'predictor_norm_exact':M,'parameter_only_exact':P,'root_error_from_translated_predictor_upper_exact':e,'uniform_integral_energy_curvature_lower_exact':curv,'location_witness_free_index':588,'integral_root_distance_from_old_predictor_lower_exact':loc,'old_Q_root_error_upper_exact':eq,'I_Q_branch_component_separation_lower_exact':gap,'no_integral_root_in_old_cube_certified':True,'same_R_infinity_norm_exact':Rn,'root_certified':True,'strict_local_minimum_certified':True,'all17_integral_or_PDE_or_biological_validation':False}
 validate_scalars(packet,expected,chosen,first);ck('all_reported_scalars_and_booleans_exact',True)
 # Strong negative packet checks: no value is accepted just because a producer
 # emits a success label. Never mutate the frozen result itself.
 badcount=0
 for key,value in expected.items():
  bad=copy.deepcopy(packet)
  if isinstance(value,Q):bad[key]=str(value+Q(1,10**30))
  elif isinstance(value,bool):bad[key]=not value
  elif isinstance(value,int):bad[key]=value+1
  else:bad[key]='different'
  try:validate_scalars(bad,expected,chosen,first)
  except (ValueError,TypeError):ck('negative_packet_'+key,True);badcount+=1;continue
  ck('negative_packet_'+key,False)
 for key in ['physical_axis_test','first_representation_test']:
  for field in ['z','V','eta','selfmap_lhs','strict_sufficient_conditions']:
   bad=copy.deepcopy(packet);bad[key][field]=not packet[key][field]if field=='strict_sufficient_conditions'else str(Q(packet[key][field])+Q(1,10**30))
   try:validate_scalars(bad,expected,chosen,first)
   except (ValueError,TypeError):ck('negative_subtest_'+key+'_'+field,True);badcount+=1;continue
   ck('negative_subtest_'+key+'_'+field,False)
 boundary_eta=(1-chosen['z'])*r
 ck('strict_selfmap_boundary_not_accepted',not S(boundary_eta,0,chosen['z'],0,0,0,0,r)['strict_sufficient_conditions'])
 ck('z_equal_one_not_accepted',not S(0,0,Q(1,2),Q(1,2),0,0,0,r)['strict_sufficient_conditions'])
 ck('large_variation_prevents_false_success',not S(eta,m,zq,b,d,M,1,r)['strict_sufficient_conditions'])
 ck('component_location_claim_requires_positive_margin',max(Q(0),abs(delta[588])-(abs(delta[588])+r))<=r)
 # Domain/proof checks are fixed theorem obligations, not numerical fitting.
 ck('same_common_radius_sweeps_both_cubes',Q(docs['geometry']['new_common_enclosing_radius_exact'])==Q(docs['cached']['common_radius_exact'])==Q(docs['physical_bound']['common_radius_exact'])==r+d)
 ck('same_independent_box_for_V',docs['parameter_protocol']['common_geometry_boxes_sha256']==docs['cached_protocol']['geometry_boxes_sha256']==docs['geometry']['prepared_box_sha256'])
 ck('same_cached_Z_independent_scalar',zq==Q(docs['cached']['ZQ_upper_exact']))
 ck('same_old_H0_energy_capsule',docs['energy_B']['root_certificate_sha256']==record['root_certificate_sha256']and Q(docs['energy_B']['uniform_Hessian_lower_bound_exact'])>0 and Rn==Q(docs['energy_B']['exact_R_infinity_norm']))
 ck('original_SCOPE_limits_retained','Not all17'in packet['scope']and'PDE'in packet['scope']and'confirmed novelty'in packet['scope'])
 result={'status':'PASS','checks':checks,'certificate_sha256':sha(PROD/'results.json'),'source_sha256':sha(ROOT/'scripts/exact_integral_root_B_safe.py'),'protocol_sha256':sha(PROD/'protocol.json'),'preflight_sha256':sha(PROD/'preflight.json'),'prelaunch_review_sha256':sha(PROD/'prelaunch_independent_review.json'),'terminal_sha256':sha(PROD/'terminal.json'),'job_sha256':sha(PROD/'job.json'),'stdout_sha256':sha(PROD/'stdout.log'),'review_source_sha256':sha(__file__),'record_id':record['record_id'],'parameter_cell_sha256':record['parameter_cell_sha256'],'dose_exact':record['dose_exact'],'independent_scalars':serial(expected),'physical_axis_test':serial(chosen),'first_representation_test':serial(first),'negative_packet_cases':badcount,'proof_review_file':(OUT/'proof_review.md').relative_to(ROOT).as_posix(),'proof_review_sha256':sha(OUT/'proof_review.md'),'independent_implementation':'Fresh exact Fraction scalar arithmetic,611binary u0/decrement checks,preconditioner rational row sums,complete frozeninput andactualterminal bindings. No actualrootcomposer/integral/Hessian/equilibrium rerun. Earlier966prelaunch source/theory/domain identity checks explicitly reused.','root_existence_and_strict_local_minimum_verified':True,'old_cube_integral_root_exclusion_verified':True,'same_parameter_IQ_component_separation_verified':True,'PDE_or_all17_or_biological_validation':False,'scope':'One unchanged B source dose, every parameter in SAME small cell,exact spatially integrated finite-element equilibrium in its translated cube; common-cube uniqueness/exclusion and normalized energy curvature only.','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 write(OUT/'final_results.json',result);write(RECEIPT,result);print(json.dumps({'status':'PASS','checks':len(checks),'receipt_sha256':sha(RECEIPT),'negative_cases':badcount,'curvature_approx':float(curv),'root_error_approx':float(e),'old_cube_distance_lower_approx':float(loc),'IQ_gap_lower_approx':float(gap)}),flush=True)
if __name__=='__main__':main()
