"""Independent scalar/identity PRELAUNCH review; no producer output or integrals."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import ast,base64,hashlib,io,json,sys,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;sys.path.insert(0,str(ROOT/'scripts'))
import exact_integral_root_B_safe as W
PROD=W.OUT
EXPECTED_SOURCE='8404a27fa339f4990c9a3fb784a84f5576cfa1fea8eaa20f1f58f8cf6ff6633c'
EXPECTED_PROTOCOL='ce2a816ba037725395cc4dc2b6a8f419edd1c645c071d542893615b228fefd41'
EXPECTED_PREFLIGHT='547f58b221ee0ef1d9857b908edd5f7d5543fefe7cf650755691badc4d8dc1df'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:serial(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [serial(v)for v in x]
 return x
checks={}
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
def reject(k,fn):
 try:fn()
 except (ValueError,AssertionError,KeyError,StopIteration):checks[k]=True;return
 ck(k,False)
def scalar(etaQ,m,zQ,b,d,s,param,r):
 etaQ,m,zQ,b,d,s,param,r=map(Q,(etaQ,m,zQ,b,d,s,param,r))
 if any(v<0 for v in(etaQ,m,zQ,b,d,s,param))or r<=0:raise ValueError('invalid scalar premises')
 z=zQ+b;v=b*s+param;eta=etaQ+m+v+z*d;lhs=eta+z*r
 return {'z':z,'V':v,'eta':eta,'selfmap_lhs':lhs,'strict_sufficient_conditions':bool(z<1 and lhs<r)}
def main():
 if(PROD/'prelaunch_independent_review.json').exists()or(OUT/'results.json').exists():raise ValueError('write-once review already exists')
 ck('frozen_source',sha(W.SOURCE)==EXPECTED_SOURCE);ck('frozen_protocol',sha(PROD/'protocol.json')==EXPECTED_PROTOCOL);ck('frozen_preflight',sha(PROD/'preflight.json')==EXPECTED_PREFLIGHT)
 p=read(PROD/'protocol.json');pre=read(PROD/'preflight.json');ck('not_executed_or_asserted',not(PROD/'results.json').exists()and pre['root_certified']is False and p['root_asserted_at_preparation']is False);ck('preflight_literal_checks_andbindings',W.passed(pre['checks'])and pre['source_sha256']==EXPECTED_SOURCE and pre['protocol_sha256']==EXPECTED_PROTOCOL)
 for k,h in p['source_and_input_hashes'].items():ck('frozen_'+k,sha(ROOT/k)==h)
 docs={k:read(ROOT/v)for k,v in p['source_files'].items()if k not in ['vectors','lemma']};record=docs['retube']['selected_record'];cache=read(ROOT/record['point_cache_file']);root=read(ROOT/record['root_certificate_file'])
 for k in ['point_cache_file','predictor_file','root_certificate_file','independent_root_capsule_file']:ck('authoritative_selected_'+k,sha(ROOT/record[k])==record[k.replace('_file','_sha256')])
 raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True);ck('embedded_binary_cache',hashlib.sha256(raw).hexdigest()==cache['frozen_binary_cache']['sha256']==record['frozen_binary_cache_sha256'])
 with np.load(io.BytesIO(raw),allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
 with np.load(ROOT/p['source_files']['vectors'],allow_pickle=False)as z:u0=z['u0'].copy();u1=z['u1'].copy()
 ck('611_finite_vectors_and_matrix',u0.shape==u1.shape==(611,)and arrays['preconditioner'].shape==(611,611)and np.isfinite(u0).all()and np.isfinite(u1).all()and np.isfinite(arrays['preconditioner']).all());ck('all611_u0_exact_binary_match',all(Q.from_float(float(a))==Q.from_float(float(b))for a,b in zip(u0,arrays['u'])))
 d=[Q.from_float(float(a))-Q.from_float(float(b))for a,b in zip(u0,u1)];dn=max(map(abs,d));point=docs['point_exact'];m=max(max(abs(Q(a)-v),abs(Q(b)-v))for(a,b),v in zip(point['RF_endpoint_sums'],d))
 ck('all611_point_endpoints_ordered',len(point['RF_endpoint_sums'])==611 and all(Q(a)<=Q(b)for a,b in point['RF_endpoint_sums']));ck('point_exact_auditbinding',sha(ROOT/p['source_files']['point_exact'])==docs['point_review']['exact_reduction_sha256'])
 # Rnorm is independently recomputed by rational binary row sums; no matrix
 # products or Hessian assembly. This tests values, not just norm metadata.
 rnorm=max(sum((abs(Q.from_float(float(v)))for v in row),Q(0))for row in arrays['preconditioner'])
 ck('actual_binary_Rnorm',rnorm==Q(point['Rnorm_exact'])==Q(docs['energy_B']['exact_R_infinity_norm']))
 param=docs['parameter'];zq=Q(docs['cached_review']['ZQ_upper_exact']);b=Q(docs['physical_bound']['norm_bounds']['b_absR_H_upper_exact']);b0=Q(docs['first_bound']['norm_bounds']['b_absR_H_upper_exact']);etaQ=max(Q(point['eta_author_exact']),Q(point['eta_independent_exact']));r=Q(p['root_cube_radius_exact']);M=Q(param['predictor_displacement_norm_upper_exact']);vp=Q(param['parameter_only_variation_upper_exact'])
 candidate=scalar(etaQ,m,zq,b,dn,M,vp,r);first=scalar(etaQ,m,zq,b0,dn,M,vp,r)
 # Evaluate the exact source's identity predicates and in-memory scalar flow,
 # removing only gate-loading and final output effects. No producer run/file.
 tree=ast.parse(W.SOURCE.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='run');body=[]
 for n in run.body:
  if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='p'for t in n.targets):continue
  if isinstance(n,ast.Expr)and isinstance(n.value,ast.Call)and isinstance(n.value.func,ast.Name)and n.value.func.id in ['gate','write','print']:continue
  body.append(n)
 env=dict(W.__dict__);env.update({'p':p});exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),'<actual_identity_flow_without_IO>','exec'),env)
 for k,v in env['checks'].items():ck('actual_composer_predicate_'+k,v)
 ck('independent_scalar_matches_actual_chosen',candidate==env['chosen']);ck('independent_scalar_matches_actual_first',first==env['first']);ck('independent_decrement_and_m',dn==env['dnorm']and m==env['m']);ck('only_exact_rational_absR_b_used',env['b']==b and env['firstb']==b0)
 e=candidate['eta']/(1-candidate['z'])if candidate['strict_sufficient_conditions']else None;curv=(1-candidate['z'])/rnorm if e is not None else None;loc=max(Q(0),abs(d[588])-e)if e is not None else None;eq=etaQ/(1-zq)if e is not None else None;gap=max(Q(0),abs(d[588])-e-eq)if e is not None else None
 ck('independent_error_curvature_location_gap',e==env['error']and curv==env['curvature']and loc==env['location_lower']and eq==env['oldQ_error']and gap==env['branch_gap']);ck('exclusion_boolean_exact',env['old_cube_excluded']==bool(e is not None and loc>r))
 # Root existence is NOT asserted by this prelaunch receipt even if the
 # scalar candidate is positive. Actual composition remains parent's task.
 fixtures=0
 for a,z,bb,ds,ss,rad in product((Q(0),Q(1,100)),(Q(0),Q(1,3)),(Q(0),Q(2,3)),(Q(0),Q(1,7)),(Q(0),Q(1,5)),(Q(1,100),Q(1))):
  args=(a,Q(1,1000),z,bb,ds,ss,Q(1,10000),rad);ck('independent_scalar_fixture_'+str(fixtures),scalar(*args)==W.sufficient(*args));fixtures+=1
 for ix in range(8):
  args=[Q(0)]*8;args[-1]=Q(1);args[ix]=Q(-1);reject('negative_scalar_rejected_'+str(ix),lambda args=args:W.sufficient(*args))
 reject('zero_radius_rejected',lambda:W.sufficient(0,0,0,0,0,0,0,0))
 valid={'status':'PASS','checks':{'one':True},'source_sha256':'s','protocol_sha256':'p','preflight_sha256':'f'};W.gate(valid,'s','p','f');ck('gate_accepts_bound_schema',True)
 for label,bad in [('empty',{}),('list',[True]),('false',{'x':False}),('truthy_integer',{'x':1}),('truthy_string',{'x':'true'}),('missing',None)]:reject('gate_rejects_'+label,lambda bad=bad:W.gate({**valid,'checks':bad},'s','p','f'))
 for key in ['source_sha256','protocol_sha256','preflight_sha256']:reject('gate_rejects_'+key,lambda key=key:W.gate({**valid,key:'bad'},'s','p','f'))
 # The theorem is verified analytically in proof_review.md. Standard sufficient
 # conditions are not manufactured by sampling the physical Hessian spectrum.
 ck('common_box_independence_V_explicit','independent'in read(ROOT/p['source_files']['parameter_protocol'])['full_b_requirement'].lower())
 ck('scope_excludes_PDE_all17_biology','Not all17'in p['scope']and'PDE'in p['scope']and'confirmed novelty'in p['scope'])
 # Adversarial identity mismatch: execute actual source predicate with a single
 # altered declared cell/dose/base/domain/cache/u0 field, never modify files.
 identity_cases=[('C_hash',"all(x['parameter_cell_sha256']==p['parameter_cell_sha256']for x in [geometry,cached,physical,spd])",'geometry','parameter_cell_sha256','wrong'),('dose',"cached['dose_exact']==p['dose_exact']==record['dose_exact']",'cached','dose_exact','1'),('geometry_domain',"physical['common_geometry_domain_id']==box['domain_id']==docs['cached_protocol']['geometry_domain_id']",'physical','common_geometry_domain_id','wrong')]
 for label,expr,name,key,value in identity_cases:
  en=dict(env);en[name]={**env[name],key:value};ck('identity_mismatch_rejected_'+label,eval(expr,en)is False)
 result={'status':'PASS','checks':checks,'source_sha256':EXPECTED_SOURCE,'protocol_sha256':EXPECTED_PROTOCOL,'preflight_sha256':EXPECTED_PREFLIGHT,'review_source_sha256':sha(__file__),'initial_findings_sha256':sha(OUT/'initial_findings.json'),'candidate_scalar_evaluation':serial(candidate),'first_representation_scalar_evaluation':serial(first),'independent_exact_decrement_norm':str(dn),'independent_point_mismatch':str(m),'independent_exact_Rnorm':str(rnorm),'conditional_error':str(e)if e is not None else None,'conditional_curvature':str(curv)if curv is not None else None,'conditional_old_cube_distance_lower':str(loc)if loc is not None else None,'conditional_I_Q_component_gap_lower':str(gap)if gap is not None else None,'actual_composer_or_integrals_run':False,'root_certified_by_this_prelaunch_receipt':False,'theory_review':'See proof_review.md. Scalar evaluation and predicate replay do not replace the pending parent composition; no matrices/integrals/equilibrium solves assembled.','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 W.write(OUT/'results.json',result);W.write(PROD/'prelaunch_independent_review.json',result);print(json.dumps({'status':'PASS','checks':len(checks),'review_sha256':sha(PROD/'prelaunch_independent_review.json'),'conditional_z':float(candidate['z']),'conditional_eta':float(candidate['eta']),'conditional_selfmap':float(candidate['selfmap_lhs']),'conditional_error':float(e)if e else None,'conditional_location':float(loc)if loc else None,'conditional_gap':float(gap)if gap else None}),flush=True)
if __name__=='__main__':main()
