"""Prepare/run the base-point parameter derivative of exact-FE/cache force defect.
Run requires a separate bound prelaunch review. Preparation performs no actual
B-field moment integration or quadrature-cache force sum. No root/V claim.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,datetime,hashlib,json
import numpy as np
from flint import arb,arb_mat,ctx
import exact_bulk_bridge_safe as B
import exact_surface_bridge_safe as S
import exact_cache_consistency as BC
import exact_surface_cache_consistency as SC
T=B.T;J=B.JET;core=B.core
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path(__file__)
OUT=ROOT/'runs/exact_integration_bridge_v1/parameter_defect_B_safe_v2'
RECORD=ROOT/B.RECORD
POINT=ROOT/'runs/exact_integration_bridge_v1/complete_point_B_safe_v2/results.json'
POINTREVIEW=ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json'
SURFACE=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2/results.json'
GEOM=ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry'
ORDER=tuple(T.ORDER)

def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):B.write(p,x)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def A(x):return x if isinstance(x,arb)else core.A(Q(x))
def unpack(p):
 lo,hi=Q(p['lower_exact']),Q(p['upper_exact'])
 if lo>hi:raise ValueError('ordered exact enclosure required')
 return A(lo).union(A(hi))
def valid_checks(c):return isinstance(c,dict)and bool(c)and all(v is True for v in c.values())
def bind(d,path,value):
 key=Path(path).relative_to(ROOT).as_posix()if Path(path).is_absolute()else Path(path).as_posix()
 if key in d and d[key]!=value:raise ValueError('conflicting immutable binding '+key)
 if sha(ROOT/key)!=value:raise ValueError('changed immutable input '+key)
 d[key]=value

def review_gate(r,s,p,f):
 if not isinstance(r,dict)or r.get('status')!='PASS'or not valid_checks(r.get('checks'))or r.get('source_sha256')!=s or r.get('protocol_sha256')!=p or r.get('preflight_sha256')!=f:raise ValueError('exact-bound nonempty literal PASS prelaunch review required')

def split_surface(base_force,poly,K,tension):
 """F=K P+(T-K)Q exactly; interval base_force gives a Q enclosure."""
 K,tension=Q(K),Q(tension)
 if K<0 or tension<0 or tension==K:raise ValueError('nonnegative base moduli and nonzero T-K required')
 if len(base_force)!=len(poly)or not base_force:raise ValueError('matching nonempty force vectors required')
 sqrt=[(f-A(K)*p)/A(tension-K)for f,p in zip(base_force,poly)]
 return {'partial_K':[p-q for p,q in zip(poly,sqrt)],'partial_T':sqrt,'sqrt_component':sqrt}

def predictor_bounds(Smat,cell,theta):
 a=np.asarray(Smat)
 if a.ndim!=2 or a.shape!=(611,5)or not np.issubdtype(a.dtype,np.number)or not np.isfinite(a).all():raise ValueError('finite611x5 predictor required')
 delta=[]
 for k in ORDER:
  lo,hi=map(Q,cell[k]);b=Q(theta[k])
  if lo>hi or not lo<=b<=hi:raise ValueError('same closed parameter cell must contain base')
  delta.append((lo-b,hi-b))
 hull=[]
 for row in a:
  terms=[(T.rational(v)*lo,T.rational(v)*hi)for v,(lo,hi)in zip(row,delta)]
  hull.append((sum((min(t)for t in terms),Q(0)),sum((max(t)for t in terms),Q(0))))
 return delta,hull,max(max(abs(lo),abs(hi))for lo,hi in hull)

def conditional_variation(b,predictor_norm,partial_norms,deltas):
 """Arithmetic ONLY: caller must supply a reviewed BOX x C full-H bound b."""
 b,m=Q(b),Q(predictor_norm)
 if b<0 or m<0 or len(partial_norms)!=5 or len(deltas)!=5:raise ValueError('five nonnegative derivative bounds and predictor bound required')
 norms=list(map(Q,partial_norms))
 if any(v<0 for v in norms):raise ValueError('nonnegative derivative bounds required')
 radii=[]
 for lo,hi in deltas:
  lo,hi=Q(lo),Q(hi)
  if lo>hi:raise ValueError('ordered parameter increments required')
  radii.append(max(abs(lo),abs(hi)))
 return b*m+sum((v*r for v,r in zip(norms,radii)),Q(0))

def poly_gradient3(q):
 a,b,c=q;v=b*b+c*c
 return [a*v,a*a*b,a*a*c]

def exact_surface_poly_cache(target,fields,face,edge,rule):
 """POINT Qstar/cache difference for K=T=1, exactly rational before weights."""
 n=len(fields['free_local']);out=[Q(0)]*n;Y=target.nodes()
 for iq in range(len(rule['surface_s'][face])):
  maps=SC.actual_maps(target,face,edge,iq)
  if maps['global_indices']!=fields['global_indices']or maps['selected_local']!=fields['free_local']:raise ValueError('same surface free projection required')
  y=[Y[node][c]for node in maps['nodes']for c in(0,1)];Bs=maps['star_full'];Bc=maps['cache_normalized_full'];qs=[sum((a*b for a,b in zip(row,y)),Q(0))for row in Bs];qc=[sum((a*b for a,b in zip(row,y)),Q(0))for row in Bc]
  gs,gc=poly_gradient3(qs),poly_gradient3(qc);w=Q(maps['auxiliary_normalized_weight_exact'])
  if w!=rule['surface_weights'][face][edge][iq]:raise ValueError('same OLD positive normalized weight once')
  for i,k in enumerate(maps['selected_local']):out[i]+=w*sum((Bs[j][k]*gs[j]-Bc[j][k]*gc[j]for j in range(3)),Q(0))
 return out

def surface_poly_defect(target,face,edge,rule):
 fields=target.surface(face,edge);jet=J.surface_squared_jacobian_jet(fields['q_physical_R'],fields['B_physical_R']);fun=S.SurfaceFunctional(target,fields,face,edge,rule);corr=exact_surface_poly_cache(target,fields,face,edge,rule)
 vals=[sum((A(v/2)*fun.defect(*key)for key,v in p.items()),arb(0))+A(corr[i])for i,p in enumerate(jet['gradient'])]
 return fields,vals,{'signed_polynomial_moment_count':len(fun.requested_defect_keys),'finite_cache_correction_exact':list(map(str,corr))}

def exact_gain_cache(target,fields,index,rule,dose):
 n=len(fields['free_local']);out=[Q(0)]*n;Y=target.nodes()
 for iq in range(len(rule['volume_barycentric'])):
  maps=BC.exact_bulk_maps(target,index,iq)
  if maps['global_indices']!=fields['global_indices']or maps['selected_local']!=fields['free_local']:raise ValueError('same bulk free projection required')
  y=[Y[node][c]for node in maps['nodes']for c in(0,1)];Bs=maps['star_full'];Bc=maps['cache_full'];qs=[sum((a*b for a,b in zip(row,y)),Q(0))for row in Bs];qc=[sum((a*b for a,b in zip(row,y)),Q(0))for row in Bc];w=Q(maps['auxiliary_normalized_weight_exact'])
  if w!=rule['volume_weights'][index][iq]:raise ValueError('same OLD normalized volume weight once')
  for i,k in enumerate(maps['selected_local']):out[i]+=dose*w*sum((Bs[j][k]*qs[j]-Bc[j][k]*qc[j]for j in(1,3)),Q(0))
 return out

def gain_defect(target,index,rule,dose):
 f=target.element(index);n=len(f['free_local'])
 if f['active']==0:return f,[arb(0)for _ in range(n)],{'inactive_exact_zero':True}
 if f['active']!=1:raise ValueError('elementconstant binary active mask required')
 fun=B.Functional(target,f,index,rule);corr=exact_gain_cache(target,f,index,rule,dose);out=[]
 for i in range(n):
  poly={}
  for j in(1,3):poly=T.add(poly,T.mul(f['q'][j],f['B'][j][i]))
  out.append(A(dose)*fun.apply(poly)+A(corr[i]))
 return f,out,{'inactive_exact_zero':False,'signed_polynomial_moment_count':len(fun.requested_defect_keys),'finite_cache_correction_exact':list(map(str,corr))}

def synthetic_tests():
 checks={}
 def ck(k,v):
  checks[k]=bool(v)
  if not v:raise AssertionError(k)
 def reject(k,fn):
  try:fn()
  except (ValueError,KeyError):checks[k]=True;return
  ck(k,False)
 import sympy as sy
 a,b,c,K,tt=sy.symbols('a b c K T',positive=True);G=a*a*(b*b+c*c);E=K*G/2+(tt-K)*sy.sqrt(G)+K/2
 for x in(a,b,c):
  ck('surface_partialK_symbolic_'+str(x),sy.simplify(sy.diff(E,x,K)-sy.diff(G/2-sy.sqrt(G),x))==0)
  ck('surface_partialT_symbolic_'+str(x),sy.simplify(sy.diff(E,x,tt)-sy.diff(sy.sqrt(G),x))==0)
 q=[Q(6,5),Q(4,5),Q(3,5)];pg=poly_gradient3(q);full,_,_=core.surface_local(list(map(A,q)),A(3),A(7),A(1));split=split_surface(full,list(map(A,pg)),3,7)
 for i,v in enumerate([1,Q(24,25),Q(18,25)]):ck('extract_sqrt_positive_den_'+str(i),split['partial_T'][i].contains(A(v)))
 fullneg,_,_=core.surface_local(list(map(A,q)),A(3),A(0),A(1));negative=split_surface(fullneg,list(map(A,pg)),3,0)
 for i,v in enumerate([1,Q(24,25),Q(18,25)]):ck('extract_sqrt_negative_den_'+str(i),negative['partial_T'][i].contains(A(v)))
 reject('zero_TminusK_rejected',lambda:split_surface(full,list(map(A,pg)),3,3));reject('mismatched_vector_rejected',lambda:split_surface(full,[A(1)],3,7))
 x,y,z,w,h,g,dose=sy.symbols('a b c d h g dose');active=(y*y+w*w)*g*dose/2
 ck('gain_dose_derivative_symbolic',sy.diff(active,g)==dose*(y*y+w*w)/2)
 for x in(y,w):ck('gain_strain_force_derivative_'+str(x),sy.diff(active,x,g)==dose*x)
 # Same simple matched maps: normalized pullbacks are identical, while an exact
 # rational geometry perturbation is retained, not rounded to zero.
 identity=[[Q(int(i==j))for j in range(3)]for i in range(3)];v=q;gg=poly_gradient3(v)
 ck('matched_KT1_maps_zero_exact',all(sum(identity[j][i]*gg[j]for j in range(3))-sum(identity[j][i]*gg[j]for j in range(3))==0 for i in range(3)))
 eps=Q(1,2**48);changed=[[a+(eps if i==j==0 else 0)for j,a in enumerate(row)]for i,row in enumerate(identity)];qq=[sum(a*b for a,b in zip(row,v))for row in changed];g2=poly_gradient3(qq)
 ck('mismatched_KT1_map_nonzero_exact',any(sum(changed[j][i]*g2[j]-identity[j][i]*gg[j]for j in range(3))!=0 for i in range(3)))
 cell={k:['-1','1']for k in ORDER};base={k:'0'for k in ORDER};Smat=np.zeros((611,5));Smat[0]=[1,-2,3,0,0];delta,hulls,norm=predictor_bounds(Smat,cell,base)
 ck('predictor_exact_signed_box_hull',hulls[0]==(-6,6)and norm==6);ck('predictor_all_other_rows_zero',all(v==(0,0)for v in hulls[1:]))
 ck('variation_exact_arithmetic',conditional_variation(Q(1,2),6,[1,2,3,4,5],delta)==18)
 reject('negative_b_rejected',lambda:conditional_variation(-1,6,[1]*5,delta));reject('negative_derivative_norm_rejected',lambda:conditional_variation(1,6,[-1]*5,delta));reject('wrong_parameter_count_rejected',lambda:conditional_variation(1,6,[1]*4,delta))
 bad=Smat.copy();bad[0,0]=float('nan');reject('nonfinite_predictor_rejected',lambda:predictor_bounds(bad,cell,base));reject('predictor_shape_rejected',lambda:predictor_bounds(Smat[:3],cell,base))
 valid={'status':'PASS','checks':{'a':True},'source_sha256':'s','protocol_sha256':'p','preflight_sha256':'f'};review_gate(valid,'s','p','f');ck('bound_review_accepted',True)
 for label,c in [('empty',{}),('list',[True]),('false',{'a':False}),('integer',{'a':1}),('string',{'a':'yes'})]:reject('review_bad_checks_'+label,lambda c=c:review_gate({**valid,'checks':c},'s','p','f'))
 for key in('source_sha256','protocol_sha256','preflight_sha256'):reject('review_wrong_'+key,lambda key=key:review_gate({**valid,key:'wrong'},'s','p','f'))
 # Synthetic polynomial integral uses analytic exact moments, not an actual
 # biological FE element: integral_R=0..1 R*(R^2) =1/4, point mass at1/2.
 ck('synthetic_surface_polynomial_defect_exact',S.M.surface_moment(Q(0),Q(1),2).constant-Q(1,2)*Q(1,2)**2==Q(1,8))
 return checks

def prepare():
 if(OUT/'protocol.json').exists():raise ValueError('write-once frozen preparation already exists')
 ctx.prec=512;core.A.cache_clear();checks=synthetic_tests();bindings={};row=read(RECORD);point=read(POINT);audit=read(POINTREVIEW)
 if audit.get('status')!='PASS_SAFE_V2_COMPLETE_POINT_REDUCTION'or not valid_checks(audit.get('checks'))or audit['source_and_input_hashes'].get(POINT.relative_to(ROOT).as_posix(),audit['source_and_input_hashes'].get(str(POINT.relative_to(ROOT))))!=sha(POINT)or audit['old_v1_claims_restored']is not False:raise ValueError('independent SAFE complete-point review must bind actual basepoint')
 for path,h in audit['source_and_input_hashes'].items():bind(bindings,path,h)
 for item in audit['all168_packets']:bind(bindings,item['file'],item['sha256'])
 independent_reduction=ROOT/'runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/exact_168_reduction.json'
 if sha(independent_reduction)!=audit['exact_reduction_sha256']:raise ValueError('independent exact reduction hash required')
 bind(bindings,independent_reduction,sha(independent_reduction))
 for obj in [SOURCE,Path(B.__file__),Path(S.__file__),Path(T.__file__),Path(J.__file__),Path(BC.__file__),Path(SC.__file__),Path(B.M.__file__),Path(core.__file__),Path(BC.target_core.__file__),RECORD,POINT,POINTREVIEW,SURFACE,ROOT/row['point_cache_file'],ROOT/row['predictor_file'],GEOM/'protocol.json',GEOM/'coordinate_boxes.json',ROOT/'evidence/exact_retube_B_geometry_independent_review.json',ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json']:
  bind(bindings,obj,sha(obj))
 # Freeze recursive exact map helper dependency; only its Fraction map factories
 # are consumed, not its old interval-polynomial operations.
 bind(bindings,ROOT/'scripts/exact_fe_target.py',sha(ROOT/'scripts/exact_fe_target.py'))
 geometry=read(GEOM/'protocol.json')
 for path,h in geometry['source_and_input_hashes'].items():bind(bindings,path,h)
 target=T.Target(read(ROOT/row['point_cache_file']));theta={k:Q(v)for k,v in point['base_parameters_exact'].items()};expected={k:T.rational(target.packet['parameters'][k])for k in ORDER if k!='common_gain'};expected['common_gain']=Q(row['base_strength_exact'])/Q(row['dose_exact'])
 if theta!=expected or row['parameter_cell_sha256']!=read(GEOM/'coordinate_boxes.json')['parameter_cell_sha256']:raise ValueError('exact common theta/base/dose/C identity required')
 with np.load(ROOT/row['predictor_file'],allow_pickle=False)as z:Smat=z['S'].copy()
 delta,hulls,m=predictor_bounds(Smat,row['parameter_cell'],theta);box=read(GEOM/'coordinate_boxes.json');centers=target.nodes()
 if len(hulls)!=len(target.a['free'])or len(box['current_centers_exact'])!=325:raise ValueError('same611/650 coordinates required')
 for k,(lo,hi)in enumerate(hulls):
  dof=int(target.a['free'][k]);node,component=divmod(dof,2);l=Q(box['current_centers_exact'][node][component])-Q(box['coordinate_radii_exact'][node][component]);h=Q(box['current_centers_exact'][node][component])+Q(box['coordinate_radii_exact'][node][component]);base=centers[node][component]
  if not l<=base<=h or not l<=base+lo<=base+hi<=h:raise ValueError('full independent common box must contain u0 and predictor endpoints')
 samehulls=all((Q(a),Q(b))==hh for(a,b),hh in zip(box['predictor_hulls_exact'],hulls))
 if not samehulls:raise ValueError('predictor theta convention differs from fresh geometry')
 surface=read(SURFACE);edges=[];known={(Path(item['file']).as_posix(),item['sha256'])for item in audit['all168_packets']if item['kind']=='surface'}
 if {(v['face'],v['edge'])for v in surface['edges']}!={(a,b)for a in(0,1)for b in range(12)}:raise ValueError('all24 basepoint surface edges required')
 for edge in surface['edges']:
  face= edge['face'];path=ROOT/edge['file'];d=read(path);f=target.surface(face,edge['edge']);K=theta['basal_area_modulus'if face==0 else'apical_area_modulus'];tt=theta['basal_surface_tension'if face==0 else'apical_surface_tension']
  if tt==K or d['global_indices']!=f['global_indices']or len(d['force_integral_minus_cache'])!=len(f['free_local'])or (Path(edge['file']).as_posix(),sha(path))not in known:raise ValueError('audited same-base surface packet and nonzero T-K required')
  edges.append({**edge,'base_K_exact':str(K),'base_T_exact':str(tt),'T_minus_K_exact':str(tt-K)})
 p={'source_and_input_hashes':bindings,'record_id':row['record_id'],'parameter_order':list(ORDER),'parameter_cell':row['parameter_cell'],'parameter_cell_sha256':row['parameter_cell_sha256'],'base_parameters_exact':{k:str(v)for k,v in theta.items()},'dose_exact':row['dose_exact'],'point_cache_file':row['point_cache_file'],'point_cache_sha256':sha(ROOT/row['point_cache_file']),'predictor_file':row['predictor_file'],'predictor_sha256':sha(ROOT/row['predictor_file']),'base_point_review_sha256':sha(POINTREVIEW),'common_geometry_boxes_file':(GEOM/'coordinate_boxes.json').relative_to(ROOT).as_posix(),'common_geometry_boxes_sha256':sha(GEOM/'coordinate_boxes.json'),'parameter_deltas_exact':[[str(a),str(b)]for a,b in delta],'predictor_displacement_norm_upper_exact':str(m),'all611_predictor_hulls_match_fresh_geometry':samehulls,'all_u0_to_predictor_segments_inside_independent_coordinate_box':True,'surface_edges':edges,'bulk_active_elements':[i for i in range(144)if T.rational(target.a['active'][i,0])==1],'precision_bits':512,'new_square_root_integrations':False,'no_new_equilibrium_optimization':True,'run_policy':'After bound review, compute only polynomial moment defects+exact same-Y finite-cache corrections, recover surface sqrt derivative from audited basepoint intervals. No full nonlinear reintegration.','full_b_requirement':'Uniform bound for ||R(HI-HQ)||infinity on INDEPENDENT full canonical nodal BOX times same C. Coupled predictor-family tube alone is insufficient. Bulk AND surface required; no full b is assumed here.','variation_formula':'V <= b_full * sup_C||S(theta-theta0)||infinity + sum_j max_C|theta_j-theta0_j| * ||R partial_theta_j(FI-FQ)(u0)||infinity.','V_established':False,'exact_integral_root_certified':False}
 write(OUT/'synthetic_tests.json',{'status':'PASS_SYNTHETICS_ONLY','checks':checks});write(OUT/'protocol.json',p);write(OUT/'preflight.json',{'status':'PASS_PREPARATION_ONLY','checks':{**checks,'same_SAFE_audited_basepoint':True,'all24_surface_packets_bound':True,'exact_predictor_hulls_match':samehulls,'independent_box_segment_containment':True},'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'synthetic_tests_sha256':sha(OUT/'synthetic_tests.json'),'actual_B_polynomial_moments_or_cache_sums_run':False,'V_established':False})
 print(json.dumps({'status':'PASS_PREPARATION_ONLY','checks':len(checks)+4,'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'predictor_norm_approx':float(m)}),flush=True)

def run():
 p=read(OUT/'protocol.json');review_gate(read(OUT/'prelaunch_independent_review.json'),sha(SOURCE),sha(OUT/'protocol.json'),sha(OUT/'preflight.json'))
 if any((OUT/k).exists()for k in ['started.json','results.json','terminal.json']):raise ValueError('single actual polynomial job only')
 for k,h in p['source_and_input_hashes'].items():
  if sha(ROOT/k)!=h:raise ValueError('changed frozen input '+k)
 ctx.prec=p['precision_bits'];core.A.cache_clear();write(OUT/'started.json',{'utc':now(),'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json')});target=T.Target(read(ROOT/p['point_cache_file']));rule=target.auxiliary_rule();dF=[[arb(0)for _ in range(611)]for _ in ORDER];rows=[]
 for edge in p['surface_edges']:
  face,index=edge['face'],edge['edge'];f,poly,info=surface_poly_defect(target,face,index,rule);packet=read(ROOT/edge['file']);split=split_surface([unpack(v)for v in packet['force_integral_minus_cache']],poly,Q(edge['base_K_exact']),Q(edge['base_T_exact']));ik=1 if face==0 else 0;it=3 if face==0 else 2
  for j,gid in enumerate(f['global_indices']):dF[ik][gid]+=split['partial_K'][j];dF[it][gid]+=split['partial_T'][j]
  dest=OUT/'surfaces'/f'face{face}_edge{index:02d}.json';write(dest,{'face':face,'edge':index,'global_indices':f['global_indices'],'source_base_packet_sha256':sha(ROOT/edge['file']),'poly_force_defect':[core.pack(v)for v in poly],'partial_K':[core.pack(v)for v in split['partial_K']],'partial_T':[core.pack(v)for v in split['partial_T']],'diagnostics':info});rows.append({'kind':'surface','file':dest.relative_to(ROOT).as_posix(),'sha256':sha(dest)});print(json.dumps({'phase':'surface_parameter_polynomial','face':face,'edge':index}),flush=True)
 for index in p['bulk_active_elements']:
  f,vals,info=gain_defect(target,index,rule,Q(p['dose_exact']))
  for j,gid in enumerate(f['global_indices']):dF[4][gid]+=vals[j]
  dest=OUT/'bulk_gain'/f'{index:03d}.json';write(dest,{'element':index,'global_indices':f['global_indices'],'gain_partial_force_defect':[core.pack(v)for v in vals],'diagnostics':info});rows.append({'kind':'bulk_gain','file':dest.relative_to(ROOT).as_posix(),'sha256':sha(dest)});print(json.dumps({'phase':'gain_polynomial','element':index}),flush=True)
 R=arb_mat(611,611,[A(T.rational(v))for v in target.a['preconditioner'].ravel()]);partials=[];norms=[]
 for key,col in zip(ORDER,dF):
  result=R*arb_mat(611,1,col);norm=core.upper(max(core.upper(result[i,0])for i in range(611)));norms.append(Q(str(norm.fmpq())));partials.append({'parameter':key,'force_partial':[core.pack(v)for v in col],'preconditioned_force_partial':[core.pack(result[i,0])for i in range(611)],'infinity_norm_upper':core.pack(norm)})
 deltas=[tuple(map(Q,v))for v in p['parameter_deltas_exact']];param_bound=sum((n*max(abs(a),abs(b))for n,(a,b)in zip(norms,deltas)),Q(0))
 for k,h in p['source_and_input_hashes'].items():
  if sha(ROOT/k)!=h:raise ValueError('input drift during polynomial job '+k)
 result={'status':'BASE_POINT_PARAMETER_FORCE_DEFECT_COMPONENT_UNREVIEWED','protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(SOURCE),'local_packets':rows,'parameter_partials':partials,'parameter_only_variation_upper_exact':str(param_bound),'predictor_displacement_norm_upper_exact':p['predictor_displacement_norm_upper_exact'],'full_b_bound_supplied':False,'V_established':False,'exact_integral_root_certified':False,'scope':'Exact-affine parameter part at SAME basepoint plus predictor norm; full independent-box bulk+surface Hessian bound and root-reuse conditions remain separate.'};write(OUT/'results.json',result);write(OUT/'terminal.json',{'status':'FINISHED','exit_code':0,'results_sha256':sha(OUT/'results.json'),'utc':now()});print(json.dumps({'status':result['status'],'parameter_bound_approx':float(param_bound)}),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);a=ap.parse_args();prepare()if a.command=='prepare'else run()
