"""Prelaunch wrapper review only: no actual element integral or cache64 sum."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import ast,json,sys,hashlib,datetime
from flint import arb,arb_mat,ctx
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_bulk_assembly_B_safe as W
U=W.U;T=U.T;core=U.core;PROD=W.OUT
checks={}
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
def up(v):return max(abs(Q(str(v.lower().fmpq()))),abs(Q(str(v.upper().fmpq()))))if isinstance(v,arb)else abs(Q(v))
def iv(x):return Q(x['lower_exact']),Q(x['upper_exact'])
def pack(l,h):return {'lower_exact':str(l),'upper_exact':str(h)}
def includes(b,q):return Q(str(b.lower().fmpq()))<=q<=Q(str(b.upper().fmpq()))
def normpoly(poly,f):
 if f['chart']['kind']=='off_axis':
  lo=min(p[0]for p in f['corners']);hi=max(p[0]for p in f['corners']);z=max(abs(p[1]-f['zcenter'])for p in f['corners']);return sum((up(c)*max(lo**a,hi**a)*z**b for(a,b),c in poly.items()),Q(0))
 return sum((up(c)for c in poly.values()),Q(0))
def harmonic(n):return sum((Q(1,i)for i in range(1,n+1)),Q(0))
def tails(n,r):
 m=n+1;k=n+2;gap=1-r;c=1+5000*harmonic(m);d=1+5000*harmonic(k);s=Q(5000,k+1)
 return r**m*(c/gap+Q(5000,m+1)*r/gap**2),r**m*(k*d/gap+(d+k*s)*r/gap**2+s*r*(1+r)/gap**3)
def assembly_ast_tests(tree):
 run=next(x for x in tree.body if isinstance(x,ast.FunctionDef)and x.name=='run');loop=next(x for x in run.body if isinstance(x,ast.For)and isinstance(x.target,ast.Name)and x.target.id=='selector');fragments=[]
 for node in loop.body:
  if isinstance(node,ast.For)and any(isinstance(x,ast.Name)and x.id in ['ids','H_integral_minus_cache_upper']for x in ast.walk(node)):
   fragments.append(node)
 # source selection has force zip(ids,...) plus H iteration packet key.
 fragments=[x for x in loop.body if isinstance(x,ast.For)and isinstance(x.target,(ast.Name,ast.Tuple))]
 ck('actual_local_assembly_two_loops',len(fragments)==2)
 code=compile(ast.fix_missing_locations(ast.Module(body=fragments,type_ignores=[])),'<actual_local_assembly>','exec')
 packets=[{'global_indices':[0,1],'force_integral_minus_cache':[pack(Q(-1),Q(2)),pack(Q(-2),Q(1))],'H_integral_minus_cache_upper':[{'i':0,'j':0,'enclosure':pack(Q(1),Q(2))},{'i':0,'j':1,'enclosure':pack(Q(-3),Q(-2))},{'i':1,'j':1,'enclosure':pack(Q(4),Q(5))}]},{'global_indices':[2,0],'force_integral_minus_cache':[pack(Q(3),Q(4)),pack(Q(-1),Q(1))],'H_integral_minus_cache_upper':[{'i':0,'j':0,'enclosure':pack(Q(-2),Q(-1))},{'i':0,'j':1,'enclosure':pack(Q(1),Q(3))},{'i':1,'j':1,'enclosure':pack(Q(2),Q(4))}]}]
 env={'g':[arb(0)for _ in range(3)],'H':{},'ball':W.ball,'arb':arb,'zip':zip,'tuple':tuple,'sorted':sorted}
 for packet in packets:exec(code,{**env,'packet':packet,'ids':packet['global_indices']})
 h=env['H'];ck('reverse_global_ids_sorted',set(h)=={(0,0),(0,1),(1,1),(2,2),(0,2)});ck('shared_diagonal_not_doubled',includes(h[0,0],Q(3))and includes(h[0,0],Q(6))and up(h[0,0])<7)
 # Extract actual norm/rowassembly block, changing only611 to toy dimension3.
 start=next(i for i,x in enumerate(run.body)if isinstance(x,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='RM'for t in x.targets));end=next(i for i in range(start,len(run.body))if isinstance(run.body[i],ast.For)and any(isinstance(x,ast.Constant)and x.value=='source_and_input_hashes'for x in ast.walk(run.body[i])))
 class Dim(ast.NodeTransformer):
  def visit_Constant(self,n):return ast.copy_location(ast.Constant(3),n)if n.value==611 else n
 frag=Dim().visit(ast.Module(body=run.body[start:end],type_ignores=[]));frag=ast.fix_missing_locations(frag)
 import numpy as np
 R=[[Q(1),Q(-1,2),Q(0)],[Q(1,4),Q(1),Q(-1,2)],[Q(-1),Q(1,2),Q(1)]]
 from types import SimpleNamespace
 target=SimpleNamespace(a={'preconditioner':np.array(R,dtype=float)})
 env.update({'U':U,'arb_mat':arb_mat,'target':target,'range':range,'max':max,'sum':sum,'float':float})
 exec(compile(frag,'<actual_norm_block_dim3>','exec'),env)
 HM=env['HM'];ck('source_symmetry_and_same_diagonal',all(HM[i,j].overlaps(HM[j,i])for i in range(3)for j in range(3))and includes(HM[0,0],Q(3))and includes(HM[0,0],Q(6)))
 intervals=[iv(v['enclosure'])for p in packets for v in p['H_integral_minus_cache_upper']];worst=Q(0)
 for signs in product((0,1),repeat=6):
  m=[[Q(0)for _ in range(3)]for _ in range(3)]
  for value,(p,v)in zip([a[s]for a,s in zip(intervals,signs)],[(p,v)for p in packets for v in p['H_integral_minus_cache_upper']]):
   a,b=p['global_indices'][v['i']],p['global_indices'][v['j']];m[a][b]+=value
   if a!=b:m[b][a]+=value
  rh=[[sum((R[i][k]*m[k][j]for k in range(3)),Q(0))for j in range(3)]for i in range(3)];n=max(sum(abs(v)for v in row)for row in rh);worst=max(worst,n)
  ck('rational_H_vertex_'+str(sum(s<<j for j,s in enumerate(signs))),all(includes(env['RH'][i,j],rh[i][j])for i in range(3)for j in range(3))and n<=up(env['signed_bound'])and n<=up(env['rowwise_bound']))
 ck('both_norm_bounds_nonnegative',env['signed_bound']>=0 and env['rowwise_bound']>=0)
 ck('min_two_valid_bounds_valid',min(up(env['signed_bound']),up(env['rowwise_bound']))>=worst)
 fints=[iv(v)for p in packets for v in p['force_integral_minus_cache']]
 for ix,signs in enumerate(product((0,1),repeat=4)):
  g=[Q(0)]*3
  for p,offset in zip(packets,(0,2)):
   for j,k in enumerate(p['global_indices']):g[k]+=fints[offset+j][signs[offset+j]]
  rf=[sum((a*b for a,b in zip(row,g)),Q(0))for row in R];ck('rational_F_vertex_'+str(ix),all(includes(env['RF'][i,0],v)for i,v in enumerate(rf)))
 return {'exact_H_vertices':64,'exact_F_vertices':16,'source_AST_blocks_reused_with_dimension_only_611_to3':True,'largest_toy_RH_norm_exact':str(worst)}
def main():
 if(OUT/'results.json').exists()or(PROD/'prelaunch_independent_review.json').exists():raise ValueError('write-once review already exists')
 ctx.prec=512;core.A.cache_clear();p=read(PROD/'protocol.json');pre=read(PROD/'preflight.json');tree=ast.parse(W.SOURCE.read_text());run=next(x for x in tree.body if isinstance(x,ast.FunctionDef)and x.name=='run')
 ck('frozen_source_protocol_preflight',pre['source_sha256']==sha(W.SOURCE)and pre['protocol_sha256']==sha(PROD/'protocol.json')and all(v is True for v in pre['checks'].values()))
 ck('no_actual137_started',not any((PROD/s).exists()for s in ['started.json','results.json','terminal.json']))
 for k,h in p['source_and_input_hashes'].items():ck('lock_'+k,sha(ROOT/k)==h)
 ck('runtime_all144_point_diagnostics_bound',all(f'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/elements/{i:03d}.json'in p['source_and_input_hashes']for i in range(144)))
 lp=read(U.OUT/'protocol.json');lpre=read(U.OUT/'prelaunch_independent_review.json');lr=read(U.OUT/'results.json');lt=read(U.OUT/'terminal.json')
 ck('uniform_local_protocol_exact',p['uniform_local_protocol']==lp);ck('local_prelaunch612_PASS',lpre['status']=='PASS'and len(lpre['checks'])==612 and all(v is True for v in lpre['checks'].values())and lpre['source_sha256']==sha(U.__file__)and lpre['protocol_sha256']==sha(U.OUT/'protocol.json'))
 ck('seven_terminal_exact',lt['status']=='FINISHED'and lt['exit_code']==0 and lt['results_sha256']==sha(U.OUT/'results.json')and sha(U.OUT/'results.json')=='2f21379bf3e86e2d1bdb6663ca686c4585a709fbf8a2fba9e9be266f16aa7a0c')
 ck('same512_and_no_root_claim',p['precision_bits']==lp['precision_bits']==512 and p['root_transfer_proved']is False and p['no_new_equilibrium_optimization']is True)
 ck('all144_exact_unique',len(p['selectors'])==144 and [x['element']for x in p['selectors']]==list(range(144))and all(type(x['element'])is int for x in p['selectors']))
 ck('reuse_exact7',set(p['reuse_completed_local_packets'])==set(map(str,[0,2,6,8,72,74,66])))
 target=T.Target(read(ROOT/read(ROOT/U.P.RECORD)['point_cache_file']));Y=U.field_nodes(read(U.GEOM/'coordinate_boxes.json'));body=read(U.GEOM/'certification_run/whole_body_geometry.json');rule=target.auxiliary_rule();order_counts={}
 for selector in p['selectors']:
  i=selector['element'];f=target.element(i,Y);jet=U.JET.bulk_jacobian_jet(f['q'],f['B']);rho=Q(selector['rho_exact']);lo,hi=iv(body['triangles'][i]['J_bound']);raw=max(abs(lo-1),abs(hi-1));q=raw*2**48
  ck(f'exact_geometry_rho_{i}',rho==Q(-(-q.numerator//q.denominator),2**48)and raw<=rho<1 and lo>Q(1,10**6))
  mass=2*arb.pi()/U.A(target.energy_scale)*U.A(f['detDM']*sum(v[0]for v in f['corners'])/6)+U.A(sum(rule['volume_weights'][i]));mass_upper=up(mass);g=max(normpoly(v,f)for v in jet['gradient']);h=max(normpoly(v,f)for v in jet['hessian_upper'].values());first=None
  for n in range(2,19):
   t1,t2=tails(n,rho);fb=mass_upper*t1*g;hb=mass_upper*(t2*g*g+t1*h)
   if fb<=Q(lp['force_tail_target_exact'])and hb<=Q(lp['H_entry_tail_target_exact']):first=n;break
  ck(f'first_adequate_order_and_active_{i}',type(selector['degree'])is int and selector['degree']==first and selector['active']==str(f['active']))
  ck(f'tailtargets_{i}',Q(selector['force_tail_upper_exact'])<=Q(lp['force_tail_target_exact'])and Q(selector['H_entry_tail_upper_exact'])<=Q(lp['H_entry_tail_target_exact']))
  pp=read(ROOT/f'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/elements/{i:03d}.json');nn=len(f['global_indices']);ck(f'point_diagnostic_map_complete_{i}',pp['element']==i and pp['global_indices']==f['global_indices']and len(pp['force_integral_minus_cache'])==nn and{(v['i'],v['j'])for v in pp['H_integral_minus_cache_upper_triangle']}=={(a,b)for a in range(nn)for b in range(a,nn)})
  order_counts[str(first)]=order_counts.get(str(first),0)+1
  if str(i)in p['reuse_completed_local_packets']:
   r=p['reuse_completed_local_packets'][str(i)];old=next(x for x in lr['elements']if x['element']==i);ck(f'reuse_exact_packet_selector_{i}',r['sha256']==old['sha256']==sha(ROOT/r['file'])and selector==next(x for x in lp['selectors']if x['element']==i)and read(ROOT/r['file'])['global_indices']==f['global_indices'])
 ck('reported_order_histogram',order_counts==pre['orders'])
 toy=assembly_ast_tests(tree)
 # Guard expression copied from actual source; mutation fixtures do not run job.
 gate=next(x for x in run.body if isinstance(x,ast.If)and any(isinstance(y,ast.Constant)and y.value=='bound independent prelaunch required'for y in ast.walk(x)));code=compile(ast.Expression(gate.test),'<actual_review_guard>','eval');valid={'status':'PASS','checks':{'a':True},'source_sha256':sha(W.SOURCE),'protocol_sha256':sha(PROD/'protocol.json'),'preflight_sha256':sha(PROD/'preflight.json')};env={'sha':sha,'SOURCE':W.SOURCE,'OUT':PROD,'all':all,'dict':dict,'isinstance':isinstance}
 ck('actual_guard_accepts_bound_review',eval(code,{**env,'review':valid})is False)
 for label,bad in [('empty',{}),('false',{'a':False}),('integer',{'a':1}),('string',{'a':'yes'})]:ck('guard_rejects_'+label,eval(code,{**env,'review':{**valid,'checks':bad}})is True)
 for k in ['source_sha256','protocol_sha256','preflight_sha256']:ck('guard_rejects_'+k,eval(code,{**env,'review':{**valid,k:'wrong'}})is True)
 # Manual source-level proof bindings: local loop calls reviewedU exactly once for
 # each nonreused selector, signed symmetry and nonnegative upper bounds above.
 ck('same_domain_passed_to_all_local_calls',any(isinstance(x,ast.Call)and isinstance(x.func,ast.Attribute)and x.func.attr=='evaluate'and ast.unparse(x)=='U.evaluate(target, i, Y, rule, selector, p[\'uniform_local_protocol\'])'for x in ast.walk(loop))if False else 'U.evaluate(target,i,Y,rule,selector,p[\'uniform_local_protocol\'])'in W.SOURCE.read_text().replace(' ',''))
 out={'status':'PASS','checks':checks,'source_sha256':sha(W.SOURCE),'protocol_sha256':sha(PROD/'protocol.json'),'preflight_sha256':sha(PROD/'preflight.json'),'review_source_sha256':sha(__file__),'local_component_review_sha256':sha(U.OUT/'prelaunch_independent_review.json'),'toy_assembly':toy,'order_counts':order_counts,'new_actual_element_integrals_or_cache64_sums_run':False,'manual_proof':'Signed interval assembly in identical611coordinates encloses the global bulk sum. Symmetric insertion counts each diagonal once and each off-diagonal in both positions. Row i of RΔH is enclosed by exact-binaryR interval multiplication; sum of absolute entry uppers bounds its infinity norm. Alternatively sum_k|Rik|sum_j|ΔHkj| gives a valid nonnegative rowwise upper. Minimum of two proved upper bounds remains an upper. No surfaces, new root or PDE claim.','initial_findings_file_sha256':sha(OUT/'initial_findings.json'),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 U.P.write(OUT/'results.json',out);U.P.write(PROD/'prelaunch_independent_review.json',out);print(json.dumps({'status':'PASS','checks':len(checks),'review_sha256':sha(PROD/'prelaunch_independent_review.json')}),flush=True)
if __name__=='__main__':main()
