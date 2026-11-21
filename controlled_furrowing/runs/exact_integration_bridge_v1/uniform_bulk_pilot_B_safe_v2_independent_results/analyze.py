"""Read-only terminal7 analysis; no actual element integrals/cache sums."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib,datetime
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_bulk_defect_B_safe as U
P=U.P;T=U.T;J=U.JET;core=U.core
PROD=ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_pilot_B_safe_v2'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def iv(x):return Q(x['lower_exact']),Q(x['upper_exact'])
def upper(x):return max(abs(Q(str(x.lower().fmpq()))),abs(Q(str(x.upper().fmpq())))) if isinstance(x,arb) else abs(Q(x))
def radius(x):return (Q(str(x.upper().fmpq()))-Q(str(x.lower().fmpq())))/2 if isinstance(x,arb) else Q(0)
def stats(xs):
 v=[iv(x) for x in xs];return {'norm_upper_exact':str(max(max(abs(l),abs(h)) for l,h in v)),'norm_lower_exact':str(max(max(l,-h,Q(0)) for l,h in v)),'max_halfwidth_exact':str(max((h-l)/2 for l,h in v))}
def overlap(a,b):return max(iv(a)[0],iv(b)[0])<=min(iv(a)[1],iv(b)[1])
def contains(a,b):return iv(a)[0]<=iv(b)[0] and iv(b)[1]<=iv(a)[1]
def harmonic(n):return sum((Q(1,j) for j in range(1,n+1)),Q(0))
def tail(n,r):
 m=n+1;k=n+2;gap=1-r;c=1+5000*harmonic(m);d=1+5000*harmonic(k);s=Q(5000,k+1)
 return r**m*(c/gap+Q(5000,m+1)*r/gap**2),r**m*(k*d/gap+(d+k*s)*r/gap**2+s*r*(1+r)/gap**3)
checks={}
def ck(k,b):
 checks[k]=bool(b)
 if not b:raise AssertionError(k)
def majorant(p,f,which=upper):
 if f['chart']['kind']=='off_axis':
  lo=min(x[0]for x in f['corners']);hi=max(x[0]for x in f['corners']);z=max(abs(x[1]-f['zcenter'])for x in f['corners'])
  return sum((which(c)*max(lo**i,hi**i)*z**j for(i,j),c in p.items()),Q(0))
 return sum((which(c)for c in p.values()),Q(0))
ctx.prec=512;core.A.cache_clear();p=read(PROD/'protocol.json');res=read(PROD/'results.json');term=read(PROD/'terminal.json');pre=read(PROD/'prelaunch_independent_review.json')
ck('results_SHA_authorized_terminal',sha(PROD/'results.json')=='2f21379bf3e86e2d1bdb6663ca686c4585a709fbf8a2fba9e9be266f16aa7a0c')
ck('terminal_FINISHED_zero',term['status']=='FINISHED'and term['exit_code']==0 and term['results_sha256']==sha(PROD/'results.json'))
ck('prelaunch612_exactbinding',sha(PROD/'prelaunch_independent_review.json')=='49f69791f98ecc4e5fc1f335194f71cd7b0a2f1c512e0b73b23d6f4641033efe'and pre['status']=='PASS'and len(pre['checks'])==612 and all(v is True for v in pre['checks'].values())and pre['protocol_sha256']==sha(PROD/'protocol.json')and pre['source_sha256']==sha(U.__file__))
for k,h in p['source_and_input_hashes'].items():ck('locked_'+k,sha(ROOT/k)==h)
ck('results_protocol',res['protocol_sha256']==sha(PROD/'protocol.json'))
ck('exact7', [x['element']for x in res['elements']]==[0,2,6,8,72,74,66])
record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule();box=read(U.GEOM/'coordinate_boxes.json');Y=U.field_nodes(box);base=target.nodes();body=read(U.GEOM/'certification_run/whole_body_geometry.json')
ck('all650_base_nodes_inside_domain',all(Q(c)-Q(r)<=x<=Q(c)+Q(r)for row,cc,rr in zip(base,box['current_centers_exact'],box['coordinate_radii_exact'])for x,c,r in zip(row,cc,rr)))
ck('same_C_basegain_contained',p['parameter_cell']==record['parameter_cell']and Q(p['parameter_cell']['common_gain'][0])<=Q(record['base_strength_exact'])/Q(record['dose_exact'])<=Q(p['parameter_cell']['common_gain'][1]))
rows=[]
for manifest,sel in zip(res['elements'],p['selectors']):
 i=manifest['element'];path=PROD/manifest['file'];d=read(path);pointfile=ROOT/f'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/elements/{i:03d}.json';pt=read(pointfile)
 ck(f'packet_SHA_{i}',sha(path)==manifest['sha256']);ck(f'packetselector_{i}',d['element']==sel['element']and d['degree']==sel['degree']and Q(d['rho_exact'])==Q(sel['rho_exact'])and d['active']==sel['active']);ck(f'allpacketchecks_{i}',all(v is True for v in d['checks'].values()))
 f=target.element(i,Y);jet=J.bulk_jacobian_jet(f['q'],f['B']);fp=target.element(i);jp=J.bulk_jacobian_jet(fp['q'],fp['B']);n=len(f['global_indices']);H={(a['i'],a['j']):a for a in d['H_integral_minus_cache_upper']};HP={(a['i'],a['j']):a for a in pt['H_integral_minus_cache_upper_triangle']};ff=d['force_integral_minus_cache'];pf=pt['force_integral_minus_cache']
 ck(f'complete_free_map_{i}',d['global_indices']==pt['global_indices']==f['global_indices']and len(ff)==n and set(H)==set(HP)=={(a,b)for a in range(n)for b in range(a,n)})
 ck(f'basecoeff_J_inclusion_{i}',all(Q(str(v.lower().fmpq()))<=jp['value'].get(k,Q(0))<=Q(str(v.upper().fmpq()))for k,v in jet['value'].items()))
 ck(f'force_overlap_{i}',all(overlap(a,b)for a,b in zip(ff,pf)));ck(f'H_overlap_{i}',all(overlap(H[k]['enclosure'],HP[k]['enclosure'])for k in H))
 rho=Q(sel['rho_exact']);lo,hi=iv(body['triangles'][i]['J_bound']);raw=max(abs(lo-1),abs(hi-1));scaled=raw*2**48
 ck(f'rho_outward_exact_{i}',rho==Q(-(-scaled.numerator//scaled.denominator),2**48)and lo>Q(1,10**6))
 t1,t2=tail(d['degree'],rho);ck(f'independent_scalar_tails_{i}',(t1,t2)==U.M.volumetric_derivative_tail_bounds(d['degree'],rho,1,5000))
 # Scalar mass only; no nonlinear integral or cache64 field evaluation.
 mass=2*arb.pi()/U.A(target.energy_scale)*U.A(f['detDM']*sum(x[0]for x in f['corners'])/6)+U.A(sum(rule['volume_weights'][i]));masshi=Q(str(mass.upper().fmpq()))
 norms=[majorant(v,f)for v in jet['gradient']]
 for a,stored in enumerate(d['force_two_mass_tail_upper']):ck(f'independent_Ftail_{i}_{a}',iv(stored)[0]<=masshi*t1*norms[a]<=iv(stored)[1])
 for (a,b),v in H.items():ck(f'independent_Htail_{i}_{a}_{b}',iv(v['two_mass_tail_upper'])[0]<=masshi*(t2*norms[a]*norms[b]+t1*majorant(jet['hessian_upper'][a,b],f))<=iv(v['two_mass_tail_upper'])[1])
 ck(f'cache_positive_domain_{i}',d['cache_consistency']['all64_convex_strain_lines_positive']is True and iv(d['cache_consistency']['minimum_cache_convex_line_J_lower'])[0]>Q(1,10**6))
 sf=stats(ff);sh=stats([v['enclosure']for v in H.values()]);sp=stats(pf);hp=stats([v['enclosure']for v in HP.values()]);ft=max(iv(x)[1]for x in d['force_two_mass_tail_upper']);ht=max(iv(x['two_mass_tail_upper'])[1]for x in H.values())
 row={'element':i,'chart':f['chart']['kind'],'active':d['active'],'degree':d['degree'],'rho_exact':str(rho),'local_free_dofs':n,'source_packet_sha256':sha(path),'point_packet_sha256':sha(pointfile),'force':sf,'H_entries':sh,'point_force':sp,'point_H_entries':hp,'force_tail_upper_exact':str(ft),'H_tail_upper_exact':str(ht),'point_force_intervals_contained_count':sum(contains(a,b)for a,b in zip(ff,pf)),'point_H_intervals_contained_count':sum(contains(H[k]['enclosure'],HP[k]['enclosure'])for k in H),'J_minus_one_raw_coefficient_majorant_exact':str(majorant(T.sub(jet['value'],T.ONE),f)),'J_minus_one_raw_coefficient_radius_majorant_exact':str(majorant(T.sub(jet['value'],T.ONE),f,radius))}
 if f['chart']['kind']=='off_axis':
  rep=J.shifted_numerator(T.sub(jet['value'],T.ONE),f['corners']);z=max(abs(x[1]-f['zcenter'])for x in f['corners']);den=rep['denominator_power'];spoly=rep['numerator']
  row['shifted_J_minus_one_majorant_exact']=str(sum((upper(v)*z**j for(i,j),v in spoly.items()),Q(0))/rep['Rlo']**den)
  row['shifted_J_minus_one_radius_majorant_exact']=str(sum((radius(v)*z**j for(i,j),v in spoly.items()),Q(0))/rep['Rlo']**den)
 for k in ['force','H_entries','point_force','point_H_entries']:row[k+'_approx']={a:float(Q(b))for a,b in row[k].items()}
 rows.append(row)
summary={'status':'PASS_LOCAL_OUTPUT_INSPECTION_ONLY','checks':checks,'producer_results_sha256':sha(PROD/'results.json'),'source_sha256':sha(U.__file__),'protocol_sha256':sha(PROD/'protocol.json'),'prelaunch_review_sha256':sha(PROD/'prelaunch_independent_review.json'),'analysis_source_sha256':sha(__file__),'rows':rows,'actual_integral_or_cache64_replays':False,'global611_or_root_claim':False,'decomposition_limit':'Uniform packets retain total F/H and tails, not separate polynomial/cache terms; attribution of non-tail width cannot be certified from outputs alone.','point_overlap_scope':'Necessary diagnostic only; reviewed domain and algebra justify uniform intervals, not overlap alone.','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
P.write(OUT/'results.json',summary)
print(json.dumps({'status':summary['status'],'checks':len(checks),'rows':[{'e':r['element'],'chart':r['chart'],'F':r['force_approx'],'H':r['H_entries_approx'],'pointF':r['point_force_approx']['norm_upper_exact'],'pointH':r['point_H_entries_approx']['norm_upper_exact'],'ft':float(Q(r['force_tail_upper_exact'])),'ht':float(Q(r['H_tail_upper_exact'])),'rawJmajorant':float(Q(r['J_minus_one_raw_coefficient_majorant_exact'])),'shiftJmajorant':float(Q(r['shifted_J_minus_one_majorant_exact']))if 'shifted_J_minus_one_majorant_exact'in r else None}for r in rows]}),flush=True)
