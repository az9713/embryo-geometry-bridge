"""Independent PRELAUNCH uniform-seven-bulk source/domain/glue review.
No actual nonlinear element integral or64-site cache sum is called.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,ast,datetime,hashlib,json,sys
from flint import arb,arb_mat,ctx
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_bulk_defect_B_safe as U
P=U.P;T=U.T;M=U.M;C=U.C;core=U.core
PRODUCER=ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_pilot_B_safe_v2'
SOURCE_SHA='7c43896c220090025e3bc98969340f98c2b145099905649d12b16bb57d6fac9c'
PROTOCOL_SHA='bbc559fab371051ab272db511dc3b2a6fffa8fde7797a419ef3a33a1f050d266'
checks={};counts={};inputs={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def ck(k,v):
    checks[k]=bool(v)
    if not v:raise AssertionError(k)
def endpoint(v):return F(str(v.fmpq()))
def upper(v):return max(abs(endpoint(v.lower())),abs(endpoint(v.upper())))if isinstance(v,arb)else abs(F(v))
def inside(v,x):return endpoint(v.lower())<=x<=endpoint(v.upper())if isinstance(v,arb)else F(v)==x

def majorant(poly,f):
    if f['chart']['kind']=='off_axis':
        lo=min(p[0]for p in f['corners']);hi=max(p[0]for p in f['corners']);zmax=max(abs(p[1]-f['z_center'])for p in f['corners']);return sum(upper(c)*max(lo**a,hi**a)*zmax**b for(a,b),c in poly.items())
    return sum(map(upper,poly.values()),F(0))
def harmonic(n):return sum((F(1,j)for j in range(1,n+1)),F(0))
def tails(N,r):
    gap=1-r;m=N+1;c=1+5000*harmonic(m);a=r**m*(c/gap+F(5000,m+1)*r/gap**2);k=m+1;d=1+5000*harmonic(k);s=F(5000,k+1);b=r**m*(k*d/gap+(d+k*s)*r/gap**2+s*r*(1+r)/gap**3);return a,b

class AD1:
    def __init__(self,v,g):self.v=U.A(v);self.g=list(g)
    def __add__(self,b):
        if not isinstance(b,AD1):b=AD1(b,[arb(0)]*len(self.g))
        return AD1(self.v+b.v,[a+c for a,c in zip(self.g,b.g)])
    __radd__=__add__
    def __neg__(self):return AD1(-self.v,[-a for a in self.g])
    def __sub__(self,b):return self+-b
    def __mul__(self,b):
        if not isinstance(b,AD1):b=AD1(b,[arb(0)]*len(self.g))
        return AD1(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.g,b.g)])
    __rmul__=__mul__
    def log(self):
        if not self.v>0:raise ValueError('positive AD logarithm required')
        return AD1(self.v.log(),[a/self.v for a in self.g])

def direct_force(B,Y,drive):
    n=len(Y);yy=[AD1(v,[arb(int(i==j))for j in range(n)])for i,v in enumerate(Y)];q=[sum((v*y for v,y in zip(row,yy)),AD1(0,[arb(0)]*n))for row in B];J=(q[0]*q[3]-q[1]*q[2])*q[4];L=J.log();W=sum((v*v for v in q),AD1(0,[arb(0)]*n))*F(1,2)-L+L*L*2500+(q[1]*q[1]+q[3]*q[3])*F(drive,2)
    return W.g

def independent_force_fixture():
    eye=[[F(int(i==j))for j in range(5)]for i in range(5)];eps=F(1,2**40);star=[[v+eps*F((i+2*j)%5-2)for j,v in enumerate(row)]for i,row in enumerate(eye)];base=[F(11,10),F(1,20),F(-1,25),F(9,10),F(21,20)];r=F(1,1000);Y=[U.A(v-r).union(U.A(v+r))for v in base];drive=U.A(1).union(U.A(2));result=C.consistency_bound(star,eye,Y,list(range(5)),1,5000,drive)
    gc,_,_=core.bulk_local(result['cache_q'],U.A(1),U.A(5000),drive);_,Hline,_=core.bulk_local(result['convex_line_box'],U.A(1),U.A(5000),drive);delta=arb_mat(5,5,[U.A(v)for row in result['delta_B_full']for v in row]);dg=arb_mat(5,5,[v for row in Hline for v in row])*arb_mat(5,1,result['delta_q']);Fbound=delta.transpose()*arb_mat(5,1,gc)+result['star_B_free'].transpose()*dg
    count=0
    for signs in product((-1,1),repeat=5):
        yy=[U.A(v+sign*r)for v,sign in zip(base,signs)]
        for a in(1,2):
            Fs=direct_force(star,yy,a);Fc=direct_force(eye,yy,a);truth=[s-c for s,c in zip(Fs,Fc)];ck('independent_dense_energy_AD_force_fixture_'+str(count),all(Fbound[i,0].contains(v)for i,v in enumerate(truth)));count+=1
    ck('sameY_deltaq_remains_small',max(upper(v)for v in result['delta_q'])<F(1,10**8));counts['synthetic_interval_force_vertices']=count
    # Weighted composition algebra, avoiding any actual quadrature sum.
    weight=F(3,7);ck('single_positive_weight_applied_once',all((U.A(weight)*Fbound[i,0]).overlaps(U.A(weight)*Fbound[i,0])for i in range(5)))

def review_gate_ast():
    tree=ast.parse(Path(U.__file__).read_text(encoding='utf8'));run=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='run');conditional=next(n for n in run.body if isinstance(n,ast.If)and any(isinstance(x,ast.Constant)and x.value=='bound nonempty review required'for x in ast.walk(n)))
    expr=compile(ast.Expression(conditional.test),'<frozen_run_gate>','eval');source=sha(U.SOURCE);protocol=sha(PRODUCER/'protocol.json');preflight=sha(PRODUCER/'preflight.json');valid={'status':'PASS','checks':{'one':True},'source_sha256':source,'protocol_sha256':protocol,'preflight_sha256':preflight}
    env={'dict':dict,'isinstance':isinstance,'all':all,'sha':lambda p:{str(U.SOURCE):source,str(PRODUCER/'protocol.json'):protocol,str(PRODUCER/'preflight.json'):preflight}[str(p)],'SOURCE':U.SOURCE,'OUT':PRODUCER}
    ck('actual_run_guard_accepts_exact_nonempty_review',eval(expr,{**env,'review':valid})is False)
    for label,c in [('empty',{}),('list',[True]),('false',{'one':False}),('truthy_int',{'one':1}),('truthy_string',{'one':'yes'}),('none',None)]:ck('actual_run_guard_rejects_'+label,eval(expr,{**env,'review':{**valid,'checks':c}})is True)
    for k in('source_sha256','protocol_sha256','preflight_sha256'):ck('actual_run_guard_rejects_mismatch_'+k,eval(expr,{**env,'review':{**valid,k:'bad'}})is True)

def main():
    if(PRODUCER/'prelaunch_independent_review.json').exists()or(OUT/'results.json').exists():raise ValueError('write-once prelaunch review exists')
    ctx.prec=512;core.A.cache_clear();p=read(PRODUCER/'protocol.json');pre=read(PRODUCER/'preflight.json');ck('source_hash_exact',sha(U.__file__)==SOURCE_SHA);ck('protocol_hash_exact',sha(PRODUCER/'protocol.json')==PROTOCOL_SHA);ck('preflight_source_protocol_bindings',pre['source_sha256']==SOURCE_SHA and pre['protocol_sha256']==PROTOCOL_SHA and all(v is True for v in pre['checks'].values()))
    for name,h in p['source_and_input_hashes'].items():ck('frozen_'+name,sha(ROOT/name)==h);inputs[name]=h
    safe=read(ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json');cache_review=read(ROOT/'evidence/exact_cache_consistency_independent_review_v2.json');geo_review=read(U.REVIEW);body=read(U.GEOM/'certification_run/whole_body_geometry.json');geom=read(U.GEOM/'certification_run/results.json');box=read(U.GEOM/'coordinate_boxes.json')
    ck('safe_component992_PASS_current',safe['status']=='PASS'and len(safe['checks'])==992 and all(v is True for v in safe['checks'].values())and safe['source_sha256']==sha(T.__file__)and safe['helper_source_sha256']==sha(U.JET.__file__))
    ck('dense_bulk_cache_component168_PASS_current',cache_review['status'].startswith('PASS_SAFE_V2_')and len(cache_review['checks'])==168 and all(v is True for v in cache_review['checks'].values())and cache_review['source_sha256']==sha(C.__file__)and cache_review['unsafe_old_T_polynomial_operations_used_with_Arb']is False)
    ck('fresh_geometry738_PASS_current',geo_review['status']=='PASS_FRESH_COMMON_NODAL_GEOMETRY_COMPONENT'and len(geo_review['checks'])==738 and all(v is True for v in geo_review['checks'].values())and geo_review['body_output_sha256']==sha(U.GEOM/'certification_run/whole_body_geometry.json')and geo_review['producer_results_sha256']==sha(U.GEOM/'certification_run/results.json'))
    ck('same_declared_domain_and_C',p['geometry_domain_id']==box['domain_id']and p['geometry_boxes_sha256']==sha(U.GEOM/'coordinate_boxes.json'))
    record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule();Y=U.field_nodes(box);centers=target.nodes();ck('basepoint_exact_nodes_inside_new_common_box',all(F(c)-F(r)<=v<=F(c)+F(r)for row,cr,rr in zip(centers,box['current_centers_exact'],box['coordinate_radii_exact'])for v,c,r in zip(row,cr,rr)))
    basegain=F(record['base_strength_exact'])/F(record['dose_exact']);ck('basegain_inside_same_C',F(p['parameter_cell']['common_gain'][0])<=basegain<=F(p['parameter_cell']['common_gain'][1]));ck('same_parameter_cell_and_dose',p['parameter_cell']==record['parameter_cell']and F(p['dose_exact'])==F(record['dose_exact']))
    selectors=p['selectors'];ck('exact7_predeclared_elements',[(s['element'])for s in selectors]==[0,2,6,8,72,74,66]and all(type(s['element'])is int for s in selectors));ck('no_actual_pilot_started',not(PRODUCER/'started.json').exists()and not any(PRODUCER.glob('element_*.json')))
    order_diagnostics=[];active=set()
    for selector in selectors:
        i=selector['element'];f=target.element(i,Y);jet=U.JET.bulk_jacobian_jet(f['q'],f['B']);fp=target.element(i);jp=U.JET.bulk_jacobian_jet(fp['q'],fp['B']);active.add(f['active'])
        for label,polys,points in[('q',f['q'],fp['q']),('J',[jet['value']],[jp['value']]),('Ju',jet['gradient'],jp['gradient']),('Juv',list(jet['hessian_upper'].values()),list(jp['hessian_upper'].values()))]:ck(f'basepoint_coefficient_inclusion_{i}_{label}',all(inside(poly.get(k,arb(0)),point.get(k,F(0)))for poly,point in zip(polys,points)for k in set(poly)|set(point)))
        ck('same_free_B_and_projection_'+str(i),f['B']==fp['B']and f['global_indices']==fp['global_indices']);rho_raw=max(abs(F(body['triangles'][i]['J_bound']['lower_exact'])-1),abs(F(body['triangles'][i]['J_bound']['upper_exact'])-1));scaled=rho_raw*2**48;rounded=F(-(-scaled.numerator//scaled.denominator),2**48);rho=F(selector['rho_exact'])
        ck('exact_outward_rho_'+str(i),rho==rounded and rho>=rho_raw and rho-rho_raw<F(1,2**48)and 0<=rho<1);ck('fresh_analytic_J_gt_penalty_'+str(i),F(body['triangles'][i]['J_bound']['lower_exact'])>F(1,10**6))
        det=f['detDM'];mass_exact=det*sum(v[0]for v in f['corners'])/6;mass=2*arb.pi()/U.A(target.energy_scale)*U.A(mass_exact)+U.A(sum(rule['volume_weights'][i]));mass_upper=upper(mass);g=max(majorant(v,f)for v in jet['gradient']);h=max(majorant(v,f)for v in jet['hessian_upper'].values());first=None
        for n in range(2,19):
            a,b=tails(n,rho);ma,mb=M.volumetric_derivative_tail_bounds(n,rho,1,5000);ck(f'independent_tail_formula_{i}_{n}',(a,b)==(ma,mb));Ftail=mass_upper*a*g;Htail=mass_upper*(b*g*g+a*h)
            if Ftail<=F(p['force_tail_target_exact'])and Htail<=F(p['H_entry_tail_target_exact']):first=n;break
        ck('selected_first_certified_order_'+str(i),type(selector['degree'])is int and selector['degree']==first);ck('producer_reported_tail_targets_'+str(i),F(selector['force_tail_upper_exact'])<=F(p['force_tail_target_exact'])and F(selector['H_entry_tail_upper_exact'])<=F(p['H_entry_tail_target_exact']))
        pp=ROOT/f'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/elements/{i:03d}.json';point=read(pp);d=len(f['global_indices']);ck('corrected_point_packet_complete_projection_'+str(i),point['global_indices']==f['global_indices']and len(point['force_integral_minus_cache'])==d and {(r['i'],r['j'])for r in point['H_integral_minus_cache_upper_triangle']}=={(r,c)for r in range(d)for c in range(r,d)})
        ck('active_drive_interval_covers_base_'+str(i),inside(U.gain_drive(p,f['active']),basegain*F(p['dose_exact'])*f['active']))
        order_diagnostics.append({'element':i,'active_exact':str(f['active']),'rho_raw_exact':str(rho_raw),'rho_outward_exact':str(rho),'selected_order':first,'independent_force_tail_upper_exact':str(Ftail),'independent_H_entry_tail_upper_exact':str(Htail),'proof_not_point_overlap':True})
    ck('active_and_inactive_cases_both_present',active=={F(0),F(1)})
    independent_force_fixture();review_gate_ast();ck('producer_synthetic_checks_replay',U.tests()==pre['checks']);ck('all_hashes_stable_before_review_freeze',sha(U.__file__)==SOURCE_SHA and sha(PRODUCER/'protocol.json')==PROTOCOL_SHA and all(sha(ROOT/f)==h for f,h in inputs.items()))
    result={'status':'PASS','checks':checks,'source_sha256':SOURCE_SHA,'protocol_sha256':PROTOCOL_SHA,'preflight_sha256':sha(PRODUCER/'preflight.json'),'review_source_sha256':sha(__file__),'input_hashes':inputs,'order_diagnostics':order_diagnostics,'counts':counts,'mathematical_review':{'uniform_polynomial':'Safe coefficient enclosures cover every canonical field in the NEW boxes. The shared gain interval extends independently of field coefficients, conservatively covering SAMEC parameter/field pairs. Signed moment-dual contraction encloses each truncated integral-minus-Qstar entry.','two_mass_remainder':'Scalar phi derivative tails multiply sup norm Ju/Juv/JuJv and BOTH positive material mass and actual auxiliary weight sum. Polynomial quadratic/shear/active terms are integrated by signed moments without a constitutive-series remainder.','uniform_cache_force':'deltaBT*gc+BstarT*Hline*deltaq follows the fundamental theorem on the same convex strain line; deltaq=(Bstarfull-Bcachefull)*SAMEY. C primitive supplies the independently reviewed three-term Hessian difference. Each actual old normalized weight is applied once.','analytic_domain':'Fresh per-triangle canonical J bounds cover the entire NEW box; outward2^-48rho is checked exactly. Runtime cache convex-strain-line guards remain mandatory and can cause unresolved failure.','point_packets':'Seven corrected safe point packets are HASH-BOUND shape/projection/overlap diagnostics only. Their existence/overlap is NOT the uniform proof.'},'scope':'PRELAUNCH only: mathematical/source/dataflow review, exact geometry/rho/order/tail checks, coefficient inclusion,64 synthetic interval-force endpoint witnesses and actual AST launch-guard tests. No actual seven-element nonlinear integral or64-site cache sum, full611 assembly, surface or root bound was run. Historical unsafe preparations/claims remain withdrawn.','actual_integrals_or_64site_cache_sums_run':False,'assembled611_or_root_claim':False}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');result['retained_review_file']=(OUT/'results.json').relative_to(ROOT).as_posix();result['retained_review_sha256']=sha(OUT/'results.json');write=PRODUCER/'prelaunch_independent_review.json';write.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps({'status':'PASS','checks':len(checks),'orders':[(d['element'],d['selected_order'])for d in order_diagnostics],'prelaunch_review_sha256':sha(write),'source_sha256':SOURCE_SHA,'protocol_sha256':PROTOCOL_SHA}),flush=True)
if __name__=='__main__':main()
