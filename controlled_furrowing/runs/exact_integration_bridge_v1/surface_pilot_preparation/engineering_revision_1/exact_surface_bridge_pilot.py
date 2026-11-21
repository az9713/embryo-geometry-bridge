"""Point surface exact-integral-minus-cache pilot: preparation first.
The engine exists, but --prepare never integrates an actual field. No old
source/cache is altered, and no uniform root or continuum result is implied.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import argparse,datetime,hashlib,json,sys,time
import numpy as np
from flint import arb,arb_mat,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runs/exact_integration_bridge_v1/bridge_design_probe'))
import exact_fe_target as T
import exact_fe_moments as M
import algebraic_root_feasibility as core
import polynomial_jet as JET
OUT=ROOT/'runs/exact_integration_bridge_v1/surface_pilot_preparation'
PREP='runs/exact_integration_bridge_v1/target_preparation/protocol.json'
A=core.A

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output already exists: '+str(p))
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def exact(v):return T.rational(v)
def merge(bindings,k,v):
    T.bind(bindings,k,v)
    if sha(ROOT/k)!=v:raise ValueError('changed immutable input: '+k)
def assert_bindings(p):
    for k,h in p['source_and_input_hashes'].items():
        if sha(ROOT/k)!=h:raise ValueError('changed frozen pilot binding: '+k)
def interval_poly(poly,lo,hi):
    """Exact rational enclosure after centering a univariate polynomial."""
    if lo>hi or any(i<0 or j!=0 for i,j in poly):raise ValueError('ordered polynomial chart interval required')
    c=(lo+hi)/2;h=(hi-lo)/2
    shifted=T.compose(poly,{(0,0):c,(1,0):h},{})
    mid=shifted.get((0,0),F(0));radius=sum(abs(v)for k,v in shifted.items()if k!=(0,0))
    return mid-radius,mid+radius

def idiv(a,b):
    if b[0]<=0<=b[1]:raise ValueError('interval denominator contains zero')
    vals=[x/y for x in a for y in b];return min(vals),max(vals)
def isquare(a):return (F(0)if a[0]<=0<=a[1]else min(x*x for x in a),max(x*x for x in a))
def positive_G_bounds(fields,max_depth=8,rho_goal=F(1,8)):
    """Whole CLOSED canonical edge, including regularized axis endpoint.
    Bounds are for the frozen point field, not a nodal/parameter tube.
    """
    if not isinstance(max_depth,int)or max_depth<0 or max_depth>16 or not 0<rho_goal<1:raise ValueError('valid finite subdivision policy required')
    history=[];selected=None
    for depth in range(max_depth+1):
        leaves=[]
        for k in range(2**depth):
            lo,hi=F(k,2**depth),F(k+1,2**depth)
            q=[interval_poly(p,lo,hi)for p in fields['q_numerators']]
            den=interval_poly(fields['hoop_denominator'],lo,hi);q[0]=idiv(q[0],den)
            aa,bb,cc=map(isquare,q);G=(aa[0]*(bb[0]+cc[0]),aa[1]*(bb[1]+cc[1]))
            leaves.append({'s_exact':[str(lo),str(hi)],'q_bounds_exact':[[str(x)for x in p]for p in q],'G_bounds_exact':list(map(str,G))})
        lower=min(F(x['G_bounds_exact'][0])for x in leaves);upper=max(F(x['G_bounds_exact'][1])for x in leaves);hoop_lower=min(F(x['q_bounds_exact'][0][0])for x in leaves)
        q0=(lower+upper)/2;rho=(upper-lower)/(upper+lower)if lower>0 else None
        history.append({'depth':depth,'leaf_count':len(leaves),'G_lower_exact':str(lower),'G_upper_exact':str(upper),'rho_exact':str(rho)if rho is not None else None})
        selected={'G_lower_exact':str(lower),'G_upper_exact':str(upper),'hoop_lower_exact':str(hoop_lower),'q0_exact':str(q0),'rho_exact':str(rho)if rho is not None else None,'depth':depth,'leaves':leaves,'history':history,'rho_goal_exact':str(rho_goal),'rho_goal_achieved':rho is not None and rho<=rho_goal}
        if lower>0 and hoop_lower>0 and rho<=rho_goal:break
    if lower<=0 or hoop_lower<=0 or q0<=0 or not 0<=rho<1:raise ValueError('positive G and convergent relative expansion not certified')
    selected['status']='CERTIFIED_POSITIVE_POINT_G_ON_FULL_EDGE'
    return selected

def binomial_coefficients(alpha,order):
    alpha=exact(alpha)
    if alpha not in(F(-1,2),F(-3,2))or not isinstance(order,int)or isinstance(order,bool)or order<0:raise ValueError('inverse half/threehalf and nonnegative integer order required')
    c=[F(1)]
    for n in range(order):c.append(c[-1]*(alpha-n)/F(n+1))
    return c

def binomial_tail(alpha,order,rho):
    """Exact absolute tail on |x|<=rho<1; no sampled error substitution."""
    alpha=exact(alpha);rho=exact(rho);c=binomial_coefficients(alpha,order+1)
    if not 0<=rho<1:raise ValueError('rho must be in [0,1)')
    if rho==0:return F(0)
    m=order+1
    if alpha==F(-1,2):return abs(c[m])*rho**m/(1-rho)
    # |binom(-3/2,n)|=(2n+1)|binom(-1/2,n)| <=2n+1.
    universal=rho**m*((2*m+1)/(1-rho)+2*rho/(1-rho)**2)
    ratio=rho*F(m*2+3,2*(m+1))
    return min(universal,abs(c[m])*rho**m/(1-ratio))if ratio<1 else universal

def normalized_series(e,q0,order,K,tension):
    """A=K/2+(T-K)/(2sqrtG), B=-(T-K)/(4G**1.5)."""
    q0=exact(q0)
    if q0<=0:raise ValueError('positive expansion center required')
    c1=binomial_coefficients(F(-1,2),order);c3=binomial_coefficients(F(-3,2),order)
    root=A(q0).sqrt();pref1=A((tension-K)/2)/root;pref3=-A((tension-K)/4)/(A(q0)*root)
    p=T.ONE;outA=T.constant(A(K/2));outB={};e={k:A(v)for k,v in e.items()}
    for n in range(order+1):
        outA=T.add(outA,T.scale(p,pref1*A(c1[n])));outB=T.add(outB,T.scale(p,pref3*A(c3[n])))
        if n<order:p=T.mul(p,e)
    return outA,outB

def pnorm(poly,lo,hi):
    if not 0<=lo<hi or any(j!=0 or (i<0 and lo==0)for i,j in poly):raise ValueError('regular physical R interval required')
    raw=sum(abs(v)*(hi**i if i>=0 else lo**i)for (i,j),v in poly.items())
    den=max(0,-min((i for i,j in poly),default=0));num=T.mul(poly,T.power(T.X,den))
    local=T.compose(num,{(0,0):(lo+hi)/2,(1,0):(hi-lo)/2},{})
    shifted=sum(map(abs,local.values()))/(lo**den if den else 1)
    return min(raw,shifted)

class SurfaceFunctional:
    def __init__(self,target,fields,face,edge,rule):
        self.target=target;self.fields=fields;self.points=[((fields['Rlo']+fields['L']*s,F(0)),w)for s,w in zip(rule['surface_s'][face],rule['surface_weights'][face][edge])]
        self.mass_integral=self.integral(0,0);self.mass_rule=A(sum(w for pt,w in self.points))
        if not self.mass_integral>0 or not self.mass_rule>0:raise ValueError('strict positive material and auxiliary measures required')
    @lru_cache(maxsize=None)
    def integral(self,a,b):
        if b!=0:raise ValueError('surface physical-R polynomial required')
        return self.target.normalization()*M.surface_moment(self.fields['Rlo'],self.fields['Rhi'],a).ball(ctx.prec)
    @lru_cache(maxsize=None)
    def defect(self,a,b):return self.integral(a,b)-A(M.auxiliary_moment(self.points,a,b))

def cache_correction(target,fields,face,edge,rule,K,tension):
    """Finite Qstar minus old binary-cache, SAME exact Xbinary+free-u nodes."""
    cache=target.a;ids=cache[f'surface_{face}_edges'][edge];Y=target.nodes();coords=[A(Y[int(j)][c])for j in ids for c in(0,1)]
    n=len(fields['free_local']);outg=[arb(0)for _ in range(n)];outH=arb_mat(n,n)
    r0=cache[f'surface_{face}_r0'];l0=cache[f'surface_{face}_l0'];Ns=cache[f'surface_{face}_N'];dNs=cache[f'surface_{face}_dN']
    for iq,(s,w)in enumerate(zip(rule['surface_s'][face],rule['surface_weights'][face][edge])):
        qs,Bs=T.evaluate_surface(fields,s);bs=arb_mat(3,n,[A(v)for row in Bs for v in row]);g,h,_=core.surface_local(list(map(A,qs)),A(K),A(tension),A(1))
        gf=bs.transpose()*arb_mat(3,1,g);HH=bs.transpose()*arb_mat(3,3,[v for row in h for v in row])*bs
        B=[[arb(0)for _ in range(6)]for _ in range(3)]
        for j in range(3):B[0][2*j]=A(float(Ns[iq,j]));B[1][2*j]=A(float(dNs[iq,j]));B[2][2*j+1]=A(float(dNs[iq,j]))
        qc=[sum((v*x for v,x in zip(row,coords)),arb(0))for row in B];bc=arb_mat(3,n,[row[k]for row in B for k in fields['free_local']])
        if not A(float(r0[edge,iq]))>0 or not A(float(l0[edge]))>0:raise ValueError('positive old cached r0/l0 required')
        scale=1/(A(float(r0[edge,iq]))*A(float(l0[edge])))
        gc,hc,_=core.surface_local(qc,A(K),A(tension),scale)
        gf-=bc.transpose()*arb_mat(3,1,gc);HH-=bc.transpose()*arb_mat(3,3,[v for row in hc for v in row])*bc
        for i in range(n):
            outg[i]+=A(w)*gf[i,0]
            for j in range(n):outH[i,j]+=A(w)*HH[i,j]
    return outg,outH

def run_edge(target,face,edge,rule,domain,order):
    """Actual integration engine, callable ONLY after separate root approval."""
    fields=target.surface(face,edge);jet=JET.surface_squared_jacobian_jet(fields['q_physical_R'],fields['B_physical_R']);n=jet['local_dofs']
    K=exact(target.packet['parameters']['basal_area_modulus'if face==0 else'apical_area_modulus']);tension=exact(target.packet['parameters']['basal_surface_tension'if face==0 else'apical_surface_tension'])
    q0=F(domain['q0_exact']);rho=F(domain['rho_exact']);e=T.scale(T.sub(jet['value'],T.constant(q0)),1/q0)
    Apoly,Bpoly=normalized_series(e,q0,order,K,tension);functional=SurfaceFunctional(target,fields,face,edge,rule)
    polynomial=JET.contract_volumetric_defect(jet,Apoly,Bpoly,functional.defect)
    mass=core.upper(functional.mass_integral+functional.mass_rule);sqrtq0=A(q0).sqrt()
    tailA=A(abs(tension-K)/2*binomial_tail(F(-1,2),order,rho))/sqrtq0
    tailB=A(abs(tension-K)/4*binomial_tail(F(-3,2),order,rho))/(A(q0)*sqrtq0)
    norms=[pnorm(p,fields['Rlo'],fields['Rhi'])for p in jet['gradient']]
    corrF,corrH=cache_correction(target,fields,face,edge,rule,K,tension)
    Ft=[mass*tailA*A(norms[i])for i in range(n)];Ht={(i,j):mass*(tailA*A(pnorm(p,fields['Rlo'],fields['Rhi']))+tailB*A(norms[i]*norms[j]))for (i,j),p in jet['hessian_upper'].items()}
    inflate=lambda v,t:v+(-core.upper(t)).union(core.upper(t))
    Fout=[inflate(v,Ft[i])+corrF[i]for i,v in enumerate(polynomial['gradient'])];Hout={k:inflate(v,Ht[k])+corrH[k[0],k[1]]for k,v in polynomial['hessian_upper'].items()}
    return {'status':'BOUNDED_POINT_SURFACE_EDGE_UNREVIEWED','face':face,'edge':edge,'global_indices':fields['global_indices'],'q0_exact':str(q0),'rho_exact':str(rho),'order':order,'precision_bits':ctx.prec,'A_support':len(Apoly),'B_support':len(Bpoly),'material_mass':core.pack(functional.mass_integral),'auxiliary_mass':core.pack(functional.mass_rule),'force_integral_minus_cache':[core.pack(v)for v in Fout],'force_auxiliary_minus_cache':[core.pack(v)for v in corrF],'force_two_mass_tail_upper':[core.pack(v)for v in Ft],'H_integral_minus_cache_upper':[{'i':i,'j':j,'enclosure':core.pack(v),'auxiliary_minus_cache':core.pack(corrH[i,j]),'two_mass_tail_upper':core.pack(Ht[i,j])}for (i,j),v in Hout.items()],'unique_moment_defects':functional.defect.cache_info().currsize,'scope':'ONE frozen point edge; no full surface/field sum, parameter/root tube or exact-integral root transfer.'}


def tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    for n in range(33):
        c1=binomial_coefficients(F(-1,2),n)[n];c3=binomial_coefficients(F(-3,2),n)[n]
        ck('exact_half_binomial_'+str(n),c1==F((-1)**n*comb(2*n,n),4**n))
        ck('exact_threehalf_binomial_'+str(n),c3==(2*n+1)*c1)
    for alpha in(F(-1,2),F(-3,2)):
        for N in(0,2,8):
            rho=F(1,5);tail=binomial_tail(alpha,N,rho)
            for x in(-rho,rho):
                P=sum(c*x**n for n,c in enumerate(binomial_coefficients(alpha,N)))
                actual=1/(1+A(x)).sqrt()
                if alpha==F(-3,2):actual/=1+A(x)
                ck(f'binomial_tail_encloses_endpoint_{alpha}_{N}_{x}',core.upper(actual-A(P))<=A(tail))
    ck('zero_rho_zero_tail',binomial_tail(F(-3,2),8,F(0))==0)
    try:binomial_tail(F(-1,2),8,F(1));ck('rho_one_rejected',False)
    except ValueError:checks['rho_one_rejected']=True
    try:binomial_coefficients(F(1,2),4);ck('unsupported_exponent_rejected',False)
    except ValueError:checks['unsupported_exponent_rejected']=True
    # Closed-form synthetic integral and a positive auxiliary point mass.
    for alpha in(F(-1,2),F(-3,2)):
        rho=F(1,5);N=4;cs=binomial_coefficients(alpha,N);weight=F(2,3);x=rho/3
        exactI=2*((1+A(rho)).sqrt()-1)/A(rho)if alpha==F(-1,2)else 2*(1-1/(1+A(rho)).sqrt())/A(rho)
        exactQ=A(weight)/(1+A(x)).sqrt()
        if alpha==F(-3,2):exactQ/=1+A(x)
        Pdelta=sum(c*rho**n/F(n+1)for n,c in enumerate(cs))-weight*sum(c*x**n for n,c in enumerate(cs))
        ck('two_mass_synthetic_integral_'+str(alpha),core.upper(exactI-exactQ-A(Pdelta))<=A((1+weight)*binomial_tail(alpha,N,rho)))
    rho=F(1,100);N=2;cs=binomial_coefficients(F(-1,2),N)
    actual=1/(1+A(rho)).sqrt()-1/(1-A(rho)).sqrt();trunc=sum(c*(rho**n-(-rho)**n)for n,c in enumerate(cs));error=core.upper(actual-A(trunc));tail=A(binomial_tail(F(-1,2),N,rho))
    ck('one_mass_remainder_can_fail_for_signed_rule_difference',error>tail and error<=2*tail)
    ck('positive_R_quadrature_toy_nonzero_exact_defect',M.surface_moment(F(1),F(2),2).constant-F(3,2)*F(3,2)**2==F(3,8))
    ck('axis_R_weight_counted_once_toy',M.surface_moment(F(0),F(1),2).constant==F(1,4))
    for lo,hi in[(F(0),F(1)),(F(1),F(3))]:
        ref=[[lo,F(0)],[(lo+hi)/2,F(0)],[hi,F(0)]];f=T.surface_fields(ref,ref);g=positive_G_bounds(f,max_depth=4)
        ck('identity_G_bound_positive_'+str(lo),F(g['G_lower_exact'])>0 and F(g['G_lower_exact'])<=1<=F(g['G_upper_exact'])and F(g['rho_exact'])<1)
    try:positive_G_bounds(T.surface_fields([[F(0),F(0)],[F(1,2),F(0)],[F(1),F(0)]],[[F(0),F(0)]]*3),max_depth=1);ck('collapsed_surface_rejected',False)
    except ValueError:checks['collapsed_surface_rejected']=True
    try:positive_G_bounds(T.surface_fields([[F(0),F(0)],[F(1,2),F(0)],[F(1),F(0)]],[[F(0),F(0)],[F(-1,2),F(0)],[F(-1),F(0)]]),max_depth=1);ck('negative_hoop_even_with_positive_G_rejected',False)
    except ValueError:checks['negative_hoop_even_with_positive_G_rejected']=True
    # Local force/Hessian identity via G=Js² at an exact rational sqrt point.
    q=[F(6,5),F(4,5),F(3,5)];K=F(3);tension=F(7);polys=list(map(T.constant,q));B=[[T.ONE if i==j else{}for j in range(3)]for i in range(3)];jet=JET.surface_squared_jacobian_jet(polys,B)
    root=F(6,5);aa=K/2+(tension-K)/(2*root);bb=-(tension-K)/(4*root**3);gg,HH,_=core.surface_local(list(map(A,q)),A(K),A(tension),A(1))
    for i in range(3):
        Gi=T.evaluate(jet['gradient'][i],0,0);ck('local_surface_gradient_G_identity_'+str(i),gg[i].contains(A(aa*Gi)))
        for j in range(i,3):
            Gj=T.evaluate(jet['gradient'][j],0,0);Gij=T.evaluate(jet['hessian_upper'][i,j],0,0);ck(f'local_surface_H_G_identity_{i}_{j}',HH[i][j].contains(A(aa*Gij+bb*Gi*Gj)))
    # Exercise the actual finite-cache engine on a perfectly matched toy rule.
    from types import SimpleNamespace
    ref=[[F(1),F(0)],[F(2),F(0)],[F(3),F(0)]];f=T.surface_fields(ref,ref);a={'surface_0_edges':np.array([[0,1,2]]),'surface_0_N':np.array([[0.,1.,0.]]),'surface_0_dN':np.array([[-1.,0.,1.]]),'surface_0_r0':np.array([[2.]]),'surface_0_l0':np.array([2.])};toy=SimpleNamespace(a=a,nodes=lambda:ref);rule={'surface_s':[[F(1,2)]],'surface_weights':[[[F(1,3)]]]}
    cg,ch=cache_correction(toy,f,0,0,rule,F(2),F(3));ck('matched_finite_cache_toy_zero_force_H',all(v.contains(0)and core.upper(v)<A(F(1,10**100))for v in cg+list(ch.entries())))
    a['surface_0_r0']=np.array([[2.+2.**-20]]);cg,ch=cache_correction(toy,f,0,0,rule,F(2),F(3));ck('perturbed_binary_reference_cache_not_silently_equal',any(not v.contains(0)for v in cg+list(ch.entries())))
    return checks

def prepare():
    if(OUT/'protocol.json').exists():raise ValueError('preparation is write-once')
    ctx.prec=512;checks=tests();prep=read(ROOT/PREP);row=next(r for r in prep['records']if r['record_id']=='optoGEF2_small_col_B');bindings={}
    for k,h in prep['source_hashes'].items():merge(bindings,k,h)
    for k,h in prep['input_hashes'].items():merge(bindings,k,h)
    paths=[Path(__file__),Path(M.__file__),Path(JET.__file__),ROOT/PREP,ROOT/'runs/exact_integration_bridge_v1/bridge_design_probe/polynomial_jet_tests.json',ROOT/'evidence/exact_fe_moments_independent_review.json',ROOT/'evidence/exact_integration_target_theory.json']
    for path in paths:merge(bindings,path.relative_to(ROOT).as_posix(),sha(path))
    momentreview=read(ROOT/'evidence/exact_fe_moments_independent_review.json')
    if momentreview['status']!='PASS_INDEPENDENT_MOMENT_COMPONENT'or any(v is not True for v in momentreview['checks'].values())or momentreview['source_sha256']!=sha(M.__file__):raise ValueError('current independently reviewed exact moment component required')
    root=read(ROOT/row['root_certificate_file'])
    if sha(ROOT/row['point_cache_file'])!=row['point_cache_sha256']or root['artifact_hashes'][row['point_cache_file']]!=row['point_cache_sha256']:raise ValueError('authoritative B cache binding mismatch')
    target=T.Target(read(ROOT/row['point_cache_file']));cases=[]
    for face in(0,1):
        for edge in(0,1,11):
            fields=target.surface(face,edge);domain=positive_G_bounds(fields,max_depth=8,rho_goal=F(1,8));jet=JET.surface_squared_jacobian_jet(fields['q_physical_R'],fields['B_physical_R'])
            if any(j!=0 for p in[jet['value']]+jet['gradient']+list(jet['hessian_upper'].values())for i,j in p):raise ValueError('univariate physical-R jets required')
            cases.append({'face':face,'name':'basal'if face==0 else'apical','edge':edge,'Rlo_exact':str(fields['Rlo']),'Rhi_exact':str(fields['Rhi']),'global_indices':fields['global_indices'],'domain':domain,'local_dofs':jet['local_dofs'],'point_field_only':True})
    checks['all6_selected_point_surface_G_positive_entire_edges']=all(F(c['domain']['G_lower_exact'])>0 and F(c['domain']['hoop_lower_exact'])>0 and F(c['domain']['rho_exact'])<1 for c in cases)
    checks['selected_old_B_cache_authoritative_hash_bound']=True;checks['no_actual_surface_field_integral_evaluated']=True
    if not all(checks.values()):raise AssertionError(checks)
    p={'version':1,'utc':stamp(),'source_and_input_hashes':bindings,'record':row,'cases':cases,'precision_bits':512,'truncation_order':24,'domain_subdivision_max_depth':8,'domain_rho_goal_exact':'1/8','constitutive_formula':'Es=K/2*G+(T-K)*sqrt(G)+K/2; A=K/2+(T-K)/(2sqrtG), B=-(T-K)/(4G^1.5); F=A*Gi,H=A*Gij+B*Gi*Gj.','expansion':'G=q0*(1+x); rational q0 from certified whole-edge G bounds; exact binomial coefficients, Arb q0 powers.','tail_policy':'Both exact positive material mass and sum of ACTUAL old-normalized positive auxiliary weights multiply derivative tails. Signed polynomial moment defects before absolute values.','finite_correction':'Exact Qstar minus frozen binary-cache surface force/H using SAME exact Xbinary+free-u coordinates and old cached N,dN,r0,l0,weights.','radius_factor':'2*pi/S * integral_R R*f(R)dR; exact surface_moment includes R exactly once.','input_kind':'Frozen B proof-center point ONLY, no parameter/root tube.','run_scope':'SIX SELECTED SURFACE EDGES ONLY; no full surface/field sum or root transfer.','run_requires_prelaunch_review':{'file':'prelaunch_independent_review.json','status':'PASS','bindings':['source_sha256','protocol_sha256','preflight_sha256'],'all_checks_literal_true':True},'actual_field_integration_run':False,'source_time_PDE_or_root_transfer_certified':False,'external_author_contact_sent':False}
    write(OUT/'protocol.json',p);write(OUT/'guard_tests.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'scope':'Exact scalar/normalization/derivative and toy finite-cache guards plus actual point G-domain preparation. No actual nonlinear surface integral.'})
    preflight={'status':'PASS_PREPARATION_ONLY','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'guard_tests_sha256':sha(OUT/'guard_tests.json'),'bound_inputs':len(bindings),'actual_field_integration_run':False,'actual_surface_cache_correction_run':False,'root_transfer_proved':False}
    write(OUT/'preflight.json',preflight);print(json.dumps({'status':preflight['status'],'checks':len(checks),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'cases':[(c['face'],c['edge'],c['domain']['depth'],float(F(c['domain']['rho_exact'])))for c in cases]}),flush=True)

def run():
    protocol=read(OUT/'protocol.json');preflight=read(OUT/'preflight.json');review=read(OUT/'prelaunch_independent_review.json')
    if review['status']!='PASS'or any(v is not True for v in review['checks'].values())or review['source_sha256']!=sha(__file__)or review['protocol_sha256']!=sha(OUT/'protocol.json')or review['preflight_sha256']!=sha(OUT/'preflight.json'):raise ValueError('exact matching independent prelaunch review required')
    assert_bindings(protocol);ctx.prec=protocol['precision_bits'];target=T.Target(read(ROOT/protocol['record']['point_cache_file']));rule=target.auxiliary_rule();folder=OUT/'integration_run';write(folder/'started.json',{'utc':stamp(),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':sha(OUT/'prelaunch_independent_review.json')});outputs=[]
    for case in protocol['cases']:
        face,edge=case['face'],case['edge'];result=run_edge(target,face,edge,rule,case['domain'],protocol['truncation_order']);name=f'face{face}_edge{edge:02d}.json';write(folder/name,result);outputs.append({'file':name,'sha256':sha(folder/name)});print(json.dumps({'phase':'edge_complete','face':face,'edge':edge}),flush=True)
    assert_bindings(protocol);write(folder/'results.json',{'status':'BOUNDED_POINT_SURFACE_PILOT_COMPLETE_UNREVIEWED','edges':outputs,'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),'all_frozen_bindings_unchanged':True,'root_transfer_proved':False,'scope':protocol['run_scope']})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run','test']);args=ap.parse_args();ctx.prec=512
    if args.command=='prepare':prepare()
    elif args.command=='run':run()
    else:print(json.dumps({'status':'PASS','checks':tests()}))
