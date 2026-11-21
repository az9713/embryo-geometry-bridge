"""SAFE V2 independent six-edge POINT replay. prepare freezes before run integrates.
Raw P2 maps and generic polynomial jets are implemented here. Reviewed scalar
constitutive formulas, rational moments and binomial tails are explicit reuse.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import datetime,hashlib,importlib.util,json,sys,time,traceback
import numpy as np
from flint import arb,ctx
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target_safe as T
import exact_fe_moments as M
import algebraic_root_feasibility as C
import exact_surface_bridge_safe as P
PREP=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2'
A=C.A
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise RuntimeError('Preserve '+str(p))
    p.write_text(json.dumps(d,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def selfchecks():
    ctx.prec=640
    checks={}
    def check(k,v):
        if v is not True:raise AssertionError(k)
        checks[k]=True
    check('safe_target_expected_hash',sha(T.__file__)=='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954')
    check('safe_surface_library_expected_hash',sha(P.__file__)=='3f8627e87b2fcf95a7109aecf10e1780a534292d2ce53640dea8ae8cd9e621b3')
    check('safe_helper_expected_hash',sha(P.JET.__file__)=='1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373')
    review=read(ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json')
    check('safe_kernel_independent_review_pass',review['status']=='PASS' and bool(review['checks']) and all(v is True for v in review['checks'].values()))
    for exponent in (1,80,300):
        radius=A(Q(1,10**exponent))
        ball=arb(0,C.upper(radius))
        check(f'zero_crossing_{exponent}_not_exact_zero',not ball.is_zero())
        check(f'zero_crossing_{exponent}_legacy_ne_false',(ball!=0) is False)
        poly=T.constant(ball)
        check(f'zero_crossing_{exponent}_retained',(0,0) in poly and not poly[(0,0)].is_zero())
        moved=T.mul(poly,{(0,0):Q(1),(1,0):Q(2)})
        check(f'zero_crossing_{exponent}_full_product_support',set(moved)=={(0,0),(1,0)})
        check(f'zero_crossing_{exponent}_evaluated_enclosure',T.evaluate(moved,Q(1,2),Q(0)).contains(2*radius) and T.evaluate(moved,Q(1,2),Q(0)).contains(-2*radius))
        near=T.sub(T.constant(arb(1,C.upper(radius))),T.ONE)
        check(f'near_cancelled_point_{exponent}_retained',(0,0) in near and not near[(0,0)].is_zero())
        jet=jet_variable(near,[T.ONE])
        square=jmul(jet,jet)
        check(f'near_cancelled_jet_gradient_{exponent}_retained',bool(square[1][0]) and not square[1][0][(0,0)].is_zero())
        value=horner({(1,0):Q(1)},[ball,Q(1)])
        check(f'horner_zero_crossing_constant_{exponent}_retained',(0,0) in value and (1,0) in value)
    check('exact_Arb_zero_removed',T.constant(arb(0))=={})
    check('exact_fraction_zero_removed',T.constant(Q(0))=={})
    check('rational_near_cancellation_exact',T.sub(T.constant(Q(1)+Q(1,10**100)),T.ONE)=={(0,0):Q(1,10**100)})
    check('no_unsafe_target_module_imported','exact_fe_target' not in sys.modules)
    check('no_unsafe_surface_module_imported','exact_surface_bridge_pilot' not in sys.modules)
    reused=P.tests()
    check('reviewed_safe_surface_synthetic_checks_pass',bool(reused) and all(v is True for v in reused.values()))
    packet={'status':'PASS_SAFE_V2_REPLAY_SELFCHECKS_PREPARATION_ONLY','source_sha256':sha(__file__),'checks':checks,'explicitly_reused_safe_surface_synthetic_checks':reused,'actual_integrals_run':False,'created_utc':stamp()}
    write(OUT/'selfchecks.json',packet)
    print(json.dumps({'status':packet['status'],'own_checks':len(checks),'reused_checks':len(reused),'source_sha256':sha(__file__),'receipt_sha256':sha(OUT/'selfchecks.json')}),flush=True)


def prepare(expected_producer_hash):
    proto=read(PREP/'protocol.json');manifest=read(PREP/'results.json');bindings={Path(__file__).relative_to(ROOT).as_posix():sha(__file__)}
    assert sha(PREP/'results.json')==expected_producer_hash
    assert manifest['status']=='SAFE_V2_ALL24_POINT_SURFACE_ASSEMBLED_UNREVIEWED' and len(manifest['edges'])==24
    assert manifest['protocol_sha256']==sha(PREP/'protocol.json') and manifest['all_frozen_bindings_unchanged'] is True
    assert sha(ROOT/'scripts/exact_surface_bridge_assembly_safe.py')=='1012f6f679d48661890310daa0a4d05064a380930a55843e22a88995b0680fb3'
    assert sha(PREP/'protocol.json')=='55e0df292e44081d58863dc8c17b1b568efd6e1ece7177f38da203de6f0d24e7'
    for k,h in proto['source_and_input_hashes'].items():assert sha(ROOT/k)==h;bindings[k]=h
    for path in [PREP/'protocol.json',PREP/'preflight.json',PREP/'prelaunch_independent_review.json',PREP/'results.json',OUT/'selfchecks.json',ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json']:
        bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    selfcheck=read(OUT/'selfchecks.json');assert selfcheck['source_sha256']==sha(__file__) and all(v is True for v in selfcheck['checks'].values())
    selected=[c for c in proto['cases'] if c['edge'] in (0,1,11)]
    assert [(c['face'],c['edge']) for c in selected]==[(f,e) for f in (0,1) for e in (0,1,11)]
    for item in manifest['edges']:
        path=ROOT/item['file'];assert sha(path)==item['sha256'];bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    packet={'version':2,'source_sha256':sha(__file__),'created_utc':stamp(),'input_hashes':bindings,'precision_bits':640,'truncation_order':28,'cases':selected,'point_cache_file':proto['point_cache_file'],'producer_results_sha256':expected_producer_hash,'domain_bounds_reuse':'Reviewed safe all24 preparation; independent raw q/B and whole-edge domain verification at freeze. POINT domains only.','independence':'Own raw P2 values/bases, own generic polynomial jet AD, Horner series and full-product monomial contractions; SAFE T polynomial arithmetic/M moments/core constitutive primitives and SAFE P binomial tails explicitly reused. Independent normalized-cache correction. Shared safe kernel is not an independent arithmetic implementation.','scope':'Only six fixed B point edges; all local force/H upper enclosures compared. No fullfield, parameter/tube/root/PDE claim; no restoration of old unsafe metadata.','no_new_equilibrium':True}
    target=T.Target(read(ROOT/packet['point_cache_file']))
    for case in selected:
        raw=raw_edge(target,case['face'],case['edge']);field=target.surface(case['face'],case['edge'])
        assert raw['q']==field['q'] and raw['B']==field['B'] and raw['global_indices']==case['global_indices']
        assert P.positive_G_bounds(field)==case['domain']
    write(OUT/'protocol.json',packet)
    print(json.dumps({'status':'PREPARED_SAFE_V2_ONLY','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'inputs':len(bindings)}),flush=True)

def coeffs(cs):return T.clean({(i,0):v for i,v in enumerate(cs)})
def physical(cs,lo,L):
    p={}
    for n,c in enumerate(cs):
        for k in range(n+1):p[k,0]=p.get((k,0),Q(0))+c*comb(n,k)*(-lo)**(n-k)/L**n
    return T.clean(p)
def deriv(c):return [(i+1)*x for i,x in enumerate(c[1:])]
def raw_edge(target,face,edge):
    a=target.a;ids=list(map(int,a[f'surface_{face}_edges'][edge]));idx={int(k):i for i,k in enumerate(a['free'])};Y=[[T.rational(x)for x in r]for r in a['reference']]
    for i,k in enumerate(a['free']):Y[int(k)//2][int(k)%2]+=T.rational(a['u'][i])
    phi=[[Q(1),Q(-3),Q(2)],[Q(0),Q(4),Q(-4)],[Q(0),Q(-1),Q(2)]];lo=T.rational(a['reference'][ids[0],0]);hi=T.rational(a['reference'][ids[2],0]);L=hi-lo
    assert L>0 and lo>=0
    R=[sum((phi[j][i]*Y[node][0]for j,node in enumerate(ids)),Q(0))for i in range(3)];Z=[sum((phi[j][i]*Y[node][1]for j,node in enumerate(ids)),Q(0))for i in range(3)]
    radius=physical(R,lo,L)
    if lo==0:
        assert R[0]==0;hoop=physical([x/L for x in R[1:]],lo,L)
    else:hoop={(i-1,j):v for (i,j),v in radius.items()}
    q=[hoop,physical([x/L for x in deriv(R)],lo,L),physical([x/L for x in deriv(Z)],lo,L)]
    dofs=[2*n+c for n in ids for c in[0,1]];free=[i for i,k in enumerate(dofs)if k in idx];B=[[],[],[]]
    for k in free:
        node,component=divmod(k,2);p=phi[node]
        if component==0:
            h=physical([x/L for x in p[1:]],lo,L)if lo==0 else{(i-1,j):v for (i,j),v in physical(p,lo,L).items()}
            B[0].append(h);B[1].append(physical([x/L for x in deriv(p)],lo,L));B[2].append({})
        else:B[0].append({});B[1].append({});B[2].append(physical([x/L for x in deriv(p)],lo,L))
    return {'q':q,'B':B,'ids':ids,'Y':Y,'lo':lo,'hi':hi,'L':L,'free':free,'global_indices':[idx[dofs[k]]for k in free]}

def jet_variable(value,grad):return(value,grad,{(i,j):{}for i in range(len(grad))for j in range(i,len(grad))})
def jadd(x,y):return(T.add(x[0],y[0]),[T.add(a,b)for a,b in zip(x[1],y[1])],{k:T.add(x[2][k],y[2][k])for k in x[2]})
def jmul(x,y):
    n=len(x[1]);g=[T.add(T.mul(x[1][i],y[0]),T.mul(x[0],y[1][i]))for i in range(n)];H={}
    for i in range(n):
        for j in range(i,n):H[i,j]=T.add(T.add(T.mul(x[2][i,j],y[0]),T.mul(x[0],y[2][i,j])),T.add(T.mul(x[1][i],y[1][j]),T.mul(x[1][j],y[1][i])))
    return(T.mul(x[0],y[0]),g,H)
def horner(e,coefficients):
    p={}
    for c in reversed(coefficients):p=T.add(T.mul(p,e),T.constant(c))
    return p
def norm(p,lo,hi):
    assert all(j==0 for i,j in p) and all(isinstance(v,Q) for v in p.values())
    raw=sum((abs(v)*(hi**i if i>=0 else lo**i)for (i,j),v in p.items()),Q(0));d=max(0,-min((i for i,j in p),default=0));n={i+d:v for (i,j),v in p.items()};c=(lo+hi)/2;h=(hi-lo)/2
    shifted=[sum((v*comb(i,k)*c**(i-k)*h**k for i,v in n.items()if i>=k),Q(0))for k in range(max(n,default=0)+1)]
    return min(raw,sum(map(abs,shifted),Q(0))/(lo**d if d else 1))
def overlap(a,b):return not(a.upper()<b.lower()or b.upper()<a.lower())

def replay_case(target,case,N):
    face,edge=case['face'],case['edge'];raw=raw_edge(target,face,edge);n=len(raw['free']);q=raw['q'];B=raw['B'];j=[jet_variable(q[i],B[i])for i in range(3)];G,Gi,Gij=jmul(jmul(j[0],j[0]),jadd(jmul(j[1],j[1]),jmul(j[2],j[2])))
    q0=Q(case['domain']['q0_exact']);rho=Q(case['domain']['rho_exact']);e=T.scale(T.sub(G,T.constant(q0)),1/q0);par=target.packet['parameters'];K=T.rational(par['basal_area_modulus'if face==0 else'apical_area_modulus']);tension=T.rational(par['basal_surface_tension'if face==0 else'apical_surface_tension']);root=A(q0).sqrt()
    AA=T.add(T.constant(A(K/2)),T.scale(horner(e,P.binomial_coefficients(Q(-1,2),N)),A((tension-K)/2)/root));BB=T.scale(horner(e,P.binomial_coefficients(Q(-3,2),N)),-A((tension-K)/4)/(A(q0)*root))
    a=target.a;samples=[T.rational(x[2])+T.rational(x[1])/2 for x in a[f'surface_{face}_N']];weights=[T.rational(x)/target.energy_scale for x in a[f'surface_{face}_weights'][edge]];rs=[raw['lo']+raw['L']*s for s in samples];memo={}
    def moment(k):
        if k not in memo:memo[k]=target.normalization()*M.surface_moment(raw['lo'],raw['hi'],k).ball(ctx.prec)-A(sum((w*r**k for w,r in zip(weights,rs)),Q(0)))
        return memo[k]
    def delta(poly):
        assert all(j==0 for i,j in poly)
        return sum((C.A(c)*moment(i)if not isinstance(c,arb)else c*moment(i)for (i,j),c in poly.items()),arb(0))
    mass=C.upper(target.normalization()*M.surface_moment(raw['lo'],raw['hi'],0).ball(ctx.prec)+A(sum(weights,Q(0))));ta=A(abs(tension-K)/2*P.binomial_tail(Q(-1,2),N,rho))/root;tb=A(abs(tension-K)/4*P.binomial_tail(Q(-3,2),N,rho))/(A(q0)*root)
    F=[delta(T.mul(AA,p))for p in Gi];H={(i,k):delta(T.add(T.mul(AA,p),T.mul(BB,T.mul(Gi[i],Gi[k]))))for(i,k),p in Gij.items()};nf=[norm(p,raw['lo'],raw['hi'])for p in Gi];tf=[mass*ta*A(v)for v in nf];th={(i,k):mass*(ta*A(norm(p,raw['lo'],raw['hi']))+tb*A(nf[i]*nf[k]))for(i,k),p in Gij.items()}
    # Exact old-cache coordinate normalization gives SAME scalarJs at scale1.
    CF=[arb(0)for _ in range(n)];CH={(i,k):arb(0)for i in range(n)for k in range(i,n)};phi=a[f'surface_{face}_N'];dphi=a[f'surface_{face}_dN'];r0=a[f'surface_{face}_r0'];l0=a[f'surface_{face}_l0']
    for iq,(r,w)in enumerate(zip(rs,weights)):
        qs=[T.evaluate(p,r,0)for p in q];bs=[[T.evaluate(p,r,0)for p in row]for row in B];gs,hs,_=C.surface_local(list(map(A,qs)),A(K),A(tension),A(1))
        qc=[sum((T.rational(phi[iq,j])*raw['Y'][node][0]for j,node in enumerate(raw['ids'])),Q(0)),sum((T.rational(dphi[iq,j])*raw['Y'][node][0]for j,node in enumerate(raw['ids'])),Q(0)),sum((T.rational(dphi[iq,j])*raw['Y'][node][1]for j,node in enumerate(raw['ids'])),Q(0))]
        D=[raw['L']/(T.rational(r0[edge,iq])*T.rational(l0[edge])),1/raw['L'],1/raw['L']];qc=[v*d for v,d in zip(qc,D)];bc=[[],[],[]]
        for k in raw['free']:
            j,comp=divmod(k,2);bc[0].append(D[0]*T.rational(phi[iq,j])if comp==0 else Q(0));bc[1].append(D[1]*T.rational(dphi[iq,j])if comp==0 else Q(0));bc[2].append(D[2]*T.rational(dphi[iq,j])if comp==1 else Q(0))
        gc,hc,_=C.surface_local(list(map(A,qc)),A(K),A(tension),A(1))
        for i in range(n):
            CF[i]+=A(w)*sum((A(bs[k][i])*gs[k]-A(bc[k][i])*gc[k]for k in range(3)),arb(0))
            for k in range(i,n):CH[i,k]+=A(w)*sum((A(bs[v][i])*hs[v][z]*A(bs[z][k])-A(bc[v][i])*hc[v][z]*A(bc[z][k])for v in range(3)for z in range(3)),arb(0))
    inflate=lambda x,t:x+arb(0,C.upper(t))
    F=[inflate(v,tf[i])+CF[i]for i,v in enumerate(F)];H={k:inflate(v,th[k])+CH[k]for k,v in H.items()}
    prod=read(PREP/'edges'/f'face{face}_edge{edge:02d}.json');unpack=lambda x:A(Q(x['lower_exact'])).union(A(Q(x['upper_exact'])));checks={'same_local_global_map':raw['global_indices']==prod['global_indices']}
    for i,v in enumerate(F):checks[f'force_interval_overlap_{i}']=overlap(v,unpack(prod['force_integral_minus_cache'][i]))
    ph={(x['i'],x['j']):x['enclosure']for x in prod['H_integral_minus_cache_upper']}
    for k,v in H.items():checks[f'Hessian_interval_overlap_{k[0]}_{k[1]}']=overlap(v,unpack(ph[k]))
    assert all(checks.values()),checks
    return {'status':'PASS_SAFE_V2_INDEPENDENT_POINT_SURFACE_EDGE_REPLAY','face':face,'edge':edge,'order':N,'precision_bits':ctx.prec,'checks':checks,'global_indices':raw['global_indices'],'force_enclosures':[C.pack(v)for v in F],'Hessian_enclosures':[{'i':i,'j':j,'enclosure':C.pack(v)}for(i,j),v in H.items()],'producer_file_sha256':sha(PREP/'edges'/f'face{face}_edge{edge:02d}.json'),'scope':'Point edge only, all force/H intervals overlap independentN28/640 enclosure; no actual-root or tube claim.'}

def run():
    proto=read(OUT/'protocol.json');assert proto['source_sha256']==sha(__file__);assert all(sha(ROOT/k)==h for k,h in proto['input_hashes'].items());ctx.prec=proto['precision_bits'];p=proto;target=T.Target(read(ROOT/proto['point_cache_file']));write(OUT/'job.json',{'status':'RUNNING','started_utc':stamp(),'protocol_sha256':sha(OUT/'protocol.json')});outputs=[];start=time.monotonic()
    try:
        for case in p['cases']:
            r=replay_case(target,case,proto['truncation_order']);f=OUT/f"face{case['face']}_edge{case['edge']:02d}.json";write(f,r);outputs.append({'file':f.name,'sha256':sha(f),'checks':len(r['checks'])});print(json.dumps({'edge_complete':[case['face'],case['edge']],'checks':len(r['checks'])}),flush=True)
        assert all(sha(ROOT/k)==h for k,h in proto['input_hashes'].items());write(OUT/'results.json',{'status':'PASS_SIX_SAFE_V2_INDEPENDENT_POINT_SURFACE_REPLAYS','outputs':outputs,'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),'all_bindings_current':True,'elapsed_seconds':time.monotonic()-start,'root_transfer_proved':False})
        (OUT/'job.json').write_text(json.dumps({'status':'FINISHED','exit_code':0,'finished_utc':stamp(),'elapsed_seconds':time.monotonic()-start},indent=2)+'\n',encoding='utf-8')
    except Exception as e:
        (OUT/'job.json').write_text(json.dumps({'status':'FAILED','error':repr(e),'traceback':traceback.format_exc(),'elapsed_seconds':time.monotonic()-start},indent=2)+'\n',encoding='utf-8');raise

if __name__=='__main__':
    if sys.argv[1:]==['selfchecks']:selfchecks()
    elif len(sys.argv)==3 and sys.argv[1]=='prepare':prepare(sys.argv[2])
    elif sys.argv[1:]==['run']:run()
    else:raise SystemExit('positional selfchecks | prepare EXPECTED_PRODUCER_SHA | run required')
