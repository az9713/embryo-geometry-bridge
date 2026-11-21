"""Small symbolic/source tests only; no 611D force/Hessian/integral assembly."""
from pathlib import Path
from fractions import Fraction as Q
import base64,hashlib,io,json,datetime
import numpy as np
import sympy as s

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
checks={}
e,mu,kappa=s.symbols('e mu kappa',real=True);J=s.symbols('J',positive=True)
phi=-mu*s.log(J)+kappa*s.log(J)**2/2
checks['bulk_phi_prime_identity']=s.simplify(s.diff(phi,J)-(kappa*s.log(J)-mu)/J)==0
checks['bulk_phi_second_identity']=s.simplify(s.diff(phi,J,2)-(kappa+mu-kappa*s.log(J))/J**2)==0
checks['bulk_phi_third_identity']=s.simplify(s.diff(phi,J,3)-(2*kappa*s.log(J)-3*kappa-2*mu)/J**3)==0
P=s.series(s.diff(phi,J).subs(J,1+e),e,0,9).removeO().expand()
H=s.series(s.diff(phi,J,2).subs(J,1+e),e,0,9).removeO().expand()
for n in range(9):
    checks[f'closed_phi_prime_coefficient_{n}']=s.simplify(P.coeff(e,n)-(-1)**n*(-mu-kappa*s.harmonic(n)))==0
    checks[f'closed_phi_second_coefficient_{n}']=s.simplify(H.coeff(e,n)-(-1)**n*(n+1)*(mu+kappa*s.harmonic(n+1)))==0
r=s.symbols('r');G=1/(1-r)
checks['geometric_k_sum_identity']=s.simplify(r*s.diff(G,r)-r/(1-r)**2)==0
checks['geometric_k_squared_sum_identity']=s.simplify(r*s.diff(r*s.diff(G,r),r)-r*(1+r)/(1-r)**3)==0
q,K,T=s.symbols('q K T',positive=True);psi=K*(s.sqrt(q)-1)**2/2+T*s.sqrt(q)
checks['surface_first_identity']=s.simplify(s.diff(psi,q)-(K/2+(T-K)/(2*s.sqrt(q))))==0
checks['surface_second_identity']=s.simplify(s.diff(psi,q,2)+(T-K)/(4*q**s.Rational(3,2)))==0
checks['surface_third_identity']=s.simplify(s.diff(psi,q,3)-3*(T-K)/(8*q**s.Rational(5,2)))==0
ac,ai,rr,ll=s.symbols('ac ai rr ll',positive=True)
checks['surface_same_primitive_metric_rescaling']=s.simplify(ac*(rr*ai/ac)*ll-ai*rr*ll)==0
for p in [s.Rational(1,2),s.Rational(3,2)]:
    series=s.series((1+e)**(-p),e,0,8).removeO().expand()
    checks[f'surface_binomial_coefficients_p{p}']=all(s.simplify(series.coeff(e,n)-(-1)**n*s.rf(p,n)/s.factorial(n))==0 for n in range(8))

# Axis chart fixtures: only small symbolic maps, no actual current fields.
x,z,ss,tt=s.symbols('x z ss tt');a,b,c,d,f=s.symbols('a b c d f')
lam=[1-ss,ss*(1-tt),ss*tt]
N=[v*(2*v-1)for v in lam]+[4*lam[0]*lam[1],4*lam[1]*lam[2],4*lam[2]*lam[0]]
radial=s.expand(sum(N[j]*v for j,v in enumerate([0,a,b,c,d,f])))
checks['axis_vertex_Duffy_R_divides_s']=s.denom(s.cancel(radial/ss))==1
checks['axis_vertex_physical_r_divisibility_not_assumed']=s.rem(x+z,x,x)==z
checks['axis_vertex_direction_dependent_example_regularized']=s.simplify(((x+z)/x).subs({x:ss,z:ss*tt})-(1+tt))==0
l0,l1,l2=s.symbols('l0 l1 l2');R=a*l0*(2*l0-1)+4*c*l0*l1+4*f*l2*l0
checks['axis_edge_R_divides_reference_r']=s.denom(s.cancel(R/l0))==1
edge=ss*(2*ss-1)*a+4*ss*(1-ss)*b
checks['axis_surface_R_divides_t_before_AD']=s.denom(s.cancel(edge/ss))==1
checks['axis_surface_free_variations_also_divide_t']=all(s.denom(s.cancel(s.diff(edge,v)/ss))==1 for v in [a,b])

# Exact radial-Laurent moments versus independent symbolic integration.
rl=s.Rational(1,2);ru=s.Rational(3,2);zl=-s.Rational(1,5);zu=zl+s.Rational(2,3)*(x-rl)/(ru-rl)
def primitive_power(k):
    return s.log(ru/rl)if k==-1 else (ru**(k+1)-rl**(k+1))/(k+1)
for aa in [-3,-2,-1,0,1,2]:
    for bb in range(4):
        poly=s.Poly(s.expand((zu**(bb+1)-zl**(bb+1))/(bb+1)),x)
        moment=sum(coeff*primitive_power(aa+degree[0])for degree,coeff in poly.terms())
        direct=s.integrate(x**aa*s.integrate(z**bb,(z,zl,zu)),(x,rl,ru))
        checks[f'Laurent_triangle_moment_{aa}_{bb}']=s.simplify(moment-direct)==0

# Binary nominal Gauss nodes do not have exactly zero rational moment defects.
node=Q.from_float(float(1/np.sqrt(3)))
checks['binary_Gauss2_quadratic_moment_defect_nonzero']=2*node**2!=Q(2,3)
checks['symbolic_parameter_cancellation_vs_independent_interval']=Q(1,6)*0**2+Q(4,6)*Q(1,2)**2+Q(1,6)*1**2==Q(1,3)

def harmonic(n):return sum((Q(1,k)for k in range(1,n+1)),Q(0))
rho=Q(21,500);tails=[]
for n in[8,10,12,16,20]:
    C=1+5000*harmonic(n+2);B=Q(5000,n+3)
    e1=rho**(n+1)*((1+5000*harmonic(n+1))/(1-rho)+Q(5000)*rho/((n+2)*(1-rho)**2))
    e2=rho**(n+1)*((n+2)*C/(1-rho)+((n+2)*B+C)*rho/(1-rho)**2+B*rho*(1+rho)/(1-rho)**3)
    partial1=sum(((1+5000*harmonic(k))*rho**k for k in range(n+1,n+51)),Q(0))
    partial2=sum(((k+1)*(1+5000*harmonic(k+1))*rho**k for k in range(n+1,n+51)),Q(0))
    checks[f'derivative_tail_majorizes_first50_terms_N{n}']=partial1<=e1 and partial2<=e2
    tails.append({'degree':n,'rho_exact':str(rho),'phi_prime_tail_upper_exact':str(e1),'phi_second_tail_upper_exact':str(e2),'phi_prime_tail_float':float(e1),'phi_second_tail_float':float(e2),'scope':'Scalar constitutive derivative tails only, BEFORE chain rules, weights, spatial moments, coefficient wrapping, cache consistency, preconditioning or root criteria.'})

# Extract only shared metadata arrays from one hash-linked actual cache.
proto=ROOT/'runs/cohort_energy_stability_v1/results.json';stab=read(proto);row=stab['records'][0]
rp=ROOT/row['file'];assert sha(rp)==row['sha256'];binding=read(rp)['root_binding'];cp=ROOT/binding['point_cache_file'];assert sha(cp)==binding['point_cache_sha256']
packet=read(cp)['frozen_binary_cache'];raw=base64.b64decode(packet['data']);assert hashlib.sha256(raw).hexdigest()==packet['sha256']
with np.load(io.BytesIO(raw))as cache:
    N=cache['shape'];vol=[]
    for i,nrow in enumerate(N):
        qrow=[Q.from_float(float(x))for x in nrow];xi=qrow[1]+(qrow[3]+qrow[4])/2;eta=qrow[2]+(qrow[4]+qrow[5])/2;zero=1-xi-eta
        assert min(xi,eta,zero)>0
        vol.append({'row':i,'lambda0_exact':str(zero),'lambda1_exact':str(xi),'lambda2_exact':str(eta)})
    surfaces=[]
    for k in[0,1]:
        for i,nrow in enumerate(cache[f'surface_{k}_N']):
            qrow=[Q.from_float(float(x))for x in nrow];t=qrow[2]+qrow[1]/2;assert 0<t<1
            surfaces.append({'surface':k,'row':i,'t_exact':str(t)})
checks['all64_intrinsic_auxiliary_nodes_strictly_in_triangle']=len(vol)==64
checks['all24_intrinsic_surface_nodes_strictly_in_edge']=len(surfaces)==24
checks['no_611D_integral_force_Hessian_or_solver_called']=True
assert all(checks.values()),[k for k,v in checks.items()if not v]
witness={'status':'PASS_INTRINSIC_AUXILIARY_NODE_METADATA','source_cache_file':binding['point_cache_file'],'source_cache_sha256':sha(cp),'frozen_cache_sha256':packet['sha256'],'volume_nodes':vol,'surface_nodes':surfaces,'scope':'Q* rule points only, reconstructed intrinsically from exact binary shape entries. Not a claim that original cached arrays form a coherent basis or exactly Gaussian rule.'}
(OUT/'auxiliary_node_witness.json').write_text(json.dumps(witness,indent=2)+'\n',encoding='utf-8')
result={'status':'PASS_SMALL_SYMBOLIC_AND_SOURCE_DESIGN_CHECKS','checks':checks,'check_count':len(checks),'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(Path(__file__)),'theory_sha256':sha(OUT/'TARGET-AND-INTEGRATION-DESIGN.md'),'auxiliary_node_witness_sha256':sha(OUT/'auxiliary_node_witness.json'),'input_hashes':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p)for p in[proto,rp,cp]},'scalar_derivative_tail_examples':tails,'scope':'Design/source and independent symbolic fixtures only; no actual integral error/root transfer or all17 integration certification.'}
(OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':result['status'],'checks':len(checks),'N16_scalar_tails':tails[3],'source_sha256':result['source_sha256']},ensure_ascii=True))
