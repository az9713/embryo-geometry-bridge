"""Independent component tests: no full-field integration, assembly or solves."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import argparse,base64,hashlib,io,itertools,json,sys
import numpy as np
from flint import arb,acb,fmpq,ctx

HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];SOURCE=CF/'scripts/exact_fe_moments.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def A(x):
    v=Q(x);return arb(fmpq(v.numerator,v.denominator))
def pack(x):return {'lower_exact':str(x.lower().fmpq()),'upper_exact':str(x.upper().fmpq())}
def poly_mul(a,b):
    out={}
    for x,v in a.items():
        for y,w in b.items():
            z=tuple(i+j for i,j in zip(x,y));out[z]=out.get(z,Q(0))+v*w
    return out
def poly_pow(p,n):
    out={(0,0,0):Q(1)}
    for _ in range(n):out=poly_mul(out,p)
    return out
def barycentric_moment(vertices,a,b,shift):
    # Independent simplex Dirichlet moments; never slices triangle radially.
    v=[[Q(x) for x in row]for row in vertices]
    det=abs((v[1][0]-v[0][0])*(v[2][1]-v[0][1])-(v[2][0]-v[0][0])*(v[1][1]-v[0][1]))
    monomials=[(1,0,0),(0,1,0),(0,0,1)]
    pr={m:p[0]for m,p in zip(monomials,v)};pz={m:p[1]-shift for m,p in zip(monomials,v)}
    polynomial=poly_mul(poly_pow(pr,a),poly_pow(pz,b))
    return det*sum((c*Q(factorial(i)*factorial(j)*factorial(k),factorial(i+j+k+2))for(i,j,k),c in polynomial.items()),Q(0))
def independent_coefficients(order,mu,kappa):
    # Formal convolution of log(1+x) and geometric inverse; no harmonics.
    log=[Q(0)]+[Q((-1)**(n+1),n)for n in range(1,order+3)]
    inv=[Q((-1)**n)for n in range(order+3)]
    first=[-mu*inv[n]+kappa*sum((log[k]*inv[n-k]for k in range(n+1)),Q(0))for n in range(order+2)]
    return first[:order+1],[(n+1)*first[n+1]for n in range(order+1)]

def main(expected):
    if sha(SOURCE)!=expected:raise ValueError('Frozen producer SHA mismatch; do not test another revision')
    sys.path.insert(0,str(CF/'scripts'));import exact_fe_moments as m
    ctx.prec=256;checks={};numeric=[];guards=[]
    def record(name,truth):checks[name]=bool(truth)
    triangles=[[(0,0),(1,0),(0,1)],[(1,-1),(2,2),(3,0)],[(1,4),(1,1),(4,-2)],[(Q(1,3),Q(-2,7)),(Q(4,3),Q(1,2)),(Q(3,2),Q(9,7))]]
    number=0
    for ti,verts in enumerate(triangles):
        for a,b in itertools.product(range(5),range(5)):
            shift=Q(2,7) if (a+b)%2 else Q(0)
            want=barycentric_moment(verts,a,b,shift)
            for pi,perm in enumerate(itertools.permutations(verts)):
                got=m.triangle_moment(perm,a,b,shift)
                record(f'exact_barycentric_T{ti}_a{a}_b{b}_permutation{pi}',got.constant==want and not got.logarithms);number+=1
    # Known analytic Laurent fixtures independent of any code's strip solver.
    l=m.triangle_moment([(1,0),(2,0),(1,1)],-1,0)
    record('analytic_log_triangle_2log2_minus1',l.constant==-1 and dict(l.logarithms)=={Q(2):Q(2)})
    l=m.triangle_moment([(1,0),(2,0),(1,1)],-2,0)
    record('analytic_laurent_triangle_1_minus_log2',l.constant==1 and dict(l.logarithms)=={Q(2):Q(-1)})
    for a in [-5,-2,-1,0,3]:
        moment=m.surface_moment(1,3,a)
        if a==-2:
            record(f'analytic_surface_moment_{a}',moment.constant==0 and dict(moment.logarithms)=={Q(3):Q(1)})
        else:
            record(f'analytic_surface_moment_{a}',moment.constant==Q(3**(a+2)-1,a+2)if a+2>=0 else moment.constant==(Q(3)**(a+2)-1)/Q(a+2))
    normalized=m.normalized_bulk_moment(triangles[0],0,0,2)
    record('normalized_bulk_known_pi_over6',normalized.overlaps(arb.pi()/6))
    # Independent rigorous numerical integrator. Callback is meromorphic (integer
    # powers only); poles automatically produce non-finite balls for Arb's probe.
    fixtures=[('A',[(1,0),(3,0),(1,2)],[(1,3,0,0,-1,3)]),('B',[(1,-1),(2,2),(3,0)],[(1,2,Q(1,2),Q(-3,2),3,-4),(2,3,Q(1,2),Q(-3,2),-2,6)]),('C',[(1,4),(1,1),(4,-2)],[(1,4,-1,2,-2,6)])]
    for name,verts,strips in fixtures:
        for a,b in [(-8,0),(-3,1),(-1,0),(-1,3),(2,2)]:
            shift=Q(2,7);quad=acb(0)
            for lo,hi,ml,cl,mu,cu in strips:
                def f(z,analytic):
                    lower=A(ml)*z+A(cl)-A(shift);upper=A(mu)*z+A(cu)-A(shift)
                    return z**a*(upper**(b+1)-lower**(b+1))/(b+1)
                quad+=acb.integral(f,lo,hi,rel_tol=A(Q(1,2**200)),abs_tol=A(Q(1,2**200)))
            target=m.triangle_moment(verts,a,b,shift).ball(256)
            finite=quad.is_finite();overlap=finite and quad.imag.contains(0) and quad.real.overlaps(target)
            narrow=finite and (quad.real.upper()-quad.real.lower())<A(Q(1,2**150)) and (target.upper()-target.lower())<A(Q(1,2**150))
            record(f'Arb_independent_Laurent_{name}_a{a}_b{b}',overlap and narrow)
            numeric.append({'fixture':name,'a':a,'b':b,'shift':str(shift),'symbolic_ball':pack(target),'independent_certified_quadrature_real':pack(quad.real),'imag_contains_zero':quad.imag.contains(0),'overlap_and_width_lt_2pow_minus150':bool(overlap and narrow)})
    # Integrate s^(a+2)*t^b or s^(a+1)*(1-s)*t^b independently.
    for kind in ['axis_vertex','axis_edge']:
        for a,b in itertools.product(range(6),range(4)):
            radial=Q(7,3);det=Q(-11,5)
            wanted=abs(det)*radial*Q(1,b+1)* (Q(1,a+3)if kind=='axis_vertex'else Q(1,a+2)-Q(1,a+3))
            record(f'axis_measure_{kind}_{a}_{b}',m.axis_chart_moment(kind,a,b,radial,det)==wanted)
    # Cross-check series against a separate formal convolution and their exact
    # analytic primitive at the worst sign x=-rho, with high precision enclosures.
    for mu,kappa in [(Q(1),Q(5000)),(Q(2,7),Q(3,2)),(Q(-2),Q(3))]:
        for n in [0,1,2,8,12,20]:
            c1,c2=m.volumetric_derivative_coefficients(n,mu,kappa)
            want=independent_coefficients(n,mu,kappa)
            record(f'convolution_coefficients_{mu}_{kappa}_{n}',(c1,c2)==want)
            for rho in [Q(1,1000),Q(21,500),Q(1,2),Q(9,10)]:
                b1,b2=m.volumetric_derivative_tail_bounds(n,rho,mu,kappa)
                for sign in [-1,1]:
                    x=sign*rho;j=1+x;logj=A(j).log()
                    actual1=(A(kappa)*logj-A(mu))/A(j)
                    actual2=(A(kappa)+A(mu)-A(kappa)*logj)/A(j)**2
                    p1=A(sum((v*x**i for i,v in enumerate(c1)),Q(0)));p2=A(sum((v*x**i for i,v in enumerate(c2)),Q(0)))
                    record(f'closed_form_tail_{mu}_{kappa}_{n}_{rho}_{sign}',(actual1-p1).abs_upper()<=A(b1) and (actual2-p2).abs_upper()<=A(b2))
    record('rho_zero_tails_exact_zero',m.volumetric_derivative_tail_bounds(20,0,1,5000)==(0,0))
    record('zero_material_tails_exact_zero',m.volumetric_derivative_tail_bounds(8,Q(9,10),0,0)==(0,0))
    # Separate verification of the author's surface-series theory; not an
    # assertion that this root module exports a surface-tail implementation.
    for p in [Q(1,2),Q(3,2)]:
        for n in [0,1,8,16]:
            coefficients=[Q(1)]
            for i in range(n+1):coefficients.append(coefficients[-1]*(p+i)/(i+1))
            for rho in [Q(21,500),Q(1,4),Q(1,2)]:
                ratio=rho if p==Q(1,2) else rho*(n+1+p)/(n+2)
                if ratio>=1:raise ValueError('independent fixture geometric ratio invalid')
                bound=coefficients[n+1]*rho**(n+1)/(1-ratio)
                actual=A(1-rho).rsqrt()
                if p==Q(3,2):actual/=A(1-rho)
                polynomial=A(sum((coefficients[i]*rho**i for i in range(n+1)),Q(0)))
                record(f'separate_surface_theory_tail_{p}_{n}_{rho}',(actual-polynomial).abs_upper()<=A(bound))
    lm=m.LogMoment(0.1,((Q(2),Q(1,3)),(2,Q(2,3))))
    record('LogMoment_direct_constructor_exact_binary_float_and_duplicate_logs',isinstance(lm.constant,Q) and lm.constant==Q.from_float(.1) and dict(lm.logarithms)=={Q(2):Q(1)})
    # These source helper boundaries must reject Boolean/fractional exponents.
    bad=[('interval_bool',lambda:m.interval_power_moment(1,2,True)),('interval_float',lambda:m.interval_power_moment(1,2,1.2)),('triangle_radial_bool',lambda:m.triangle_moment(triangles[1],True,0)),('triangle_axial_negative',lambda:m.triangle_moment(triangles[1],1,-1)),('normalized_radial_bool',lambda:m.normalized_bulk_moment(triangles[1],True,0,1)),('normalized_radial_float',lambda:m.normalized_bulk_moment(triangles[1],0.5,0,1)),('surface_bool',lambda:m.surface_moment(1,2,True)),('surface_float',lambda:m.surface_moment(1,2,0.5)),('auxiliary_bool',lambda:m.auxiliary_moment([((1,1),1)],True)),('auxiliary_float',lambda:m.auxiliary_moment([((1,1),1)],0.5)),('auxiliary_negative_axial',lambda:m.auxiliary_moment([((1,1),1)],0,-1)),('axis_bool',lambda:m.axis_chart_moment('axis_vertex',True,0,1,1)),('axis_unknown',lambda:m.axis_chart_moment('guess',0,0,1,1)),('axis_divergent_unregularized',lambda:m.triangle_moment(triangles[0],-1,0)),('zero_normalization',lambda:m.normalized_bulk_moment(triangles[1],1,0,0)),('bad_rho',lambda:m.volumetric_derivative_tail_bounds(2,1,1,2)),('negative_order',lambda:m.volumetric_derivative_coefficients(-1,1,2)),('LogMoment_nonpositive_base',lambda:m.LogMoment(0,((0,1),))),('LogMoment_bool_constant',lambda:m.LogMoment(True)),('LogMoment_bool_coefficient',lambda:m.LogMoment(0,((2,True),)))]
    for name,call in bad:
        rejected=False;detail='accepted'
        try:call()
        except (ValueError,TypeError)as error:rejected=True;detail=type(error).__name__
        record('guard_'+name,rejected);guards.append({'case':name,'rejected':rejected,'outcome':detail})
    # Context restoration and actual-rule defect retention.
    ctx.prec=181;m.triangle_moment(fixtures[0][1],-1,0).ball(256)
    record('LogMoment_ball_restores_precision',ctx.prec==181);ctx.prec=256
    defect=m.bulk_moment_defect([(1,0),(2,0),(1,1)],[((Q(4,3),Q(1,3)),Q(1))],0,0,1)
    record('actual_rule_nonzero_constant_defect_retained',not defect.contains(0))
    record('source_hash_unchanged_during_component_tests',sha(SOURCE)==expected)
    # Actual geometry/constants only: independently integrate affine reference r
    # by simplex moments, then retain rational cache-weight errors with exact pi.
    joinpath=CF/'evidence/shared_cell_cohort_inference.json';joined=json.loads(joinpath.read_text());p=joined['proof_provenance'][0]
    rootpath=CF/p['root_certificate_file']
    if sha(rootpath)!=p['root_certificate_sha256']:raise ValueError('selected root hash differs')
    root=json.loads(rootpath.read_text());cachepath=CF/root['point_cache_file'];cache=json.loads(cachepath.read_text())
    raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('actual cache checksum')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
    sigma=Q.from_float(float(arrays['energy_scale'][0]));actual=[]
    for index,t in enumerate(arrays['triangles']):
        verts=[[Q.from_float(float(v))for v in row]for row in arrays['reference'][t[:3]]]
        exact_r_moment=barycentric_moment(verts,1,0,Q(0))
        independent_integral=2*arb.pi()*A(exact_r_moment)/A(sigma)
        rational_rule=sum((Q.from_float(float(w))/sigma for w in arrays['volume_weights'][index]),Q(0))
        independent_defect=independent_integral-A(rational_rule)
        target=m.bulk_moment_defect(verts,[((1,0),rational_rule)],0,0,sigma)
        ok=target.overlaps(independent_defect)and not target.contains(0)and(target.upper()-target.lower())<A(Q(1,2**150))
        record(f'actual144_weighted_constant_defect_{index}',ok)
        actual.append({'element':index,'independent_defect':pack(independent_defect),'source_defect':pack(target),'nonzero_and_agreement':bool(ok)})
    record('all144_actual_reference_elements_examined',len(actual)==144)
    out={'status':'PASS_INDEPENDENT_MOMENT_COMPONENT'if all(checks.values())else'FAIL_COMPONENT_REVIEW','checks':checks,'source_sha256':expected,'review_source_sha256':sha(Path(__file__)),'precision_bits':256,'exact_polynomial_fixture_count':number,'certified_Laurent_quadrature_cases':numeric,'actual_element_weighted_constant_defects':actual,'actual_input_sha256':{str(v.relative_to(CF)).replace('\\','/'):sha(v)for v in[joinpath,rootpath,cachepath]},'guard_cases':guards,'scope':'Different-algorithm exact simplex polynomial moments, two analytic logarithmic identities, rigorous acb meromorphic quadrature, independent formal power coefficients and closed-form derivative-tail checks, all144 actual affine-r weighted constant defects. Surface half-power geometric-tail theory is separately checked; it is not exported by the source. No whole-field integration or root transfer.','not_covered':['All17FEgradient/Hessian integration-error bounds','An exported surfacehalf-powercoefficients/tails implementation (absent from reviewed module)','General unequal-nonaxis-radius chart (axis helper assumes constant reference_radius)','Continuumerror/sourcevalidation/novelty']}
    out['review_doc_sha256']=sha(HERE/'REVIEW.md')
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');(CF/'evidence/exact_fe_moments_independent_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(CF/'evidence/exact_fe_moments_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--source-sha',required=True);args=parser.parse_args();main(args.source_sha)
