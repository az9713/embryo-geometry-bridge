"""Reference-moment component verification; no force/root/PDE certification."""
from pathlib import Path
import sys, json, hashlib, datetime, math, itertools
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_moments as m
from flint import arb, ctx
OUT=Path(__file__).resolve().parent
checks={}
def ck(k,v):
    checks[k]=bool(v)
    if not v: raise AssertionError(k)
def reject(k,fn):
    try: fn()
    except (ValueError,ZeroDivisionError): checks[k]=True;return
    raise AssertionError(k)
def equal(k,x,y):
    d=x.ball(512)-(y.ball(512) if isinstance(y,m.LogMoment) else m.A(y))
    ck(k,d.contains(0) and d.rad()<arb('1e-100'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

ctx.prec=512
tri=((0,0),(1,0),(0,1))
for a,b in itertools.product(range(6),repeat=2):
    z=m.triangle_moment(tri,a,b)
    ck(f'simplex_{a}_{b}',not z.logarithms and z.constant==F(math.factorial(a)*math.factorial(b),math.factorial(a+b+2)))
rect=[((1,-2),(3,-2),(3,1)),((1,-2),(3,1),(1,1))]
for a,b in itertools.product((-5,-2,-1,0,3),(0,1,4)):
    axial=((F(1)-F(1,3))**(b+1)-(F(-2)-F(1,3))**(b+1))/F(b+1)
    observed=sum((m.triangle_moment(t,a,b,F(1,3)) for t in rect),m.LogMoment())
    equal(f'rectangle_{a}_{b}',observed,m.interval_power_moment(1,3,a)*axial)
gold=((1,0),(2,0),(2,1))
equal('golden_inverse',m.triangle_moment(gold,-1,0),m.LogMoment(F(1),((F(2),F(-1)),)))
equal('golden_inverse_square',m.triangle_moment(gold,-2,0),m.LogMoment(F(-1,2),((F(2),F(1)),)))
equal('three_radial_cuts',m.triangle_moment(((1,0),(2,1),(3,0)),-1,0),m.LogMoment(F(0),((F(3),F(3)),(F(2),F(-4)))))
for a,b in ((-3,3),(0,5),(4,2)):
    base=m.triangle_moment(gold,a,b,F(1,7))
    for i,p in enumerate(itertools.permutations(gold)):
        equal(f'permutation_{a}_{b}_{i}',m.triangle_moment(p,a,b,F(1,7)),base)
    translated=tuple((x,F(z)+F(5,9)) for x,z in gold)
    equal(f'axial_translation_{a}_{b}',m.triangle_moment(translated,a,b,F(1,7)+F(5,9)),base)
    mid=(F(3,2),F(0))
    equal(f'subdivision_{a}_{b}',m.triangle_moment((gold[0],mid,gold[2]),a,b,F(1,7))+m.triangle_moment((mid,gold[1],gold[2]),a,b,F(1,7)),base)
equal('axis_vertex_weighted',m.triangle_moment(((0,0),(1,0),(1,2)),1,0),F(2,3))
equal('axis_edge_weighted',m.triangle_moment(((0,0),(0,2),(1,0)),1,0),F(1,3))
ck('axis_vertex_chart_measure',m.axis_chart_moment('axis_vertex',0,0,1,2)==F(2,3))
ck('axis_edge_chart_measure',m.axis_chart_moment('axis_edge',0,0,1,2)==F(1,3))

# Independently derive the two Taylor coefficient lists with symbolic derivatives.
import sympy as sp
x=sp.Symbol('x'); phi=-sp.log(x)+2500*sp.log(x)**2
for order in (0,4,12,16):
    c1,c2=m.volumetric_derivative_coefficients(order,1,5000)
    for derivative,coeffs in ((1,c1),(2,c2)):
        expr=sp.diff(phi,x,derivative)
        expected=[]
        for n in range(order+1):
            expected.append(F(sp.diff(expr,x,n).subs(x,1))/math.factorial(n))
        ck(f'symbolic_derivative_{derivative}_N{order}',coeffs==expected)
    rho=F(21,500);t1,t2=m.volumetric_derivative_tail_bounds(order,rho,1,5000)
    for j,e in enumerate((-rho,-rho/2,F(0),rho/3,rho)):
        z=1+m.A(e); first=(5000*z.log()-1)/z;second=(5001-5000*z.log())/z**2
        for d,(truth,cs,tail) in enumerate(((first,c1,t1),(second,c2,t2)),1):
            p=sum((m.A(v)*m.A(e)**n for n,v in enumerate(cs)),arb(0))
            ck(f'tail_sample_d{d}_N{order}_{j}',abs(truth-p)<m.A(tail) if e else (truth-p).contains(0))

reject('bool_normalized_exponent',lambda:m.normalized_bulk_moment(tri,True,0,1))
reject('bool_surface_exponent',lambda:m.surface_moment(1,2,True))
reject('invalid_auxiliary_axial_exponent',lambda:m.auxiliary_moment([((1,0),1)],0,-1))
reject('axis_negative_Laurent',lambda:m.triangle_moment(tri,-1,0))
reject('degenerate',lambda:m.triangle_moment(((1,0),(2,0),(3,0)),0,0))
reject('rho_one',lambda:m.volumetric_derivative_tail_bounds(4,1,1,5000))
reject('bool_order',lambda:m.volumetric_derivative_coefficients(True,1,5000))
ck('LogMoment_exact_float_constructor',m.LogMoment(0.1).constant==F.from_float(0.1))
ck('LogMoment_duplicate_logs_combined',m.LogMoment(0,((2,1),(2,-1),(1,9))).logarithms==())
reject('LogMoment_nonpositive_base',lambda:m.LogMoment(0,((0,1),)))
reject('LogMoment_bool_constructor',lambda:m.LogMoment(True))

# Actual cached geometry: a reference-measure check, not a field integral bound.
from exact_fe_target import Target,read,rational
protocol=read(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json')
cachepath=next(p for p in protocol['input_hashes'] if '/00_optoGEF2_small_col_B/' in p and p.endswith('/point_cache.json'))
target=Target(read(ROOT/cachepath)); rule=target.auxiliary_rule(); a=target.a
defects=[]
for k,ids in enumerate(a['triangles']):
    vertices=tuple(tuple(rational(v) for v in row) for row in a['reference'][ids[:3]])
    points=[]
    for bary,w in zip(rule['volume_barycentric'],rule['volume_weights'][k]):
        pt=tuple(sum(lam*v[d] for lam,v in zip(bary,vertices)) for d in (0,1))
        points.append((pt,w))
    d=m.bulk_moment_defect(vertices,points,0,0,target.energy_scale,precision_bits=512)
    defects.append({'element':k,'exact_integral_minus_auxiliary_constant_moment':str(d),'nonzero_certified':not d.contains(0)})
ck('actual_144_positive_measure_moments',len(defects)==144)
ck('actual_binary_rule_not_assumed_exact',any(d['nonzero_certified'] for d in defects))
ck('precision_restored',ctx.prec==512)
sources={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts/exact_fe_moments.py',ROOT/'scripts/exact_fe_target.py',Path(__file__),ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json',ROOT/cachepath)}
receipt={'status':'PASS_COMPONENT_ONLY','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'count':len(checks),'source_hashes':sources,'actual_measure_defects':defects,'limits':['Reference moments and scalar series only.','No integrated force or Hessian bound, root persistence, continuum error, source-data fit or novelty established.']}
out=OUT/'verification_v2.json'
if out.exists():raise ValueError('write-once receipt exists')
out.write_text(json.dumps(receipt,indent=2,ensure_ascii=True)+'\n',encoding='utf8')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'nonzero_binary_measure_defects':sum(d['nonzero_certified'] for d in defects),'receipt_sha256':sha(out)}))
