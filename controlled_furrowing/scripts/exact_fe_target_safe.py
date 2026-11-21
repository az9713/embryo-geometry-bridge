"""Canonical exact-P2/affine-material target, distinct from frozen binary Q-cache.
SAFE V2: exact-zero-only interval cleanup. Preparation/local algebra only; no spatial integration, root proof or solve.
Polynomial keys are exponents of (R,zeta) off axis or (s,t) in an axis chart.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import argparse,base64,datetime,hashlib,io,json,math
import numpy as np
from flint import arb,ctx
import algebraic_root_feasibility as core
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/target_preparation_safe_v2'
ORDER=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')
def rational(v):
    if isinstance(v,(bool,np.bool_)): raise ValueError('Boolean is not a coefficient')
    if isinstance(v,(float,np.floating)):
        if not math.isfinite(float(v)): raise ValueError('finite binary required')
        return Q.from_float(float(v))
    return Q(v)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')
def bind(d,k,v):
    if k in d and d[k]!=v:raise ValueError('conflicting immutable binding: '+k)
    d[k]=v

def cadd(a,b):return (a if isinstance(a,arb)else core.A(a))+(b if isinstance(b,arb)else core.A(b))if isinstance(a,arb)or isinstance(b,arb)else a+b
def cmul(a,b):return (a if isinstance(a,arb)else core.A(a))*(b if isinstance(b,arb)else core.A(b))if isinstance(a,arb)or isinstance(b,arb)else a*b
def cdiv(a,b):return (a if isinstance(a,arb)else core.A(a))/(b if isinstance(b,arb)else core.A(b))if isinstance(a,arb)or isinstance(b,arb)else a/b
def is_exact_zero(value):
    if isinstance(value,arb):
        if not value.is_finite():raise ValueError('finite interval coefficient required')
        return value.is_zero()
    return value==0

def clean(p):
    # Arb !=0 means PROVABLY nonzero, not logical negation of exact zero.
    # Retain every finite nonzero-width zero-crossing ball.
    return {tuple(k):v for k,v in p.items() if not is_exact_zero(v)}
def constant(c):return clean({(0,0):c})
def add(a,b):
    p=dict(a)
    for k,v in b.items():p[k]=cadd(p.get(k,0),v)
    return clean(p)
def scale(a,c):return clean({k:cmul(v,c) for k,v in a.items()})
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
    p={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():p[(i+k,j+l)]=cadd(p.get((i+k,j+l),0),cmul(v,w))
    return clean(p)
def power(a,n):
    if not isinstance(n,int)or isinstance(n,bool):raise ValueError('integer non-Boolean polynomial power required')
    if n<0:
        if len(a)!=1:raise ValueError('negative power requires a nonzero monomial')
        (i,j),v=next(iter(a.items()));return {(i*n,j*n):v**n}
    p=constant(Q(1))
    for _ in range(n):p=mul(p,a)
    return p
def diff(a,axis):
    if axis not in (0,1):raise ValueError('axis 0 or 1 required')
    p={}
    for key,v in a.items():
        e=key[axis]
        if e:
            k=list(key);k[axis]-=1;p[tuple(k)]=v*e
    return clean(p)
def compose(a,x,y):
    out={}
    for (i,j),v in a.items():out=add(out,scale(mul(power(x,i),power(y,j)),v))
    return out
def shift(a,x,y):return compose(a,{(1,0):Q(1),(0,0):x},{(0,1):Q(1),(0,0):y})
def evaluate(a,x,y):
    if any((i<0 and x==0)or(j<0 and y==0) for i,j in a):raise ValueError('Laurent singular evaluation')
    out=0
    for (i,j),v in a.items():out=cadd(out,cmul(cmul(v,x**i),y**j))
    return out
def divide_s(a):
    if any(i<1 for i,j in a):raise ValueError('axis factor not exactly divisible by s')
    return {(i-1,j):v for (i,j),v in a.items()}
ONE=constant(Q(1));X={(1,0):Q(1)};Y={(0,1):Q(1)}
def p2(lam):
    return [mul(l,sub(scale(l,2),ONE)) for l in lam]+[scale(mul(lam[i],lam[j]),4)for i,j in ((0,1),(1,2),(2,0))]
def linear_combination(polys,coeff):
    p={}
    for a,c in zip(polys,coeff):p=add(p,scale(a,c))
    return p

def material_basis(corners):
    x=[[rational(v)for v in row]for row in corners]
    zcenter=sum((p[1]for p in x),Q(0))/3
    a=x[1][0]-x[0][0];b=x[2][0]-x[0][0];c=x[1][1]-x[0][1];d=x[2][1]-x[0][1];det=a*d-b*c
    if det<=0 or min(p[0]for p in x)<0:raise ValueError('positive oriented nonnegative material triangle required')
    rr=sub(X,constant(x[0][0]));zz=add(Y,constant(zcenter-x[0][1]))
    xi=scale(sub(scale(rr,d),scale(zz,b)),1/det)
    eta=scale(add(scale(rr,-c),scale(zz,a)),1/det)
    lam=[sub(sub(ONE,xi),eta),xi,eta]
    return {'corners':x,'zcenter':zcenter,'z_center':zcenter,'det':det,'detDM':det,'lambda':lam,'phi':p2(lam)}

def axis_chart(corners):
    zero=[i for i,p in enumerate(corners)if p[0]==0];off=[i for i in range(3)if i not in zero]
    if len(zero)not in(1,2):raise ValueError('axis chart requires exactly one or two axis vertices')
    lam=[None]*3
    if len(zero)==2:
        lam[off[0]]=X;lam[zero[0]]=mul(sub(ONE,X),sub(ONE,Y));lam[zero[1]]=mul(sub(ONE,X),Y);jac=sub(ONE,X);kind='axis_edge'
    else:
        lam[zero[0]]=sub(ONE,X);lam[off[0]]=mul(X,sub(ONE,Y));lam[off[1]]=mul(X,Y);jac=X;kind='axis_vertex'
    radii={corners[i][0]for i in off}
    if len(radii)!=1 or next(iter(radii))<=0:raise ValueError('declared R=s*R1 chart requires equal positive off-axis radii')
    rr=linear_combination(lam,[p[0]for p in corners]);zz=linear_combination(lam,[p[1]for p in corners])
    return {'kind':kind,'lambda':lam,'R':rr,'Z':zz,'R1':next(iter(radii)),'jacobian':jac}

def element_fields(corners,current,free_local=None):
    """q and local-free B polynomials; current includes every fixed coordinate.
    No physical value at the isolated chart corner is postulated: t may persist.
    """
    g=material_basis(corners);phi=g['phi'];grads=[[diff(p,k)for k in(0,1)]for p in phi]
    radial=linear_combination(phi,[p[0]for p in current]);vertical=linear_combination(phi,[p[1]for p in current])
    qs=[diff(radial,0),diff(radial,1),diff(vertical,0),diff(vertical,1)]
    canonical_ref=g['corners']+[[sum(g['corners'][j][c] for j in pair)/2 for c in(0,1)]for pair in((0,1),(1,2),(2,0))]
    free_local=[k for k in range(12)if k%2 or canonical_ref[k//2][0]>0]if free_local is None else list(free_local)
    if len(set(free_local))!=len(free_local) or any(k not in range(12)for k in free_local):raise ValueError('unique valid local DOFs required')
    Bs=[]
    for k in free_local:
        j,component=divmod(k,2);Bs.append([grads[j][0]if component==0 else{},grads[j][1]if component==0 else{},grads[j][0]if component==1 else{},grads[j][1]if component==1 else{}])
    on_axis=any(p[0]==0 for p in g['corners'])
    if on_axis:
        chart=axis_chart(g['corners']);r,z=chart['R'],sub(chart['Z'],constant(g['zcenter']))
        transform=lambda p:compose(p,r,z)
        qs=[transform(p)for p in qs]
        regular_radial={}
        for p,node in zip(phi,current):
            if not is_exact_zero(node[0]):regular_radial=add(regular_radial,scale(divide_s(transform(p)),node[0]))
        qs.append(scale(regular_radial,1/chart['R1']))
        for col,k in zip(Bs,free_local):col[:]=[transform(p)for p in col]+[scale(divide_s(transform(phi[k//2])),1/chart['R1'])if k%2==0 else{}]
        measure=scale(mul(chart['R'],chart['jacobian']),g['det'])
        variables=('s','t')
    else:
        chart={'kind':'off_axis','R':X,'jacobian':ONE};qs.append(mul(radial,{(-1,0):Q(1)}))
        for col,k in zip(Bs,free_local):col.append(mul(phi[k//2],{(-1,0):Q(1)})if k%2==0 else{})
        measure=X;variables=('R','zeta')
    return {**g,'q':qs,'B':[[c[i]for c in Bs]for i in range(5)],'free_local':free_local,'chart':chart,'variables':variables,'measure_polynomial':measure}

def surface_fields(reference,current,free_local=None):
    """s in [0,1], edge node order [left,mid,right]; canonical horizontal R."""
    ref=[[rational(v)for v in row]for row in reference]
    L=ref[2][0]-ref[0][0]
    if L<=0 or ref[0][0]<0:raise ValueError('outward horizontal canonical edge required')
    rref=add(constant(ref[0][0]),scale(X,L))
    phi=[mul(sub(ONE,X),sub(ONE,scale(X,2))),scale(mul(X,sub(ONE,X)),4),mul(X,sub(scale(X,2),ONE))]
    radial=linear_combination(phi,[p[0]for p in current]);vertical=linear_combination(phi,[p[1]for p in current])
    free_local=[k for k in range(6)if k%2 or ref[k//2][0]>0]if free_local is None else list(free_local)
    if ref[0][0]==0:
        a=linear_combination([scale(divide_s(p),1/L)if not is_exact_zero(row[0]) else {}for p,row in zip(phi,current)],[row[0]for row in current]);basis_a=[scale(divide_s(phi[k//2]),1/L)if k%2==0 else{}for k in free_local]
    else:
        # Rational denominator stays explicit; polynomial numerators available for series.
        a=radial;basis_a=[phi[k//2]if k%2==0 else{}for k in free_local]
    q=[a,scale(diff(radial,0),1/L),scale(diff(vertical,0),1/L)]
    B=[basis_a,[scale(diff(phi[k//2],0),1/L)if k%2==0 else{}for k in free_local],[scale(diff(phi[k//2],0),1/L)if k%2 else{}for k in free_local]]
    physical_s=scale(sub(X,constant(ref[0][0])),1/L)
    physical_q=[compose(p,physical_s,{})for p in q];physical_B=[[compose(p,physical_s,{})for p in row]for row in B]
    if ref[0][0]!=0:
        physical_q[0]=mul(physical_q[0],{(-1,0):Q(1)});physical_B[0]=[mul(p,{(-1,0):Q(1)})for p in physical_B[0]]
    return {'q_physical_R':physical_q,'B_physical_R':physical_B,'Rlo':ref[0][0],'Rhi':ref[2][0],'q_numerators':q,'B_numerators':B,'hoop_denominator':ONE if ref[0][0]==0 else rref,'phi':phi,'Rref':rref,'L':L,'measure_polynomial':scale(rref,L),'variables':('s','unused'),'free_local':free_local}

def evaluate_bulk(fields,x,y):return [evaluate(p,x,y)for p in fields['q']]
def evaluate_surface(fields,s):
    den=evaluate(fields['hoop_denominator'],s,0)
    if not den>0:raise ValueError('positive regularized surface hoop denominator required')
    q=[evaluate(p,s,0)for p in fields['q_numerators']];q[0]=cdiv(q[0],den)
    B=[[cdiv(evaluate(p,s,0),(den if i==0 else 1))for p in row]for i,row in enumerate(fields['B_numerators'])]
    return q,B

def decode(packet):
    c=packet['frozen_binary_cache'];raw=base64.b64decode(c['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=c['sha256']:raise ValueError('cache embedded bytes mismatch')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:a={k:z[k].copy()for k in z.files}
    if any(not np.isfinite(v).all()for v in a.values()):raise ValueError('finite cache arrays required')
    return a
class Target:
    def __init__(self,packet):
        self.packet=packet;self.a=decode(packet);a=self.a;n=len(a['u'])
        if n!=611 or a['reference'].shape!=(325,2)or a['triangles'].shape!=(144,6)or a['free'].shape!=(n,)or len(set(a['free']))!=n or min(a['free'])<0 or max(a['free'])>=650:raise ValueError('selected 12x6 P2/free geometry required')
        if not np.issubdtype(a['free'].dtype,np.integer)or not np.issubdtype(a['triangles'].dtype,np.integer)or min(a['triangles'].ravel())<0 or max(a['triangles'].ravel())>=325 or any(len(set(row))!=6 for row in a['triangles']):raise ValueError('integer valid connectivity/free map required')
        if any(int(k)%2==0 and a['reference'][int(k)//2,0]==0 for k in a['free']):raise ValueError('radial axis coordinates must be fixed')
        self.index={int(k):i for i,k in enumerate(a['free'])};self.energy_scale=rational(a['energy_scale'][0])
        if self.energy_scale<=0:raise ValueError('positive exact binary normalization required')
        if a['active'].shape!=(144,1)or set(np.unique(a['active']))!={0.,1.}or np.count_nonzero(a['active'])!=24:raise ValueError('elementconstant24active mask required')
    def nodes(self,u=None):
        u=self.a['u']if u is None else u
        if len(u)!=len(self.a['free']):raise ValueError('matching free vector required')
        y=[[rational(v)for v in row]for row in self.a['reference']]
        for i,k in enumerate(self.a['free']):y[int(k)//2][int(k)%2]=cadd(y[int(k)//2][int(k)%2],u[i]if isinstance(u[i],arb)else rational(u[i]))
        return y
    def node_tube(self,S,cell,base,radius):
        if np.asarray(S).shape!=(len(self.a['u']),5)or not np.isfinite(S).all():raise ValueError('finite n by5 predictor required')
        r=rational(radius)
        if r<=0:raise ValueError('positive normal corrector radius required')
        delta=[(rational(cell[k][0])-rational(base[k]),rational(cell[k][1])-rational(base[k]))for k in ORDER]
        if any(lo>hi for lo,hi in delta):raise ValueError('ordered cell required')
        lo=[];hi=[]
        for i,row in enumerate(S):
            terms=[(rational(v)*l,rational(v)*h)for v,(l,h)in zip(row,delta)]
            lo.append(rational(self.a['u'][i])+sum(min(t)for t in terms)-r);hi.append(rational(self.a['u'][i])+sum(max(t)for t in terms)+r)
        u=[core.A(l).union(core.A(h))for l,h in zip(lo,hi)]
        return self.nodes(u)
    def element(self,index,nodes=None):
        ids=self.a['triangles'][index];ref=self.a['reference'][ids[:3]];y=self.nodes()if nodes is None else nodes
        dofs=[2*int(j)+c for j in ids for c in(0,1)];free=[i for i,k in enumerate(dofs)if k in self.index]
        f=element_fields(ref,[y[int(j)]for j in ids],free);f['global_indices']=[self.index[dofs[i]]for i in free];f['active']=rational(self.a['active'][index,0]);return f
    def surface(self,k,e,nodes=None):
        ids=self.a[f'surface_{k}_edges'][e];y=self.nodes()if nodes is None else nodes
        dofs=[2*int(j)+c for j in ids for c in(0,1)];free=[i for i,v in enumerate(dofs)if v in self.index]
        f=surface_fields(self.a['reference'][ids],[y[int(j)]for j in ids],free);f['global_indices']=[self.index[dofs[i]]for i in free];return f
    def normalization(self):return 2*arb.pi()/core.A(self.energy_scale)
    def auxiliary_rule(self):
        a=self.a;volume=[]
        if np.any(a['volume_weights']<=0)or any(np.any(a[f'surface_{k}_weights']<=0)for k in(0,1)):raise ValueError('strictly positive auxiliary binary weights required')
        for N in a['shape']:
            N=list(map(rational,N));xi=N[1]+(N[3]+N[4])/2;eta=N[2]+(N[4]+N[5])/2
            if min(xi,eta,1-xi-eta)<=0:raise ValueError('strictly interior Qstar point required')
            volume.append([1-xi-eta,xi,eta])
        surfaces=[]
        for k in(0,1):
            s=[rational(N[2])+rational(N[1])/2 for N in a[f'surface_{k}_N']]
            if any(not 0<x<1 for x in s):raise ValueError('strictly interior Qstar surface point required')
            surfaces.append(s)
        return {'volume_barycentric':volume,'surface_s':surfaces,'volume_weights':[[rational(w)/self.energy_scale for w in row]for row in a['volume_weights']],'surface_weights':[[[rational(w)/self.energy_scale for w in row]for row in a[f'surface_{k}_weights']]for k in(0,1)],'weight_semantics':'Exact OLD normalized binary weights for an auxiliary finite rule only; no Gaussian exactness claim.'}


def tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    def reject(k,fn):
        try:fn();ck(k,False)
        except (ValueError,ZeroDivisionError):checks[k]=True
    z=arb(0,1)
    ck('Arb_uncertain_constant_retained',constant(z).get((0,0))is z)
    ck('Arb_uncertain_addition_retained',(0,0)in add(constant(arb(1)),constant(arb(-1,1))))
    ck('Arb_uncertain_multiplication_retained',(0,0)in mul(constant(z),constant(arb(2))))
    ck('Arb_exact_zero_removed',constant(arb(0))=={})
    ck('Arb_zero_overlap_not_zero',not is_exact_zero(z))
    ck('Arb_uncertain_derivative_retained',(0,0)in diff({(1,0):z},0))
    ck('Arb_uncertain_composition_retained',(0,0)in compose({(1,0):z},ONE,Y))
    ck('Arb_uncertain_power_retained',(0,0)in power(constant(z),2))
    ck('Arb_uncertain_evaluation_contains_nonzero_witness',evaluate(constant(z),Q(0),Q(0)).contains(core.A(Q(1))))
    reject('Boolean_power_rejected',lambda:power(X,True))
    ck('polynomial_product_derivative',diff(mul(add(X,Y),sub(X,Y)),0)==scale(X,2))
    ck('polynomial_shift_composition',shift(mul(X,Y),Q(2),Q(-3))==mul(add(X,constant(Q(2))),sub(Y,constant(Q(3)))))
    ck('Laurent_integer_operations',mul(power(X,-1),power(X,2))==X)
    reject('nonmonomial_negative_power_rejected',lambda:power(add(X,Y),-1))
    reject('Laurent_axis_zero_rejected',lambda:evaluate(power(X,-1),Q(0),Q(1)))
    reject('Boolean_binary_rejected',lambda:rational(True))
    reject('nonfinite_binary_rejected',lambda:rational(float('nan')))
    bindings={'a':'a'};reject('conflicting_binding_rejected',lambda:bind(bindings,'a','b'))
    triangles=[[[Q(1),Q(0)],[Q(2),Q(0)],[Q(1),Q(1)]],[[Q(0),Q(0)],[Q(1),Q(0)],[Q(0),Q(1)]],[[Q(0),Q(0)],[Q(1),Q(0)],[Q(1),Q(1)]]]
    for n,corners in enumerate(triangles):
        ref=corners+[[sum(corners[j][c]for j in pair)/2 for c in(0,1)]for pair in((0,1),(1,2),(2,0))]
        g=material_basis(corners);ck(f'partition_{n}',linear_combination(g['phi'],[Q(1)]*6)==ONE)
        ck(f'exact_R_reproduction_{n}',linear_combination(g['phi'],[x[0]for x in ref])==X)
        ck(f'exact_Z_reproduction_{n}',linear_combination(g['phi'],[x[1]for x in ref])==add(Y,constant(g['zcenter'])))
        field=element_fields(corners,ref)
        for point in [(Q(1,3),Q(2,5)),(Q(0),Q(0)),(Q(0),Q(1))] if n else [(Q(4,3),Q(1,7))]:
            ck(f'canonical_identity_{n}_{point}',evaluate_bulk(field,*point)==[1,0,0,1,1])
        stretched=[[Q(2)*x,Q(1,4)*z]for x,z in ref];sf=element_fields(corners,stretched)
        q=evaluate_bulk(sf,Q(1,3)if n else Q(4,3),Q(2,5));ck(f'affine_isochoric_J1_{n}',(q[0]*q[3]-q[1]*q[2])*q[4]==1)
        # Differentiate a polynomial nodal perturbation, not a solved equilibrium.
        k=field['free_local'][0];h=Q(1,101);plus=[row[:]for row in ref];minus=[row[:]for row in ref];plus[k//2][k%2]+=h;minus[k//2][k%2]-=h
        fp=element_fields(corners,plus,field['free_local']);fm=element_fields(corners,minus,field['free_local'])
        ck(f'B_exact_nodal_derivative_{n}',all(scale(sub(a,b),1/(2*h))==row[0]for a,b,row in zip(fp['q'],fm['q'],field['B'])))
    bad=[[Q(0),Q(0)],[Q(1),Q(0)],[Q(2),Q(1)]]
    reject('unsupported_vertex_axis_R1_rejected',lambda:axis_chart(bad))
    for i,r0 in enumerate([Q(0),Q(1)]):
        ref=[[r0,Q(0)],[r0+Q(1,2),Q(0)],[r0+1,Q(0)]];sf=surface_fields(ref,ref)
        for t in [Q(0),Q(1,3),Q(1)]:
            q,B=evaluate_surface(sf,t);ck(f'surface_canonical_identity_{i}_{t}',q==[1,1,0])
            R=r0+t;ck(f'surface_R_conversion_{i}_{t}',q==[evaluate(p,R,0)for p in sf['q_physical_R']])
        k=sf['free_local'][0];h=Q(1,113);plus=[row[:]for row in ref];minus=[row[:]for row in ref];plus[k//2][k%2]+=h;minus[k//2][k%2]-=h
        qp,_=evaluate_surface(surface_fields(ref,plus,sf['free_local']),Q(1,3));qm,_=evaluate_surface(surface_fields(ref,minus,sf['free_local']),Q(1,3));_,B=evaluate_surface(sf,Q(1,3))
        ck(f'surface_B_exact_nodal_derivative_{i}',all((a-b)/(2*h)==row[0]for a,b,row in zip(qp,qm,B)))
    # A regular axis field can have a direction-dependent chart-corner limit.
    corners=triangles[2];ref=corners+[[sum(corners[j][c]for j in pair)/2 for c in(0,1)]for pair in((0,1),(1,2),(2,0))]
    cur=[[r+z/10,z]for r,z in ref];field=element_fields(corners,cur)
    ck('axis_regularization_retains_directional_corner_dependence',evaluate(field['q'][4],0,0)!=evaluate(field['q'][4],0,1))
    # Arb coefficients remain valid after exact basis cancellation, including s=0.
    tube=[[core.A(r)+arb(0,core.A(Q(1,100000)))if r else Q(0),core.A(z)]for r,z in ref]
    tf=element_fields(corners,tube)
    ck('Arb_axis_coefficients_finite',all(v.is_finite()if isinstance(v,arb)else True for v in evaluate_bulk(tf,core.A(0),core.A(Q(1,2)))))
    # Independent symbolic scalar energy chain rule along a polynomial direction.
    import sympy as sy
    t=sy.symbols('t');v=[sy.Rational(6,5)+t/7,sy.Rational(1,10)-t/9,sy.Rational(-1,20)+t/11,sy.Rational(9,10)+t/13,sy.Rational(11,10)+t/17]
    J=(v[0]*v[3]-v[1]*v[2])*v[4];W=sy.Rational(3,2)*(sum(x*x for x in v)-3)-3*sy.log(J)+sy.Rational(7,2)*sy.log(J)**2+sy.Rational(5,2)*(v[1]**2+v[3]**2)
    q=[core.A(Q(str(x.subs(t,0))))for x in v];d=[core.A(Q(str(sy.diff(x,t))))for x in v];g,H,_=core.bulk_local(q,core.A(3),core.A(7),core.A(5))
    def arb_symbolic(expr):
        L=sy.log(sy.Rational(str(Q(1)))) # actual expression replaced below
        # Evaluate exact rationals and log atoms by interval operations.
        if expr.is_Rational:return core.A(Q(int(expr.p),int(expr.q)))
        if expr.func==sy.log:return arb_symbolic(expr.args[0]).log()
        if expr.is_Add:return sum((arb_symbolic(x)for x in expr.args),arb(0))
        if expr.is_Mul:
            out=arb(1)
            for x in expr.args:out*=arb_symbolic(x)
            return out
        if expr.is_Pow and expr.args[1].is_Integer:return arb_symbolic(expr.args[0])**int(expr.args[1])
        raise ValueError('unsupported symbolic check expression')
    first=sum(g[i]*d[i]for i in range(5));second=sum(H[i][j]*d[i]*d[j]for i in range(5)for j in range(5))
    ck('local_energy_directional_gradient_symbolic',first.overlaps(arb_symbolic(sy.diff(W,t).subs(t,0))))
    ck('local_energy_directional_Hessian_symbolic',second.overlaps(arb_symbolic(sy.diff(W,t,2).subs(t,0))))
    return checks

def cache_defects(target):
    """ALL 9216 points; exact Qstar vs cached shapes/physical gradients/qr.
    Merely computes interpolation defects, not an energy or force integral.
    """
    a=target.a;rule=target.auxiliary_rule();shape_max=Q(0);grad_max=Q(0);qr_max=Q(0);shape_changed=grad_changed=0;fixed_offsets={}
    for lam,N in zip(rule['volume_barycentric'],a['shape']):
        exact=[l*(2*l-1)for l in lam]+[4*lam[i]*lam[j]for i,j in((0,1),(1,2),(2,0))]
        for x,y in zip(exact,N):
            defect=abs(x-rational(y));shape_max=max(shape_max,defect);shape_changed+=defect!=0
    classes=Counter()
    for ti,ids in enumerate(a['triangles']):
        g=material_basis(a['reference'][ids[:3]]);co=g['corners'];axis=sum(x[0]==0 for x in co);classes[{0:'off_axis',1:'axis_vertex',2:'axis_edge'}[axis]]+=1
        A=co[1][0]-co[0][0];B=co[2][0]-co[0][0];C=co[1][1]-co[0][1];D=co[2][1]-co[0][1];det=g['det']
        dl=[(-(D-C)/det,(B-A)/det),(D/det,-B/det),(-C/det,A/det)]
        for iq,lam in enumerate(rule['volume_barycentric']):
            gr=[[(4*lam[k]-1)*dl[k][j]for j in(0,1)]for k in range(3)]+[[4*(lam[k]*dl[l][j]+lam[l]*dl[k][j])for j in(0,1)]for k,l in((0,1),(1,2),(2,0))]
            for k,row in enumerate(gr):
                for j in(0,1):
                    defect=abs(row[j]-rational(a['shape_grad'][ti,iq,k,j]));grad_max=max(grad_max,defect);grad_changed+=defect!=0
            R=sum(l*x[0]for l,x in zip(lam,co));qr_max=max(qr_max,abs(R-rational(a['qr'][ti,iq])))
        mids=[[(co[i][j]+co[k][j])/2 for j in(0,1)]for i,k in((0,1),(1,2),(2,0))]
        for node,canonical in zip(ids[3:],mids):
            for c in(0,1):
                dof=2*int(node)+c;gap=rational(a['reference'][node,c])-canonical[c]
                if dof not in target.index and gap:fixed_offsets[str(dof)]=str(gap)
    return {'all_volume_points':144*64,'shape_coefficients_checked':64*6,'physical_gradient_coefficients_checked':144*64*6*2,'shape_max_abs_defect_exact':str(shape_max),'shape_coefficients_changed':int(shape_changed),'physical_gradient_max_abs_defect_exact':str(grad_max),'physical_gradient_coefficients_changed':int(grad_changed),'material_radius_max_abs_defect_exact':str(qr_max),'fixed_binary_reference_midpoint_offsets_exact':fixed_offsets,'element_classes':dict(classes),'Qstar_is_exact_Gaussian_rule':False}

def select_records():
    joinedfile='evidence/shared_cell_cohort_inference.json';joined=read(ROOT/joinedfile);bindings={joinedfile:sha(ROOT/joinedfile)};records=[]
    for row in joined['proof_provenance']:
        for name in('root_certificate','root_review','observer_certificate','application_review'):
            f=row[name+'_file'];expected=row[name+'_sha256'];bind(bindings,f,expected)
            if sha(ROOT/f)!=expected:raise ValueError('changed selected proof packet: '+f)
        root=read(ROOT/row['root_certificate_file']);review=read(ROOT/row['root_review_file'])
        matches=[r for r in review['record_audits']if r['record_id']==row['record_id']]
        if len(matches)!=1:raise ValueError('one independent root capsule required')
        capsule=matches[0]['audit']
        if capsule['status']!='PASS'or any(v is not True for v in capsule['checks'].values())or capsule['certificate_sha256']!=row['root_certificate_sha256']or capsule['uniform_quantifier_verified']is not True:raise ValueError('matching independently audited selected root required')
        if root['status']!='CERTIFIED_UNIFORM_LOCAL_ROOT'or root['record_id']!=row['record_id']or rational(root['dose_exact'])!=rational(row['source_dose_exact']):raise ValueError('selected source record mismatch')
        for f,h in root['artifact_hashes'].items():
            bind(bindings,f,h)
            if sha(ROOT/f)!=h:raise ValueError('changed root input: '+f)
        cachefile=root['point_cache_file'];predfile=root['predictor_file']
        if cachefile not in root['artifact_hashes']or predfile not in root['artifact_hashes']:raise ValueError('cache/predictor must have authoritative root binding')
        packet=read(ROOT/cachefile);target=Target(packet);state=packet['input_state'];bind(bindings,state,packet['input_state_sha256'])
        if sha(ROOT/state)!=packet['input_state_sha256']:raise ValueError('proof center state changed')
        with np.load(ROOT/predfile,allow_pickle=False)as z:S=z['S'].copy()
        radius=rational(root['uniform_cell_result']['selected_radius_exact'])
        base={k:str(rational(packet['parameters'][k]))for k in ORDER[:4]};base['common_gain']=str(rational(root['base_strength_exact'])/rational(root['dose_exact']))
        cell=root['parameter_cell']
        if any(capsule[k]!=root[k]for k in('parameter_cell_sha256','model_observation_contract_sha256','branch_domain_id')):raise ValueError('root family capsule context mismatch')
        # Exercise whole nodal-tube construction; not an exact-target root assertion.
        tube=target.node_tube(S,cell,base,radius)
        if any(not(v.is_finite()if isinstance(v,arb)else True)for r in tube for v in r):raise ValueError('nonfinite nodal tube')
        records.append({'record_id':row['record_id'],'root_certificate_file':row['root_certificate_file'],'root_certificate_sha256':row['root_certificate_sha256'],'independent_root_capsule_file':row['root_review_file'],'independent_root_capsule_sha256':row['root_review_sha256'],'point_cache_file':cachefile,'point_cache_sha256':root['artifact_hashes'][cachefile],'frozen_binary_cache_sha256':packet['frozen_binary_cache']['sha256'],'predictor_file':predfile,'predictor_sha256':root['artifact_hashes'][predfile],'proof_center_state_file':state,'proof_center_state_sha256':packet['input_state_sha256'],'dose_exact':str(rational(root['dose_exact'])),'base_parameters_exact':base,'parameter_cell':cell,'parameter_cell_sha256':root['parameter_cell_sha256'],'model_observation_contract_sha256':root['model_observation_contract_sha256'],'selected_corrector_radius_exact':str(radius),'branch_domain_id':root['branch_domain_id'],'new_exact_target_root_certified':False})
    if len(records)!=17 or len({r['record_id']for r in records})!=17 or len({r['parameter_cell_sha256']for r in records})!=1:raise ValueError('all17 same-cell records required')
    return records,bindings

def prepare():
    ctx.prec=160
    if(OUT/'protocol.json').exists():raise ValueError('write-once protocol exists')
    checks=tests();records,bindings=select_records();first=Target(read(ROOT/records[0]['point_cache_file']));defects=cache_defects(first)
    if not rational(defects['shape_max_abs_defect_exact'])>0 or not rational(defects['physical_gradient_max_abs_defect_exact'])>0:raise ValueError('actual exact-vs-cache defects unexpectedly erased')
    checks['all17_selected_cache_and_independent_capsules_bound']=True;checks['all17_predictor_full_corrector_tubes_constructed']=True
    checks['actual_shape_and_gradient_defects_retained']=True;checks['axis_classes6plus6']=defects['element_classes']=={'axis_edge':6,'axis_vertex':6,'off_axis':132}
    checks['two_fixed_binary_Z_offsets_retained']=len(defects['fixed_binary_reference_midpoint_offsets_exact'])==2
    # Shared geometry cache equality is byte-value equality, not a new interpolation claim.
    geometry=[k for k in first.a if k not in('u','preconditioner')]
    for record in records[1:]:
        target=Target(read(ROOT/record['point_cache_file']))
        if any(not np.array_equal(target.a[k],first.a[k])for k in geometry):raise ValueError('selected17 geometry/integration caches differ')
    checks['geometry_auxiliary_rules_identical_across17']=True
    target_examples={}
    for kind in('off_axis','axis_edge','axis_vertex'):
        ti=next(i for i,ids in enumerate(first.a['triangles'])if {0:'off_axis',1:'axis_vertex',2:'axis_edge'}[sum(first.a['reference'][ids[:3],0]==0)]==kind)
        f=first.element(ti);x,y=(Q(1,3),Q(2,5))if kind!='off_axis'else (sum(p[0]for p in f['corners'])/3,Q(0))
        q=evaluate_bulk(f,x,y);core.bulk_local(list(map(core.A,q)),core.A(1),core.A(5000),core.A(first.packet['activation_strength_binary_exact'])*core.A(f['active']))
        target_examples[kind]={'element':ti,'variables':f['variables'],'z_center_exact':str(f['z_center']),'detDM_exact':str(f['detDM']),'R1_exact':str(f['chart'].get('R1')),'q_at_declared_sample_exact':list(map(str,q))}
    checks['actual_center_three_chart_positive_local_kernel_samples']=True
    if not all(checks.values()):raise AssertionError(checks)
    protocol={'version':1,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':{Path(__file__).relative_to(ROOT).as_posix():sha(__file__),'scripts/algebraic_root_feasibility.py':sha(core.__file__)},'input_hashes':bindings,'records':records,'target_definition':{'material_map':'P1 affine triangle from exact binary vertex corners; exact mathematical P2 basis [v0,v1,v2,e01,e12,e20].','current_nodes':'Y_j = exactbinary X_j + free injection u, including fixed binary coordinates unchanged; u=0 is not silently idealized to identity.','volume_energy':'(2*pi()/S_binary) integral_material R [mu/2*(q.q-3)-mu*logJ+kappa/2*logJ^2+chi*g*dose/2*(q_b^2+q_d^2)] dR dZ','volume_q':'(r_R,r_Z,z_R,z_Z,r/R); q and B=dq/du exactly polynomial/Laurent before evaluation.','surface_energy':'(2*pi()/S_binary) integral_edge Rref [K/2*(Js-1)^2+T*Js] dRref, Js=(r/Rref)*sqrt((dr/dRref)^2+(dz/dRref)^2).','surface_reference':'Canonical horizontal R(s)=Rlo+L*s with exact endpoint L, not cached P2 r0 or l0.','normalization':'Same positive exact binary energy_scale S; pi is mathematical pi, not floating pi.','activation':'Exactly24 elementconstant active triangles from old cache.','axis':'6 axis-edge+6 axis-vertex Duffy charts, referenceR=s*R1; all free radial basis/current radial polynomial factors s cancelled symbolically before q/B evaluation. No unique physical limit at isolated chart corner postulated.','auxiliary_rule':'Qstar from exact xi=N1+(N01+N12)/2, eta=N2+(N12+N20)/2; surface s=Nright+Nmid/2. Weights exactly old normalized binary weights. Not claimed Gaussian-exact.','bridge_not_yet_proved':'Exact-target minus Qstar moment defects plus finite Qstar-minus-cache correction must be validated later.'},'no_integration_or_new_solves':True,'no_exact_target_root_certificate':True,'no_PDE_continuum_or_source_time_claim':True,'external_author_contact_sent':False}
    write(OUT/'protocol.json',protocol);write(OUT/'guard_tests.json',{'status':'PASS','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'scope':'Exact local algebra/selected-target preparation only; no integration/root bridge.'})
    write(OUT/'cache_interpolation_defects.json',{'status':'EXACT_DEFECTS_RETAINED','defects':defects,'chart_samples':target_examples,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    preflight={'status':'PASS_TARGET_PREPARATION_ONLY','checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'guard_tests_sha256':sha(OUT/'guard_tests.json'),'cache_defects_sha256':sha(OUT/'cache_interpolation_defects.json'),'all17_prepared':True,'exact_FE_root_bridge_proved':False,'full_integration_run':False,'new_equilibrium_solve':False,'runtime_scope':'Only sparse polynomial preparation, exact local tests, binary/proof binding and cache interpolation defect calculations.'}
    write(OUT/'preflight.json',preflight)
    print(json.dumps({'status':preflight['status'],'checks':len(checks),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'records':17}),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['test','prepare']);args=parser.parse_args();ctx.prec=160
    if args.command=='test':print(json.dumps({'status':'PASS','checks':tests()}))
    else:prepare()
