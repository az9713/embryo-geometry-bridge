"""SAFE V2 independent local point library; no automatic execution or old-output premises."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,datetime,hashlib,json,sys,time
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2]
sys.path.insert(0,str(CF/'scripts'));sys.path.insert(0,str(CF/'runs/exact_integration_bridge_v1/safe_kernel_v2'))
import exact_fe_target_safe as T
import exact_fe_moments as M
import polynomial_jet_safe as DUAL
import algebraic_root_feasibility as C
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
    if Path(p).exists():raise ValueError('immutable output exists: '+str(p))
    Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(x,indent=2)+'\n')
def rational(x):return F.from_float(float(x))if isinstance(x,(float,np.floating))else F(x)
def product(items):
    p=T.ONE
    for x in items:p=T.mul(p,x)
    return p
def independent_jet(q,B):
    # Different construction: derivatives of two signed cubic monomials,
    # rather than differentiating a meridional determinant and multiplying h.
    terms=[(1,(0,3,4)),(-1,(1,2,4))];n=len(B[0]);value={};g=[{}for _ in range(n)];H={(i,j):{}for i in range(n)for j in range(i,n)}
    for sign,indices in terms:
        value=T.add(value,T.scale(product(q[k]for k in indices),sign))
        for i in range(n):
            for a in range(3):g[i]=T.add(g[i],T.scale(product([B[indices[a]][i]]+[q[indices[b]]for b in range(3)if b!=a]),sign))
        for i,j in H:
            for a in range(3):
                for b in range(3):
                    if a==b:continue
                    H[i,j]=T.add(H[i,j],T.scale(product([B[indices[a]][i],B[indices[b]][j]]+[q[indices[c]]for c in range(3)if c not in(a,b)]),sign))
    return {'value':value,'gradient':g,'hessian_upper':H,'local_dofs':n}
def coefficients(n):
    log=[F(0)]+[F((-1)**(k+1),k)for k in range(1,n+2)]
    first=[F(-(-1)**k)+5000*sum((log[j]*(-1)**(k-j)for j in range(k+1)),F(0))for k in range(n+2)]
    return first[:n+1],[(k+1)*first[k+1]for k in range(n+1)]
def series(e,n):
    first,second=coefficients(n);power=T.ONE;a={};b={};ee={k:C.A(v)for k,v in e.items()}
    for k in range(n+1):
        a=T.add(a,T.scale(power,C.A(first[k])));b=T.add(b,T.scale(power,C.A(second[k])))
        if k<n:power=T.mul(power,ee)
    return a,b
def majorant(p,f):
    if f['chart']['kind']!='off_axis':return sum(map(abs,p.values()),F(0))
    rlo=min(v[0]for v in f['corners']);rhi=max(v[0]for v in f['corners']);z=max(abs(v[1]-f['zcenter'])for v in f['corners'])
    return sum((abs(c)*max(rlo**i,rhi**i)*z**j for(i,j),c in p.items()),F(0))
class Functional:
    def __init__(self,target,f,index,rule):
        self.target=target;self.f=f;self.points=[];self.lambdas=rule['volume_barycentric'];verts=f['corners'];zero=[i for i,x in enumerate(verts)if x[0]==0];off=[i for i in range(3)if i not in zero]
        for lam,w in zip(self.lambdas,rule['volume_weights'][index]):
            if f['chart']['kind']=='axis_vertex':
                ss=lam[off[0]]+lam[off[1]];xy=(ss,lam[off[1]]/ss)
            elif f['chart']['kind']=='axis_edge':xy=(lam[off[0]],lam[zero[1]]/(1-lam[off[0]]))
            else:xy=(sum(l*v[0]for l,v in zip(lam,verts)),sum(l*v[1]for l,v in zip(lam,verts))-f['zcenter'])
            self.points.append((xy,w))
        self.balls=[(C.A(x),C.A(y),C.A(w))for(x,y),w in self.points]
        self.massI=self.integral(0,0);self.massQ=C.A(sum(w for _,w in self.points))
        if not self.massI>0 or not self.massQ>0:raise ValueError('positive two-mass measures required')
    @lru_cache(None)
    def integral(self,i,j):
        if self.f['chart']['kind']=='off_axis':return M.normalized_bulk_moment(self.f['corners'],i,j,self.target.energy_scale,self.f['zcenter'],640)
        return self.target.normalization()*C.A(M.axis_chart_moment(self.f['chart']['kind'],i,j,self.f['chart']['R1'],self.f['det']))
    @lru_cache(None)
    def node_power(self,k,d,e):return self.balls[k][d]**e
    @lru_cache(None)
    def defect(self,i,j):
        discrete=sum((w*self.node_power(k,0,i)*self.node_power(k,1,j)for k,(_,_,w)in enumerate(self.balls)),arb(0))
        return self.integral(i,j)-discrete
    def apply(self,p):return sum((T.cmul(c,self.defect(*e))for e,c in p.items()),arb(0))
def direct_basis(corners,lam):
    a=corners[1][0]-corners[0][0];b=corners[2][0]-corners[0][0];c=corners[1][1]-corners[0][1];d=corners[2][1]-corners[0][1];det=a*d-b*c
    dl=[((c-d)/det,(b-a)/det),(d/det,-b/det),(-c/det,a/det)]
    shape=[x*(2*x-1)for x in lam]+[4*lam[i]*lam[j]for i,j in[(0,1),(1,2),(2,0)]]
    grad=[[(4*lam[i]-1)*dl[i][k]for k in range(2)]for i in range(3)]+[[4*(lam[i]*dl[j][k]+lam[j]*dl[i][k])for k in range(2)]for i,j in[(0,1),(1,2),(2,0)]]
    return shape,grad,sum(x*y[0]for x,y in zip(lam,corners))
def correction(target,f,index,functional,drive):
    arrays=target.a;nodes=target.nodes();ys=[nodes[int(k)]for k in arrays['triangles'][index]];n=len(f['free_local']);g=[arb(0)for _ in range(n)];H={(i,j):arb(0)for i in range(n)for j in range(i,n)};basis_checks=0
    for iq,(lam,(xy,w))in enumerate(zip(functional.lambdas,functional.points)):
        N,D,rr=direct_basis(f['corners'],lam);star=[[F(0)for _ in range(12)]for _ in range(5)];old=[[F(0)for _ in range(12)]for _ in range(5)]
        for k in range(6):
            star[0][2*k]=D[k][0];star[1][2*k]=D[k][1];star[2][2*k+1]=D[k][0];star[3][2*k+1]=D[k][1];star[4][2*k]=N[k]/rr
            x=rational(arrays['shape_grad'][index,iq,k,0]);z=rational(arrays['shape_grad'][index,iq,k,1]);old[0][2*k]=x;old[1][2*k]=z;old[2][2*k+1]=x;old[3][2*k+1]=z;old[4][2*k]=rational(arrays['shape'][iq,k])/rational(arrays['qr'][index,iq])
        flat=[v for pair in ys for v in pair];qs=[sum(x*y for x,y in zip(row,flat))for row in star]
        if qs!=[T.evaluate(p,*xy)for p in f['q']]:raise ValueError('independent canonical-point field mismatch')
        for k in range(5):
            for i,local in enumerate(f['free_local']):
                if star[k][local]!=T.evaluate(f['B'][k][i],*xy):raise ValueError('independent canonical-point free basis mismatch')
                basis_checks+=1
        qc=[sum(x*y for x,y in zip(row,flat))for row in old]
        gs,hs,_=C.bulk_local(list(map(C.A,qs)),C.A(1),C.A(5000),C.A(drive));gc,hc,_=C.bulk_local(list(map(C.A,qc)),C.A(1),C.A(5000),C.A(drive))
        for i,li in enumerate(f['free_local']):
            g[i]+=C.A(w)*sum((C.A(star[k][li])*gs[k]-C.A(old[k][li])*gc[k]for k in range(5)),arb(0))
            for j in range(i,n):
                lj=f['free_local'][j];H[i,j]+=C.A(w)*sum((C.A(star[k][li])*hs[k][l]*C.A(star[l][lj])-C.A(old[k][li])*hc[k][l]*C.A(old[l][lj])for k in range(5)for l in range(5)),arb(0))
    return g,H,basis_checks
def bounds(x):return F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
def overlaps(x,p):
    lo,hi=bounds(x);return lo<=F(p['upper_exact'])and F(p['lower_exact'])<=hi
def selftests():
    checks={};q=[T.constant(F(2)),T.constant(F(1,3)),T.constant(F(-1,7)),T.constant(F(3,2)),T.constant(F(6,5))];B=[[T.constant(F(i+1,j+3))for j in range(2)]for i in range(5)]
    ours=independent_jet(q,B);reference=DUAL.bulk_jacobian_jet(q,B)
    checks['cubic_monomial_jet_matches_reviewed_determinant_jet']=ours==reference
    checks['convolution_coefficients_order2']=coefficients(2)==([F(-1),F(5001),F(-7501)],[F(5001),F(-15002),F(27503)])
    corners=[[F(1),F(0)],[F(2),F(0)],[F(1),F(1)]];N,D,r=direct_basis(corners,[F(1,3)]*3)
    checks['direct_canonical_basis_partition']=sum(N)==1 and all(sum(x[k]for x in D)==0 for k in range(2))and r==F(4,3)
    checks['component_source_current']=sha(CF/'scripts/exact_fe_moments.py')=='52bd1e4b54f4ae03ca66592bf4ab295fef30713e5fe171e6fc3e955acf2b9925'
    if not all(checks.values()):raise ValueError(checks)
    return checks
