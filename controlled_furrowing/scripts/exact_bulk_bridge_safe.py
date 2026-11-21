"""Validated elementwise exact-integral minus binary-cache force/H pilot.

SAFE V2: exact-zero-only sparse interval arithmetic. POINT operator only; no standalone preparation/run or uniform root bridge.
Positive-J scope is supplied by the separately reviewed exact body geometry.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import sys,json,hashlib,datetime,argparse,time
import numpy as np
from flint import arb,arb_mat,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'runs/exact_integration_bridge_v1/safe_kernel_v2'))
import exact_fe_target_safe as T
import exact_fe_moments as M
import polynomial_jet_safe as JET
import algebraic_root_feasibility as core
OUT=ROOT/'runs/exact_integration_bridge_v1/bulk_operator_safe_v2'
SOURCE=Path(__file__)
RECORD='runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json'
PREP='runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json'
GEOMETRY='runs/exact_integration_bridge_v1/independent_design/geometry_alignment.json'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def fp(p):return {k:core.A(v) for k,v in p.items()}
def inflated(value,tail):return value+(-tail).union(tail)
def pnorm(p,fields):
    """Exact coefficient majorant on a containing chart rectangle."""
    if fields['chart']['kind']=='off_axis':
        lo=min(v[0] for v in fields['corners']);hi=max(v[0] for v in fields['corners'])
        z=max(abs(v[1]-fields['zcenter']) for v in fields['corners'])
        return sum((abs(v)*(hi**i if i>=0 else lo**i)*z**j for (i,j),v in p.items()),F(0))
    return sum((abs(v) for v in p.values()),F(0))

class Functional:
    def __init__(self,target,fields,index,rule):
        self.t=target;self.f=fields;self.index=index
        self.vertices=tuple(tuple(v) for v in fields['corners'])
        self.points=[];self.requested_defect_keys=set()
        zero=[i for i,v in enumerate(self.vertices) if v[0]==0];off=[i for i in range(3) if i not in zero]
        for lam,w in zip(rule['volume_barycentric'],rule['volume_weights'][index]):
            if fields['chart']['kind']=='axis_vertex':
                s=lam[off[0]]+lam[off[1]];pt=(s,lam[off[1]]/s)
            elif fields['chart']['kind']=='axis_edge':
                s=lam[off[0]];pt=(s,lam[zero[1]]/(1-s))
            else:
                pt=(sum(l*v[0] for l,v in zip(lam,self.vertices)),sum(l*v[1] for l,v in zip(lam,self.vertices))-fields['zcenter'])
            self.points.append((pt,w))
        self.mass_integral=self.integral(0,0);self.mass_rule=core.A(sum(w for pt,w in self.points))
        if not self.mass_integral>0 or not self.mass_rule>0:raise ValueError('positive measures required')
        self.point_balls=[(core.A(pt[0]),core.A(pt[1]),core.A(w)) for pt,w in self.points]
    @lru_cache(maxsize=None)
    def integral(self,a,b):
        f=self.f;c=f['chart']
        if c['kind']=='off_axis':return M.normalized_bulk_moment(self.vertices,a,b,self.t.energy_scale,f['zcenter'],ctx.prec)
        return self.t.normalization()*core.A(M.axis_chart_moment(c['kind'],a,b,c['R1'],f['det']))
    @lru_cache(maxsize=None)
    def defect(self,a,b):
        # Outward balls avoid an enormous rational LCM for powers of distinct
        # positive-radius nodes. The actual rational nodes/weights are enclosed,
        # never replaced by regenerated Gaussian abscissae.
        self.requested_defect_keys.add((a,b))
        q=sum((w*self.node_power(i,0,a)*self.node_power(i,1,b) for i,(_,_,w) in enumerate(self.point_balls)),arb(0))
        return self.integral(a,b)-q
    @lru_cache(maxsize=None)
    def node_power(self,index,axis,exponent):
        return self.point_balls[index][axis]**exponent
    def apply(self,p):
        return sum((T.cmul(v,self.defect(*k)) for k,v in p.items()),arb(0))

def series(e,order):
    c1,c2=M.volumetric_derivative_coefficients(order,1,5000)
    p=T.ONE;A={};B={};e=fp(e)
    for n,(a,b) in enumerate(zip(c1,c2)):
        A=T.add(A,T.scale(p,core.A(a)));B=T.add(B,T.scale(p,core.A(b)))
        if n<order:p=T.mul(p,e)
    return A,B

def cache_correction(target,f,index,functional,drive):
    """Exact coherent auxiliary minus old binary strain-cache, at same POINT."""
    a=target.a;ids=a['triangles'][index];y=target.nodes();coords=[core.A(y[int(j)][c]) for j in ids for c in (0,1)]
    n=len(f['free_local']);outg=[arb(0) for _ in range(n)];outH=arb_mat(n,n)
    for iq,(pt,w) in enumerate(functional.points):
        Bstar=arb_mat(5,n,[core.A(T.evaluate(p,*pt)) for row in f['B'] for p in row])
        qstar=[core.A(T.evaluate(p,*pt)) for p in f['q']]
        g,h,_=core.bulk_local(qstar,core.A(1),core.A(5000),core.A(drive))
        gg=Bstar.transpose()*arb_mat(5,1,g);HH=Bstar.transpose()*arb_mat(5,5,[v for row in h for v in row])*Bstar
        cacheB=[[arb(0) for _ in range(12)] for _ in range(5)]
        for j in range(6):
            gx=core.A(float(a['shape_grad'][index,iq,j,0]));gz=core.A(float(a['shape_grad'][index,iq,j,1]))
            cacheB[0][2*j]=gx;cacheB[1][2*j]=gz;cacheB[2][2*j+1]=gx;cacheB[3][2*j+1]=gz
            cacheB[4][2*j]=core.A(float(a['shape'][iq,j]))/core.A(float(a['qr'][index,iq]))
        qc=[sum((v*z for v,z in zip(row,coords)),arb(0)) for row in cacheB]
        Bs=arb_mat(5,n,[row[j] for row in cacheB for j in f['free_local']])
        gc,hc,_=core.bulk_local(qc,core.A(1),core.A(5000),core.A(drive))
        gg-=Bs.transpose()*arb_mat(5,1,gc);HH-=Bs.transpose()*arb_mat(5,5,[v for row in hc for v in row])*Bs
        for i in range(n):
            outg[i]+=core.A(w)*gg[i,0]
            for j in range(n):outH[i,j]+=core.A(w)*HH[i,j]
    return outg,outH

def run_element(target,index,rule,order,rho,drive):
    start=time.monotonic();f=target.element(index);jet=JET.bulk_jacobian_jet(f['q'],f['B']);n=jet['local_dofs']
    fun=Functional(target,f,index,rule);e=T.sub(jet['value'],T.ONE)
    A,B=series(e,order)
    print(json.dumps({'phase':'series','element':index,'A_terms':len(A),'B_terms':len(B),'seconds':time.monotonic()-start}),flush=True)
    vol=JET.contract_volumetric_defect(jet,A,B,fun.defect)
    tail1,tail2=M.volumetric_derivative_tail_bounds(order,rho,1,5000)
    mass=core.upper(fun.mass_integral+fun.mass_rule)
    normJu=[pnorm(p,f) for p in jet['gradient']]
    diag=[F(1),F(1)+drive,F(1),F(1)+drive,F(1)]
    Fpoly=[];Hpoly={};Ftail=[];Htail={}
    for i in range(n):
        p={}
        for k in range(5):p=T.add(p,T.scale(T.mul(f['q'][k],f['B'][k][i]),diag[k]))
        Fpoly.append(fun.apply(p)+vol['gradient'][i]);Ftail.append(mass*core.A(tail1*normJu[i]))
    for i,j in jet['hessian_upper']:
        p={}
        for k in range(5):p=T.add(p,T.scale(T.mul(f['B'][k][i],f['B'][k][j]),diag[k]))
        Hpoly[i,j]=fun.apply(p)+vol['hessian_upper'][i,j]
        Htail[i,j]=mass*core.A(tail2*normJu[i]*normJu[j]+tail1*pnorm(jet['hessian_upper'][i,j],f))
    corrF,corrH=cache_correction(target,f,index,fun,drive)
    fullF=[inflated(v,Ftail[i])+corrF[i] for i,v in enumerate(Fpoly)]
    fullH={(i,j):inflated(v,Htail[i,j])+corrH[i,j] for (i,j),v in Hpoly.items()}
    return {'element':index,'chart':f['chart']['kind'],'active':str(f['active']),'global_indices':f['global_indices'],'local_dofs':n,
        'truncation_order':order,'rho_exact':str(rho),'A_support':len(A),'B_support':len(B),
        'dual_statistics':{k:v for k,v in vol.items() if k.startswith('factor_') or k.startswith('dual_')},
        'unique_moment_defects':len(fun.requested_defect_keys),
        'force_integral_minus_auxiliary_polynomial':[core.pack(v) for v in Fpoly],
        'force_integral_minus_auxiliary_tail_upper':[core.pack(v) for v in Ftail],
        'force_auxiliary_minus_cache':[core.pack(v) for v in corrF],
        'force_integral_minus_cache':[core.pack(v) for v in fullF],
        'H_integral_minus_cache_upper_triangle':[{'i':i,'j':j,'enclosure':core.pack(v),'tail_upper':core.pack(Htail[i,j]),'cache_correction':core.pack(corrH[i,j])} for (i,j),v in fullH.items()],
        'max_force_defect_upper':core.pack(max(core.upper(v) for v in fullF)),
        'max_H_entry_defect_upper':core.pack(max(core.upper(v) for v in fullH.values())),
        'force_max_tail_upper':core.pack(max(core.upper(v) for v in Ftail)),
        'H_max_tail_upper':core.pack(max(core.upper(v) for v in Htail.values())),
        'seconds':time.monotonic()-start,
        'scope':'Pointwise element contribution ONLY; exact integral minus frozen binary Q-cache includes two-mass series remainders and finite strain-cache correction. No full assembled force, surface contribution, tube-uniform bound or root transfer.'}
