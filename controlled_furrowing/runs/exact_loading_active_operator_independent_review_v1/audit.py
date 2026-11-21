"""Cold independent source and toy algebra readback of the exact loading active operator."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import itertools
import json
from math import comb
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PROD=HERE.parent/"exact_loading_active_operator_preparation_v1"
PINS={"source_sha256":"440c5ee10e1eba29565c664a64d80af88f7a43f2c09e2c99c74848f8d2d6469a",
      "protocol_sha256":"e897ab039e5109eb454922f41d1318df9023109bba5adaa1bf71ed55da3354a2",
      "preflight_sha256":"ffac861a067dc8e877ab8e6c741a729cc3ccb5c60aef54a83917acfd730f027f"}
KL=(ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/negative_KL_active_energy_v1/active_energy.py")
KL_SHA="072ed83a30565be34bd5fe7dde559526e089a070df0ad95f0e61611db9e9e76d"


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def module(p,h,name):
    assert sha(p)==h
    spec=spec_from_file_location(name,p)
    m=module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def docs(p):
    return json.loads(p.read_bytes())


def reject(f,fail=(KeyError,ValueError,TypeError)):
    try:
        f()
    except fail:
        return
    raise AssertionError("bad fixture admitted")


def own_gram(points):
    """Independent Cartesian P2 interpolation, power-basis iterated integral."""
    import sympy as S
    r,z,s,t=S.symbols("R Z s t")
    x=[S.Matrix([S.Rational(a.numerator,a.denominator) for a in q]) for q in points]
    middle=[(x[i]+x[j])/2 for i,j in ((0,1),(1,2),(2,0))]
    basis=(1,r,z,r*r,r*z,z*z)
    M=S.Matrix([[S.sympify(m).subs({r:a,z:b}) for m in basis] for a,b in x+middle]).inv()
    grad=S.Matrix([[S.diff(m,z) for m in basis]])*M
    D=S.Matrix.hstack(x[1]-x[0],x[2]-x[0])
    d=D.det()
    pos=x[0]+D*S.Matrix([s,t])
    p=grad.subs({r:pos[0],z:pos[1]},simultaneous=True)
    def integrate(expr):
        value=Q(0)
        for (i,j),coeff in S.Poly(expr,s,t).terms():
            # int_0^1 s^i * ((1-s)^(j+1)/(j+1)) ds, by binomial expansion.
            integral=sum(Q((-1)**k*comb(j+1,k),(j+1)*(i+k+1)) for k in range(j+2))
            value+=Q(str(coeff))*integral
        return value
    return tuple(tuple(integrate(d*pos[0]*u*v) for v in p) for u in p)


def fixture_assembly(A):
    r=[(0,0),(2,0),(2,3),(0,3)]
    edges=[(0,1),(1,2),(0,2),(2,3),(0,3)]
    ref=[tuple(map(Q,t)) for t in r]+[tuple(Q(r[a][k]+r[b][k],2) for k in (0,1)) for a,b in edges]
    tri=[[0,1,2,4,5,6],[0,2,3,6,7,8]]
    # Make one stored binary-like fixed midpoint deliberately noncanonical.
    ref[5]=(ref[5][0],ref[5][1]+Q(1,128))
    free=[1,2,3,6,7,8,9,12,13,14,15,17]
    own=[own_gram([ref[j] for j in cells[:3]]) for cells in tri]
    prod=[A.active_gram([ref[j] for j in cells[:3]]) for cells in tri]
    assert prod==own and all(own[k][i][j]==own[k][j][i] for k,i,j in itertools.product(range(2),range(6),range(6)))
    # Full all-node two-component reference Gram before projecting.
    full=[[Q(0) for unused in range(18)] for _ in range(18)]
    for cells,G in zip(tri,own,strict=True):
        for a,b in itertools.product(range(6),repeat=2):
            for c in (0,1):
                full[2*cells[a]+c][2*cells[b]+c]+=G[a][b]
    A0=[[full[i][j] for j in free] for i in free]
    X=[val for row in ref for val in row]
    off=[sum(x*y for x,y in zip(full[i],X,strict=True)) for i in free]
    rows,b,leaves=A.scatter_active(ref,tri,free,[0,1])
    assert b==off and [dict((j,x) for j,x in enumerate(row) if x) for row in A0]==rows
    assert b!=[sum(a*X[j] for a,j in zip(row,free,strict=True)) for row in A0]
    assert len(leaves)==2 and all(leaf["axial_gram_exact"]==own[i] for i,leaf in enumerate(leaves))
    n=len(free)
    R=[[Q(((i+3)*(j+7))%15-7,8) for j in range(n)] for i in range(n)]
    u=[Q(i-5,16) for i in range(n)]
    v=[Q((i*7)%23-10,13) for i in range(n)]
    RA=[[sum(R[i][k]*A0[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Rb=[sum(c*d for c,d in zip(rr,b,strict=True)) for rr in R]
    dot=lambda p,q:sum(a*b for a,b in zip(p,q,strict=True))
    exact_result=([dot(w,u)+h for w,h in zip(RA,Rb,strict=True)],
                  [dot(w,v) for w in RA],[sum(map(abs,w)) for w in RA],
                  [sum(map(abs,w)) for w in R])
    assert A.multiply_signed(rows,b,R,u,v)==exact_result
    assert any(sum(abs(R[i][k]*A0[k][j]) for k in range(n)) > abs(RA[i][j])
               for i,j in itertools.product(range(n),repeat=2))
    return n


def machinarctan(x, n=180):
    z=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(n))
    rest=(-1)**n*x**(2*n+1)/(2*n+1)
    return min(z,z+rest),max(z,z+rest)


def scale_extension(A):
    a,b=machinarctan(Q(1,5))
    c,d=machinarctan(Q(1,239))
    lo,hi=16*a-4*d,16*b-4*c
    g,s=Q(1234567,891011),Q(8101,32768)
    I=A.directed_factor(g,s)
    assert 0<I[0]<2*g/s*lo<=2*g/s*hi<I[1]
    for q in (Q(-3,8),Q(0),Q(7,11)):
        assert A.scale_signed(q,I)==(min(q*I[0],q*I[1]),max(q*I[0],q*I[1]))
    for B,M in ((Q(100),Q(60)),(Q(1,100),Q(0)),(Q(0),Q(0))):
        eta,z,r,Rnorm,gap=Q(1,40),Q(1,4),Q(1,10),Q(3),Q(2,3)
        v=A.old_K_extension(eta,z,r,Rnorm,gap,B,M)
        h=v["positive_dose_increment_exact"]
        e,m=eta+h*B,z+h*M
        margin=(1-z)*r-eta
        assert h>0 and h<=gap/2 and h*(B+r*M)<=margin/2
        assert v["strict_selfmap_margin_lower_exact"]==r-e-m*r>=margin/2>0
        assert v["eta_upper_exact"]==e and v["Z_upper_exact"]==m
        assert v["posterior_upper_exact"]==e/(1-m)<r and v["normalized_H_curvature_lower_exact"]==(1-m)/Rnorm
        assert v["reaches_L"] is v["whole_KL_path_certified"] is False
    reject(lambda:A.old_K_extension(Q(1,2),Q(1,2),Q(1),Q(1),Q(1),Q(0),Q(0)))

