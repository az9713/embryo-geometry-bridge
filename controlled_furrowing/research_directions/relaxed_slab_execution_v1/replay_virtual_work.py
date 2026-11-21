"""Separate integer fixed-point interval replay, without the Arb integration code."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as Q
from math import isqrt
import argparse
import importlib.util
import json

HERE=Path(__file__).resolve().parent
MECH=HERE.parents[1]
BITS=96
S=1<<BITS


def ceil_div(a,b):return -((-a)//b)


class I:
    def __init__(self,lo,hi=None):
        self.lo=int(lo);self.hi=int(lo if hi is None else hi);assert self.lo<=self.hi
    @staticmethod
    def q(x):
        x=Q(x);return I((x.numerator*S)//x.denominator,ceil_div(x.numerator*S,x.denominator))
    def __add__(a,b):
        if not isinstance(b,I):b=I.q(b)
        return I(a.lo+b.lo,a.hi+b.hi)
    __radd__=__add__
    def __neg__(a):return I(-a.hi,-a.lo)
    def __sub__(a,b):return a+-b if isinstance(b,I) else a+-I.q(b)
    def __rsub__(a,b):return I.q(b)+-a
    def __mul__(a,b):
        if not isinstance(b,I):b=I.q(b)
        c=[x*y for x in [a.lo,a.hi] for y in [b.lo,b.hi]]
        return I(min(c)//S,ceil_div(max(c),S))
    __rmul__=__mul__
    def __truediv__(a,b):
        if not isinstance(b,I):b=I.q(b)
        assert b.lo>0
        c=[(x*S,y) for x in [a.lo,a.hi] for y in [b.lo,b.hi]]
        return I(min(n//d for n,d in c),max(ceil_div(n,d) for n,d in c))
    def sqrt(a):
        assert a.lo>=0
        return I(isqrt(a.lo*S),isqrt(a.hi*S)+1)
    def positive(a):return I(max(a.lo,0),max(a.hi,0))
    def exact(a):return [str(Q(a.lo,S)),str(Q(a.hi,S))]


def read_cases():
    p=MECH/'runs/nonaxisymmetric_adequacy_audit_v1/surface_stress.py'
    assert sha256(p.read_bytes()).hexdigest()=='f33aaf277eea10fd364a529f95d369551c99b73747f7856cee3af4318b60b3a3'
    spec=importlib.util.spec_from_file_location('_same_archived_input_identity',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cases,_,pins=module.inputs();return cases,pins


def evaluate(c,n):
    ids=list(map(int,c['a']['surface_1_edges'][10]));node=ids[1]
    assert 2*node in set(map(int,c['a']['free']))
    assert sum(node in list(map(int,e)) for side in [0,1] for e in c['a'][f'surface_{side}_edges'])==1
    R0,Rm,R1=[c['reference'][k][0] for k in ids];length=R1-R0
    assert R0>0 and Rm==(R0+R1)/2
    free=set(map(int,c['a']['free']));e=c['posterior'];coords=[]
    for dim in [0,1]:
        row=[]
        for node_id in ids:
            k=2*node_id+dim;delta=e if k in free else Q(0)
            a=I.q(c['center'][k]-delta);b=I.q(c['center'][k]+delta);row.append(I(a.lo,b.hi))
        coords.append(row)
    K=I.q(c['effective_material']['apical_area_modulus']);T=I.q(c['effective_material']['apical_surface_tension'])
    out=I.q(0)
    for j in range(n):
        a=I.q(Q(j,n));b=I.q(Q(j+1,n));t=I(a.lo,b.hi)
        # Recenter the nodal polynomial before interval evaluation. This is the
        # same P2 interpolant, but does not repeatedly multiply full coordinate
        # magnitudes by dependent interval basis functions.
        nmid=4*t*(1-t);dnmid=(4-8*t)/length
        a,m,b=coords[0];delta=m-(a+b)*Q(1,2)
        rv=a+(b-a)*t+delta*nmid
        rr=(b-a)/length+delta*dnmid
        a,m,b=coords[1];delta=m-(a+b)*Q(1,2)
        zr=(b-a)/length+delta*dnmid
        ell=(rr*rr+zr*zr).sqrt();assert ell.lo>0
        J=rv/(I.q(R0)+t*length)*ell
        compression=(K-T-K*J).positive()
        derivative=nmid*ell+rv*rr*dnmid/ell
        out=out+compression*derivative*(length/Q(n))
    return {'case':c['case_id'],'reduced_virtual_work_exact':out.exact(),
            'strict_negative':out.hi<0,'decimal_display_only':[out.lo/S,out.hi/S],
            'subdivisions':n,'fixed_point_bits':BITS,'node':node,'component':'radial'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--subdivisions',type=int,default=4096);a=p.parse_args();assert not a.output.exists()
    cases,pins=read_cases();rows=[evaluate(c,a.subdivisions) for c in cases]
    passed=all(x['strict_negative'] for x in rows)
    r={'status':'PASS_SEPARATE_FIXED_POINT_STRICT_NEGATIVE_REPLAY' if passed else 'UNRESOLVED_SIGN',
       'method':'Integer-directed dyadic rounding and integer square roots; endpoint-linear plus midpoint-bubble P2 representation. No Arb or producer integration import.',
       'shared_input_loader_only':True,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
       'rows':rows,'input_sha256':pins,'scope':'Independent arithmetic implementation, not independent external proof review.'}
    a.output.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':r['status'],'rows':[{k:x[k] for k in ['case','strict_negative','decimal_display_only','subdivisions']} for x in rows]}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
