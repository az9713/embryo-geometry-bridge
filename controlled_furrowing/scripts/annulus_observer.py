"""Exact-rational necessary tests for an unknown horizontal apical-void slice.

No images or experimental outcomes are read. Passing necessary conditions is
NOT validation. Fractions preserve binary floats exactly; strings are decimal.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isfinite,isqrt
from pathlib import Path
import hashlib,json


def rational(value):
    if isinstance(value,float) and not isfinite(value):raise ValueError('finite coordinates/errors required')
    return value if isinstance(value,Q) else Q(value)


def serial(value):
    if isinstance(value,Q):return {'numerator':value.numerator,'denominator':value.denominator,'display':float(value)}
    if isinstance(value,Interval):return {'lower':serial(value.lo),'upper':serial(value.hi),'lower_closed':value.lc,'upper_closed':value.rc}
    if isinstance(value,dict):return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


@dataclass(frozen=True)
class Interval:
    lo:Q
    hi:Q
    lc:bool=True
    rc:bool=True

    def __post_init__(self):
        object.__setattr__(self,'lo',rational(self.lo));object.__setattr__(self,'hi',rational(self.hi))
        if self.lo>self.hi or (self.lo==self.hi and not(self.lc and self.rc)):raise ValueError('empty interval')

    def contains(self,x):
        x=rational(x)
        return (x>self.lo or (self.lc and x==self.lo)) and (x<self.hi or (self.rc and x==self.hi))


def intersection(a,b):
    lo=max(a.lo,b.lo);hi=min(a.hi,b.hi)
    lc=(a.lc if lo==a.lo else True) and (b.lc if lo==b.lo else True)
    rc=(a.rc if hi==a.hi else True) and (b.rc if hi==b.hi else True)
    if lo>hi or (lo==hi and not(lc and rc)):return None
    return Interval(lo,hi,lc,rc)


def _lower(i,bound,closed):
    lo=max(i.lo,bound);lc=i.lc if lo==i.lo else True
    if lo==bound:lc=lc and closed
    if lo>i.hi or (lo==i.hi and not(lc and i.rc)):return None
    return Interval(lo,i.hi,lc,i.rc)


def _upper(i,bound,closed):
    hi=min(i.hi,bound);rc=i.rc if hi==i.hi else True
    if hi==bound:rc=rc and closed
    if hi<i.lo or (hi==i.lo and not(i.lc and rc)):return None
    return Interval(i.lo,hi,i.lc,rc)


def _profile(profile):
    p=[(rational(r),rational(z)) for r,z in profile]
    if len(p)<2 or p[0][0]<0 or any(p[k][0]>=p[k+1][0] for k in range(len(p)-1)):raise ValueError('nonnegative strictly increasing radial grid required')
    return p


def p1_value(profile,r):
    p=_profile(profile);r=rational(r)
    if not p[0][0]<=r<=p[-1][0]:raise ValueError('radius outside profile domain')
    for (r0,z0),(r1,z1) in zip(p[:-1],p[1:]):
        if r0<=r<=r1:return z0+(r-r0)*(z1-z0)/(r1-r0)
    raise AssertionError('unreachable profile evaluation')


def p1_height_range(profile,radius_interval):
    p=_profile(profile);i=Interval(*radius_interval)
    if i.lo<0 or i.lo<p[0][0] or i.hi>p[-1][0]:raise ValueError('radius uncertainty interval not covered by profile')
    values=[p1_value(p,i.lo),p1_value(p,i.hi)]+[z for r,z in p if i.lo<=r<=i.hi]
    return Interval(min(values),max(values))


def _necessary_from_ranges(inner,outer,surface_error,slice_interval=None,occupancy=None):
    error=rational(surface_error)
    if error<0:raise ValueError('nonnegative assumed uniform surface error required')
    common=intersection(Interval(inner.lo-error,inner.hi+error),Interval(outer.lo-error,outer.hi+error))
    if common is not None and slice_interval is not None:common=intersection(common,Interval(*slice_interval))
    # Occupancy entries are conservative bounds: b_true >= lower, a_true <= upper.
    for lower_b,upper_a in occupancy or []:
        if common is None:break
        common=_lower(common,lower_b,False)
        if common is not None:common=_upper(common,upper_a,False)
    return {'status':'NOT_EXCLUDED_BY_NECESSARY_TEST' if common is not None else 'EXCLUDED_UNDER_STATED_ASSUMPTIONS',
            'inner_height_enclosure':inner,'outer_height_enclosure':outer,
            'assumed_uniform_surface_error':error,'possible_shared_height_enclosure':common,
            'logical_scope':'Necessary horizontal apical-void level/optional occupancy constraints only; not sufficiency, biological validation, optical calibration or a PDE error certificate.'}


def p1_necessary_test(apical,inner_radius_interval,outer_radius_interval,surface_error=0,
                      basal=None,occupancy_radii=None,slice_interval=None):
    a=_profile(apical);inner=Interval(*inner_radius_interval);outer=Interval(*outer_radius_interval)
    if inner.lo>=outer.hi:return {'status':'EXCLUDED_UNDER_STATED_ASSUMPTIONS','reason':'No positive-width ordered annulus is possible in supplied radius intervals.'}
    occupancy=[]
    if occupancy_radii is not None:
        if basal is None:raise ValueError('basal profile required for tissue occupancy witnesses')
        b=_profile(basal);error=rational(surface_error)
        for r in occupancy_radii:occupancy.append((p1_value(b,r)-error,p1_value(a,r)+error))
    return _necessary_from_ranges(p1_height_range(a,(inner.lo,inner.hi)),p1_height_range(a,(outer.lo,outer.hi)),surface_error,slice_interval,occupancy)


def p1_pair_branches(apical,basal=None,inner_radius_interval=None,outer_radius_interval=None):
    """Enumerate exact affine-radius branches for a section with one void annulus.

    Central disk and outer tissue ring are required; the annulus is absent
    because the apical graph is at/below the plane. A basal profile, if supplied,
    must lie below the apical profile and below the entire admissible slice.
    """
    a=_profile(apical)
    if a[0][0]!=0:raise ValueError('full central-disk enumeration requires radius0')
    b=_profile(basal) if basal is not None else None
    if b is not None:
        if b[0][0]!=0 or b[-1][0]!=a[-1][0]:raise ValueError('matching full basal/apical radial domains required')
        for r in sorted({r for r,_ in a}|{r for r,_ in b}):
            if p1_value(a,r)<=p1_value(b,r):raise ValueError('positive graph thickness required')
    branches=[]
    for i,((r0,z0),(r1,z1)) in enumerate(zip(a[:-1],a[1:])):
        if z1>=z0:continue
        for j in range(i+1,len(a)-1):
            s0,w0=a[j];s1,w1=a[j+1]
            if w1<=w0:continue
            ri=intersection(Interval(r0,r1),Interval(*(inner_radius_interval or (r0,r1))))
            ro=intersection(Interval(s0,s1),Interval(*(outer_radius_interval or (s0,s1))))
            if ri is None or ro is None:continue
            ih=p1_height_range(a,(ri.lo,ri.hi));oh=p1_height_range(a,(ro.lo,ro.hi));level=intersection(ih,oh)
            if level is None:continue
            level=_lower(level,max(z for _,z in a[i+1:j+1]),True)
            if level is None:continue
            level=_upper(level,min(z for _,z in a[:i+1]),False)
            if level is None:continue
            level=_upper(level,min(z for _,z in a[j+1:]),False)
            if level is None:continue
            if b is not None:
                level=_lower(level,max(z for _,z in b),False)
                if level is None:continue
            ai=(r1-r0)/(z1-z0);bi=r0-ai*z0
            ao=(s1-s0)/(w1-w0);bo=s0-ao*w0
            slope=ao-ai;offset=bo-bi
            if slope>0:level=_lower(level,-offset/slope,False)
            elif slope<0:level=_upper(level,-offset/slope,False)
            elif offset<=0:level=None
            if level is None:continue
            branches.append({'inner_edge':i,'outer_edge':j,'slice_height':level,
                             'inner_radius_affine':(ai,bi),'outer_radius_affine':(ao,bo),
                             'area_inner_polynomial':(ai*ai,2*ai*bi,bi*bi),
                             'area_outer_polynomial':(ao*ao,2*ao*bo,bo*bo)})
    return {'branches':branches,'scope':'Exact rational P1 one-apical-void-annulus graph enumeration; areas are pi times the reported quadratic polynomials.',
            'basal_occupancy_verified':b is not None,'duplicates_at_shared_edge_heights_allowed':True}


def _quadratic_controls(controls):
    curves=[];previous=None
    for segment in controls:
        if len(segment)!=3:raise ValueError('quadratic controls must be [start,mid,end]')
        q=[tuple(map(rational,p)) for p in segment]
        if previous is not None and q[0]!=previous:raise ValueError('surface segments must share exact endpoints')
        c0=q[0];c1=tuple(-3*q[0][k]+4*q[1][k]-q[2][k] for k in range(2));c2=tuple(2*q[0][k]-4*q[1][k]+2*q[2][k] for k in range(2))
        if min(c1[0],c1[0]+2*c2[0])<0 or q[2][0]<=q[0][0] or q[0][0]<0:raise ValueError('radially monotone quadratic surface required')
        curves.append((q,c0,c1,c2));previous=q[-1]
    if not curves:raise ValueError('nonempty quadratic surface required')
    return curves


def _poly(c0,c1,c2,t):return c0+c1*t+c2*t*t


def _root_bracket(c0,c1,c2,radius,steps):
    left,right=Q(0),Q(1)
    if radius==c0:return left,left
    if radius==c0+c1+c2:return right,right
    for _ in range(steps):
        mid=(left+right)/2;value=_poly(c0,c1,c2,mid)
        if value==radius:return mid,mid
        if value<radius:left=mid
        else:right=mid
    return left,right


def p2_height_range(controls,radius_interval,bisection_steps=64):
    """Rigorous rational enclosure on a physical-radius interval, no float roots."""
    if not isinstance(bisection_steps,int) or bisection_steps<1:raise ValueError('positive integer root-bisection steps required')
    curves=_quadratic_controls(controls);r=Interval(*radius_interval)
    if r.lo<curves[0][0][0][0] or r.hi>curves[-1][0][2][0]:raise ValueError('radius uncertainty interval not covered by quadratic surface')
    lows=[];highs=[];brackets=[]
    for q,c0,c1,c2 in curves:
        cut=intersection(r,Interval(q[0][0],q[2][0]))
        if cut is None:continue
        lower=_root_bracket(c0[0],c1[0],c2[0],cut.lo,bisection_steps)
        upper=_root_bracket(c0[0],c1[0],c2[0],cut.hi,bisection_steps)
        lo,hi=lower[0],upper[1]
        values=[_poly(c0[1],c1[1],c2[1],lo),_poly(c0[1],c1[1],c2[1],hi)]
        if c2[1]!=0:
            critical=-c1[1]/(2*c2[1])
            if lo<=critical<=hi:values.append(_poly(c0[1],c1[1],c2[1],critical))
        lows.append(min(values));highs.append(max(values));brackets.append({'radius':cut,'lower_root':lower,'upper_root':upper,'enlarged_parameter_interval':Interval(lo,hi)})
    return {'height':Interval(min(lows),max(highs)),'root_brackets':brackets,'scope':'Exact rational enclosure of supplied binary/decimal control geometry, not model-to-biology/PDE/optical error.'}


def p2_necessary_test(apical_controls,inner_radius_interval,outer_radius_interval,surface_error=0,
                      basal_controls=None,occupancy_radii=None,slice_interval=None,bisection_steps=64):
    inner=Interval(*inner_radius_interval);outer=Interval(*outer_radius_interval)
    if inner.lo>=outer.hi:return {'status':'EXCLUDED_UNDER_STATED_ASSUMPTIONS','reason':'No positive-width ordered annulus possible.'}
    ir=p2_height_range(apical_controls,(inner.lo,inner.hi),bisection_steps)
    outer_result=p2_height_range(apical_controls,(outer.lo,outer.hi),bisection_steps)
    occupancy=[]
    if occupancy_radii is not None:
        if basal_controls is None:raise ValueError('basal controls required for occupancy witnesses')
        error=rational(surface_error)
        for r in occupancy_radii:
            ah=p2_height_range(apical_controls,(r,r),bisection_steps)['height'];bh=p2_height_range(basal_controls,(r,r),bisection_steps)['height']
            occupancy.append((bh.lo-error,ah.hi+error))
    result=_necessary_from_ranges(ir['height'],outer_result['height'],surface_error,slice_interval,occupancy)
    result.update({'inner_root_brackets':ir['root_brackets'],'outer_root_brackets':outer_result['root_brackets']})
    return result


def pi_enclosure():
    def atan_bounds(den,n):
        x=Q(1,den);s=sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)),Q(0));term=(-1)**n*x**(2*n+1)/Q(2*n+1)
        return min(s,s+term),max(s,s+term)
    a,b=atan_bounds(5,40),atan_bounds(239,12)
    return Interval(16*a[0]-4*b[1],16*a[1]-4*b[0])


def _sqrt_enclosure(q,bits):
    if q<0:raise ValueError('nonnegative square-root argument required')
    scale=1<<bits;n=isqrt(q.numerator*scale*scale//q.denominator)
    lower=Q(n,scale);upper=lower if lower*lower==q else Q(n+1,scale)
    return lower,upper


def area_equivalent_radius_interval(area_interval,bits=80):
    area=Interval(*area_interval)
    if area.lo<0 or not isinstance(bits,int) or bits<1:raise ValueError('nonnegative area and positive precision required')
    pi=pi_enclosure();lower=_sqrt_enclosure(area.lo/pi.hi,bits)[0];upper=_sqrt_enclosure(area.hi/pi.lo,bits)[1]
    return {'radius':Interval(lower,upper),'pi_enclosure':pi,'area':area,
            'interpretation':'Translation-invariant equal-area circle radius; scalar shape relaxation, not evidence the measured contour is circular or an anatomical occupancy contour.'}


def _self_test():
    top=[(0,2),(1,1),(2,0),(3,1),(4,2)];basal=[(0,-1),(4,-1)];checks={}
    ok=p1_necessary_test(top,('0.9','1.1'),('2.9','3.1'),basal=basal,occupancy_radii=[0,4])
    bad=p1_necessary_test(top,('0.4','0.6'),('2.9','3.1'))
    rescued=p1_necessary_test(top,('0.4','0.6'),('2.9','3.1'),surface_error='0.2')
    checks['P1_known_pair_not_excluded']=ok['status'].startswith('NOT_EXCLUDED')
    checks['P1_wrong_pair_excluded']=bad['status'].startswith('EXCLUDED')
    checks['larger_assumed_error_relaxes_exclusion']=rescued['status'].startswith('NOT_EXCLUDED')
    checks['occupancy_strict_equality_excluded']=p1_necessary_test(top,(0,0),(4,4),basal=basal,occupancy_radii=[0,4])['status'].startswith('EXCLUDED')
    branches=p1_pair_branches(top,basal)['branches']
    for sigma in [Q(1,4),Q(3,4),Q(1),Q(3,2)]:
        possible=[b for b in branches if b['slice_height'].contains(sigma)]
        checks['P1_branch_covers_'+str(sigma)]=bool(possible)
        for b in possible:
            ri=b['inner_radius_affine'][0]*sigma+b['inner_radius_affine'][1];ro=b['outer_radius_affine'][0]*sigma+b['outer_radius_affine'][1]
            checks['P1_equal_levels_'+str(sigma)]=p1_value(top,ri)==sigma==p1_value(top,ro) and ri<ro
    checks['P1_zero_width_and_empty_centre_excluded']=not any(b['slice_height'].contains(Q(0)) or b['slice_height'].contains(Q(2)) for b in branches)
    checks['P1_plateau_critical_height_retained']=any(b['slice_height'].contains(Q(1)) for b in p1_pair_branches([(0,2),(1,1),(2,1),(3,2)],[(0,0),(3,0)])['branches'])
    checks['P1_two_annuli_not_mistaken_for_one']=not p1_pair_branches([(0,2),(1,1),(2,2),(3,1),(4,2)],basal)['branches']
    controls=[[(0,4),(2,0),(4,4)]]
    checks['P2_known_pair_not_excluded']=p2_necessary_test(controls,(1,1),(3,3))['status'].startswith('NOT_EXCLUDED')
    checks['P2_wrong_pair_excluded']=p2_necessary_test(controls,(1,1),('3.5','3.5'))['status'].startswith('EXCLUDED')
    curved=[[(0,1),(Q(3,2),0),(4,1)]]
    checks['P2_nonlinear_radial_map_equal_levels']=p2_necessary_test(curved,(Q(5,8),Q(5,8)),(Q(21,8),Q(21,8)))['status'].startswith('NOT_EXCLUDED')
    for k in range(1,40):
        t=Q(k,40);radius=2*t+2*t*t;height=(2*t-1)**2
        bound=p2_height_range(curved,(radius,radius))
        checks['P2_rational_root_enclosure_'+str(k)]=bound['height'].contains(height)
    checks['area_zero_radius_exact']=area_equivalent_radius_interval((0,0))['radius']==Interval(0,0)
    radius=area_equivalent_radius_interval((3,4))['radius'];pi=pi_enclosure()
    checks['area_radius_bound_algebra']=radius.lo**2*pi.hi<=3 and radius.hi**2*pi.lo>=4
    checks['pi_Machin_known_decimal_enclosed']=pi.lo<Q('3.14159265358979323846264338327950288419716939937510582097495')<pi.hi
    shifted=[(r,z+7) for r,z in top]
    checks['vertical_translation_preserves_decision']=p1_necessary_test(shifted,('0.4','0.6'),('2.9','3.1'))['status']==bad['status']
    rejected=False
    try:p2_height_range([[(0,1),(3,0),(2,1)]],(0,1))
    except ValueError:rejected=True
    checks['nonmonotone_surface_rejected']=rejected
    root=Path(__file__).resolve().parents[1]
    receipt={'status':'PASS_SYNTHETIC_COMPONENT' if all(checks.values()) else 'FAIL','checks':checks,'count':len(checks),'examples':{'P1_correct':ok,'P1_wrong':bad,'P1_error_relaxation':rescued,'P1_pair_branches':p1_pair_branches(top,basal),'area_radius':area_equivalent_radius_interval((3,4))},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'evidence_limit':'No measured Movie3 anatomical radius intervals are supplied. Tests are synthetic geometry only; no target-outcome fit or biological exclusion.'}
    (root/'evidence/annulus_observer_verification.json').write_text(json.dumps(serial(receipt),indent=2)+'\n');print(json.dumps({'status':receipt['status'],'count':len(checks),'failures':[k for k,v in checks.items() if not v]},indent=2));return all(checks.values())


if __name__=='__main__':raise SystemExit(0 if _self_test() else 1)
