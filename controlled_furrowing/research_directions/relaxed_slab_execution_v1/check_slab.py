"""Exact identities and interval angle brackets supporting the written slab proof."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import math
import csv
import numpy as np
from flint import arb,ctx

HERE=Path(__file__).resolve().parent


def A(q):return arb(str(q))


def root(q):
    delta=q/4;lo=Q(0);hi=Q(str((A(q)/2).acos().lower().fmpq()))
    def f(x):return (1-A(delta))*A(x).cos()+A(delta)*A(x).sin()/A(x)-A(q) if x else 1-A(q)
    assert f(lo).lower()>0 and f(hi).upper()<0
    for _ in range(150):
        mid=(lo+hi)/2;z=f(mid)
        if z.lower()>0:lo=mid
        elif z.upper()<0:hi=mid
        else:raise ArithmeticError('undecided bisection sign; do not branch on an uncertain ball')
    assert f(lo).lower()>0 and f(hi).upper()<0
    assert A(hi).cos().lower()>A(q/2).upper()
    return lo,hi


def profile(t,alpha,q):
    delta=q/4;a=(1-delta)/2
    pieces=[(a,alpha,0.),(delta/2,alpha,-4*alpha/delta),
            (a,-alpha,0.),(delta/2,-alpha,4*alpha/delta)]
    p=h=0.;remaining=t;theta=alpha
    for length,start,slope in pieces:
        dt=min(max(remaining,0.),length)
        if slope:
            end=start+slope*dt
            p+=(math.sin(end)-math.sin(start))/slope-q*dt
            h+=(-math.cos(end)+math.cos(start))/slope
        else:p+=(math.cos(start)-q)*dt;h+=math.sin(start)*dt
        theta=start+slope*dt;remaining-=length
        if remaining<=0:break
    return p,h,theta


def main():
    ctx.prec=256
    import sympy as s
    l1,l3,c,ep,d,cp,p,h,pp,hp=s.symbols('l1 l3 c ep d cp p h pp hp')
    matrix=s.Matrix([[l1+c*pp,ep/d*cp*p],[c*hp,l3+ep/d*cp*h]])
    rhs=l3*(l1+c*pp)+ep/d*cp*((l1+c*pp)*h-c*hp*p)
    assert s.expand(matrix.det()-rhs)==0
    theta,v,dt=s.symbols('theta v dt',real=True)
    numerator=(v*s.cos(theta))*(v*s.cos(theta)*dt)- (v*s.sin(theta))*(-v*s.sin(theta)*dt)
    assert s.trigsimp(numerator-v*v*dt)==0
    alpha,delta=s.symbols('alpha delta',positive=True)
    w=s.symbols('w',real=True)
    assert s.simplify(delta/(2*alpha)*s.integrate(s.cos(w),(w,-alpha,alpha))-delta*s.sin(alpha)/alpha)==0
    rows=[]
    for q in [Q(1,5),Q(1,2),Q(4,5),Q(19,20)]:
        lo,hi=root(q);C=max(Q(1),18/q)
        rows.append({'compression_ratio_exact':str(q),'delta_exact':str(q/4),
                     'angle_bracket_exact':[str(lo),str(hi)],
                     'angle_decimal_display_only':[float(lo),float(hi)],
                     'mean_equation_lower_positive':True,'mean_equation_upper_negative':True,
                     'cosine_lower_above_q_over_2':True,
                     'collar_C_exact':str(C),'body_J_lower_exact':str(q/4)})
    q=Q(4,5);lo,hi=root(q);C=18/q;M=Q(3,2)
    Fmax2=2+1+(2*M/C)**2+(1+M/C)**2
    jmin=q/4;jmax=1+q/4
    logbound=max((-A(jmin).log()).upper(),A(jmax).log().upper())
    # Illustration uses mu=kappa=a_max=1, v=lambda2=lambda3=1.
    MW=(A(Fmax2)+3+2*logbound)/2+logbound**2/2+A(Fmax2)/2
    Cbulk=2*A(C)*MW
    Cbend=8*A(hi)**2/A(q/4)
    records=[]
    for N in [64,128,256,512,1024,2048,4096]:
        eps=Q(1,N);b=eps**3;depth=C*eps;assert 2*depth<1
        energy=(Cbulk*A(eps)+Cbend*A(b/eps**2)).upper()
        uniform=(A(5).sqrt()*A(eps)).upper()
        records.append({'N':N,'epsilon_exact':str(eps),'b_exact':str(b),'collar_depth_exact':str(depth),
                        'energy_error_upper_exact':str(energy.fmpq()),
                        'uniform_position_error_upper_exact':str(uniform.fmpq()),
                        'energy_display':float(energy),'position_display':float(uniform)})
    samples=[(t,*profile(t,(float(lo)+float(hi))/2,float(q))) for t in np.linspace(0,1,1001)]
    assert abs(samples[-1][1])<1e-13 and abs(samples[-1][2])<1e-13
    with (HERE/'slab_profile.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.writer(f);writer.writerow(['phase','horizontal_correction','vertical_correction','angle']);writer.writerows(samples)
    r={'status':'PASS_SLAB_IDENTITIES_AND_INTERVAL_CONSTANTS','source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
       'written_theorem_sha256':sha256((HERE/'THEOREM.md').read_bytes()).hexdigest(),
       'symbolic_checks':{'collar_determinant_identity':True,'constant_speed_curve_curvature_numerator':True,'angle_transition_mean':True},
       'interval_angle_examples':rows,'illustrative_bound_rows':records,
       'rate_status':'Constructive O(b^(1/3)) upper estimate; no matching lower bound or optimality claim.',
       'scope':'Checks support the explicit written periodic affine-slab proof; no full cylinder or biological conclusion.'}
    (HERE/'slab_verification.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':r['status'],'angle_brackets':len(rows),'symbolic_identities':3,'bound_rows':len(records),
                      'example_angle_display':[float(lo),float(hi)]}))


if __name__=='__main__':main()
