"""Sparse polynomial jets and moment-dual contractions; no integration policy.
All constitutive series, moment defects, tail bounds and admissibility proofs
must be supplied by the caller. This helper never loads or solves a field.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target_safe as T
TARGET_SHA256='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954'
if hashlib.sha256(Path(T.__file__).read_bytes()).hexdigest()!=TARGET_SHA256:raise ValueError('frozen target kernel changed')

def bulk_jacobian_jet(q,B):
    """J, J_u and J_uv; B has5rows and local-free columns, q linear in u."""
    if len(q)!=5 or len(B)!=5 or len({len(row)for row in B})!=1:raise ValueError('matching five-dimensional q/B required')
    a,b,c,d,h=q;n=len(B[0]);det=T.sub(T.mul(a,d),T.mul(b,c));detu=[]
    for k in range(n):detu.append(T.sub(T.add(T.mul(B[0][k],d),T.mul(a,B[3][k])),T.add(T.mul(B[1][k],c),T.mul(b,B[2][k]))))
    gradient=[T.add(T.mul(v,h),T.mul(det,B[4][k]))for k,v in enumerate(detu)];H={}
    for k in range(n):
        for l in range(k,n):
            d2=T.sub(T.add(T.mul(B[0][k],B[3][l]),T.mul(B[0][l],B[3][k])),T.add(T.mul(B[1][k],B[2][l]),T.mul(B[1][l],B[2][k])))
            H[k,l]=T.add(T.mul(d2,h),T.add(T.mul(detu[k],B[4][l]),T.mul(detu[l],B[4][k])))
    return {'value':T.mul(det,h),'gradient':gradient,'hessian_upper':H,'local_dofs':n}

def surface_squared_jacobian_jet(q,B):
    """G=Js², Gu,Guv for q=(r/R,r'/L,z'/L), no sqrt branch selected."""
    if len(q)!=3 or len(B)!=3 or len({len(row)for row in B})!=1:raise ValueError('matching three-dimensional q/B required')
    a,b,c=q;n=len(B[0]);v=T.add(T.mul(b,b),T.mul(c,c));a2=T.mul(a,a);gradient=[];H={}
    for k in range(n):
        vk=T.scale(T.add(T.mul(b,B[1][k]),T.mul(c,B[2][k])),2);ak=T.scale(T.mul(a,B[0][k]),2)
        gradient.append(T.add(T.mul(ak,v),T.mul(a2,vk)))
    for k in range(n):
        for l in range(k,n):
            ak=T.scale(T.mul(a,B[0][k]),2);al=T.scale(T.mul(a,B[0][l]),2)
            vk=T.scale(T.add(T.mul(b,B[1][k]),T.mul(c,B[2][k])),2);vl=T.scale(T.add(T.mul(b,B[1][l]),T.mul(c,B[2][l])),2)
            vv=T.scale(T.add(T.mul(B[1][k],B[1][l]),T.mul(B[2][k],B[2][l])),2);aa=T.scale(T.mul(B[0][k],B[0][l]),2)
            H[k,l]=T.add(T.add(T.mul(aa,v),T.mul(ak,vl)),T.add(T.mul(al,vk),T.mul(a2,vv)))
    return {'value':T.mul(a2,v),'gradient':gradient,'hessian_upper':H,'local_dofs':n}

def shifted_numerator(poly,corners,scaled_radial=True):
    """Exact p=P(x,zeta)/(Rc+h*x)^d, no finite approximation of 1/R.
    zeta is unchanged. The containing radial interval is a numerical conditioning
    coordinate, not a new physical material map or domain enlargement claim.
    """
    radii=[T.rational(v[0])for v in corners];lo,hi=min(radii),max(radii)
    if lo<=0 or hi<=lo:raise ValueError('positive-width off-axis radial interval required')
    rc=(lo+hi)/2;h=(hi-lo)/2 if scaled_radial else F(1);den=max(0,-min((i for i,j in poly),default=0))
    numerator=T.mul(poly,T.power(T.X,den))
    if any(i<0 for i,j in numerator):raise ValueError('numerator failed polynomial conversion')
    shifted=T.compose(numerator,{(0,0):rc,(1,0):h},T.Y)
    return {'numerator':shifted,'denominator_power':den,'Rcenter':rc,'x_scale':h,'Rlo':lo,'Rhi':hi,'variables':('x','zeta'),'identity':'R=Rcenter+x_scale*x; original= numerator(x,zeta)/R**denominator_power'}

def recover_physical_R(representation):
    r=representation;x=T.scale(T.sub(T.X,T.constant(r['Rcenter'])),1/r['x_scale'])
    return T.mul(T.compose(r['numerator'],x,T.Y),T.power(T.X,-r['denominator_power']))

def moment_dual(coefficients,factor_support,moment_defect):
    """d[k]=sum_l coefficients[l]*DeltaMoment(k+l), signed before abs.
    The caller may memoize moment_defect across17 identical geometries/rules.
    """
    result={}
    for i,j in factor_support:
        v=0
        for (k,l),c in coefficients.items():v=T.cadd(v,T.cmul(c,moment_defect(i+k,j+l)))
        result[i,j]=v
    return result

def contract_factor(poly,dual):
    if any(k not in dual for k in poly):raise ValueError('dual lacks factor support')
    result=0
    for k,c in poly.items():result=T.cadd(result,T.cmul(c,dual[k]))
    return result

def contract_volumetric_defect(jet,phi_prime_poly,phi_second_poly,moment_defect):
    """All local gradient/Hessian truncation defects without full products.
    This returns POLYNOMIAL/TRUNCATED contributions only: not remainder bounds,
    full energy, quadrature-cache correction or a certified force/Hessian error.
    """
    g=jet['gradient'];H=jet['hessian_upper'];products={(i,j):T.mul(g[i],g[j])for i,j in H}
    supportA={k for p in list(g)+list(H.values())for k in p};supportB={k for p in products.values()for k in p}
    A=moment_dual(phi_prime_poly,supportA,moment_defect);B=moment_dual(phi_second_poly,supportB,moment_defect)
    return {'gradient':[contract_factor(p,A)for p in g],'hessian_upper':{k:T.cadd(contract_factor(H[k],A),contract_factor(products[k],B))for k in H},'factor_support_A':len(supportA),'factor_support_B':len(supportB),'dual_coefficient_multiply_adds':len(phi_prime_poly)*len(supportA)+len(phi_second_poly)*len(supportB),'scope':'Truncated constitutive polynomial moment defects only; no tail or finite-cache correction included.'}
