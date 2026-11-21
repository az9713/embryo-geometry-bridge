"""Whole-body exact geometry for a P2 nodal interval box, including the axis.

This extends fixed-coefficient geometric checks to every field in a supplied
coordinate box. It does not certify that a root lies in that box; an actual
independently reviewed root receipt is required for the production application.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import argparse,base64,hashlib,io,json
import numpy as np
from mixed_geometry_certificate_v2 import strong_monotonicity_certificate

ROOT=Path(__file__).resolve().parents[1]
INDICES=((2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1))
VERTICES=((Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1)))


def q(v):return v if isinstance(v,Q) else Q(v)


@dataclass(frozen=True)
class I:
    lo:Q
    hi:Q
    def __post_init__(self):
        object.__setattr__(self,'lo',q(self.lo));object.__setattr__(self,'hi',q(self.hi))
        if self.lo>self.hi:raise ValueError('ordered exact interval required')
    def __add__(self,b):
        b=interval(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-interval(b)
    def __rsub__(self,b):return interval(b)+-self
    def __mul__(self,b):
        b=interval(b);v=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)];return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=interval(b)
        if b.lo<=0<=b.hi:raise ValueError('division interval contains zero')
        return self*I(1/b.hi,1/b.lo)
    def square(self):
        low=0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2)
        return I(low,max(self.lo**2,self.hi**2))


def interval(v):return v if isinstance(v,I) else I(q(v),q(v))
def packed(v):
    v=interval(v);return {'lower_exact':str(v.lo),'upper_exact':str(v.hi)}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def multinomial(alpha):
    value=factorial(sum(alpha))
    for a in alpha:value//=factorial(a)
    return value


def nodal_derivatives(xi,eta):
    return ((-3+4*xi+4*eta,-3+4*xi+4*eta),(4*xi-1,0),(0,4*eta-1),
            (4-8*xi-4*eta,-4*xi),(4*eta,4*xi),(-4*eta,4-4*xi-8*eta))


def triangle_box(reference,current):
    x=[[q(v) for v in row] for row in reference]
    y=[[interval(v) for v in row] for row in current]
    if len(x)!=3 or len(y)!=6 or any(len(row)!=2 for row in x+y):raise ValueError('one affine reference/P2 current triangle required')
    a=x[1][0]-x[0][0];b=x[2][0]-x[0][0];c=x[1][1]-x[0][1];d=x[2][1]-x[0][1]
    area=a*d-b*c
    if area<=0 or min(row[0] for row in x)<0:raise ValueError('positive reference orientation and nonnegative radii required')
    inverse=((d/area,-b/area),(-c/area,a/area));G=[];monotone=[]
    for xi,eta in VERTICES:
        derivatives=nodal_derivatives(xi,eta)
        matrix=[[sum((y[k][i]*derivatives[k][j] for k in range(6)),I(0,0)) for j in (0,1)] for i in (0,1)]
        G.append(matrix)
        F=[[sum((matrix[i][k]*inverse[k][j] for k in (0,1)),I(0,0)) for j in (0,1)] for i in (0,1)]
        aa=F[0][0];bb=(F[0][1]+F[1][0])/2;dd=F[1][1];det=aa*dd-bb.square()
        good=aa.lo>0 and det.lo>0 and aa.hi+dd.hi>0
        monotone.append({'symmetric_gradient':[packed(aa),packed(bb),packed(dd)],
                         'positive_definite_for_entire_box':good,
                         'lambda_lower_exact':str(det.lo/(aa.hi+dd.hi)) if good else None})
    # Affine derivative fields have degree-one Bernstein controls at vertices.
    # Their determinant is quadratic; form degree-two coefficients directly.
    determinant=[]
    for i,j in ((0,0),(1,1),(2,2),(0,1),(1,2),(0,2)):
        value=(G[i][0][0]*G[j][1][1]+G[j][0][0]*G[i][1][1]
               -G[i][0][1]*G[j][1][0]-G[j][0][1]*G[i][1][0])/(2*area)
        determinant.append(value)
    # Direct nodal-to-Bernstein conversion preserves exactly fixed axis zeros.
    radial=[y[i][0] for i in range(3)]+[2*y[m][0]-(y[i][0]+y[j][0])/2 for i,j,m in ((0,1,3),(1,2,4),(2,0,5))]
    numerator={}
    for alpha,A in zip(INDICES,determinant):
        for beta,B in zip(INDICES,radial):
            gamma=tuple(i+j for i,j in zip(alpha,beta));weight=Q(multinomial(alpha)*multinomial(beta),multinomial(gamma))
            numerator[gamma]=numerator.get(gamma,I(0,0))+weight*A*B
    coefficients=[];ratios=[];axis=True
    for gamma,value in sorted(numerator.items()):
        denominator=sum((Q(gamma[i],4)*x[i][0] for i in range(3)),Q(0))
        if denominator==0:axis &= value.lo==value.hi==0
        else:ratios.append(value/denominator)
        coefficients.append({'index':list(gamma),'numerator':packed(value),'denominator_exact':str(denominator),
                             'ratio':packed(value/denominator) if denominator else None})
    valid_axis=axis and bool(ratios)
    J=I(min(v.lo for v in ratios),max(v.hi for v in ratios)) if valid_axis else None
    return {'positive_J_for_entire_box':bool(J and J.lo>0),'axis_zero_coefficients_exact_for_entire_box':bool(axis),
            'J_bound':packed(J) if J else None,'correlated_degree4_coefficients':coefficients,
            'strong_monotonicity_for_entire_box':all(v['positive_definite_for_entire_box'] for v in monotone),
            'strong_monotonicity_lower_exact':str(min(Q(v['lambda_lower_exact']) for v in monotone)) if all(v['positive_definite_for_entire_box'] for v in monotone) else None,
            'vertex_symmetric_gradient_bounds':monotone}


def certify_box(reference_nodes,current_centers,triangles,coordinate_radii):
    if not (len(reference_nodes)==len(current_centers)==len(coordinate_radii)):raise ValueError('matching node counts required')
    x=[[q(v) for v in row] for row in reference_nodes];y=[[q(v) for v in row] for row in current_centers]
    radii=[[q(v) for v in row] for row in coordinate_radii]
    if any(len(row)!=2 for row in x+y+radii) or any(v<0 for row in radii for v in row):raise ValueError('two coordinates and nonnegative radii required')
    n=max(i for t in triangles for i in t[:3])+1
    topology=strong_monotonicity_certificate(x[:n],x,triangles)
    if not topology.get('convex_reference_rectangle') or not topology.get('conforming_P2_traces'):
        raise ValueError('convex rectangular conforming reference not certified')
    Y=[[I(v-e,v+e) for v,e in zip(row,errors)] for row,errors in zip(y,radii)]
    records=[triangle_box([x[i] for i in t[:3]],[Y[i] for i in t]) for t in triangles]
    axis_preserved=all(Y[i][0].lo==Y[i][0].hi==0 for i,row in enumerate(x) if row[0]==0)
    local=all(r['positive_J_for_entire_box'] for r in records)
    mono=all(r['strong_monotonicity_for_entire_box'] for r in records)
    bounds=[I(Q(r['J_bound']['lower_exact']),Q(r['J_bound']['upper_exact'])) for r in records if r['J_bound']]
    return {'status':'CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING' if local and mono and axis_preserved else 'UNRESOLVED_GEOMETRY_BOX',
            'all_element_J_positive':local,'uniform_strong_monotonicity':mono,'axis_preserved_for_entire_box':axis_preserved,
            'global_injectivity_for_entire_box':local and mono and axis_preserved,
            'J_bound':packed(I(min(v.lo for v in bounds),max(v.hi for v in bounds))) if len(bounds)==len(records) else None,
            'strong_monotonicity_lower_exact':str(min(Q(r['strong_monotonicity_lower_exact']) for r in records)) if mono else None,
            'reference_topology':topology,'triangles':records,
            'scope':'Every conforming P2 field in the supplied nodal-coordinate box. No root existence, physical parameter coverage, continuum error or biological validation supplied by this function.'}


def synthetic():
    x=[[Q(0),Q(0)],[Q(1),Q(0)],[Q(1),Q(1)],[Q(0),Q(1)],
       [Q(1,2),Q(0)],[Q(1),Q(1,2)],[Q(1,2),Q(1,2)],[Q(1,2),Q(1)],[Q(0),Q(1,2)]]
    t=[[0,1,2,4,5,6],[0,2,3,6,7,8]];zero=[[Q(0),Q(0)] for row in x];checks={}
    identity=certify_box(x,x,t,zero)
    checks['identity_exact_J_one']=identity['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING' and identity['J_bound']==packed(1)
    affine=[[2*r,z/4] for r,z in x];iso=certify_box(x,affine,t,zero)
    checks['affine_isochoric_J_one']=iso['J_bound']==packed(1)
    small=[[Q(0) if r==0 else Q(1,10000),Q(1,10000)] for r,z in x]
    perturbed=certify_box(x,x,t,small)
    checks['nonzero_entire_box_certified']=perturbed['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING'
    checks['axis_uncertainty_not_certified']=certify_box(x,x,t,[[Q(1,10000)]*2 for _ in x])['status']=='UNRESOLVED_GEOMETRY_BOX'
    checks['large_box_unresolved_not_inversion_claim']=certify_box(x,x,t,[[0 if r==0 else 10,10] for r,z in x])['status']=='UNRESOLVED_GEOMETRY_BOX'
    checks['reflected_map_not_certified']=certify_box(x,[[r,-z] for r,z in x],t,zero)['status']=='UNRESOLVED_GEOMETRY_BOX'
    for a,b in ((-2,3),(-1,-Q(1,2)),(Q(1,3),2)):
        sq=I(a,b).square();checks[f'interval_square_{a}_{b}']=all(sq.lo<=v*v<=sq.hi for v in (q(a),q(b),(q(a)+q(b))/2))
    checks['interval_product_with_negative']=I(-2,1)*I(-3,2)==I(-4,6)
    for i,call in enumerate([lambda:I(1,0),lambda:I(1,2)/I(-1,1),lambda:certify_box(x,x,t,[[-1,0] for _ in x])]):
        try:call();checks[f'invalid_input_{i}']=False
        except ValueError:checks[f'invalid_input_{i}']=True
    out={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'identity':identity,'nonzero_box':perturbed,'source_sha256':sha(__file__),
         'scope':'Exact analytical/synthetic interval geometry; actual root application separate.'}
    (ROOT/'evidence/root_body_geometry_synthetics.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'checks':len(checks)}))
    if not all(checks.values()):raise RuntimeError('exact geometry synthetics failed')


def apply():
    path=ROOT/'evidence/root_body_geometry_verification.json'
    if path.exists():raise RuntimeError('existing root-body receipt; do not overwrite or rerun')
    cert_path=ROOT/'evidence/algebraic_root_certification_feasibility.json';review_path=ROOT/'evidence/algebraic_root_independent_review.json'
    c=json.loads(cert_path.read_text(encoding='utf-8'));r=json.loads(review_path.read_text(encoding='utf-8'))
    if r.get('status')!='PASS' or not r.get('checks') or not all(v is True for v in r['checks'].values()):raise RuntimeError('independent actual root review required')
    if r['certificate_sha256']!=sha(cert_path) or r['source_sha256']!=c['script_sha256']:raise RuntimeError('root review bindings differ')
    source=ROOT/c['input_state']
    if sha(source)!=c['input_state_sha256']:raise RuntimeError('root input changed')
    state=json.loads(source.read_text(encoding='utf-8'));blob=base64.b64decode(c['frozen_binary_cache']['data'])
    if hashlib.sha256(blob).hexdigest()!=c['frozen_binary_cache']['sha256']:raise RuntimeError('cache digest differs')
    with np.load(io.BytesIO(blob),allow_pickle=False) as z:
        reference=z['reference'];free=z['free'];u=z['u']
        if not np.array_equal(reference,np.asarray(state['reference_nodes'])) or not np.array_equal(u,np.asarray(state['free_displacement'])):
            raise RuntimeError('cache field differs')
    producer=min(Q(row['root_distance_upper']['upper_exact']) for row in c['radius_attempts'] if row['strict_selfmap_and_contraction'])
    radius=max(producer,Q(r['root_distance_upper_exact']));indices={int(k):i for i,k in enumerate(free)}
    centers=[];radii=[]
    for node,xy in enumerate(reference):
        centers.append([Q.from_float(float(v))+(Q.from_float(float(u[indices[2*node+j]])) if 2*node+j in indices else Q(0)) for j,v in enumerate(xy)])
        radii.append([radius if 2*node+j in indices else Q(0) for j in (0,1)])
    geometry=certify_box(state['reference_nodes'],centers,state['triangles'],radii)
    out={'status':'PASS_ACTUAL_ROOT_ENTIRE_BODY_GEOMETRY' if geometry['global_injectivity_for_entire_box'] else 'UNRESOLVED_ACTUAL_ROOT_GEOMETRY',
         'root_certificate_sha256':sha(cert_path),'root_review_sha256':sha(review_path),'input_sha256':sha(source),
         'node_radius_dimensionless_exact':str(radius),'centers':'Exact frozen reference plus archived free vector; fixed coordinates have zero radius',
         'source_sha256':sha(__file__),'dependency_hashes':{n:sha(Path(__file__).with_name(n)) for n in ('mixed_geometry_certificate_v2.py','mixed_geometry_certificate.py','certify_discrete_geometry.py')},
         'geometry':geometry,'scope':'The one independently certified old weak-load root lies in a whole-body positive-J, axis-preserving, globally injective nodal box. No continuum, stability, parameter inference or new17-field claim.'}
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'triangles':len(geometry['triangles']),'J_bound':geometry['J_bound'],'monotonicity_lower_exact':geometry['strong_monotonicity_lower_exact']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','apply']);args=parser.parse_args()
    {'test':synthetic,'apply':apply}[args.stage]()
