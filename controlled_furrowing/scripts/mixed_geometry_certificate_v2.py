"""Additional exact P2 certificates; the original verifier is unchanged.

Combines det(Fm)*R before Bernstein division by r, and checks sufficient
strong monotonicity on a convex rectangular reference meridian.
"""
from fractions import Fraction as Q
from math import factorial
from mixed_geometry_certificate import (points,packed,pair,cross,certify_p2_triangle,
                                        power_coefficients,derivative,evaluate)

INDICES2=((2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1))

def multinomial(index):
    result=factorial(sum(index))
    for v in index:result//=factorial(v)
    return result

def correlated_J_certificate(reference_vertices,current_nodes):
    """Whole-element rational Bernstein enclosure of J=det(Fm)*R/r."""
    base=certify_p2_triangle(reference_vertices,current_nodes)
    x=points(reference_vertices);det=list(map(Q,base['meridional_determinant_bernstein']))
    R=list(map(Q,base['current_radius_bernstein']))
    numerator={}
    for alpha,a in zip(INDICES2,det):
        for beta,b in zip(INDICES2,R):
            gamma=tuple(i+j for i,j in zip(alpha,beta))
            numerator[gamma]=numerator.get(gamma,Q(0))+Q(multinomial(alpha)*multinomial(beta),multinomial(gamma))*a*b
    coeff=[];ratios=[];axis=True
    for gamma in sorted(numerator):
        denominator=sum((Q(gamma[i],4)*x[i][0] for i in range(3)),Q(0))
        value=numerator[gamma]
        if denominator==0:axis &= value==0
        elif denominator>0:ratios.append(value/denominator)
        else:raise ValueError('nonnegative reference radius required')
        coeff.append({'index':list(gamma),'numerator':str(value),'reference_radius':str(denominator),
                      'ratio':str(value/denominator) if denominator else None})
    bounded=bool(axis and ratios)
    return {'status':'CERTIFIED_DISCRETE_J_INTERVAL' if bounded else 'UNCERTIFIED',
            'axis_zero_coefficients_exact':bool(axis),'coefficients':coeff,
            'J_bound':pair(ratios) if bounded else None,
            'positive_J_certified':bool(bounded and min(ratios)>0),
            'legacy_product_J_bound':base['J_bound'],
            'scope':'Exact degree4 rational Bernstein enclosure of the stored P2 map; not necessarily attained extrema or a continuum error bound.'}

def strong_monotonicity_certificate(reference_vertices,current_nodes,p2_triangles):
    from certify_discrete_geometry import disk_boundary,simple_polygon
    x=points(reference_vertices);y=points(current_nodes);n=len(x);tri=[list(map(int,t)) for t in p2_triangles]
    try:
        cycle=disk_boundary([t[:3] for t in tri],n);simple_polygon(x,cycle)
        r0=min(p[0] for p in x);r1=max(p[0] for p in x);z0=min(p[1] for p in x);z1=max(p[1] for p in x)
        if not r0<r1 or not z0<z1:raise ValueError('nondegenerate rectangle required')
        if not {(r0,z0),(r1,z0),(r1,z1),(r0,z1)}<=set(x[i] for i in cycle):raise ValueError('missing rectangle corners')
        for i,j in zip(cycle,cycle[1:]+cycle[:1]):
            if not ((x[i][0]==x[j][0] and x[i][0] in (r0,r1)) or (x[i][1]==x[j][1] and x[i][1] in (z0,z1))):raise ValueError('reference boundary not rectangular')
        mids={}
        for t in tri:
            for i,j,m in ((t[0],t[1],t[3]),(t[1],t[2],t[4]),(t[2],t[0],t[5])):
                key=tuple(sorted((i,j)))
                if m<n or key in mids and mids[key]!=m:raise ValueError('nonconforming P2 midpoint')
                mids[key]=m
        if len(set(mids.values()))!=len(mids):raise ValueError('midpoint reused across edges')
        bounds=[];records=[]
        for t in tri:
            v=[x[i] for i in t[:3]];a=v[1][0]-v[0][0];b=v[2][0]-v[0][0];c=v[1][1]-v[0][1];d=v[2][1]-v[0][1];det=a*d-b*c
            if det<=0:raise ValueError('reference orientation not positive')
            inv=((d/det,-b/det),(-c/det,a/det))
            polys=[power_coefficients([y[i][j] for i in t]) for j in (0,1)]
            for xi,eta in ((Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))):
                G=[[evaluate(derivative(poly,j),xi,eta) for j in (0,1)] for poly in polys]
                F=[[sum(G[i][k]*inv[k][j] for k in (0,1)) for j in (0,1)] for i in (0,1)]
                aa=F[0][0];bb=(F[0][1]+F[1][0])/2;dd=F[1][1];sd=aa*dd-bb*bb
                if aa<=0 or sd<=0:raise ValueError('vertex symmetric gradient not positive definite')
                lower=sd/(aa+dd);bounds.append(lower)
                records.append({'symmetric_F':[str(aa),str(bb),str(dd)],'lambda_min_lower':str(lower)})
        axis=all(y[i][0]==0 for i in range(n) if x[i][0]==0)
        axis &= all(y[m][0]==0 for (i,j),m in mids.items() if x[i][0]==x[j][0]==0)
        return {'status':'CERTIFIED_STRONGLY_MONOTONE_DISCRETE_MAP','global_injectivity_certified':True,
                'convex_reference_rectangle':True,'conforming_P2_traces':True,
                'strong_monotonicity_constant_lower':packed(min(bounds)),
                'axis_preserved':bool(axis),'vertex_symmetric_gradient_certificates':records,
                'scope':'Exact sufficient discrete meridional injectivity on convex reference; positive axisymmetric J is checked separately.'}
    except (AssertionError,ValueError,IndexError,KeyError) as exc:
        return {'status':'UNCERTIFIED','global_injectivity_certified':False,'reason':str(exc) or 'topology assertion failed'}

def _self_test():
    from pathlib import Path
    import hashlib,json,random
    checks={};examples=[];rng=random.Random(917)
    ref=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))]
    sample=ref+[(Q(1,2),Q(0)),(Q(1,2),Q(1,2)),(Q(0),Q(1,2))]
    for n in range(40):
        a,b,c,d=[Q(rng.randint(-10,10),1000) for _ in range(4)]
        cur=[(r*(1+a*r+b*z),z+c*r*r+d*z*z) for r,z in sample]
        cert=correlated_J_certificate(ref,cur);examples.append(cert)
        checks[f'analytic_correlated_positive_{n}']=cert['positive_J_certified']
        for i in range(1,9):
            for j in range(9-i):
                r,z=Q(i,8),Q(j,8)
                J=((1+2*a*r+b*z)*(1+2*d*z)-b*r*2*c*r)*(1+a*r+b*z)
                bary=(1-r-z,r,z);N=Q(0)
                for row in cert['coefficients']:
                    gamma=row['index'];N+=Q(row['numerator'])*multinomial(gamma)*bary[0]**gamma[0]*bary[1]**gamma[1]*bary[2]**gamma[2]
                assert N/r==J
                assert Q(cert['J_bound']['lower']['exact'])<=J<=Q(cert['J_bound']['upper']['exact'])
    checks['independent_exact_analytic_product_identity_and_bounds']=True
    affine=correlated_J_certificate(ref,[(2*r,z/4) for r,z in sample])
    checks['constant_isochoric_affine_exact_bounds']=Q(affine['J_bound']['lower']['exact'])==Q(affine['J_bound']['upper']['exact'])==1
    rectangle=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1)),(Q(0),Q(1))]
    cur=rectangle+[(Q(1,2),Q(0)),(Q(1),Q(1,2)),(Q(1,2),Q(1,2)),(Q(1,2),Q(1)),(Q(0),Q(1,2))]
    tri=[[0,1,2,4,5,6],[0,2,3,6,7,8]]
    mono=strong_monotonicity_certificate(rectangle,cur,tri)
    checks['identity_strong_monotonicity']=mono['global_injectivity_certified'] and mono['axis_preserved']
    flipped=[(r,-z) for r,z in cur]
    checks['inverted_map_not_strongly_monotone']=not strong_monotonicity_certificate(rectangle,flipped,tri)['global_injectivity_certified']
    broken=[t[:] for t in tri];broken[1][3]=4
    checks['nonconforming_map_not_certified']=not strong_monotonicity_certificate(rectangle,cur,broken)['global_injectivity_certified']
    root=Path(__file__).resolve().parents[1]
    receipt={'status':'PASS_EXACT_V2_SYNTHETICS' if all(checks.values()) else 'FAIL','checks':checks,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'v1_source_sha256':hashlib.sha256(Path(__file__).with_name('mixed_geometry_certificate.py').read_bytes()).hexdigest(),
             'example':examples[0],'affine':affine,'identity_strong_monotonicity':mono}
    (root/'evidence/mixed_geometry_v2_synthetic.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'all_pass':all(checks.values())}))
    return all(checks.values())

if __name__=='__main__':raise SystemExit(0 if _self_test() else 1)
