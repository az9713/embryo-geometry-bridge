"""Exact rational sufficient certificates for serialized P2 meridional maps.

Node order: v0,v1,v2,m01,m12,m20. Reference geometry is affine (P1).
Input floats denote their exact binary64 rational values, not uncertain data.
Uncertified is not synonymous with physically inverted.
"""
from fractions import Fraction as Q
from collections import Counter
from pathlib import Path
import hashlib, json, math, random

ORDER='vertices_then_01_12_20'

def rational(value):
    if isinstance(value,Q):return value
    if isinstance(value,int):return Q(value)
    if isinstance(value,str):return Q(value)
    value=float(value)
    if not math.isfinite(value):raise ValueError('finite coordinates required')
    return Q.from_float(value)

def points(values):
    return [tuple(rational(x) for x in row) for row in values]

def packed(value):
    return {'exact':str(value),'float':float(value)}

def pair(values):
    return {'lower':packed(min(values)),'upper':packed(max(values))}

def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])

def power_coefficients(values):
    """Quadratic Lagrange nodal values -> power basis in xi=L1,eta=L2."""
    a,b,c,d,e,f=values
    return {(0,0):a,(1,0):-3*a-b+4*d,(0,1):-3*a-c+4*f,
            (2,0):2*a+2*b-4*d,(0,2):2*a+2*c-4*f,
            (1,1):4*a-4*d+4*e-4*f}

def derivative(poly,axis):
    result={}
    for power,value in poly.items():
        if power[axis]:
            key=list(power);key[axis]-=1;result[tuple(key)]=value*power[axis]
    return result

def product(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():
            key=(i[0]+j[0],i[1]+j[1]);result[key]=result.get(key,Q(0))+x*y
    return result

def evaluate(poly,xi,eta):
    return sum((v*xi**i*eta**j for (i,j),v in poly.items()),Q(0))

def bernstein(poly):
    """Degree-two triangular Bernstein order v0,v1,v2,e01,e12,e20."""
    a=evaluate(poly,Q(0),Q(0));b=evaluate(poly,Q(1),Q(0));c=evaluate(poly,Q(0),Q(1))
    return [a,b,c,2*evaluate(poly,Q(1,2),Q(0))-(a+b)/2,
            2*evaluate(poly,Q(1,2),Q(1,2))-(b+c)/2,
            2*evaluate(poly,Q(0),Q(1,2))-(a+c)/2]

def certify_p2_triangle(reference_vertices,current_nodes,node_order=ORDER):
    if node_order!=ORDER:raise ValueError('unsupported node order')
    x=points(reference_vertices);y=points(current_nodes)
    if len(x)!=3 or len(y)!=6 or any(len(v)!=2 for v in x+y):raise ValueError('3 reference vertices and 6 current planar nodes required')
    area=cross(*x)
    if area<=0:raise ValueError('positive reference orientation required')
    if min(p[0] for p in x)<0:raise ValueError('nonnegative reference radii required')
    polys=[power_coefficients([v[i] for v in y]) for i in (0,1)]
    ad=product(derivative(polys[0],0),derivative(polys[1],1))
    bc=product(derivative(polys[0],1),derivative(polys[1],0))
    determinant={k:(ad.get(k,Q(0))-bc.get(k,Q(0)))/area for k in set(ad)|set(bc)}
    detcoeff=bernstein(determinant)
    radialcoeff=bernstein(polys[0])
    r0,r1,r2=[v[0] for v in x]
    refcoeff=[r0,r1,r2,(r0+r1)/2,(r1+r2)/2,(r2+r0)/2]
    zeros=all(v==0 for v,r in zip(radialcoeff,refcoeff) if r==0)
    ratios=[v/r for v,r in zip(radialcoeff,refcoeff) if r>0]
    hoop_positive=bool(zeros and ratios and min(ratios)>0)
    det_positive=min(detcoeff)>0
    certified=det_positive and hoop_positive
    answer={'status':'CERTIFIED_POSITIVE_DISCRETE_J' if certified else 'UNCERTIFIED',
            'node_order':ORDER,'arithmetic':'exact rational of supplied binary64 coordinates',
            'reference_determinant':packed(area),'meridional_determinant_bernstein':[str(v) for v in detcoeff],
            'current_radius_bernstein':[str(v) for v in radialcoeff],
            'reference_radius_bernstein':[str(v) for v in refcoeff],
            'meridional_determinant_bound':pair(detcoeff),
            'axis_zero_coefficients_exact':zeros,'hoop_ratio_coefficients':[str(v) for v in ratios],
            'hoop_stretch_bound':pair(ratios) if zeros and ratios else None,
            'positive_meridional_determinant_certified':det_positive,
            'positive_hoop_stretch_certified':hoop_positive,
            'scope':'Sufficient whole-element bounds, including axis limit. Failure to certify does not prove inversion; no equilibrium/continuum error bound.'}
    if certified:
        answer['J_bound']={'lower':packed(min(detcoeff)*min(ratios)),
                           'upper':packed(max(detcoeff)*max(ratios))}
    else:answer['J_bound']=None
    return answer

def edge_controls(a,m,b):
    return [a,tuple(2*m[i]-(a[i]+b[i])/2 for i in (0,1)),b]

def certify_flat_rectangular_boundary(reference_vertices,current_nodes,p2_triangles):
    """Conservative exact Jordan certificate for a flat rectangular reference.

    Verifies disk topology, shared mids, four sides and strict coordinate
    monotonicity. Top/bottom separation uses disjoint global z hulls, stronger
    than pointwise graph separation. An overlapping hull returns UNCERTIFIED.
    """
    from certify_discrete_geometry import disk_boundary, simple_polygon
    x=points(reference_vertices);y=points(current_nodes);tri=[list(map(int,t)) for t in p2_triangles]
    n=len(x)
    try:
        if any(len(t)!=6 or min(t)<0 or max(t)>=len(y) for t in tri):raise ValueError('invalid P2 connectivity')
        if set(i for t in tri for i in t[:3])!=set(range(n)):raise ValueError('vertex indices must precede mids')
        cycle=disk_boundary([t[:3] for t in tri],n);simple_polygon(x,cycle)
        edge_mid={}
        for t in tri:
            for a,b,m in ((t[0],t[1],t[3]),(t[1],t[2],t[4]),(t[2],t[0],t[5])):
                if m<n:raise ValueError('P2 midpoint indices must follow all vertex indices')
                key=tuple(sorted((a,b)))
                if key in edge_mid and edge_mid[key]!=m:raise ValueError('nonconforming shared P2 midpoint')
                edge_mid[key]=m
        if len(set(edge_mid.values()))!=len(edge_mid):raise ValueError('midpoint reused on different edges')
        R=max(v[0] for v in x);z0=min(v[1] for v in x);z1=max(v[1] for v in x)
        if min(v[0] for v in x)!=0 or R<=0 or z0>=z1:raise ValueError('nondegenerate axis rectangle required')
        sides={'bottom':[],'top':[],'axis':[],'rim':[]}
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            if x[a][1]==x[b][1]==z0:side='bottom'
            elif x[a][1]==x[b][1]==z1:side='top'
            elif x[a][0]==x[b][0]==0:side='axis'
            elif x[a][0]==x[b][0]==R:side='rim'
            else:raise ValueError('boundary is not flat rectangular reference')
            axis=0 if side in ('top','bottom') else 1
            if x[a][axis]>x[b][axis]:a,b=b,a
            m=edge_mid[tuple(sorted((a,b)))];c=edge_controls(y[a],y[m],y[b]);sides[side].append(c)
            if not all(c[i+1][axis]-c[i][axis]>0 for i in (0,1)):raise ValueError(side+' coordinate monotonicity uncertified')
        if any(not v for v in sides.values()):raise ValueError('four nonempty sides required')
        if any(c[i][0]!=0 for c in sides['axis'] for i in range(3)):raise ValueError('axis not preserved')
        rim_radius=sides['rim'][0][0][0]
        if rim_radius<=0 or any(c[i][0]!=rim_radius for c in sides['rim'] for i in range(3)):raise ValueError('rim not vertical')
        for side in ('top','bottom'):
            if min(c[i][0] for c in sides[side] for i in range(3))!=0 or max(c[i][0] for c in sides[side] for i in range(3))!=rim_radius:raise ValueError('graph side endpoints do not span axis to rim')
        gap=min(c[i][1] for c in sides['top'] for i in range(3))-max(c[i][1] for c in sides['bottom'] for i in range(3))
        if gap<=0:raise ValueError('global Bernstein z hulls overlap; pointwise graph order may still hold')
        return {'status':'CERTIFIED_JORDAN_BOUNDARY','disk_topology':True,'conforming_P2_edges':True,
                'strict_coordinate_monotonicity':True,'top_bottom_hull_gap':packed(gap),'axis_preserved':True,
                'scope':'Exact sufficient certificate for serialized flat-reference P2 boundary; no continuum approximation claim.'}
    except (AssertionError,ValueError,IndexError,KeyError) as exc:
        return {'status':'UNCERTIFIED','reason':str(exc) or 'topology/boundary assertion failed',
                'scope':'Not a proof of self-intersection or inversion.'}

def certify_mesh(reference_vertices,current_nodes,p2_triangles):
    records=[certify_p2_triangle([reference_vertices[i] for i in t[:3]],[current_nodes[i] for i in t]) for t in p2_triangles]
    boundary=certify_flat_rectangular_boundary(reference_vertices,current_nodes,p2_triangles)
    local=all(t['status']=='CERTIFIED_POSITIVE_DISCRETE_J' for t in records)
    global_ok=local and boundary['status']=='CERTIFIED_JORDAN_BOUNDARY'
    return {'status':'CERTIFIED_DISCRETE_EMBEDDING' if global_ok else 'UNCERTIFIED',
            'all_local_J_certified':local,'boundary':boundary,'triangles':records,
            'global_injectivity_certified':global_ok,
            'scope':'Exact P2 meridian embedding and axis-preserving revolution only, conditional on supplied P1 reference disk; not continuum error or biology validation.'}

def certify_saved(path):
    """Replay one mixed result JSON, with optional outer 'result' wrapper."""
    path=Path(path);root=Path(__file__).resolve().parents[1]
    raw=json.loads(path.read_text(encoding='utf-8'));data=raw.get('result',raw)
    n=max(i for t in data['triangles'] for i in t[:3])+1
    cert=certify_mesh(data['reference_nodes'][:n],data['current_nodes'],data['triangles'])
    try:input_name=path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:input_name=path.name
    cert.update({'input_path':input_name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                 'certificate_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 'reported_stationarity':data.get('stationary_discrete_equations',data.get('converged')),
                 'reported_J_bounds_are_float_diagnostics_only':[data.get('minimum_J'),data.get('maximum_J')]})
    cert['exact_surface_extrema']=surface_extrema(data['reference_nodes'][:n],data['current_nodes'],data['triangles'],data['parameters']['height'],data['parameters']['patch_radius'])
    destination=root/'evidence'/('mixed_geometry_'+path.stem+'.json')
    destination.write_text(json.dumps(cert,indent=2)+'\n')
    print(json.dumps({'input':str(path),'status':cert['status'],'local_J_pass':cert['all_local_J_certified'],
                      'boundary_status':cert['boundary']['status'],'triangles':len(cert['triangles']),
                      'receipt':str(destination)}))
    return cert

def quadratic_extrema(values):
    """Exact attained extrema of q(t), t∈[0,1], from endpoint/midpoint nodes."""
    a,m,b=values;c0=a;c1=-3*a+4*m-b;c2=2*a-4*m+2*b
    samples=[(Q(0),a),(Q(1),b)]
    if c2:
        t=-c1/(2*c2)
        if 0<t<1:samples.append((t,c0+c1*t+c2*t*t))
    return min(samples,key=lambda p:p[1]),max(samples,key=lambda p:p[1])

def surface_extrema(reference_vertices,current_nodes,p2_triangles,height,patch_radius):
    x=points(reference_vertices);y=points(current_nodes);zlo=min(p[1] for p in x);zhi=max(p[1] for p in x)
    edges={}
    for a,b,c,d,e,f in p2_triangles:
        for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):
            if x[i][0]>x[j][0]:i,j=j,i
            edges[i,j]=m
    sides=[]
    for z in (zlo,zhi):
        rows=[]
        for (i,j),m in edges.items():
            if x[i][1]==x[j][1]==z:rows.append((x[i][0],x[j][0],[y[i][1],y[m][1],y[j][1]]))
        sides.append(sorted(rows))
    if len(sides[0])!=len(sides[1]):raise ValueError('unaligned top/bottom reference edges')
    stats={}
    for name,rows in zip(('basal','apical'),sides):
        extrema=[quadratic_extrema(row[2]) for row in rows]
        stats[name]={'minimum_Z':packed(min(v[0][1] for v in extrema)),
                     'maximum_Z':packed(max(v[1][1] for v in extrema))}
    heights=[];patch=rational(patch_radius)
    for basal,apical in zip(*sides):
        if basal[:2]!=apical[:2]:raise ValueError('unaligned top/bottom material intervals')
        if basal[1]<=patch:
            vals=[a-b for a,b in zip(apical[2],basal[2])]
            minimum,_=quadratic_extrema(vals)
            heights.append((minimum[1],basal[0]+minimum[0]*(basal[1]-basal[0])))
    if heights:
        h,r=min(heights)
        stats.update({'maximum_activated_material_column_shortening':packed(rational(height)-h),
                      'argmax_reference_radius':packed(r),
                      'scope':'Exact extrema of quadratic Z curves and aligned material-column heights, not Eulerian thickness minima.'})
    return stats

def _self_test():
    checks={};records=[]
    def nodes(x,fn):
        allx=x+[(tuple((a[k]+b[k])/2 for k in (0,1))) for a,b in ((x[0],x[1]),(x[1],x[2]),(x[2],x[0]))]
        return [fn(*p) for p in allx]
    x=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))]
    a=certify_p2_triangle(x,nodes(x,lambda r,z:(2*r,z/2)))
    checks['isochoric_meridional_affine_axis']=a['status']=='CERTIFIED_POSITIVE_DISCRETE_J' and Q(a['J_bound']['lower']['exact'])==2 and Q(a['J_bound']['upper']['exact'])==2
    # Three-dimensional volume-preserving axisymmetric affine map.
    b=certify_p2_triangle(x,nodes(x,lambda r,z:(2*r,z/4)))
    checks['isochoric_3D_affine_axis']=Q(b['J_bound']['lower']['exact'])==1 and Q(b['J_bound']['upper']['exact'])==1
    rng=random.Random(1849)
    for k in range(40):
        a,b,c,d=[Q(rng.randint(-10,10),1000) for _ in range(4)]
        fn=lambda r,z:(r*(1+a*r+b*z),z+c*r*r+d*z*z)
        cert=certify_p2_triangle(x,nodes(x,fn));records.append(cert)
        checks[f'analytic_positive_map_{k}']=cert['status']=='CERTIFIED_POSITIVE_DISCRETE_J'
        # Independent analytic derivatives at a rational barycentric lattice.
        for i in range(9):
            for j in range(9-i):
                r,z=Q(i,8),Q(j,8)
                det=(1+2*a*r+b*z)*(1+2*d*z)-b*r*2*c*r
                hoop=1+a*r+b*z
                assert Q(cert['meridional_determinant_bound']['lower']['exact'])<=det<=Q(cert['meridional_determinant_bound']['upper']['exact'])
                assert Q(cert['hoop_stretch_bound']['lower']['exact'])<=hoop<=Q(cert['hoop_stretch_bound']['upper']['exact'])
                assert Q(cert['J_bound']['lower']['exact'])<=det*hoop<=Q(cert['J_bound']['upper']['exact'])
    checks['analytic_lattice_bounds_exact']=True
    off=[(Q(1),Q(0)),(Q(2),Q(0)),(Q(1),Q(1))]
    ambiguous=certify_p2_triangle(off,nodes(off,lambda r,z:(r+Q(9,10)*z*z,z+Q(9,10)*(r-1)**2)))
    checks['positive_map_may_be_uncertified']=ambiguous['status']=='UNCERTIFIED' and Q(1)-Q(9,10)**2>0
    inverted=certify_p2_triangle(off,nodes(off,lambda r,z:(r,-z)))
    checks['inverted_map_not_certified']=inverted['status']=='UNCERTIFIED'
    # Two-triangle P2 rectangular disk, globally shared mids.
    ref=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1)),(Q(0),Q(1))]
    tri=[[0,1,2,4,5,6],[0,2,3,6,7,8]]
    cur=ref+[(Q(1,2),Q(0)),(Q(1),Q(1,2)),(Q(1,2),Q(1,2)),(Q(1,2),Q(1)),(Q(0),Q(1,2))]
    mesh=certify_mesh(ref,cur,tri)
    checks['global_identity_disk']=mesh['status']=='CERTIFIED_DISCRETE_EMBEDDING'
    broken=[row[:] for row in tri];broken[1][3]=4
    checks['nonconforming_midpoint_not_certified']=certify_mesh(ref,cur,broken)['status']=='UNCERTIFIED'
    tiny=Q(1,10**40)
    qmin,qmax=quadratic_extrema([Q(0),-tiny/4,Q(0)])
    checks['arbitrarily_small_quadratic_extremum_retained']=qmin==(Q(1,2),-tiny/4) and qmax[1]==0
    root=Path(__file__).resolve().parents[1]
    receipt={'status':'PASS_EXACT_RATIONAL_SYNTHETIC_CERTIFICATE' if all(checks.values()) else 'FAIL',
             'checks':checks,'positive_maps':len(records),'first_positive_map':records[0],
             'positive_but_uncertified_example':ambiguous,'identity_mesh':mesh,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'scope':'Synthetic exact-rational verification; no mixed solver states, experimental outcomes or continuum error enclosures.'}
    (root/'evidence/mixed_geometry_synthetic.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'all_pass':all(checks.values())}))
    return all(checks.values())

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path)
    args=parser.parse_args()
    if args.input:raise SystemExit(0 if certify_saved(args.input)['status']=='CERTIFIED_DISCRETE_EMBEDDING' else 2)
    raise SystemExit(0 if _self_test() else 1)
