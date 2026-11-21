"""Exact rational orientation/Jacobian/embedding checks for saved P1 bodies.

Certificates concern the serialized binary64-node piecewise-linear meridian,
not closeness to a continuum solution or statistical validity.
"""
from collections import Counter,defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def sign(x):
    return (x>0)-(x<0)


def on_segment(a,b,p):
    return cross(a,b,p)==0 and min(a[0],b[0])<=p[0]<=max(a[0],b[0]) and min(a[1],b[1])<=p[1]<=max(a[1],b[1])


def intersects(a,b,c,d):
    ab1,ab2=cross(a,b,c),cross(a,b,d)
    cd1,cd2=cross(c,d,a),cross(c,d,b)
    if sign(ab1)*sign(ab2)<0 and sign(cd1)*sign(cd2)<0:return True
    return (ab1==0 and on_segment(a,b,c)) or (ab2==0 and on_segment(a,b,d)) or (cd1==0 and on_segment(c,d,a)) or (cd2==0 and on_segment(c,d,b))


def disk_boundary(triangles,n):
    directed=Counter();edges=Counter();links=defaultdict(list)
    for a,b,c in triangles:
        for i,j in [(a,b),(b,c),(c,a)]:directed[i,j]+=1;edges[tuple(sorted((i,j)))]+=1
        for i,j,k in [(a,b,c),(b,c,a),(c,a,b)]:links[i].append((j,k))
    assert all(v in (1,2) for v in edges.values())
    assert n-len(edges)+len(triangles)==1
    boundary_edges=[(i,j) for (i,j) in directed if edges[tuple(sorted((i,j)))]==1]
    for (i,j),count in edges.items():
        if count==2:assert directed[i,j]==1 and directed[j,i]==1
    successor=dict(boundary_edges)
    assert len(successor)==len(boundary_edges)
    first=boundary_edges[0][0];cycle=[first];current=successor[first]
    while current!=first:
        assert current not in cycle;cycle.append(current);current=successor[current]
    assert len(cycle)==len(boundary_edges)
    boundary=set(cycle)
    # Every vertex link must be one path (boundary) or one cycle (interior).
    for vertex in range(n):
        adjacency=defaultdict(set)
        for i,j in links[vertex]:adjacency[i].add(j);adjacency[j].add(i)
        pending=[next(iter(adjacency))];seen=set()
        while pending:
            i=pending.pop()
            if i in seen:continue
            seen.add(i);pending.extend(adjacency[i]-seen)
        assert len(seen)==len(adjacency)
        degrees=sorted(len(a) for a in adjacency.values())
        if vertex in boundary:assert degrees.count(1)==2 and all(d in (1,2) for d in degrees)
        else:assert all(d==2 for d in degrees)
    return cycle


def simple_polygon(nodes,cycle):
    points=[nodes[i] for i in cycle]
    assert len(set(points))==len(points)
    length=len(points)
    for i in range(length):
        a,b=points[i],points[(i+1)%length]
        c=points[(i+2)%length]
        if cross(a,b,c)==0:
            assert sum((b[k]-a[k])*(c[k]-b[k]) for k in (0,1))>0
        for j in range(i+1,length):
            if j==i+1 or (i==0 and j==length-1):continue
            assert not intersects(a,b,points[j],points[(j+1)%length])


def certify(result):
    reference=[tuple(Q.from_float(float(x)) for x in p) for p in result['reference_nodes']]
    current=[tuple(Q.from_float(float(x)) for x in p) for p in result['current_nodes']]
    triangles=result['triangles'];cycle=disk_boundary(triangles,len(reference))
    simple_polygon(reference,cycle);simple_polygon(current,cycle)
    for old,new in zip(reference,current):
        assert old[0]>=0 and new[0]>=0
        assert (old[0]==0)==(new[0]==0)
    minimum=None;maximum=None;v0=Q(0);v1=Q(0)
    for face in triangles:
        x=[reference[i] for i in face];y=[current[i] for i in face]
        a0=cross(*x);a1=cross(*y)
        assert a0>0 and a1>0
        ratios=[new[0]/old[0] for old,new in zip(x,y) if old[0]>0]
        assert ratios and min(ratios)>0
        lo=a1/a0*min(ratios);hi=a1/a0*max(ratios)
        minimum=lo if minimum is None else min(minimum,lo)
        maximum=hi if maximum is None else max(maximum,hi)
        v0+=a0*sum(p[0] for p in x)
        v1+=a1*sum(p[0] for p in y)
    # Volume = pi/3 sum(det(Dm)*sum(vertex radii)), so pi cancels.
    ratio=v1/v0
    return {'positive_meridional_orientation_exact':True,'reference_and_current_boundary_simple_exact':True,
            'disk_vertex_links_verified':True,'axis_and_positive_radial_interior_exact':True,
            'minimum_J_exact':str(minimum),'maximum_J_exact':str(maximum),
            'minimum_J_float':float(minimum),'maximum_J_float':float(maximum),
            'volume_ratio_exact':str(ratio),'relative_volume_change_float':float(ratio-1),
            'reported_J_agreement':abs(float(minimum)-result['minimum_J'])<1e-12 and abs(float(maximum)-result['maximum_J'])<1e-12,
            'reported_volume_agreement':abs(float(ratio-1)-result['relative_volume_change'])<1e-12,
            'vertices':len(reference),'triangles':len(triangles),'boundary_vertices':len(cycle)}


def main():
    paths=sorted((ROOT/'runs/spatial_transfer').glob('basal_*/calibration_direct_*.json'))+sorted((ROOT/'runs/spatial_transfer').glob('basal_*/prediction_*.json'))
    assert len(paths)==144
    records=[];start=time.monotonic()
    for path in paths:
        source=json.loads(path.read_text(encoding='utf-8'))
        result=certify(source['result'])
        assert result['reported_J_agreement'] and result['reported_volume_agreement']
        records.append({'path':path.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),**result})
    receipt={'status':'PASS_EXACT_DISCRETE_GEOMETRY','states':len(records),'records':records,'seconds':time.monotonic()-start,
             'scope':'Exact rational verification of saved binary64-node P1 meridians. Disk topology + positive triangle orientation + Jordan boundary implies global PL homeomorphism by degree; axis preservation gives an injective revolved body.',
             'not_proved':['equilibrium error enclosure','continuum approximation error','incompressible cell dynamics','biological model validity','novelty'],
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'evidence/discrete_geometry_certificate.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:receipt[k] for k in ['status','states','seconds','scope']}))


if __name__=='__main__':
    main()
