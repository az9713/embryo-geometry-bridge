"""Fresh exact geometry replay plus independent Bernstein/SPD bound assembly."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from collections import Counter
import hashlib,json,sys
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];OUT=CF/'runs/exact_integration_bridge_v1/retube_B_geometry';RUN=OUT/'certification_run'
sys.path.insert(0,str(CF/'scripts'));import root_body_geometry as body
from nodal_observer_error_transfer import curve_bound
from quadratic_endpoint_observer import Curve,positive_thickness
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return-a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):v=[x*y for x in a for y in b];return min(v),max(v)
def scalar(a,x):return(a[0]*x,a[1]*x)if x>=0 else(a[1]*x,a[0]*x)
def square(a):return(F(0)if a[0]<=0<=a[1]else min(x*x for x in a),max(x*x for x in a))
def interval(p):return F(p['lower_exact']),F(p['upper_exact'])
def packed(a):return{'lower_exact':str(a[0]),'upper_exact':str(a[1])}
def multinomial(a):
    n=factorial(sum(a))
    for x in a:n//=factorial(x)
    return n
def shape_derivatives(k,x,y):
    # Construct each Kronecker nodal basis via independent power coefficients.
    v=[F(int(j==k))for j in range(6)];c10=-3*v[0]-v[1]+4*v[3];c01=-3*v[0]-v[2]+4*v[5];c20=2*v[0]+2*v[1]-4*v[3];c02=2*v[0]+2*v[2]-4*v[5];c11=4*(v[0]-v[3]+v[4]-v[5]);return c10+2*c20*x+c11*y,c01+c11*x+2*c02*y
def own_triangle(x,y):
    a=x[1][0]-x[0][0];b=x[2][0]-x[0][0];c=x[1][1]-x[0][1];d=x[2][1]-x[0][1];det=a*d-b*c
    if det<=0:raise ValueError('nonpositive material triangle')
    inverse=((d/det,-b/det),(-c/det,a/det));G=[];spd=[]
    for xi,eta in[(F(0),F(0)),(F(1),F(0)),(F(0),F(1))]:
        D=[shape_derivatives(k,xi,eta)for k in range(6)];g=[[sum_intervals(scalar(y[k][i],D[k][j])for k in range(6))for j in range(2)]for i in range(2)];G.append(g)
        Fm=[[sum_intervals(scalar(g[i][k],inverse[k][j])for k in range(2))for j in range(2)]for i in range(2)];aa=Fm[0][0];bb=scalar(add(Fm[0][1],Fm[1][0]),F(1,2));dd=Fm[1][1];ds=sub(mul(aa,dd),square(bb))
        if not aa[0]>0 or not ds[0]>0:raise ValueError('vertexSPD not certified')
        spd.append({'aa':aa,'bb':bb,'dd':dd,'m':ds[0]/(aa[1]+dd[1])})
    indices=[(2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1)];D=[]
    for i,j in[(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]:D.append(scalar(sub(sub(add(mul(G[i][0][0],G[j][1][1]),mul(G[j][0][0],G[i][1][1])),mul(G[i][0][1],G[j][1][0])),mul(G[j][0][1],G[i][1][0])),1/(2*det)))
    radial=[y[k][0]for k in range(3)]+[sub(scalar(y[m][0],2),scalar(add(y[i][0],y[j][0]),F(1,2)))for i,j,m in[(0,1,3),(1,2,4),(2,0,5)]];poly={}
    for alpha,A in zip(indices,D):
        for beta,B in zip(indices,radial):
            gamma=tuple(i+j for i,j in zip(alpha,beta));w=F(multinomial(alpha)*multinomial(beta),multinomial(gamma));poly[gamma]=add(poly.get(gamma,(F(0),F(0))),scalar(mul(A,B),w))
    coeff=[];ratios=[]
    for gamma,N in sorted(poly.items()):
        denominator=sum(F(gamma[i],4)*x[i][0]for i in range(3))
        if denominator==0:
            if N!=(0,0):raise ValueError('axiszero coefficient failed')
            ratio=None
        else:ratio=scalar(N,1/denominator);ratios.append(ratio)
        coeff.append((gamma,N,denominator,ratio))
    return coeff,spd,(min(v[0]for v in ratios),max(v[1]for v in ratios))
def sum_intervals(values):
    total=(F(0),F(0))
    for value in values:total=add(total,value)
    return total
def z_range(controls):
    values=[]
    for edge in controls:
        p0,pm,p1=[[F(v)for v in row]for row in edge];c0=p0[1];c1=-3*p0[1]+4*pm[1]-p1[1];c2=2*p0[1]-4*pm[1]+2*p1[1];values.extend([c0,c0+c1+c2])
        if c2:
            t=-c1/(2*c2)
            if 0<t<1:values.append(c0+c1*t+c2*t*t)
    return min(values),max(values)
def main():
    if sha(RUN/'results.json')!='788e1bc7a87fa7363bc5811b824fb5d44e1c5420428c3234283fc51831ddb9a3':raise ValueError('authoritative result changed')
    p=read(OUT/'protocol.json');result=read(RUN/'results.json');box=read(OUT/'coordinate_boxes.json');boundary=read(OUT/'boundary_preparation.json');geometry=read(RUN/'whole_body_geometry.json');checks={};bindings={}
    for k,v in p['source_and_input_hashes'].items():
        if sha(CF/k)!=v:raise ValueError('frozen source/input changed '+k)
        bindings[k]=v
    checks['prelaunch668_and_prepared_boxes_bound']=sha(OUT/'prelaunch_independent_review.json')=='97b1c8a8d4c7a3db00b51d35ef0d7396b5ee3df09972d5fc04fc68b966e31577'and sha(OUT/'coordinate_boxes.json')==p['prepared_box_sha256']and sha(OUT/'boundary_preparation.json')==p['boundary_preparation_sha256']
    checks['actual_body_output_bound']=sha(RUN/'whole_body_geometry.json')==result['body_geometry_sha256']
    # Explicitly reviewed component replay; not mislabelled a second kernel.
    replay=body.certify_box(box['reference_nodes_exact'],box['current_centers_exact'],box['triangles'],box['coordinate_radii_exact']);checks['reviewed_body_component_all144_replay_identical']=replay==geometry
    x=[[F(v)for v in row]for row in box['reference_nodes_exact']];centers=[[F(v)for v in row]for row in box['current_centers_exact']];radii=[[F(v)for v in row]for row in box['coordinate_radii_exact']];Y=[[(v-e,v+e)for v,e in zip(c,r)]for c,r in zip(centers,radii)];tri=box['triangles'];allJ=[];allm=[];coefficient_count=0;vertex_count=0
    for index,t in enumerate(tri):
        coeff,spd,jbound=own_triangle([x[k]for k in t[:3]],[Y[k]for k in t]);old=geometry['triangles'][index];good=True
        for (gamma,N,den,ratio),saved in zip(coeff,old['correlated_degree4_coefficients']):good &= list(gamma)==saved['index']and N==interval(saved['numerator'])and den==F(saved['denominator_exact'])and(ratio is None and saved['ratio']is None or ratio is not None and ratio==interval(saved['ratio']));coefficient_count+=1
        checks[f'element{index:03d}_all15_independent_Bernstein_coefficients']=bool(good and len(coeff)==len(old['correlated_degree4_coefficients'])==15)
        checks[f'element{index:03d}_fresh_J_bound_independent']=jbound==interval(old['J_bound'])and jbound[0]>F(1,1000000);allJ.append(jbound)
        for k,(v,saved)in enumerate(zip(spd,old['vertex_symmetric_gradient_bounds'])):
            checks[f'element{index:03d}_vertex{k}_SPD_independent']=v['aa']==interval(saved['symmetric_gradient'][0])and v['bb']==interval(saved['symmetric_gradient'][1])and v['dd']==interval(saved['symmetric_gradient'][2])and v['m']==F(saved['lambda_lower_exact'])and v['m']>0;allm.append(v['m']);vertex_count+=1
    J=(min(v[0]for v in allJ),max(v[1]for v in allJ));rho=max(abs(J[0]-1),abs(J[1]-1));m=min(allm)
    checks['all2160_J_coefficients_and432_SPD_vertices_checked']=len(tri)==144 and coefficient_count==2160 and vertex_count==432
    checks['fresh_global_J_rho_and_m_recomputed']=J==interval(result['new_J_interval'])and rho==F(result['fresh_rho_about1_exact'])and m==F(geometry['strong_monotonicity_lower_exact'])and J[0]>F(1,1000000)and rho<1
    checks['axis_radial_boxes_exactzero']=all(Y[i][0]==(0,0)for i,ref in enumerate(x)if ref[0]==0)
    # Constructive structured affine triangulation and shared P2 traces.
    n=max(k for t in tri for k in t[:3])+1;rv=sorted({p[0]for p in x[:n]});zv=sorted({p[1]for p in x[:n]});grid={tuple(p):i for i,p in enumerate(x[:n])};expected=set()
    for a,b in zip(rv,rv[1:]):
        for c,d in zip(zv,zv[1:]):
            p00,p10,p11,p01=[grid[k]for k in[(a,c),(b,c),(b,d),(a,d)]];expected.update([frozenset([p00,p10,p11]),frozenset([p00,p11,p01])])
    checks['material_reference_exact_rectangular_grid_triangulation']=len(grid)==len(rv)*len(zv)and len(expected)==len(tri)==144 and expected=={frozenset(t[:3])for t in tri}
    midpoint={};conforming=True
    for t in tri:
        for i,j,k in[(t[0],t[1],t[3]),(t[1],t[2],t[4]),(t[2],t[0],t[5])]:
            edge=tuple(sorted((i,j)));conforming &= k>=n and(edge not in midpoint or midpoint[edge]==k);midpoint[edge]=k
    checks['P2_shared_edge_labels_conforming']=conforming and len(set(midpoint.values()))==len(midpoint)
    checks['global_embedding_follows_from_uniform_SPD_on_convex_reference']=geometry['global_injectivity_for_entire_box']is True and m>0 and checks['material_reference_exact_rectangular_grid_triangulation']and checks['P2_shared_edge_labels_conforming']
    transfers={};errors={}
    for name in['top','bottom']:
        controls=boundary['baseline_controls_um'][name];eps=F(boundary['boundary_specific_nodal_error_um_exact'][name]);fresh=curve_bound(controls,eps,fixed_radial_endpoints=True,roi_um=boundary['same_current_radius_ROI_um_exact']);checks[name+'_reviewed_graph_component_replay']=fresh==result['boundary_graph_transfer'][name];slope=F(0);radial=None
        for edge in controls:
            a,b,c=[[F(v)for v in q]for q in edge];c1=[-3*a[k]+4*b[k]-c[k]for k in range(2)];c2=[2*a[k]-4*b[k]+2*c[k]for k in range(2)];rmin=min(c1[0],c1[0]+2*c2[0]);zmax=max(abs(c1[1]),abs(c1[1]+2*c2[1]));radial=rmin if radial is None else min(radial,rmin);slope=max(slope,zmax/rmin)
        error=F(5,4)*eps*(1+slope);errors[name]=error;checks[name+'_independent_derivative_error_formula']=radial-8*eps==F(fresh['radial_derivative_margin_exact'])and slope==F(fresh['base_graph_Lipschitz_upper_exact'])and error==F(fresh['Eulerian_vertical_error_um_exact'])and radial>8*eps
    thick=positive_thickness(Curve(boundary['baseline_controls_um']['top']),Curve(boundary['baseline_controls_um']['bottom']),list(map(F,boundary['same_current_radius_ROI_um_exact'])),F(1,10000000),8192);sep=thick-errors['top']-errors['bottom'];checks['reviewed_thickness_component_and_independent_error_subtraction']=thick==F(result['baseline_thickness_lower_um_exact'])and sep==F(result['fresh_uniform_Eulerian_thickness_lower_um_exact'])and sep>0
    topz=z_range(boundary['baseline_controls_um']['top']);bottomz=z_range(boundary['baseline_controls_um']['bottom']);independent_thickness=topz[0]-bottomz[1]-errors['top']-errors['bottom'];checks['independent_global_quadratic_extrema_thickness_positive']=independent_thickness>0
    cor=result['uniform_surface_positive_lower_corollary'];checks['surface_lower_only_corollary_multipliers']=F(cor['meridional_m_exact'])==F(cor['hoop_lower_exact'])==F(cor['canonical_radial_derivative_lower_exact'])==m and F(cor['Js_lower_exact'])==m*m and F(cor['G_lower_exact'])==m**4
    checks['no_root_or_surface_upperG_claim']=result['root_membership_or_new_integral_root_certified']is False and result['uniform_surface_G_upper_or_binomial_rho_established']is False
    checks['all_frozen_inputs_still_current']=all(sha(CF/k)==v for k,v in bindings.items())
    details={'fresh_J_interval_exact':packed(J),'fresh_rho_exact':str(rho),'uniform_symmetric_gradient_m_exact':str(m),'fresh_Eulerian_thickness_lower_um_exact':str(sep),'independent_global_extrema_thickness_lower_um_exact':str(independent_thickness),'body_elements':144,'Bernstein_J_coefficients':coefficient_count,'SPD_vertices':vertex_count}
    out={'status':'PASS_FRESH_COMMON_NODAL_GEOMETRY_COMPONENT'if all(checks.values())else'FAIL_GEOMETRY_REVIEW','checks':checks,'review_source_sha256':sha(Path(__file__)),'producer_results_sha256':sha(RUN/'results.json'),'body_output_sha256':sha(RUN/'whole_body_geometry.json'),'protocol_sha256':sha(OUT/'protocol.json'),'source_and_input_hashes':bindings,'details':details,'reuse_disclosure':'Replayed reviewed body.certify_box, curve_bound and positive_thickness. Independently reimplemented exact interval arithmetic, P2power-basis derivative construction, all144Bernstein/SPD assemblies, structuredreference/conformingtracechecks, grapherrorformula andglobalquadraticextrema bound. Not a wholly separate foundational geometry theorem.','scope':'Every canonical field in newlydeclared containingboxes; membership as anyrootisNOTassumed orproved. No H/cache/integration/root/PDE/biologicalvalidation, no surfaceGupper.'}
    HERE.mkdir(parents=True,exist_ok=True);(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');(CF/'evidence/exact_retube_B_geometry_independent_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'J_display':list(map(float,J)),'rho':float(rho),'m':float(m),'thickness_um':float(sep),'independent_global_extrema_thickness_um':float(independent_thickness),'receipt_sha256':sha(CF/'evidence/exact_retube_B_geometry_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
