"""Independent physical-moment/selector and adapter prelaunch checks."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
from flint import ctx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_axis_edge_bridge_B_safe as P
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
checks={}
def ck(k,v):
    checks[k]=bool(v)
    assert v,k
source=ROOT/'scripts/exact_uniform_axis_edge_bridge_B_safe.py'
ck('source_literal',sha(source)=='a3128a70e56507b287532f896ef4bc5cc2457603baa87173e22476ee412b6c60')
ck('protocol_literal',sha(HERE/'protocol.json')=='6fcfb2c20d32885f1b493c1be601d8ba469d9089fb4c50d2a25cf1d0c26d777a')
ck('preflight_literal',sha(HERE/'preflight.json')=='1480daa99fb626157b556c0a4244d6a860bf675476bbda09db2418727ec64096')
p=read(HERE/'protocol.json');pre=read(HERE/'preflight.json')
ck('preflight_nonempty_literal_pass',pre['status']=='PASS_PHYSICAL_AXIS_EDGE_BRIDGE_PREPARATION_ONLY'and pre['checks']and all(v is True for v in pre['checks'].values()))
for path,h in p['source_and_input_hashes'].items():ck('binding_'+path,sha(ROOT/path)==h)
ctx.prec=512;P.core.A.cache_clear()
record=read(ROOT/P.P.RECORD);target=P.T.Target(read(ROOT/record['point_cache_file']));Y=P.U.field_nodes(read(P.U.GEOM/'coordinate_boxes.json'));rule=target.auxiliary_rule();body=read(P.U.GEOM/'certification_run/whole_body_geometry.json')
ck('exact_six_axis_edge_indices',[s['element']for s in p['selectors']]==list(range(72,78)))
ck('same_C_and_dose',p['parameter_cell']==record['parameter_cell']and p['dose_exact']==record['dose_exact'])
for selector in p['selectors']:
    i=selector['element'];fields=P.AX.physical_element(target,i,nodes=Y);point=target.element(i);chart=point['chart'];corners=fields['corners'];zcenter=fields['zcenter']
    ck(f'no_negative_physical_powers_{i}',all(a>=0 and b>=0 for poly in fields['q']+[v for row in fields['B']for v in row]for a,b in poly))
    ck(f'same_free_map_{i}',fields['global_indices']==selector['global_indices']==point['global_indices'])
    # Two exact geometrical integrations of ordinary physical monomials:
    # barycentric triangle moments vs transformed Duffy beta-product moments.
    R=chart['R'];Z=P.T.sub(chart['Z'],P.T.constant(zcenter))
    for a in range(5):
        for b in range(5):
            transformed=P.T.mul(P.T.power(R,a),P.T.power(Z,b))
            duffy=sum((v*P.M.axis_chart_moment('axis_edge',j,k,chart['R1'],fields['det'])for(j,k),v in transformed.items()),F(0))
            ordinary=P.M.triangle_moment(corners,a+1,b,zcenter)
            ck(f'exact_geometry_moment_{i}_{a}_{b}',not ordinary.logarithms and ordinary.constant==duffy)
    jb=body['triangles'][i]['J_bound'];raw=max(abs(F(jb['lower_exact'])-1),abs(F(jb['upper_exact'])-1));v=raw*2**48;rho=F((v.numerator+v.denominator-1)//v.denominator,2**48)
    ck(f'exact_rho_{i}',rho==F(selector['rho_exact'])and raw<=rho<1)
    # Recompute the first acceptable order using independent exact scalar
    # derivative majorants; no field-defect moments or64-site sums are run.
    jet=P.JET.bulk_jacobian_jet(fields['q'],fields['B']);rmax=max(v[0]for v in corners);zmax=max(abs(v[1]-zcenter)for v in corners)
    def norm(poly):
        return sum(P.abs_coefficient(c)*rmax**a*zmax**b for(a,b),c in poly.items())
    g=max(map(norm,jet['gradient']));h=max(map(norm,jet['hessian_upper'].values()))
    mass=F(selector['geometry_two_mass_upper_exact'])
    possible=[]
    for n in range(2,19):
        t1,t2=P.M.volumetric_derivative_tail_bounds(n,rho,1,5000)
        if mass*t1*g<=F(1,10**22)and mass*(t2*g*g+t1*h)<=F(1,10**16):possible.append(n)
    ck(f'first_acceptable_order_{i}',possible and possible[0]==selector['degree'])
    ck(f'actual_axis_constraints_exact_{i}',all(P.T.is_exact_zero(Y[fields['nodes'][j]][0])for j in fields['axis_local_nodes']))
base={'status':'PASS','checks':{'a':True},'source_sha256':sha(source),'protocol_sha256':sha(HERE/'protocol.json'),'preflight_sha256':sha(HERE/'preflight.json')}
P.validate_review(base,base['source_sha256'],base['protocol_sha256'],base['preflight_sha256'])
for label,changes in [('empty',{}),('int',{'a':1}),('false',{'a':False}),('list',[True])]:
    try:P.validate_review({**base,'checks':changes},base['source_sha256'],base['protocol_sha256'],base['preflight_sha256']);rejected=False
    except ValueError:rejected=True
    ck('guard_'+label,rejected)
ck('no_root_claim',p['root_certified']is False)
out={**base,'checks':checks,'review_source_sha256':sha(Path(__file__)),
     'scope':'Independent adapter, exact physical/Duffy geometry moments, selectors and guard review; safe generic constitutive/moment/cache kernels and physical representation algebra reused. No actual six force/H integrals or64-site sums.'}
(HERE/'prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':'PASS','checks':len(checks),'review_sha256':sha(HERE/'prelaunch_independent_review.json')}))
