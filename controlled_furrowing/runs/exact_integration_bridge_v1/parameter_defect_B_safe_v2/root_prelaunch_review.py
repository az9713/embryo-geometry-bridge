"""Independent algebra/guard review before polynomial-only evaluation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_parameter_defect_B_safe as P
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
checks={}
def ck(k,v):
    checks[k]=bool(v)
    assert v,k
source=ROOT/'scripts/exact_parameter_defect_B_safe.py'
ck('source_hash_literal',sha(source)=='82abb15cf7546ccd9f9b9819f21bc39ca939ec68d86c48efc388179cfd805a1b')
ck('protocol_hash_literal',sha(HERE/'protocol.json')=='729dd66d1549e0faad25333f8d726d5a7ddce153ccc7acd9cb6cde7e6c5a982c')
ck('preflight_hash_literal',sha(HERE/'preflight.json')=='1b3363d6db822d923282086d6232f32523b7c4dde5dae4cb5ba57c63eebc17e9')
p=read(HERE/'protocol.json');pre=read(HERE/'preflight.json')
ck('preflight_nonempty_literal_pass',pre['status']=='PASS_PREPARATION_ONLY'and pre['checks']and all(v is True for v in pre['checks'].values()))
for name,h in p['source_and_input_hashes'].items():ck('input_'+name,sha(ROOT/name)==h)
order=['apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain']
ck('parameter_scatter_order',p['parameter_order']==order==list(P.ORDER))
ck('same_C_and_exact_base_dose',p['base_parameters_exact']==read(P.POINT)['base_parameters_exact']and p['parameter_cell']==read(P.RECORD)['parameter_cell']and p['dose_exact']==read(P.RECORD)['dose_exact'])
ck('exactly24_distinct_surface_edges',len(p['surface_edges'])==24 and {(v['face'],v['edge'])for v in p['surface_edges']}=={(a,b)for a in [0,1]for b in range(12)})
ck('surface_divisors_nonzero',all(Q(v['T_minus_K_exact'])==Q(v['base_T_exact'])-Q(v['base_K_exact'])!=0 for v in p['surface_edges']))
ctx.prec=640;P.core.A.cache_clear()
# Independent exact algebra at three different rational unit-tangent strains,
# tested with both signs of T-K. No biological element is integrated.
for case,(a,b,c)in enumerate([(Q(5,4),Q(3,5),Q(4,5)),(Q(7,6),Q(5,13),Q(12,13)),(Q(9,8),Q(-8,17),Q(15,17))]):
    poly=[a*(b*b+c*c),a*a*b,a*a*c];root=[Q(1),a*b,a*c]
    for K,T in [(Q(3),Q(7)),(Q(3),Q(0))]:
        force=[K*x+(T-K)*y for x,y in zip(poly,root)]
        split=P.split_surface(list(map(P.A,force)),list(map(P.A,poly)),K,T)
        ck(f'exact_surface_T_derivative_{case}_{T}',all(v.contains(P.A(q))for v,q in zip(split['partial_T'],root)))
        ck(f'exact_surface_K_derivative_{case}_{T}',all(v.contains(P.A(q))for v,q in zip(split['partial_K'],[x-y for x,y in zip(poly,root)])))
with np.load(ROOT/p['predictor_file'],allow_pickle=False)as z:S=z['S']
deltas=[(Q(p['parameter_cell'][k][0])-Q(p['base_parameters_exact'][k]),Q(p['parameter_cell'][k][1])-Q(p['base_parameters_exact'][k]))for k in order]
independent_hulls=[]
for row in S:
    pairs=[sorted([Q.from_float(float(s))*lo,Q.from_float(float(s))*hi])for s,(lo,hi)in zip(row,deltas)]
    independent_hulls.append((sum(v[0]for v in pairs),sum(v[1]for v in pairs)))
m=max(max(abs(a),abs(b))for a,b in independent_hulls)
ck('independent_all611_exact_predictor_hulls',len(independent_hulls)==611 and m==Q(p['predictor_displacement_norm_upper_exact']))
box=read(ROOT/p['common_geometry_boxes_file'])
ck('independent_hulls_match_common_geometry',independent_hulls==[tuple(map(Q,v))for v in box['predictor_hulls_exact']])
ck('independent_parameter_deltas',deltas==[tuple(map(Q,v))for v in p['parameter_deltas_exact']])
target=P.T.Target(read(ROOT/p['point_cache_file']))
ck('same24_binary_active_elements',p['bulk_active_elements']==[i for i in range(144)if Q.from_float(float(target.a['active'][i,0]))==1]and len(p['bulk_active_elements'])==24)
base={'status':'PASS','checks':{'a':True},'source_sha256':sha(source),'protocol_sha256':sha(HERE/'protocol.json'),'preflight_sha256':sha(HERE/'preflight.json')}
P.review_gate(base,base['source_sha256'],base['protocol_sha256'],base['preflight_sha256'])
for label,alter in [('empty',{}),('list',[True]),('int',{'a':1}),('false',{'a':False})]:
    try:P.review_gate({**base,'checks':alter},base['source_sha256'],base['protocol_sha256'],base['preflight_sha256']);rejected=False
    except ValueError:rejected=True
    ck('guard_'+label,rejected)
ck('no_root_or_V_claim',p['V_established']is False and p['exact_integral_root_certified']is False)
out={**base,'checks':checks,'review_source_sha256':sha(Path(__file__)),
     'scope':'Independent algebra, exact predictor/delta and source/guard prelaunch review; previously reviewed moment/maps and SAFE point packets reused. No actual B-field moments or cache sums executed.'}
(HERE/'prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':'PASS','checks':len(checks),'review_sha256':sha(HERE/'prelaunch_independent_review.json')}))
