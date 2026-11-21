"""Parent independent preparation/guard checks; no 611 assembly."""
from pathlib import Path
from fractions import Fraction as Q
import ast, hashlib, importlib.util, json, sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_cached_contraction_B_common as producer
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
checks={}
def ck(k,v):
    checks[k]=bool(v)
    assert v,k
source=ROOT/'scripts/exact_cached_contraction_B_common.py'
ck('literal_source_hash',sha(source)=='3d509ae81ec3e3349b9560a0c0431ed767cc2981b361f9730817083760dc2629')
ck('literal_protocol_hash',sha(HERE/'protocol.json')=='f14bdde08b79b76fe214d02d743b1fde3991cb9d7daf7c093ce6501ba65fd998')
ck('literal_preflight_hash',sha(HERE/'preflight.json')=='be6bb5888211fecef2a43eb58b8fc9be96ff32690438d49c1abe0363ee78a91c')
p=read(HERE/'protocol.json');pre=read(HERE/'preflight.json')
ck('preflight_literal_pass',pre['status']=='PASS_PREPARATION_ONLY'and pre['checks']and all(v is True for v in pre['checks'].values()))
for path,h in p['source_and_input_hashes'].items():ck('binding_'+path,sha(ROOT/path)==h)
ledger=read(HERE/'exact_containment.json')
ck('containment_ledger_hash',sha(HERE/'exact_containment.json')==p['exact_containment_sha256'])
ck('all650_global_coordinates',len(ledger['rows'])==650 and [v['global_dof']for v in ledger['rows']]==list(range(650)))
ck('all3250_declared_containments_literal',all(set(v['checks'])=={'old_cube','candidate_cube','translation_path','swept_cube','canonical_geometry_box'}and all(x is True for x in v['checks'].values())for v in ledger['rows']))
retube=ROOT/'runs/exact_integration_bridge_v1/retube_B_preparation'
row=read(retube/'protocol.json')['selected_record']
with np.load(retube/'candidate_vectors.npz',allow_pickle=False)as z:
    decrement=[Q.from_float(float(a))-Q.from_float(float(b))for a,b in zip(z['u0'],z['u1'])]
ck('independent_actual_binary_decrement_norm',max(map(abs,decrement))==Q(p['actual_applied_decrement_norm_exact']))
ck('exact_common_radius',Q(p['common_radius_exact'])==Q(p['old_radius_exact'])+max(map(abs,decrement)))
ck('same_parameter_cell_and_physical_input',p['parameter_cell']==row['parameter_cell']and p['base_parameters_exact']==row['base_parameters_exact']and p['dose_exact']==row['dose_exact'])
ck('exact_gain_times_dose',Q(p['base_strength_exact'])==Q(p['base_parameters_exact']['common_gain'])*Q(p['dose_exact']))
ck('build_point_absolute_predictor_safe',Path(str((ROOT/p['build_point_parameters']['predictor_file']).resolve())).is_absolute())
ck('build_point_cache_original',p['build_point_parameters']['point_certificate_file']==row['point_cache_file'])
# Guard mutations must be rejected before heavy work. No guard result is used
# as a scientific contraction premise.
base={'status':'PASS','checks':{'one':True},'source_sha256':sha(source),'protocol_sha256':sha(HERE/'protocol.json'),'preflight_sha256':sha(HERE/'preflight.json')}
producer.validate_review(base)
for name,alter in [('empty',{'checks':{}}),('nonboolean',{'checks':{'one':1}}),('false',{'checks':{'one':False}}),('wrong_source',{'source_sha256':'0'*64}),('wrong_protocol',{'protocol_sha256':'0'*64}),('wrong_preflight',{'preflight_sha256':'0'*64}),('wrong_status',{'status':'PASS_UNREVIEWED'})]:
    try:producer.validate_review({**base,**alter});rejected=False
    except ValueError:rejected=True
    ck('guard_rejects_'+name,rejected)
tree=ast.parse(source.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='run')
names=[n.func.id if isinstance(n.func,ast.Name)else n.func.attr for n in ast.walk(run)if isinstance(n,ast.Call)and isinstance(n.func,(ast.Name,ast.Attribute))]
ck('run_has_gate_before_point_assembly',names.index('validate_review')<names.index('build_point'))
ck('no_integral_root_claim',p['actual_root_or_selfmap_claim']is False and p['integral_H_or_PDE_bound']is False)
# Independent signed 3x3 row-majorant example.
B=[[Q(2),Q(-1),Q(0)],[Q(-3),Q(1),Q(2)]]
D=[[Q(-4),Q(2)],[Q(2),Q(1)]]
H=[[sum(B[a][i]*D[a][b]*B[b][j]for a in range(2)for b in range(2))for j in range(3)]for i in range(3)]
bound=[sum(abs(B[a][i])*sum(abs(D[a][b])*sum(map(abs,B[b]))for b in range(2))for a in range(2))for i in range(3)]
ck('signed_local_matrix_row_majorant',all(sum(map(abs,H[i]))<=bound[i]for i in range(3)))
out={**base,'checks':checks,'scope':'Independent input, exact-radius, preparation and launch-guard review. Previously reviewed dense P/core analytic kernels reused; no 611 matrix or variation run and no root/integral claim.',
     'review_source_sha256':sha(Path(__file__))}
(HERE/'prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':'PASS','checks':len(checks),'review_sha256':sha(HERE/'prelaunch_independent_review.json')}))
