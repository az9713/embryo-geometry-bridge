"""Versioned exact isolation of actual certified M cell inside broader seed."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];BASE=Path(__file__).resolve().parent;OUT=BASE/'actual_M_cover_v2'
sys.path.insert(0,str(ROOT/'scripts'))
import validated_identified_set as e
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
if (OUT/'protocol.json').exists():raise RuntimeError('Versioned broad cover exists; retain rather than overwrite')
oldfile=ROOT/'evidence/validated_identified_set_actual_M.json';old=read(oldfile)
if old['status']!='PASS'or not all(v is True for v in old['checks'].values()):raise RuntimeError('Verified actualM provider application required')
decision=old['primary_decision'];target=e.serial_cell(decision['parameter_cell']);context=decision['inference_context']
foundation=ROOT/'evidence/validated_identified_set_actual.json';f=read(foundation)
seed=f['cover_ledger']['root_cell'];ambient=f['cover_ledger']['ambient_domain']
if not e.subset(target,seed):raise RuntimeError('Actual certified cell must lie inside declared broad seed')
protocol={'scope':'Isolate one actual certified cell in the existing broader finite seed; preserve all other finite leaves and unbounded complement',
    'engine_source_sha256':sha(ROOT/'scripts/validated_identified_set.py'),'generator_sha256':sha(__file__),
    'actualM_application_sha256':sha(oldfile),'foundation_seed_receipt_sha256':sha(foundation),
    'root_seed_cell':seed,'target_cell':target,'target_cell_sha256':e.cell_hash(target),
    'ambient_domain':ambient,'split_rule':'Sorted coordinates; cut exact target lower then upper endpoints, following the child containing target. Closed boundaries retained.',
    'common_delta_exact':'51/10','other_branches_unresolved':True,'full_physical_identified_set_complete':False,'new_equilibrium_solves':0}
save(OUT/'protocol.json',protocol)
ledger=e.CoverLedger(seed,ambient,uncovered_physics=f['cover_ledger']['uncovered_physics'])
leaf='0';splits=[]
for coordinate,(lower,upper)in e.cell(target).items():
    a,b=e.cell(ledger.nodes[leaf]['cell'])[coordinate]
    if lower>a:
        previous=leaf;left,right=ledger.subdivide(leaf,coordinate,lower);leaf=right
        splits.append({'parent':previous,'coordinate':coordinate,'at_exact':str(lower),'target_child':leaf})
    a,b=e.cell(ledger.nodes[leaf]['cell'])[coordinate]
    if upper<b:
        previous=leaf;left,right=ledger.subdivide(leaf,coordinate,upper);leaf=left
        splits.append({'parent':previous,'coordinate':coordinate,'at_exact':str(upper),'target_child':leaf})
if e.cell_hash(ledger.nodes[leaf]['cell'])!=e.cell_hash(target):raise RuntimeError('Isolation did not exactly recover certified cell')
leaves=[name for name,node in ledger.nodes.items()if not node['children']]
for name in leaves:
    output=e.classify_cell(context['observations'],context['discrepancies'],context['doses'],ledger.nodes[name]['cell'],
        contract_hash=context['model_observation_contract_sha256'],package=decision['uniform_forward_package'],operator=context['operator'])
    ledger.attach(name,output)
snapshot=ledger.export();replay=e.verify_cover_snapshot(snapshot)
checks={'exact_target_isolated':e.cell_hash(ledger.nodes[leaf]['cell'])==decision['parameter_cell_sha256'],
    'nine_exact_axis_cuts':len(splits)==9,'ten_finite_leaves_retained':len(leaves)==10,
    'only_certified_target_selected_OUT':ledger.nodes[leaf]['classification']['status']=='OUT'and sum(ledger.nodes[n]['classification']['status']=='OUT'for n in leaves)==1,
    'allother9_finite_leaves_UNRESOLVED':all(ledger.nodes[n]['classification']['status']=='UNRESOLVED'for n in leaves if n!=leaf),
    'all68_rows_every_leaf':all(ledger.nodes[n]['classification']['constraint_count']==68 for n in leaves),
    'no_tiny_certificate_reused_for_broader_leaf':all(not ledger.nodes[n]['classification']['uniform_constraints_covered']for n in leaves if n!=leaf),
    'certified_target_still67missing':sum(r['status']=='UNRESOLVED'for r in ledger.nodes[leaf]['classification']['constraints'])==67,
    'all_other_branches_still_UNRESOLVED':all(ledger.nodes[n]['classification']['all_branches_existential_parameter_status']=='UNRESOLVED'for n in leaves),
    'full_exact_leaf_union_equals_broad_seed':replay['exact_closed_leaf_union_equals_root_box'],
    'unbounded_and_zero_face_remainder_UNCOVERED':replay['ambient_remainder_status']=='UNCOVERED'and all(not d['lower_strict']for d in ambient.values()),
    'standalone_snapshot_replay_equal':replay==snapshot['verification'],
    'source_inputs_unchanged':sha(ROOT/'scripts/validated_identified_set.py')==protocol['engine_source_sha256']and sha(oldfile)==protocol['actualM_application_sha256']and sha(foundation)==protocol['foundation_seed_receipt_sha256']}
result={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'protocol_sha256':sha(OUT/'protocol.json'),
    'actualM_application_sha256':sha(oldfile),'target_leaf':leaf,'exact_splits':splits,'cover_ledger':snapshot,
    'scope':'One tiny actualM selected stationary-family cell excluded at DECLAREDdelta5.1; nineotherfiniteleaves,67otherconstraints,otherbranches,zero-faces/unboundedtails and remainingphysics remain unresolved. Not fullTheta or a noise estimate.'}
save(OUT/'results.json',result);save(ROOT/'evidence/validated_identified_set_actual_M_cover.json',result)
print(json.dumps({'status':result['status'],'checks':len(checks),'splits':len(splits),'leaves':len(leaves),'target_leaf':leaf,
    'outside_finite_leaves_unresolved':len(replay['unresolved_leaves']),'ambient_remainder':replay['ambient_remainder_status'],
    'receipt_sha256':sha(ROOT/'evidence/validated_identified_set_actual_M_cover.json')}),flush=True)
if not all(checks.values()):raise RuntimeError('Broad seed cover verification failed')
