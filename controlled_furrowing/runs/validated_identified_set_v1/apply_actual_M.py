"""Consume the frozen independently reviewed ONE-channel actual-M provider."""
from pathlib import Path
from fractions import Fraction as Q
import copy,csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'actual_M'
sys.path.insert(0,str(ROOT/'scripts'))
import validated_identified_set as engine
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
source=ROOT/'scripts/validated_identified_set.py'
combined_path=ROOT/'evidence/parametric_actual_M_independent_review.json'
root_path=ROOT/'runs/parametric_root_cells_v1/actual_M/root_certificate.json'
observer_path=root_path.with_name('observer_certificate.json')
contract_path=root_path.with_name('model_observation_contract.json')
if not contract_path.exists():
    # The producer may use the shorter filename; choose only by its reviewed hash.
    wanted=read(root_path)['model_observation_contract_sha256']
    matches=[p for p in root_path.parent.glob('*.json')if sha(p)==wanted]
    if len(matches)!=1:raise RuntimeError('Unique actual reviewed contract file required')
    contract_path=matches[0]
if (OUT/'protocol.json').exists():raise RuntimeError('Actual provider application already recorded; preserve outputs')
combined=read(combined_path);r=read(root_path);o=read(observer_path)
if combined['status']!='PASS':raise RuntimeError('Final independent combined review required')
for section,path in [('root_audit',root_path),('observer_audit',observer_path)]:
    audit=combined[section]
    if audit['status']!='PASS'or not audit['checks']or not all(v is True for v in audit['checks'].values())or audit['certificate_sha256']!=sha(path):raise RuntimeError('Actual independent audit binding does not pass')
if sha(combined_path)!='4b03db70583be266a8d4ae581f19c8db0099bc028ddd6af82a5c990a91d8e600':raise RuntimeError('Expected frozen final combined review version')
test=read(ROOT/'evidence/validated_identified_set_synthetic.json')
if test['status']!='PASS'or test['source_sha256']!=sha(source):raise RuntimeError('Current engine synthetic verification required')
with (ROOT/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig')as f:rows=list(csv.DictReader(f))
observed={row['record_id']:{k:row[k]for k in engine.OBSERVABLES}for row in rows};doses={row['record_id']:row['activation_magnitude']for row in rows}
rid='optoGEF2_small_col_M';channel='apical_max_shortening_um';cell=r['parameter_cell']
if r['record_id']!=rid or o['record_id']!=rid or Q(r['dose_exact'])!=Q(doses[rid])or Q(o['dose_exact'])!=Q(doses[rid]):raise RuntimeError('Actual sourceM dose, never diagnosticI1, required')
if set(o['endpoint_intervals'])!={channel}:raise RuntimeError('Application scope is precisely one reviewed apical channel')
if r['model_observation_contract_sha256']!=sha(contract_path)or o['model_observation_contract_sha256']!=sha(contract_path):raise RuntimeError('Actual static model/observation contract must bind')
protocol={'stage':'Consume one independently reviewed actualM uniform source-channel provider against ALL68 constraints',
    'engine_source_sha256':sha(source),'generator_sha256':sha(__file__),'combined_independent_review_sha256':sha(combined_path),
    'root_certificate_sha256':sha(root_path),'observer_certificate_sha256':sha(observer_path),'contract_sha256':sha(contract_path),
    'parameter_cell_sha256':engine.cell_hash(cell),'parameter_cell':cell,
    'calibration_csv_sha256':sha(ROOT/'data/derived/calibration_small_18um.csv'),
    'primary_common_delta_exact':'51/10','delta_family_rule':'Use0,5,5.1,cutoff_lower−1e−9,lower,(lower+upper)/2,upper,upper+1e−9,5.2,6; exact rationals; these are declared tolerances, NOT noise estimates',
    'branch_law_semantics':engine.STATIONARY_CLOSURE,'same_source_observations':68,'only_source_channel':[rid,channel],
    'other67_constraints_missing':True,'all_other_branches_unresolved':True,'full_physical_identified_set_complete':False,
    'no_new_equilibrium_solves':True}
save(OUT/'protocol.json',protocol)
def ref(path):return {'file':path.relative_to(ROOT).as_posix(),'sha256':sha(path)}
adapters={}
for section,name in [('root_audit','root_review_adapter.json'),('observer_audit','observer_review_adapter.json')]:
    value=copy.deepcopy(combined[section]);value['transparent_derivation']={'combined_review':ref(combined_path),'section':section,'new_mathematical_claims':False}
    path=OUT/name;save(path,value);adapters[section]=path
package={'kind':'AUDITED_UNIFORM_CELL_FORWARD_V1','model_observation_contract_sha256':sha(contract_path),
    'branch_law_semantics':engine.STATIONARY_CLOSURE,'records':[{'record_id':rid,'root_certificate':ref(root_path),
        'root_review':ref(adapters['root_audit']),'observer_certificate':ref(observer_path),'observer_review':ref(adapters['observer_audit'])}],
    'scope':'One actualM apical endpoint only; root/observer audit adapters copy frozen independent review sections without new proof assertions'}
save(OUT/'forward_package.json',package)
delta={k:Q('5.1')for k in engine.OBSERVABLES}
primary=engine.classify_cell(observed,delta,doses,cell,contract_hash=sha(contract_path),package=package)
save(OUT/'delta_5p1_decision.json',primary)
lo,hi=map(Q,o['endpoint_intervals'][channel]);y=Q(observed[rid][channel]);lower,upper=engine.interval_residual(y,(lo,hi))
step=Q(1,10**9);tolerances=sorted(set([Q(0),Q(5),Q('5.1'),max(Q(0),lower-step),lower,(lower+upper)/2,upper,upper+step,Q('5.2'),Q(6)]))
profiles=[]
for value in tolerances:
    decision=engine.classify_cell(observed,{k:value for k in engine.OBSERVABLES},doses,cell,contract_hash=sha(contract_path),package=package)
    profiles.append({'common_delta_exact':str(value),'common_delta_display':float(value),'decision':decision})
domain={k:{'lower':0,'lower_strict':False,'upper':None}for k in cell}
ledger=engine.CoverLedger(cell,domain,uncovered_physics=['Other root branches','Other material and geometry cells including zero faces/unbounded tails','Compression-ratio variation beyondfixed5000','Common stable/quasistatic path and source temporal-peak observation law'])
ledger.attach('0',primary);export=ledger.export();engine.verify_cover_snapshot(export)
checks={'all17_records_all68_constraints_retained':primary['record_count']==17 and primary['constraint_count']==68,
    'exactly_one_reviewed_actual_source_channel':primary['uniform_constraints_covered']==[{'record_id':rid,'observable':channel}],
    'delta_5p1_selected_local_branch_OUT':primary['status']=='OUT'and primary['selected_branch_status']=='OUT',
    'other_branches_UNRESOLVED':primary['all_branches_existential_parameter_status']=='UNRESOLVED',
    'other67_rows_still_UNRESOLVED':sum(row['status']=='UNRESOLVED'for row in primary['constraints'])==67,
    'no_quasistatic_or_stable_path_claim':primary['common_path_or_energetic_stability_certified']is False,
    'parameter_cell_exact_binding':primary['parameter_cell_sha256']==r['parameter_cell_sha256']==o['parameter_cell_sha256'],
    'whole_delta_family_has68_rows':all(p['decision']['constraint_count']==68 for p in profiles),
    'strict_cutoff_decision_correct':all(p['decision']['status']==('OUT'if Q(p['common_delta_exact'])<lower else'UNRESOLVED')for p in profiles),
    'never_full_IN_from_one_row':all(p['decision']['status']!='IN'for p in profiles),
    'unbounded_and_zero_faces_retained':export['verification']['ambient_remainder_status']=='UNCOVERED'and all(not d['lower_strict']for d in export['ambient_domain'].values()),
    'standalone_cover_replay_matches':engine.verify_cover_snapshot(export)==export['verification'],
    'root_review_adapter_exact_section':all(read(adapters['root_audit'])[k]==v for k,v in combined['root_audit'].items()),
    'observer_review_adapter_exact_section':all(read(adapters['observer_audit'])[k]==v for k,v in combined['observer_audit'].items()),
    'all_frozen_inputs_unchanged':sha(source)==protocol['engine_source_sha256']and sha(combined_path)==protocol['combined_independent_review_sha256']and sha(root_path)==protocol['root_certificate_sha256']and sha(observer_path)==protocol['observer_certificate_sha256']and sha(contract_path)==protocol['contract_sha256']}
result={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'protocol_sha256':sha(OUT/'protocol.json'),
    'forward_package_sha256':sha(OUT/'forward_package.json'),'primary_decision':primary,'conditional_delta_family':profiles,
    'exact_one_channel_cutoffs':{'necessary_residual_lower_exact':str(lower),'sufficient_row_residual_upper_exact':str(upper),'display':[float(lower),float(upper)],
        'symbolic_family':'For arbitrary nonnegative other-channel allowances: delta_apical<lower excludes this local stationary family. delta_apical>=upper satisfies ONLY this row; other67rows remain unresolved.'},
    'cover_ledger':export,'source_channel_provenance':{'record_id':rid,'dose_exact':str(Q(doses[rid])),'observed_apical_um_exact':str(y),'contract_file':contract_path.relative_to(ROOT).as_posix()},
    'engine_source_sha256':sha(source),'generator_sha256':sha(__file__),'combined_independent_review_sha256':sha(combined_path),
    'scope':'First actual one-channel uniform parameter-cell selected-root-family exclusion; other67constraints and allotherbranches/unboundedphysicaldomain remain unresolved. Declareddelta is not measurednoise; no novelty or fullidentifiedset claimed.'}
save(OUT/'results.json',result);save(ROOT/'evidence/validated_identified_set_actual_M.json',result)
print(json.dumps({'status':result['status'],'checks':len(checks),'selected_status':primary['status'],'covered':len(primary['uniform_constraints_covered']),
    'missing':sum(row['status']=='UNRESOLVED'for row in primary['constraints']),'cutoffs':[float(lower),float(upper)],'receipt_sha256':sha(ROOT/'evidence/validated_identified_set_actual_M.json')}),flush=True)
if not all(checks.values()):raise RuntimeError('Actual uniform provider application failed')
