"""Independently review the fresh F/K replay's target and center semantics."""
from pathlib import Path
import ast,hashlib,importlib,json

ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs/cohort_recentered_root_review_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))

def main():
    out=RUN/'prelaunch_review.json'
    if out.exists():raise RuntimeError('preserve concrete independent replay review')
    m=importlib.import_module('review_recentered_roots');p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json')
    producer=read(ROOT/'runs/cohort_recentered_roots_v1/protocol.json');guards=m.guards()
    source=Path(m.__file__);text=source.read_text(encoding='utf-8')
    calls=[n.func.attr for n in ast.walk(ast.parse(text)) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)]
    checks={'review_source_protocol_preflight_bound':sha(source)==pre['review_source_sha256']==p['review_source_sha256'] and pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'preflight_all_boolean_PASS':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'all_actual_input_and_source_hashes_current':all(sha(ROOT/k)==h for k,h in p['file_hashes'].items()),
        'all29_lightweight_guards_reexecuted':bool(guards) and all(v is True for v in guards.values()),
        'guards_and_checks_native_JSON_serializable':bool(json.dumps(guards)),
        'fresh192bit_declared':p['precision_bits']==192,
        'sameC_and_physical_contract':p['parameter_cell_sha256']=='64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf' and p['model_observation_contract_sha256']=='2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82',
        'only_two_failed_seed_records_replayed':{r['record_id'] for r in p['records']}=={'optoGEF2_small_col_F','optoGEF2_small_col_K'},
        'no_producer_point_or_cell_routines_called':'build_point' not in calls and 'evaluate_cell' not in calls,
        'old_RF_proposal_not_misrepresented_as_independently_recomputed':p['old_RF_assembly_repeated'] is False,
        'old_u0_archive_and_new_u1_cache_checked_separately':"checks['old_matrix_archive_belongs_to_seed_u0_not_new_u1']" in text and "checks['new_cache_center_u1_matches_auxiliary_state']" in text,
        'strict_new_center_root_proof_not_inferred_from_recipe':"'uniform_strict_selfmap':bool(lhs<r)" in text and "'strict_Banach_acceptance':strict_acceptance(eta,Z,r)" in text,
        'no_replay_already_started':not (RUN/'records').exists() and not (ROOT/'evidence/cohort_recentered_roots_independent_review.json').exists()}
    for original,entry in zip(p['producer_records'],p['records']):
        root=read(ROOT/entry['file']);report=m.lineage(original,root,producer)
        checks[original['record_id']+'_all_lineage_checks_replayed']=report['status']=='PASS' and all(v is True for v in report['checks'].values())
        checks[original['record_id']+'_actual_root_payload_and_protocol_bound']=sha(ROOT/entry['file'])==entry['sha256'] and root['protocol_sha256']==p['producer_protocol_sha256']
    m.bindings(p);checks['frozen_environment_guard_reexecuted']=True
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'review_source_sha256':sha(source),
        'protocol_sha256':sha(RUN/'protocol.json'),'preflight_sha256':sha(RUN/'preflight.json'),'parent_review_source_sha256':sha(__file__),
        'mathematical_review':'The independent force/Hessian/parameter-column/Taylor-row assembly is retained from the already reviewed separate implementation, evaluated freshly at u1. Its source-context guards correctly bind u0 to the old float archive and u1 to the auxiliary cache, while bytewise physical-potential/R/S/C identity is checked independently. The Banach theorem only needs a finite chosen proof center; it does not assume the midpoint proposal is an exact Newton iterate or a root. Normal corrector radius remains separate from translation.',
        'scope':'Prelaunch source/target/lineage review. Fresh numerical replay and separate observer offset/geometry application remain unexecuted by this receipt.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':len(checks),'guards':len(guards),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
