"""Bind the two preparation-only reviews; does not launch point computation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,datetime
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    checks={};bindings={};receipts={}
    for folder,proto,pre in [('bulk_assembly_B_safe_v2','protocol_v2.json','preflight_v2.json'),('surface_assembly_B_safe_v2','protocol.json','preflight.json')]:
        b=ROOT/'runs/exact_integration_bridge_v1'/folder;rp=b/'prelaunch_independent_review.json';v=read(rp);p=read(b/proto)
        checks[folder+'_allchecks_true']=v['status']=='PASS'and bool(v['checks'])and all(x is True for x in v['checks'].values())
        checks[folder+'_launchbindings_current']=v['source_sha256']==sha(ROOT/('scripts/exact_bulk_bridge_assembly_safe.py'if folder.startswith('bulk')else'scripts/exact_surface_bridge_assembly_safe.py'))and v['protocol_sha256']==sha(b/proto)and v['preflight_sha256']==sha(b/pre)
        checks[folder+'_allinputhashes_current']=all(sha(ROOT/k)==x for k,x in p['source_and_input_hashes'].items())
        receipts[folder]={'sha256':sha(rp),'checks':len(v['checks'])};bindings[rp.relative_to(ROOT).as_posix()]=sha(rp)
    record=read(ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json')
    safe=read(ROOT/'runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json');row=next(x for x in safe['records']if x['record_id']==record['record_id'])
    checks['original_cache_proofcenter_predictor_C_dose_unchanged']=all(row[k]==record[k]for k in ['point_cache_file','predictor_file','dose_exact','parameter_cell','parameter_cell_sha256'])and all(row[k.removesuffix('_file')+'_sha256']==record['artifact_hashes'][record[k]]==sha(ROOT/record[k])for k in ['point_cache_file','predictor_file'])
    checks['original_exact_target_strength_unchanged']=Q(row['dose_exact'])*Q(row['base_parameters_exact']['common_gain'])==Q(record['base_strength_exact'])
    assert all(v is True for v in checks.values()),checks
    out={'status':'PASS_SAFE_POINT_REPLAY_PRELAUNCH_DELIVERY','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'launch_receipts':receipts,'input_hashes':bindings,'common_review_sha256':sha(HERE/'common_review.json'),'report_sha256':sha(HERE/'REVIEW.md'),'review_source_sha256':sha(__file__),'old_point_claims_restored':False,'reviewer_actual_field_integrations_run':False,'scope':'Fresh point computation may launch; subsequent independent local replay and full reduction remain required. No uniform integration bound, new root or original-five-gate completion.'}
    p=HERE/'delivery.json';assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'sha256':sha(p)}))
if __name__=='__main__':main()
