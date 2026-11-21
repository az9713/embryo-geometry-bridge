"""Read-only verification and write-once design-probe manifest."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target as T

def main():
    dest=OUT/'delivery.json'
    if dest.exists():raise ValueError('write-once delivery already exists')
    sha=T.sha;checks={};p=T.read(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json')
    checks['frozen_target_and_local_core_unchanged']=all(sha(ROOT/k)==v for k,v in p['source_hashes'].items())
    checks['all147_target_inputs_unchanged']=len(p['input_hashes'])==147 and all(sha(ROOT/k)==v for k,v in p['input_hashes'].items())
    j=T.read(OUT/'polynomial_jet_tests.json');checks['all25_helper_checks_literal_true']=len(j['checks'])==25 and all(v is True for v in j['checks'].values())
    checks['helper_and_test_sources_match']=sha(OUT/'polynomial_jet.py')==j['helper_source_sha256'] and sha(OUT/'test_polynomial_jet.py')==j['test_source_sha256']
    checks['generic_AD_verifier_matches']=sha(OUT/'verify_derivatives.py')==j['generic_AD_verifier_sha256']
    v=T.read(OUT/'derivative_verification.json');checks['all12_independent_AD_checks_literal_true']=len(v['checks'])==12 and all(x is True for x in v['checks'].values())
    probe=T.read(OUT/'results.json');checks['probe_hash_bindings_match']=sha(OUT/'results.json')==v['probe_results_sha256']==j['probe_results_sha256'] and sha(OUT/'probe.py')==v['probe_source_sha256']==probe['source_sha256']
    target=T.Target(T.read(ROOT/p['records'][0]['point_cache_file']));pairs=set();slots=0
    for ids in list(target.a['triangles'])+list(target.a['surface_0_edges'])+list(target.a['surface_1_edges']):
        fs=[target.index[2*int(n)+c]for n in ids for c in(0,1)if 2*int(n)+c in target.index];slots+=len(fs)*(len(fs)+1)//2;pairs.update((min(i,j),max(i,j))for i in fs for j in fs)
    sparse=T.read(OUT/'sparsity.json');checks['connectivity_sparsity_independently_replayed']=slots==sparse['element_and_surface_local_uppertriangle_slots']==10952 and len(pairs)==sparse['assembled_global_uppertriangle_pattern']==6676 and 2*len(pairs)-sum(i==j for i,j in pairs)==sparse['assembled_global_full_symmetric_pattern']==12741
    checks['no_new_full_integral_or_equilibrium_result_claimed']=probe['checks']['no_field_integration_or_solve']is True
    if not all(checks.values()):raise AssertionError(checks)
    files=[x for x in OUT.iterdir()if x.is_file()and x.name!='delivery.json']
    output={'status':'PASS_BOUNDED_BRIDGE_DESIGN_COMPONENT','checks':checks,'file_hashes':{x.name:sha(x)for x in sorted(files)},'target_preparation_protocol_sha256':sha(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json'),'frozen_target_source_sha256':sha(ROOT/'scripts/exact_fe_target.py'),'scope':'Actual B-center local support/conditioning probe, exact low-order polynomial jet and dual-moment algebra, executable integration design. No full-field integration, highorder power16 coefficient expansion, equilibrium solve or exact-integral root certificate. Source of next pilot and all constitutive/moment/remainder inputs still require their own review.'}
    T.write(dest,output);print(json.dumps({'status':output['status'],'checks':len(checks),'receipt_sha256':sha(dest),'files':len(files)}))
if __name__=='__main__':main()
