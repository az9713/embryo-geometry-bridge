"""Finalize only the integration checkpoint, never the original five gates."""
from pathlib import Path
import hashlib,json,datetime
ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'runs/exact_integration_bridge_v1'
OUT=ROOT/'evidence/exact_integration_bridge_final_delivery.json'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    checks={};bindings={};reviews={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    def bind(path,expected=None):
        p=ROOT/path;h=sha(p)
        if expected is not None and h!=expected:raise AssertionError('binding '+path)
        if path in bindings and bindings[path]!=h:raise AssertionError('conflicting binding '+path)
        bindings[path]=h;return p
    for path in [
        'evidence/exact_fe_moments_independent_review.json',
        'evidence/exact_bulk_pilot_independent_review.json',
        'evidence/exact_bulk_assembly_B_independent_review.json',
        'evidence/exact_complete_point_B_independent_review.json',
        'evidence/exact_cache_consistency_independent_review.json',
        'evidence/exact_integration_documentation_review.json']:
        p=bind(path);d=read(p);ck(path+'_PASS',d['status'].startswith('PASS'))
        if isinstance(d.get('checks'),dict):ck(path+'_all_checks',bool(d['checks']) and all(v is True for v in d['checks'].values()))
        reviews[path]={'status':d['status'],'sha256':sha(p),'checks':len(d.get('checks',{}))}
    protocols=[
        RUN/'target_preparation/protocol.json',RUN/'bulk_pilot/protocol_v2.json',
        RUN/'bulk_assembly_B/protocol_v2.json',RUN/'surface_pilot_preparation/protocol.json',
        RUN/'surface_assembly_B/protocol.json',RUN/'cache_consistency_preparation/protocol.json']
    for path in protocols:
        bind(path.relative_to(ROOT).as_posix());p=read(path)
        for key in ('source_hashes','input_hashes','source_and_input_hashes','target_preparation_input_hashes'):
            for name,h in p.get(key,{}).items():bind(name,h)
    ck('all_numerical_protocol_bindings_still_current',True)
    point=read(bind('runs/exact_integration_bridge_v1/complete_point_B/results.json'))
    ck('all_144_bulk_and_24_surface_at_same_point',point['all144bulk_and24surface_contributions'] and point['free_dofs']==611 and point['H_upper_slots']==6676)
    ck('point_negative_selfmap_certificate',point['old_map_center_test']['center_leaves_old_cube_certified'] and point['old_map_center_test']['not_equilibrium_nonexistence'])
    ck('no_uniform_root_claim',not point['uniform_a_or_b_established'] and not point['new_integral_root_certified'])
    bind('runs/exact_integration_bridge_v1/complete_point_B_independent_review/coordinate_witness.json')
    surfaces=read(bind('runs/exact_integration_bridge_v1/surface_pilot_independent_replay/results.json'))
    ck('six_surface_independent_replays_PASS',surfaces['status']=='PASS_SIX_INDEPENDENT_POINT_SURFACE_REPLAYS' and len(surfaces['outputs'])==6)
    browser=read(bind('evidence/readiness_snapshot_verification.json'))
    ck('browser_105_pages_210_viewports',browser['status']=='PASS' and browser['html_pages']==105 and browser['page_viewport_checks']==210 and all(browser['checks'].values()))
    for path,h in browser['snapshot_sha256'].items():bind(path,h)
    bind('evidence/report_browser_verification.json',browser['browser_receipt_sha256'])
    raw=read(bind('evidence/raw_prose_math_scan.json'))
    ck('raw_prose_scan_PASS_105',raw['status']=='PASS' and raw.get('pages',raw.get('html_pages'))==105)
    parent=read(bind('evidence/parent_index_checkpoint_verification.json'))
    ck('parent_current_landing_two_viewports_PASS',parent['status']=='PASS' and len(parent['checks'])==2 and all(v['pass'] for v in parent['checks']) and sha(ROOT.parent/'docs/index.html')==parent['page_sha256'] and sha(ROOT.parent/'scripts/build_dossier.py')==parent['builder_sha256'])
    space=read(bind('evidence/exact_integration_space_update.json'))
    ck('Space_exact_118_blocks_and_native_images',space['status']=='PASS_SPACE_EXACT_READBACK' and space['blocks']==118 and space['sequence']==16 and all(space['checks'].values()))
    for name in ['root_visual_review.json','point_map_root_visual_review.json']:
        d=read(bind('runs/exact_integration_bridge_v1/visuals/'+name));ck('visual_'+name,d['status']=='PASS_VISUAL_INSPECTION' and all(d['checks'].values()))
    for name,h in [('bulk-element-integration-defects.png','baf3aa804034ab89d0f0b55fb6710ab19a7f02d6d57d235d8b2fcc6ed469a453'),('complete-point-integral-map-failure.png','93e0a6aadfe26809a4e38dca40b3bc0e3fc96c8d638dbbdbd08948b05f37e02c')]:bind('runs/exact_integration_bridge_v1/visuals/'+name,h)
    ck('native_upload_bytes_match_local_PNGs',set(v['sha256'] for v in space['uploads'])=={'baf3aa804034ab89d0f0b55fb6710ab19a7f02d6d57d235d8b2fcc6ed469a453','93e0a6aadfe26809a4e38dca40b3bc0e3fc96c8d638dbbdbd08948b05f37e02c'})
    ck('original_PLAN_unchanged',sha(bind('PLAN.md'))=='993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500')
    gates=read(bind('evidence/gates.json'));ck('original_five_gate_goal_explicitly_incomplete',gates['goal_achieved'] is False)
    bind(Path(__file__).relative_to(ROOT).as_posix())
    result={'status':'PASS_CHECKPOINT_DELIVERY_NOT_GOAL_COMPLETE','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'reviews':reviews,'bound_files':bindings,'independent_reintegration_coverage':{'bulk_elements':6,'bulk_total':144,'surface_edges':6,'surface_total':24},'scope':'Complete point B integral/cache discrepancy and negative direct self-map certificate, component preparations, documentation/Space/browser delivery. This is not a uniform integral-root or continuum theorem, biological validation, novelty or completion of the original five gates.','remaining_original_requirements':['Revised map/preconditioner/center/domain proof with uniform parameter/predictor/corrector errors.','Full physical material/geometry identified set and stable loading paths.','Matched source geometry/reference/ROI/time and defensible discrepancy information.','Revised frozen larger-region predictions, continuum accuracy and contribution/novelty assessment.'],'external_contact_sent':False,'GitHub_publication_or_new_remote':False,'goal_achieved':False}
    if OUT.exists():raise ValueError('write-once final delivery exists')
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':result['status'],'checks':len(checks),'bound_files':len(bindings),'sha256':sha(OUT)}),flush=True)
if __name__=='__main__':main()
